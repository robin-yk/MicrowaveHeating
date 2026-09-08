"""Full vector curl-curl Maxwell solve with a passive lumped coax-probe port.

Peak phasors, exp(-i omega t). S11 is a complex amplitude ratio. The ideal
50-ohm port is a model boundary; it is not a resolved coaxial cable or antenna.
"""
from dataclasses import dataclass
import warnings
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, MatrixRankWarning
from .materials import EPS0, MU0


@dataclass
class EMResult:
    voltage: np.ndarray
    e_vector: np.ndarray
    e2: np.ndarray
    q_sample: np.ndarray
    q_quartz: np.ndarray
    incident_w: float
    reflected_w: float
    wall_w: float
    aperture_w: float
    sample_w: float
    quartz_w: float
    s11: complex
    power_residual: float
    linear_residual: float
    frequency_hz: float


class Maxwell:
    def __init__(self, grid, cfg):
        self.grid, self.cfg = grid, cfg
        self.K = (grid.C.T @ sparse.diags(grid.face_metric/MU0) @ grid.C).tocsr()
        self.wall_model = cfg["em"]["wall_model"]
        self.aperture_mass = np.zeros(grid.ne)
        if cfg["em"].get("tube_apertures",False):
            ro=cfg["geometry"]["tube_inner_diameter_m"]/2+cfg["geometry"]["tube_wall_m"]
            cx,cy=cfg["geometry"].get("tube_offset_m",[0,0])
            for cell,axis,_,area,_,face in grid.boundary:
                if int(axis)!=2:
                    continue
                cell=int(cell); fraction=0.
                for a in (np.arange(12)+.5)/12-.5:
                    for b in (np.arange(12)+.5)/12-.5:
                        xx,yy=grid.xyz[cell,:2]+np.array([a,b])*grid.cell_lengths[cell,:2]
                        fraction += ((xx-cx)**2+(yy-cy)**2<ro*ro)/144
                for edge in grid.face_edges[int(face)]:
                    direction=0 if edge<=grid.edges[0].max() else 1
                    self.aperture_mass[edge] += fraction*area/(2*grid.cell_lengths[cell,direction]**2)
        self.wall_mass = np.maximum(grid.wall_mass-self.aperture_mass,0)
        self.free = np.flatnonzero(grid.used_edges & (self.wall_mass<1e-15)) if self.wall_model == "pec" else np.flatnonzero(grid.used_edges)
        self.Kfree = self.K[self.free][:,self.free].tocsc()
        self.port = self._port()
        self.portfree = self.port[self.free]
        self.port_matrix = sparse.csc_matrix(self.portfree[:,None]) @ sparse.csc_matrix(self.portfree[None,:])

    def _port(self):
        grid, p = self.grid, self.cfg["port"]
        # A radial x-directed probe starts at the +x cavity wall at y≈0.
        # Its voltage is the integral along a configurable number of interior edges.
        j = int(np.argmin(abs(grid.axes[1]-p.get("y_m",0))))
        k = int(np.argmin(abs(grid.axes[2]-p.get("z_m",0))))
        available = [int(e) for e in grid.edges[0][:,j,k] if e in set(self.free)]
        if not available:
            raise ValueError("Port position has no interior x-directed edges")
        requested = p["probe_length_m"]
        b = np.zeros(grid.ne); accumulated = 0.
        for e in available[::-1]:
            i = int(np.argwhere(grid.edges[0] == e)[0,0])
            included = min(grid.dx[0][i], requested-accumulated)
            b[e] = included/grid.dx[0][i]
            accumulated += included
            if accumulated >= requested:
                break
        self.effective_probe_length = accumulated
        return b

    def solve(self, properties, frequency_hz, incident_w):
        if incident_w < 0 or frequency_hz <= 0:
            raise ValueError("Frequency must be positive and incident power nonnegative")
        grid = self.grid; w = 2*np.pi*frequency_hz
        mr = EPS0*grid.mass(properties["ep"])
        bed_sigma = properties.get("bed_sigma",np.zeros(grid.n))
        mi = EPS0*grid.mass(properties["epp"])+grid.mass(bed_sigma)/w
        admittance = 0j
        if self.wall_model == "impedance":
            sigma = self.cfg["em"]["wall_conductivity_s_m"]
            zs = (1-1j)*np.sqrt(w*MU0/(2*sigma))
            admittance = 1/zs
        g = 1/self.cfg["port"]["resistance_ohm"]
        A = self.Kfree - sparse.diags(w*w*(mr[self.free]+1j*mi[self.free])
            +1j*w*(admittance*self.wall_mass[self.free]+np.sqrt(EPS0/MU0)*self.aperture_mass[self.free])) - 1j*w*g*self.port_matrix
        # S11 is a linear scattering property, including at zero requested drive.
        # Compute it with a unit probe at zero power, then zero the physical field.
        vinc = np.sqrt(2*(incident_w if incident_w else 1.)/g)
        rhs = -2j*w*g*vinc*self.portfree
        # Symmetric diagonal equilibration limits the impedance-wall scale ratio.
        scale = 1/np.sqrt(np.maximum(abs(A.diagonal()),1e-30))
        S = sparse.diags(scale)
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            efree = scale*spsolve((S@A@S).tocsc(),scale*rhs)
        residual = float(np.linalg.norm(A@efree-rhs)/max(np.linalg.norm(rhs),1e-30))
        if not np.all(np.isfinite(efree)) or residual > 1e-7:
            raise RuntimeError(f"Maxwell linear solve failed: residual={residual:.3g}")
        e = np.zeros(grid.ne,complex); e[self.free] = efree
        s11 = complex(self.port@e/vinc-1)
        if incident_w == 0:
            e *= 0
        e2 = grid.cell_e2(e)*grid.active
        qs = .5*(w*EPS0*properties["bed_epp"]+bed_sigma)*e2
        qq = .5*w*EPS0*properties["quartz_epp"]*e2
        ps, pq = float(qs@grid.volume), float(qq@grid.volume)
        pw = float(.5*admittance.real*np.dot(self.wall_mass,abs(e)**2))
        pa = float(.5*np.sqrt(EPS0/MU0)*np.dot(self.aperture_mass,abs(e)**2))
        reflected = float(incident_w*abs(s11)**2)
        balance = float((incident_w-reflected-ps-pq-pw-pa)/max(incident_w,1e-12))
        if abs(balance) > 1e-6 or reflected > incident_w*(1+1e-7):
            raise RuntimeError(f"Electromagnetic power balance failed: {balance:.3g}")
        return EMResult(e,grid.cell_e(e),e2,qs,qq,float(incident_w),reflected,pw,pa,ps,pq,
                        complex(s11),balance,residual,float(frequency_hz))

    def spectrum(self, properties, frequencies, incident_w, weights=None):
        f = np.asarray(frequencies,float)
        weights = np.ones(len(f)) if weights is None else np.asarray(weights,float)
        if len(f) != len(weights) or len(f)<1 or np.any(weights<0) or weights.sum()<=0:
            raise ValueError("Dwell weights must match frequencies and have a positive sum")
        weights = weights/weights.sum()
        runs = [self.solve(properties,float(freq),incident_w) for freq in f]
        return runs, weights
