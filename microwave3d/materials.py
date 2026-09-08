"""Explicit material data and subcell volume fractions; no automatic fitting."""
import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
SIGMA = 5.670374419e-8
RGAS = 8.314462618


class Materials:
    def __init__(self, grid, cfg):
        self.grid, self.cfg = grid, cfg
        geo = cfg["geometry"]
        r = geo["tube_inner_diameter_m"]/2
        ro = r + geo["tube_wall_m"]
        cx,cy = geo.get("tube_offset_m", [0,0])
        # Cross-sectional integration preserves thin quartz instead of dropping it
        # when its centre falls between cells. Geometry convergence is still needed.
        n = cfg["mesh"].get("subcell_samples", 8)
        f_inner = np.zeros(grid.n); f_outer = np.zeros(grid.n)
        for a in (np.arange(n)+.5)/n-.5:
            for b in (np.arange(n)+.5)/n-.5:
                rr = (grid.xyz[:,0]+a*grid.cell_lengths[:,0]-cx)**2 + (grid.xyz[:,1]+b*grid.cell_lengths[:,1]-cy)**2
                f_inner += rr < r*r; f_outer += rr < ro*ro
        f_inner /= n*n; f_outer /= n*n
        in_bed = abs(grid.xyz[:,2]-geo.get("sample_z_m",0)) < geo["sample_length_m"]/2
        self.bed = f_inner * in_bed * grid.active
        self.quartz = (f_outer-f_inner) * grid.active
        self.gas = f_inner * ~in_bed * grid.active
        self.air = grid.active.astype(float)-self.bed-self.quartz-self.gas
        self.fluid_fraction = self.bed.copy()
        if cfg.get("channels"):
            from .channels import fraction
            channels = fraction(grid,cfg)
            if np.any(channels > self.bed+1e-8):
                raise ValueError("Channel volume exceeds the monolith volume; refine geometry quadrature")
            self.bed -= channels
            self.gas += channels
            self.fluid_fraction = channels
        self.tube_area = np.zeros(grid.n)
        # Exact total cylinder area distributed among surface-containing cells.
        theta = 2*np.pi*(np.arange(4096)+.5)/4096
        ii = np.searchsorted(grid.axes[0], cx+ro*np.cos(theta))-1
        jj = np.searchsorted(grid.axes[1], cy+ro*np.sin(theta))-1
        ii = np.clip(ii,0,grid.shape[0]-1); jj = np.clip(jj,0,grid.shape[1]-1)
        for k,dz in enumerate(grid.dx[2]):
            np.add.at(self.tube_area, grid.ids[ii,jj,k], 2*np.pi*ro*dz/len(theta))
        if np.any((self.tube_area>0)&~grid.active):
            raise ValueError("Quartz tube intersects the cavity wall")
        self.cavity_area = 2*np.pi*(geo["cavity_diameter_m"]/2)*(geo["cavity_height_m"]+geo["cavity_diameter_m"]/2)
        self.sample_volume = float(np.dot(self.bed,grid.volume))
        self.sample_volume_exact = np.pi*r*r*geo["sample_length_m"]
        if cfg.get("channels"):
            self.sample_volume_exact -= len(cfg["channels"]["centres_m"])*cfg["channels"]["width_m"]**2*geo["sample_length_m"]

    def at(self, temperature):
        cfg = self.cfg
        data = np.array(cfg["material"]["table"],float)
        t = np.asarray(temperature)
        if np.any(t[self.grid.active] < data[0,0]-1e-5) or np.any(t[self.grid.active] > data[-1,0]+1e-5):
            raise ValueError("Temperature outside material table; extend with justified data, not silent extrapolation")
        ep,epp,k = [np.interp(t,data[:,0],data[:,c]) for c in (1,2,3)]
        quartz = cfg["quartz"]
        dielectric_real = self.air+self.gas+self.bed*ep+self.quartz*quartz["epsilon_real"]
        dielectric_loss = self.bed*epp+self.quartz*quartz["epsilon_loss"]
        kg = .1513*((t+273.15)/293.15)**.72 if cfg["gas"]["species"] == "He" else .0258*((t+273.15)/293.15)**.72
        ka = cfg["thermal"]["air_conductivity_factor"]*.0263*((t+273.15)/293.15)**.76
        conductivity = self.air*ka + self.gas*kg + self.bed*k + self.quartz*quartz["conductivity_w_mk"]
        return {"ep":dielectric_real, "epp":dielectric_loss, "k":conductivity,
                "bed_epp":self.bed*epp, "quartz_epp":self.quartz*quartz["epsilon_loss"],
                "bed_sigma":self.bed*cfg["material"].get("electrical_conductivity_s_m",0)}


def gas_properties(temperature, cfg):
    tk = np.asarray(temperature)+273.15
    he = cfg["species"] == "He"
    molar_mass, cp, mu0, s = (0.004002602,5193.,1.96e-5,79.4) if he else (.0280134,1039.7,1.76e-5,111.)
    rho = cfg["pressure_pa"]*molar_mass/(RGAS*tk)
    mu = mu0*(tk/293.15)**1.5*(293.15+s)/(tk+s)
    mdot = cfg["flow_sccm"]*1e-6/(60*.022414)*molar_mass
    return rho,mu,cp,mdot
