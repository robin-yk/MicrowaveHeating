"""Orthogonal 3D grid and integral Yee operators, all lengths in metres.

Electric unknowns are oriented edge voltages. C maps edge voltages to face
circulations. Material masses and heat deposition use the same cell quadrature.
The cylinder is a stair-step cell union, with its geometry error reported.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse


@dataclass
class Grid:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    active: np.ndarray | None = None
    electric: bool = True

    def __post_init__(self):
        self.axes = tuple(np.asarray(a, dtype=float) for a in (self.x, self.y, self.z))
        if any(a.ndim != 1 or len(a) < 3 or not np.all(np.diff(a) > 0) for a in self.axes):
            raise ValueError("Grid axes must contain at least three increasing nodes")
        self.dx = tuple(np.diff(a) for a in self.axes)
        self.centres = tuple((a[:-1] + a[1:]) / 2 for a in self.axes)
        self.shape = tuple(len(a) - 1 for a in self.axes)
        self.n = int(np.prod(self.shape))
        self.ids = np.arange(self.n).reshape(self.shape)
        self.volume = np.einsum("i,j,k->ijk", *self.dx).ravel()
        self.xyz = np.stack(np.meshgrid(*self.centres, indexing="ij"), axis=-1).reshape(-1, 3)
        self.active = np.ones(self.n, bool) if self.active is None else np.asarray(self.active, bool).ravel()
        if len(self.active) != self.n or not np.any(self.active):
            raise ValueError("Empty or mismatched cavity mask")
        self.cell_lengths = np.stack(np.meshgrid(*self.dx, indexing="ij"), axis=-1).reshape(-1, 3)
        if self.electric:
            self._topology()
        else:
            self._thermal_topology()

    def _thermal_topology(self):
        links, boundary = [], []
        mask = self.active.reshape(self.shape)
        for idx in np.ndindex(self.shape):
            if not mask[idx]:
                continue
            i = self.ids[idx]
            for axis in range(3):
                area = self.volume[i]/self.cell_lengths[i,axis]
                for sign in (-1,1):
                    jdx = list(idx); jdx[axis] += sign
                    if 0 <= jdx[axis] < self.shape[axis] and mask[tuple(jdx)]:
                        if sign > 0:
                            j = self.ids[tuple(jdx)]
                            links.append((i,j,axis,area,self.cell_lengths[i,axis]/2,self.cell_lengths[j,axis]/2))
                    else:
                        boundary.append((i,axis,sign,area,self.cell_lengths[i,axis]/2,-1))
        self.links = np.asarray(links).reshape(-1,6)
        self.boundary = np.asarray(boundary).reshape(-1,6)

    def _topology(self):
        nx, ny, nz = self.shape
        edge_shapes = [(nx, ny+1, nz+1), (nx+1, ny, nz+1), (nx+1, ny+1, nz)]
        offset = 0
        self.edges = []
        for shape in edge_shapes:
            size = int(np.prod(shape))
            self.edges.append(np.arange(offset, offset+size).reshape(shape))
            offset += size
        self.ne = offset
        ex, ey, ez = self.edges
        self.cell_edges = np.stack([
            ex[:, :-1, :-1], ex[:, 1:, :-1], ex[:, :-1, 1:], ex[:, 1:, 1:],
            ey[:-1, :, :-1], ey[1:, :, :-1], ey[:-1, :, 1:], ey[1:, :, 1:],
            ez[:-1, :-1, :], ez[1:, :-1, :], ez[:-1, 1:, :], ez[1:, 1:, :],
        ], axis=-1).reshape(self.n, 12)
        lengths = np.stack(np.meshgrid(*self.dx, indexing="ij"), axis=-1).reshape(-1, 3)
        self.cell_lengths = lengths
        self.edge_weights = self.volume[:, None] / (4 * np.repeat(lengths**2, 4, axis=1))
        rows, cols, values, metric, face_edges = [], [], [], [], []
        face_id = 0
        self.face_lookup = []
        # Oriented circulations for faces normal to x, y, z.
        for axis, shape in enumerate([(nx+1,ny,nz), (nx,ny+1,nz), (nx,ny,nz+1)]):
            lookup = np.empty(shape, int)
            for i,j,k in np.ndindex(shape):
                idx = (i,j,k)
                if axis == 0:
                    edges = [ey[i,j,k], ez[i,j+1,k], ey[i,j,k+1], ez[i,j,k]]
                elif axis == 1:
                    edges = [ez[i,j,k], ex[i,j,k+1], ez[i+1,j,k], ex[i,j,k]]
                else:
                    edges = [ex[i,j,k], ey[i+1,j,k], ex[i,j+1,k], ey[i,j,k]]
                a = idx[axis]
                width = self.dx[axis]
                lower = list(idx); lower[axis] -= 1
                upper = list(idx)
                active = self.active.reshape(self.shape)
                dual = ((width[a-1] if a > 0 and active[tuple(lower)] else 0)
                        + (width[a] if a < len(width) and active[tuple(upper)] else 0)) / 2
                area = np.prod([self.dx[d][idx[d]] for d in range(3) if d != axis])
                rows.extend([face_id]*4); cols.extend(edges); values.extend([1,1,-1,-1])
                metric.append(dual/area); face_edges.append(edges)
                lookup[idx] = face_id; face_id += 1
            self.face_lookup.append(lookup)
        self.C = sparse.coo_matrix((values, (rows, cols)), shape=(face_id,self.ne)).tocsr()
        self.face_metric = np.array(metric)
        self.face_edges = np.array(face_edges)
        used = np.zeros(self.ne, bool)
        used[self.cell_edges[self.active].ravel()] = True
        self.used_edges = used
        # Interior cell links and cavity boundary faces, each counted once.
        links, boundary = [], []
        mask = self.active.reshape(self.shape)
        for idx in np.ndindex(self.shape):
            if not mask[idx]:
                continue
            cell = self.ids[idx]
            for axis in range(3):
                area = self.volume[cell] / self.dx[axis][idx[axis]]
                for sign in (-1, 1):
                    other_idx = list(idx); other_idx[axis] += sign
                    inside = 0 <= other_idx[axis] < self.shape[axis]
                    other_idx = tuple(other_idx)
                    if inside and mask[other_idx]:
                        if sign > 0:
                            other = self.ids[other_idx]
                            links.append((cell, other, axis, area,
                                          self.cell_lengths[cell,axis]/2,
                                          self.cell_lengths[other,axis]/2))
                    else:
                        fidx = list(idx); fidx[axis] += sign > 0
                        face = self.face_lookup[axis][tuple(fidx)]
                        boundary.append((cell, axis, sign, area, self.cell_lengths[cell,axis]/2, face))
        self.links = np.asarray(links)
        self.boundary = np.asarray(boundary)
        self.pec_edges = np.unique(self.face_edges[self.boundary[:,5].astype(int)])
        self.free_pec = np.flatnonzero(used & ~np.isin(np.arange(self.ne), self.pec_edges))
        # Surface mass integral of |E_t|² for the impedance wall.
        self.wall_mass = np.zeros(self.ne)
        for cell, axis, _, area, _, face in self.boundary:
            edges = self.face_edges[int(face)]
            for edge in edges:
                direction = next(d for d, ids in enumerate(self.edges) if ids.min() <= edge <= ids.max())
                self.wall_mass[edge] += area / (2*self.cell_lengths[int(cell),direction]**2)

    def mass(self, cell_property):
        return np.bincount(self.cell_edges.ravel(),
            weights=(self.edge_weights*np.asarray(cell_property)[:,None]*self.active[:,None]).ravel(), minlength=self.ne)

    def cell_e2(self, voltage):
        """Cell average |E|² using the same positive quadrature as mass()."""
        return np.sum(self.edge_weights * abs(voltage[self.cell_edges])**2, axis=1) / self.volume

    def cell_e(self, voltage):
        v = voltage[self.cell_edges].reshape(self.n,3,4)
        return v.mean(axis=2) / self.cell_lengths

    def closest(self, xyz, mask=None):
        valid = self.active if mask is None else self.active & mask
        ids = np.flatnonzero(valid)
        if not len(ids):
            raise ValueError("No cells at requested material/sensor location")
        return ids[np.argmin(np.sum((self.xyz[ids] - xyz)**2, axis=1))]


def clustered_axis(half, core_half, core_cells, outer_cells):
    if not (half > core_half > 0 and core_cells >= 2 and outer_cells >= 1):
        raise ValueError("Invalid clustered grid")
    core = np.linspace(-core_half, core_half, core_cells+1)
    outer = np.linspace(core_half, half, outer_cells+1)[1:]
    return np.r_[-outer[::-1], core, outer]


def reactor_grid(cfg):
    geo, mesh = cfg["geometry"], cfg["mesh"]
    r, h = geo["cavity_diameter_m"]/2, geo["cavity_height_m"]
    ro = geo["tube_inner_diameter_m"]/2 + geo["tube_wall_m"]
    x = clustered_axis(r, ro, mesh["core_xy"], mesh["outer_xy"])
    y = x.copy()
    z = clustered_axis(h/2, geo["sample_length_m"]/2, mesh["core_z"], mesh["outer_z"])
    xx,yy,zz = np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2,(z[:-1]+z[1:])/2,indexing="ij")
    mask = xx**2+yy**2 < r*r
    return Grid(x,y,z,mask)
