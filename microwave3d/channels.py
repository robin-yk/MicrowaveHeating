"""Square-channel geometry and conservative nested EM/thermal mesh transfer."""
import numpy as np
from .grid import Grid


def geometry(cfg):
    c = cfg["channels"]
    centres = np.asarray(c["centres_m"],float)
    width = float(c["width_m"])
    radius = cfg["geometry"]["tube_inner_diameter_m"]/2
    if centres.ndim != 2 or centres.shape[1] != 2 or len(centres) == 0 or not np.all(np.isfinite(centres)):
        raise ValueError("Channel centres must be a finite nonempty N x 2 array")
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Channel width must be positive")
    if np.any(np.linalg.norm(abs(centres)+width/2,axis=1) >= radius):
        raise ValueError("Channel corners must fit inside the monolith")
    for i,p in enumerate(centres):
        if np.any(np.all(abs(centres[:i]-p) <= width,axis=1)):
            raise ValueError("Square channels overlap or leave zero solid wall thickness")
    for key in ("thermal_cells_across","thermal_cells_axial","flow_cells_across","flow_cells_axial"):
        if int(c[key]) != c[key] or c[key] < 2:
            raise ValueError(f"{key} must be an integer >= 2")
    return centres,width


def fraction(grid,cfg):
    centres,width = geometry(cfg)
    offset = np.asarray(cfg["geometry"].get("tube_offset_m",[0,0]))
    result = np.zeros(grid.n)
    for centre in centres+offset:
        overlap = np.ones(grid.n)
        for d in range(2):
            lo=grid.xyz[:,d]-grid.cell_lengths[:,d]/2
            hi=lo+grid.cell_lengths[:,d]
            overlap *= np.maximum(0,np.minimum(hi,centre[d]+width/2)-np.maximum(lo,centre[d]-width/2))/grid.cell_lengths[:,d]
        result += overlap
    result*=((abs(grid.xyz[:,2]) < cfg["geometry"]["sample_length_m"]/2)*grid.active)
    # Aligned faces can leave 1e-14 numerical "solid" inside a pure gas cell.
    # Remove roundoff only, so extrema do not classify gas as sample material.
    result[np.isclose(result,1,rtol=0,atol=1e-12)]=1.
    result[abs(result)<1e-14]=0.
    return result


class NestedMesh:
    """Each thermal cell belongs to exactly one EM cell; heat transfer sums watts."""
    def __init__(self,em,cfg):
        centres,width=geometry(cfg); c=cfg["channels"]
        offset=np.asarray(cfg["geometry"].get("tube_offset_m",[0,0]))
        axes=[]
        for d in range(2):
            points=[em.axes[d]]
            for centre in centres+offset:
                points.append(np.linspace(centre[d]-width/2,centre[d]+width/2,c["thermal_cells_across"]+1))
            # Merge roundoff-equivalent nodes without creating tiny sliver cells.
            axes.append(np.unique(np.round(np.concatenate(points),14)))
        length=cfg["geometry"]["sample_length_m"]
        axes.append(np.unique(np.round(np.r_[em.axes[2],np.linspace(-length/2,length/2,c["thermal_cells_axial"]+1)],14)))
        xyz=np.stack(np.meshgrid(*[(a[1:]+a[:-1])/2 for a in axes],indexing="ij"),axis=-1).reshape(-1,3)
        indices=[np.searchsorted(em.axes[d],xyz[:,d])-1 for d in range(3)]
        self.parent=em.ids[tuple(indices)]
        self.grid=Grid(*axes,active=em.active[self.parent],electric=False)
        self.em=em
        if not np.allclose(np.bincount(self.parent,weights=self.grid.volume,minlength=em.n),em.volume,rtol=1e-9,atol=1e-20):
            raise RuntimeError("Nested mesh volumes do not close")

    def average(self,values):
        return np.bincount(self.parent,weights=np.asarray(values)*self.grid.volume,minlength=self.em.n)/self.em.volume

    def deposit(self,source,phase):
        volume=np.bincount(self.parent,weights=phase*self.grid.volume,minlength=self.em.n)
        watts=source*self.em.volume
        if np.any((volume == 0)&(abs(watts)>1e-14)):
            raise RuntimeError("EM absorption has no corresponding thermal phase")
        density=np.divide(watts,volume,out=np.zeros_like(watts),where=volume>0)
        return density[self.parent]*phase
