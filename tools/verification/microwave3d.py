"""Reproducible empty cylindrical TE111 cavity grid study (no material fitting)."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from microwave3d.grid import Grid
from microwave3d.electromagnetics import Maxwell
from microwave3d.materials import EPS0,MU0
from microwave3d.solver import load_config


def cylinder_modes(n,cfg):
    radius=cfg["geometry"]["cavity_diameter_m"]/2
    height=cfg["geometry"]["cavity_height_m"]
    x=np.linspace(-radius,radius,n+1); z=np.linspace(-height/2,height/2,n+1)
    xx,yy,zz=np.meshgrid((x[:-1]+x[1:])/2,(x[:-1]+x[1:])/2,(z[:-1]+z[1:])/2,indexing="ij")
    grid=Grid(x,x.copy(),z,xx**2+yy**2<radius**2)
    local=json.loads(json.dumps(cfg)); local["em"].update(wall_model="pec",tube_apertures=False)
    model=Maxwell(grid,local)
    mass=EPS0*grid.mass(np.ones(grid.n))[model.free]
    inv=sparse.diags(1/np.sqrt(mass))
    matrix=inv@model.Kfree@inv
    # First zero of J1'(x), TE111 cylindrical cavity analytical eigenfrequency.
    exact=np.sqrt((1.8411837813406593/radius)**2+(np.pi/height)**2)/(2*np.pi*np.sqrt(EPS0*MU0))
    eigen=eigsh(matrix.tocsc(),k=4,sigma=(2*np.pi*exact)**2,which="LM",return_eigenvectors=False,tol=1e-9)
    freqs=np.sqrt(np.maximum(eigen,0))/(2*np.pi)
    freq=float(freqs[np.argmin(abs(freqs-exact))])
    return {"cells_per_axis":n,"unknowns":len(model.free),"analytic_te111_hz":exact,"calculated_hz":freq,"relative_error":freq/exact-1,
            "nearby_eigenfrequencies_hz":sorted(freqs.tolist()),"cavity_volume_relative_error":float(grid.volume[grid.active].sum()/(np.pi*radius**2*height)-1)}


def main():
    p=argparse.ArgumentParser(); p.add_argument("--grids",default="8,12,16")
    p.add_argument("--output",default="results3d/verification/cylinder.json")
    args=p.parse_args(); cfg=load_config("examples3d/itaca-cylinder.json")
    rows=[]
    for n in map(int,args.grids.split(",")):
        row=cylinder_modes(n,cfg); rows.append(row); print(json.dumps(row),flush=True)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps({"study":"Empty PEC cylinder, TE111 eigenfrequency, stair-step boundary", "rows":rows},indent=2)+"\n")


if __name__=="__main__":main()
