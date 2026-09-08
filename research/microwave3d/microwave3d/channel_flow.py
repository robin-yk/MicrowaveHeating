"""Conservative transfer of resolved FEM channel velocities to thermal faces."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from .channels import geometry
from .navier_stokes import solve_channel
from .materials import gas_properties
from .flow import FlowResult


def channel_flow(grid,cfg):
    centres,width=geometry(cfg); c=cfg['channels']; gas=cfg['gas']
    length=cfg['geometry']['sample_length_m']
    centres=centres+np.asarray(cfg['geometry'].get('tube_offset_m',[0,0]))
    rho,mu,cp,mdot=gas_properties(gas.get('reference_temperature_c',gas['inlet_c']),gas)
    rate=mdot/len(centres)
    solution=solve_channel(width,length,rate/(rho*width**2),float(rho),float(mu),
        c['flow_cells_across'],c['flow_cells_axial'])
    if abs(solution.diagnostics['pressure_drop_pa']) > .01*gas['pressure_pa']:
        raise ValueError('Channel pressure drop exceeds the 1% constant-density pressure limit')
    pressure=np.zeros(grid.n)
    velocity=np.zeros((grid.n,3))
    all_links=[]; all_flux=[]; all_in=[]; all_out=[]; all_fin=[]; all_fout=[]
    correction_norm=raw_norm=0.; residual=0.
    for centre in centres:
        mask=(abs(grid.xyz[:,0]-centre[0])<width/2-1e-13)&(abs(grid.xyz[:,1]-centre[1])<width/2-1e-13)&(abs(grid.xyz[:,2])<length/2)
        ids=np.flatnonzero(mask)
        if not len(ids):
            raise ValueError('No thermal cells inside a channel')
        def local(x):
            p=np.array(x,copy=True)
            p[:,:2]-=centre-width/2; p[:,2]=length/2-p[:,2]
            return p
        uv,pp=solution.sample(local(grid.xyz[ids])); pressure[ids]=pp
        uv[:,2]*=-1; velocity[ids]=uv
        links=grid.links
        left,right=links[:,:2].astype(int).T
        links=links[mask[left]&mask[right]]
        left,right=links[:,:2].astype(int).T
        axes=links[:,2].astype(int)
        faces=grid.xyz[left].copy()
        faces[np.arange(len(links)),axes]+=links[:,4]
        # Gaussian face integration followed by a discrete H(div) projection. The
        # correction magnitude is reported, not hidden as FEM conservation.
        def integrated_speed(points,axes,lengths):
            result=np.zeros(len(points))
            other=np.array([[d for d in range(3) if d!=a] for a in axes])
            rows=np.arange(len(points))
            for a in (-1/np.sqrt(3),1/np.sqrt(3)):
                for b in (-1/np.sqrt(3),1/np.sqrt(3)):
                    p=points.copy()
                    p[rows,other[:,0]]+=a*lengths[rows,other[:,0]]/2
                    p[rows,other[:,1]]+=b*lengths[rows,other[:,1]]/2
                    u,_=solution.sample(local(p)); u[:,2]*=-1
                    result+=u[rows,axes]/4
            return result
        raw=rho*integrated_speed(faces,axes,grid.cell_lengths[left])*links[:,3]
        top=ids[np.isclose(grid.xyz[ids,2],grid.xyz[ids,2].max(),rtol=0,atol=1e-12)]
        bottom=ids[np.isclose(grid.xyz[ids,2],grid.xyz[ids,2].min(),rtol=0,atol=1e-12)]
        def end_mass(cells,z):
            points=grid.xyz[cells].copy(); points[:,2]=z
            speed=-integrated_speed(points,np.full(len(cells),2),grid.cell_lengths[cells])
            f=rho*speed*grid.volume[cells]/grid.cell_lengths[cells,2]
            if np.any(f < -1e-14):
                raise RuntimeError('Channel outlet backflow needs a different thermal boundary condition')
            return f*(rate/f.sum()) if rate>0 else np.zeros(len(cells))
        fin=end_mass(top,length/2); fout=end_mass(bottom,-length/2)
        numbering=np.full(grid.n,-1); numbering[ids]=np.arange(len(ids))
        il,ir=numbering[left],numbering[right]
        D=sparse.coo_matrix((np.r_[np.ones(len(left)),-np.ones(len(left))],
            (np.r_[il,ir],np.r_[np.arange(len(left)),np.arange(len(left))])),shape=(len(ids),len(left))).tocsr()
        target=np.zeros(len(ids)); np.add.at(target,numbering[top],fin); np.add.at(target,numbering[bottom],-fout)
        conductance=links[:,3]/(links[:,4]+links[:,5])
        L=D@sparse.diags(conductance)@D.T
        error=D@raw-target
        phi=np.zeros(len(ids)); phi[1:]=spsolve(L[1:,1:].tocsc(),error[1:])
        corrected=raw-conductance*(D.T@phi)
        err=float(np.max(abs(D@corrected-target))/max(rate,1e-30))
        residual=max(residual,err)
        correction_norm+=float(np.sum((corrected-raw)**2)); raw_norm+=float(np.sum(raw**2))
        all_links.append(links); all_flux.append(corrected); all_in.append(top); all_out.append(bottom); all_fin.append(fin); all_fout.append(fout)
    if residual>1e-8:
        raise RuntimeError('Channel finite-volume mass projection failed')
    diag=solution.diagnostics.copy()
    diag.update({'channels':len(centres),'equal_flow_split':True,
        'reference_temperature_c':gas.get('reference_temperature_c',gas['inlet_c']),
        'face_flux_projection_relative_l2':float(np.sqrt(correction_norm/max(raw_norm,1e-60))),
        'finite_volume_mass_relative':residual,
        'pressure_drop_definition':'FEM area-average pressure at inlet minus outlet boundary planes',
        'thermal_expansion_in_flow':False})
    return FlowResult(pressure,np.concatenate(all_links),np.concatenate(all_flux),
        np.concatenate(all_in),np.concatenate(all_fin),np.concatenate(all_out),np.concatenate(all_fout),
        cp,diag['pressure_drop_pa'],residual,0.,diag,velocity)
