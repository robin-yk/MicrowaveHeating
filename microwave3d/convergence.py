"""Independent channel-thermal, flow, or cavity-grid sensitivity runs.

Success means the LAST tested pair meets requested tolerances, not a rigorous
error bound or experimental validation. Nonmonotone sequences remain visible.
"""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from .solver import Coupled,load_config
from .cli import save


def study(cfg,levels,kind='thermal',temperature_tolerance=1.,relative_tolerance=.02,output=None):
    if len(levels)<3 or any(int(n)!=n or n<2 for n in levels) or np.any(np.diff(levels)<=0):
        raise ValueError('At least three increasing integer refinement levels >= 2 required')
    if kind not in ('thermal','flow','em'):
        raise ValueError('kind must be thermal, flow, or em')
    if kind != 'em' and not cfg.get('channels'):
        raise ValueError('Channel thermal/flow study requires channels configuration')
    rows=[]
    for level in levels:
        run=copy.deepcopy(cfg)
        if kind=='thermal':
            run['channels'].update(thermal_cells_across=level,thermal_cells_axial=6*level)
        elif kind=='flow':
            run['channels'].update(flow_cells_across=level,flow_cells_axial=4*level)
        else:
            run['mesh'].update(core_xy=level,core_z=level,outer_xy=max(2,level//2),outer_z=max(2,level//2))
        model=Coupled(run)
        summary,fields=model.solve(progress=lambda r: print(f"{kind} level {level} | iteration {r['iteration']} | step {r['undamped_step_k']:.3g} K",flush=True) if r['iteration']==1 or r['iteration']%10==0 else None)
        if output:
            save(Path(output)/f'{kind}-{level}',summary,fields)
        row={'level':level,'converged':summary['converged'],'temperature_max_c':summary['temperature_c']['sample_max'],
            'sample_absorbed_w':summary['power_w']['sample'],
            'pressure_drop_pa':summary['flow']['pressure_drop_pa'] if summary['flow'] else None,
            'gas_outlet_c':summary['temperature_c'].get('channel_gas_outlet_mass_average'),
            'thermal_shape':summary['mesh']['shape'],'em_shape':summary['mesh']['em_shape']}
        if rows:
            previous=rows[-1]
            row['temperature_change_k']=abs(row['temperature_max_c']-previous['temperature_max_c'])
            for key in ('sample_absorbed_w','pressure_drop_pa'):
                row[key+'_relative_change']=abs(row[key]-previous[key])/max(abs(row[key]),1e-12) if row[key] is not None else None
        rows.append(row)
        if not summary['converged']:
            break
    last=rows[-1]
    passed=(len(rows)==len(levels) and all(r['converged'] for r in rows)
        and last['temperature_change_k']<temperature_tolerance
        and last['sample_absorbed_w_relative_change']<relative_tolerance
        and (kind!='flow' or last['pressure_drop_pa_relative_change']<relative_tolerance))
    result={'kind':kind,'last_pair_passed':bool(passed),'temperature_tolerance_k':temperature_tolerance,
        'relative_tolerance':relative_tolerance,'rows':rows,'experimentally_validated':False,
        'scope':'Pairwise sensitivity only; not a Richardson error bound. Thermal study fixes EM mesh and FEM flow; flow study fixes thermal/EM meshes. EM refinement also subdivides the nested thermal mesh. Port and material uncertainty are independent.'}
    if output:
        path=Path(output);path.mkdir(parents=True,exist_ok=True)
        (path/'convergence.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='examples3d/itaca-monolith-ns.json')
    p.add_argument('--kind',choices=['thermal','flow','em'],default='thermal')
    p.add_argument('--levels',default='2,3,4')
    p.add_argument('--output',default='results3d/convergence')
    args=p.parse_args(argv)
    result=study(load_config(args.config),[int(x) for x in args.levels.split(',')],args.kind,output=args.output)
    print(json.dumps(result,indent=2))
    return 0 if result['last_pair_passed'] else 2


if __name__=='__main__':
    raise SystemExit(main())
