"""Command line runner, sweeps, portable field export, and offline viewer."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from .solver import Coupled, load_config
from .report import power_chart


def save(path, summary, fields):
    path=Path(path); path.mkdir(parents=True,exist_ok=True)
    (path/"result.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    np.savez_compressed(path/"fields.npz",**fields)
    # Legacy rectilinear VTK, cell ordering x fastest (Fortran array order).
    shape=tuple(summary["mesh"]["shape"])
    with (path/"fields.vtk").open("w",encoding="ascii") as f:
        f.write("# vtk DataFile Version 3.0\nMicrowave3D research fields\nASCII\nDATASET RECTILINEAR_GRID\n")
        f.write("DIMENSIONS "+" ".join(str(n+1) for n in shape)+"\n")
        for axis in "xyz":
            coords=fields[axis+"_m"]
            f.write(f"{axis.upper()}_COORDINATES {len(coords)} double\n")
            f.write(" ".join(map(str,coords))+"\n")
        f.write(f"CELL_DATA {int(np.prod(shape))}\n")
        for key,values in fields.items():
            if np.asarray(values).shape == (int(np.prod(shape)),):
                f.write(f"SCALARS {key} double 1\nLOOKUP_TABLE default\n")
                f.write("\n".join(map(str,np.asarray(values,float).reshape(shape).ravel(order="F")))+"\n")
            elif np.asarray(values).shape == (int(np.prod(shape)),3):
                f.write(f"VECTORS {key} double\n")
                for row in np.asarray(values).reshape((*shape,3)).transpose(2,1,0,3).reshape(-1,3):
                    f.write(" ".join(map(str,row))+"\n")
    viewer=Path(__file__).with_name("viewer.html")
    if viewer.exists():
        public_fields={k:np.asarray(v).tolist() for k,v in fields.items() if not k.startswith("electric_field_")}
        data=json.dumps({"summary":summary,"fields":public_fields},allow_nan=False).replace("<","\\u003c")
        (path/"viewer.html").write_text(viewer.read_text(encoding="utf-8").replace("/*__RESULT_DATA__*/", "window.RESULT="+data+";"),encoding="utf-8")


def main(argv=None):
    parser=argparse.ArgumentParser(description="3D Maxwell / Darcy or channel Navier-Stokes / thermal research model")
    parser.add_argument("--config",default="examples/itaca-cylinder.json")
    parser.add_argument("--output",default="results3d/run")
    parser.add_argument("--power",type=float,help="Incident port power, W")
    parser.add_argument("--powers",help="Comma-separated incident powers, evaluated in the given continuation order")
    parser.add_argument("--frequencies-mhz",help="Comma-separated frequency dwell points, equally weighted")
    parser.add_argument("--spectrum-only",action="store_true",help="Cold-cavity EM scan; no temperature prediction")
    args=parser.parse_args(argv)
    try:
        cfg=load_config(args.config)
        if args.power is not None:
            cfg["solver"]["power_w"]=args.power
        if args.frequencies_mhz:
            cfg["em"]["frequencies_hz"]=[float(x)*1e6 for x in args.frequencies_mhz.split(",")]
            cfg["em"]["dwell_weights"]=[1]*len(cfg["em"]["frequencies_hz"])
        model=Coupled(cfg)
        out=Path(args.output)
        if args.spectrum_only:
            t=np.full(model.grid.n,cfg["thermal"]["ambient_c"])
            runs,weights=model.maxwell.spectrum(model.em_properties(t),cfg["em"]["frequencies_hz"],cfg["solver"]["power_w"],cfg["em"]["dwell_weights"])
            rows=[{"frequency_hz":e.frequency_hz,"s11_real":e.s11.real,"s11_imag":e.s11.imag,"reflected_fraction":abs(e.s11)**2,"sample_w":e.sample_w,"quartz_w":e.quartz_w,"wall_w":e.wall_w,"aperture_w":e.aperture_w,"balance":e.power_residual} for e in runs]
            out.mkdir(parents=True,exist_ok=True)
            (out/"spectrum.json").write_text(json.dumps({"temperature_c":cfg["thermal"]["ambient_c"],"config":cfg,"rows":rows},indent=2)+"\n")
            print(json.dumps(rows,indent=2)); return 0
        powers=[float(p) for p in args.powers.split(",")] if args.powers else [cfg["solver"]["power_w"]]
        initial=None; rows=[]; all_converged=True
        for index,power in enumerate(powers):
            def progress(row):
                if row["iteration"]==1 or row["iteration"]%5==0:
                    print(f"P_in={power:g} W | iteration {row['iteration']} | T_max={row['max_temperature_c']:.2f} C | step={row['undamped_step_k']:.3g} K | thermal residual={row['thermal_power_residual']:.3g}",flush=True)
            summary,fields=model.solve(power,initial,progress)
            target=out if len(powers)==1 else out/f"{index:03d}-{power:g}w"
            save(target,summary,fields)
            all_converged &= summary["converged"]
            rows.append({"incident_w":power,"converged":summary["converged"],**summary["temperature_c"],"sample_absorbed_w":summary["power_w"]["sample"],"reflected_w":summary["power_w"]["reflected"]})
            print(json.dumps({"output":str(target),"converged":summary["converged"],"temperature_c":summary["temperature_c"],"power_w":summary["power_w"]},indent=2),flush=True)
            if not summary["converged"]:
                print("Stopped continuation: unconverged state is saved but not used as a valid branch point.",file=sys.stderr)
                break
            initial=fields["temperature_c"]
        if len(powers)>1:
            with (out/"power-sweep.csv").open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            power_chart(out/"power-sweep.html",rows)
        return 0 if all_converged else 2
    except (ValueError,RuntimeError) as error:
        print(f"microwave3d: {error}",file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
