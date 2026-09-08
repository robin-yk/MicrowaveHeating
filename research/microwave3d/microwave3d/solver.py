"""Coupled Maxwell, porous or resolved-channel flow, and steady heat transfer."""
import copy
import json
import time
import numpy as np
from .grid import reactor_grid
from .materials import Materials, EPS0
from .electromagnetics import Maxwell
from .flow import darcy, cell_velocity
from .thermal import Thermal
from .channels import NestedMesh, geometry as channel_geometry


def validate(cfg):
    def positive(v, label):
        if not np.isfinite(v) or v <= 0:
            raise ValueError(f"{label} must be positive and finite")
    if cfg.get("channels"):
        channel_geometry(cfg)
        if cfg["gas"].get("model") != "navier_stokes":
            raise ValueError("Resolved channels require gas.model = navier_stokes")
        reference=cfg["gas"].get("reference_temperature_c",cfg["gas"]["inlet_c"])
        if not np.isfinite(reference) or reference<=-273.15:
            raise ValueError("Flow reference temperature must be finite and above absolute zero")
    elif cfg["gas"].get("model","darcy") != "darcy":
        raise ValueError("Navier-Stokes requires explicit channels geometry")
    for k,v in cfg["geometry"].items():
        if k not in ("tube_offset_m","sample_z_m"):
            positive(v,k)
    geo = cfg["geometry"]
    if geo["sample_length_m"] >= geo["cavity_height_m"]:
        raise ValueError("Sample must fit within the cavity")
    if geo.get("sample_z_m",0) != 0:
        raise ValueError("This grid aligns centred bed end faces; sample_z_m must be zero")
    if np.linalg.norm(geo.get("tube_offset_m",[0,0]))+geo["tube_inner_diameter_m"]/2+geo["tube_wall_m"] >= geo["cavity_diameter_m"]/2:
        raise ValueError("Tube must fit inside the cavity")
    for k in ("core_xy","core_z","outer_xy","outer_z","subcell_samples"):
        positive(cfg["mesh"][k],k)
        if int(cfg["mesh"][k]) != cfg["mesh"][k]:
            raise ValueError("Mesh counts must be integers")
    table = np.asarray(cfg["material"]["table"],float)
    if table.ndim != 2 or table.shape[1] != 4 or len(table)<2 or not np.all(np.isfinite(table)):
        raise ValueError("Material table must contain finite [T, epsilon', epsilon'', k] rows")
    if np.any(np.diff(table[:,0])<=0) or np.any(table[:,1]<=0) or np.any(table[:,2]<0) or np.any(table[:,3]<=0):
        raise ValueError("Material table needs increasing T, positive epsilon'/k, nonnegative loss")
    if not cfg.get("channels"):
        positive(cfg["material"]["particle_diameter_m"],"particle diameter")
    sigma=cfg["material"].get("electrical_conductivity_s_m",0)
    if not np.isfinite(sigma) or sigma<0:
        raise ValueError("Electrical conductivity must be finite and nonnegative")
    if cfg["material"].get("loss_definition","dielectric_only") not in ("dielectric_only","total_measured"):
        raise ValueError("loss_definition must be dielectric_only or total_measured")
    if sigma>0 and cfg["material"].get("loss_definition")=="total_measured":
        raise ValueError("Do not add conductivity to a dielectric-loss table already containing total measured loss")
    if not cfg.get("channels") and not 0 < cfg["material"]["porosity"] < 1:
        raise ValueError("Packed-bed porosity must lie strictly between zero and one")
    for value in cfg["material"].get("permeability_ratios",[1,1,1]):
        positive(value,"permeability ratio")
    if len(cfg["material"].get("permeability_ratios",[1,1,1])) != 3:
        raise ValueError("Three permeability ratios required")
    if cfg["em"]["wall_model"] not in ("pec","impedance"):
        raise ValueError("wall_model must be pec or impedance")
    positive(cfg["em"]["wall_conductivity_s_m"],"wall conductivity")
    for k in ("epsilon_real","conductivity_w_mk"):
        positive(cfg["quartz"][k],k)
    if cfg["quartz"]["epsilon_loss"] < 0:
        raise ValueError("Quartz loss must be nonnegative")
    if not 0 <= cfg["quartz"]["thermal_emissivity"] <= 1 or not 0 < cfg["thermal"]["wall_emissivity"] <= 1:
        raise ValueError("Tube emissivity must be in [0,1], wall emissivity in (0,1]")
    for key in ("resistance_ohm","probe_length_m"):
        positive(cfg["port"][key],key)
    for key in ("temperature_tolerance_k","power_tolerance","max_iterations"):
        positive(cfg["solver"][key],key)
    if not 0 < cfg["solver"]["relaxation"] <= 1:
        raise ValueError("Relaxation must lie in (0,1]")
    if cfg["solver"].get("thermal_linear_solver","direct") not in ("direct","krylov"):
        raise ValueError("thermal_linear_solver must be direct or krylov")
    if not np.isfinite(cfg["solver"]["power_w"]) or cfg["solver"]["power_w"] < 0:
        raise ValueError("Incident power must be finite and nonnegative")
    f = np.asarray(cfg["em"]["frequencies_hz"],float)
    weights = np.asarray(cfg["em"]["dwell_weights"],float)
    if f.ndim!=1 or len(f)<1 or np.any(~np.isfinite(f)) or np.any(f<=0) or weights.shape!=f.shape or np.any(~np.isfinite(weights)) or np.any(weights<0) or weights.sum()<=0:
        raise ValueError("Provide positive frequencies and matching nonnegative dwell weights")
    if cfg["gas"]["species"] not in ("He","N2") or cfg["gas"]["flow_sccm"] < 0:
        raise ValueError("Supported gas: He or N2 at nonnegative flow")
    positive(cfg["gas"]["pressure_pa"],"gas pressure")
    positive(cfg["thermal"]["air_conductivity_factor"],"air conductivity factor")
    contact = cfg["thermal"].get("contact_h_w_m2k")
    if contact is not None:
        positive(contact,"contact conductance")
    for t in (cfg["thermal"]["ambient_c"],cfg["gas"]["inlet_c"]):
        if not table[0,0] <= t <= table[-1,0] or t <= -273.15:
            raise ValueError("Boundary temperature must lie inside the material table and above absolute zero")
    return cfg


def load_config(path):
    with open(path,encoding="utf-8") as f:
        return validate(json.load(f))


class Coupled:
    def __init__(self, cfg):
        self.cfg = validate(copy.deepcopy(cfg))
        self.em_grid = reactor_grid(self.cfg)
        self.transfer = NestedMesh(self.em_grid,self.cfg) if self.cfg.get("channels") else None
        self.grid = self.transfer.grid if self.transfer else self.em_grid
        self.materials = Materials(self.grid,self.cfg)
        self.maxwell = Maxwell(self.em_grid,self.cfg)
        self.thermal = Thermal(self.grid,self.materials,self.cfg)
        self.channel_flow = None
        self._em_cache = None
        if self.transfer and self.cfg["gas"]["enabled"]:
            from .channel_flow import channel_flow
            self.channel_flow = channel_flow(self.grid,self.cfg)

    def em_properties(self,temperature):
        p=self.materials.at(temperature)
        return {k:self.transfer.average(v) for k,v in p.items()} if self.transfer else p

    def display_em(self,values):
        return values[self.transfer.parent] if self.transfer else values

    def thermal_source(self,em,weights,phase,properties=None):
        q=sum(w*getattr(e,"q_"+phase) for w,e in zip(weights,em))
        if not self.transfer:
            return q
        if properties is None:
            raise ValueError("Fine-grid material properties required for channel heat transfer")
        result=np.zeros(self.grid.n)
        for w,e in zip(weights,em):
            loss=2*np.pi*e.frequency_hz*EPS0*properties["bed_epp" if phase=="sample" else "quartz_epp"]
            if phase=="sample":
                loss=loss+properties["bed_sigma"]
            result+=w*self.transfer.deposit(getattr(e,"q_"+phase),loss)
        return result

    def evaluate(self, temperature, power):
        p = self.materials.at(temperature)
        ep={k:self.transfer.average(v) for k,v in p.items()} if self.transfer else p
        key=(power,tuple(self.cfg["em"]["frequencies_hz"]),tuple(self.cfg["em"]["dwell_weights"]))
        arrays=[ep[k] for k in ("ep","epp","bed_epp","quartz_epp","bed_sigma")]
        if self._em_cache is not None and key == self._em_cache[0] and all(np.array_equal(a,b) for a,b in zip(arrays,self._em_cache[1])):
            em,weights=self._em_cache[2:]
        else:
            em,weights = self.maxwell.spectrum(ep,self.cfg["em"]["frequencies_hz"],power,self.cfg["em"]["dwell_weights"])
            self._em_cache=(key,[a.copy() for a in arrays],em,weights)
        source = self.thermal_source(em,weights,"sample",p)+self.thermal_source(em,weights,"quartz",p)
        flow = self.channel_flow if self.transfer else (darcy(self.grid,self.materials.fluid_fraction,temperature,self.cfg) if self.cfg["gas"]["enabled"] else None)
        audit = self.thermal.audit(temperature,p,flow,source)
        return p,em,weights,source,flow,audit

    def solve(self, power=None, initial=None, progress=None):
        cfg=self.cfg; g=self.grid; m=self.materials
        power=cfg["solver"]["power_w"] if power is None else float(power)
        if power < 0 or not np.isfinite(power):
            raise ValueError("Incident power must be finite and nonnegative")
        temperature=np.full(g.n,cfg["thermal"]["ambient_c"],dtype=float) if initial is None else np.asarray(initial,float).copy()
        if temperature.shape != (g.n,):
            raise ValueError("Initial temperature must match the mesh")
        history=[]; converged=False; started=time.monotonic()
        # Audit and convergence use EM and flow freshly solved at the current T.
        for it in range(int(cfg["solver"]["max_iterations"])):
            p,em,weights,source,flow,audit=self.evaluate(temperature,power)
            candidate=self.thermal.step(temperature,p,flow,source)
            delta=float(np.max(abs(candidate[g.active]-temperature[g.active])))
            row={"iteration":it+1,"undamped_step_k":delta,"thermal_power_residual":audit.residual,
                 "sample_absorbed_w":sum(w*e.sample_w for w,e in zip(weights,em)),
                 "max_temperature_c":float(temperature[m.bed>0].max())}
            history.append(row)
            if progress:
                progress(row)
            cell_scale=max(float(source@g.volume),1e-6)
            if delta<cfg["solver"]["temperature_tolerance_k"] and abs(audit.residual)<cfg["solver"]["power_tolerance"] and audit.max_cell_residual_w/cell_scale<cfg["solver"]["power_tolerance"]:
                converged=True; break
            if it+1<int(cfg["solver"]["max_iterations"]):
                temperature += cfg["solver"]["relaxation"]*(candidate-temperature)
        # No post-convergence field refresh is hidden from the thermal residual.
        mean=lambda name: float(sum(w*getattr(e,name) for w,e in zip(weights,em)))
        bed_weights=m.bed*g.volume
        centre=g.closest([*cfg["geometry"].get("tube_offset_m",[0,0]),0],m.bed>0)
        wallmask=(m.tube_area>0)&(abs(g.xyz[:,2])==np.min(abs(g.xyz[:,2])))
        wallt=float(np.average(temperature[wallmask],weights=m.tube_area[wallmask]))
        reflected,wall_em,aperture,sample,quartz=map(mean,["reflected_w","wall_w","aperture_w","sample_w","quartz_w"])
        exact_cavity=np.pi*(cfg["geometry"]["cavity_diameter_m"]/2)**2*cfg["geometry"]["cavity_height_m"]
        warnings=["Research prototype: not experimentally validated; material table and probe geometry are illustrative.",
                  "Stair-step cylinder and subcell mixture require geometry/grid convergence before quantitative use.",
                  "Coax cable, antenna metal, and viewing windows are not resolved. Tube apertures use a first-order absorbing boundary when enabled.",
                  "Thermal quartz is opaque diffuse-gray; no spectral transmission or external natural-convection flow is solved.",
                  "Numerical fixed-point convergence does not establish physical stability or uniqueness."]
        if len(em)>1:
            warnings.append("Frequency dwell is averaged at frozen T; valid only when thermal changes per sweep are small. Four repeated identical sweeps do not multiply average power.")
        if self.transfer:
            warnings.extend(["Square channels are resolved thermally and use a 3D P2/P1 Navier-Stokes velocity field. Identical channels have equal imposed flow; inlet/outlet manifolds are not resolved.",
                "Flow uses constant density and viscosity at the specified reference temperature; thermal expansion does not feed back into momentum. Large gas temperature rises require a variable-density extension.",
                "Channel wall heat flux uses solid/gas harmonic conduction, without an imposed heat-transfer coefficient. Coarse EM cells homogenize subcell dielectric geometry; this is not channel-resolved electromagnetics."])
        elif cfg["gas"]["enabled"]:
            warnings.append("Darcy/LTE describes a homogenized packed bed with prescribed end temperatures/pressures; it does not resolve millimetre monolith channels or open-tube flow.")
        if cfg["em"]["wall_model"]=="pec":
            warnings.append("PEC model omits metal-wall microwave loss.")
        summary={"name":cfg["name"],"converged":converged,"experimentally_validated":False,
            "iterations":len(history),"elapsed_seconds":time.monotonic()-started,
            "power_w":{"incident":power,"reflected":reflected,"sample":sample,"quartz":quartz,"metal_wall":wall_em,"aperture_escape":aperture,
                "thermal_boundary":audit.boundary_w,"thermal_radiation":audit.radiation_w,"gas_enthalpy":audit.gas_w},
            "temperature_c":{"sample_centre_cell":float(temperature[centre]),"sample_volume_average":float(temperature@bed_weights/bed_weights.sum()),
                "sample_max":float(temperature[m.bed>0].max()),"sample_min":float(temperature[m.bed>0].min()),"quartz_midplane_average":wallt},
            "residuals":{"em_relative":max(abs(e.power_residual) for e in em),"thermal_relative":audit.residual,
                "thermal_max_cell_w":audit.max_cell_residual_w,"global_relative":float((power-reflected-wall_em-aperture-audit.boundary_w-audit.radiation_w-audit.gas_w)/max(power,1e-9)),
                "flow_mass_relative":flow.mass_residual if flow else None,
                "darcy_mass_relative":flow.mass_residual if flow and not self.transfer else None},
            "flow":{"pressure_drop_pa":flow.pressure_drop_pa,"particle_reynolds_max":flow.reynolds_particle_max if not self.transfer else None,
                "diagnostics":flow.diagnostics} if flow else None,
            "mesh":{"shape":list(g.shape),"active_cells":int(g.active.sum()),"electric_unknowns":len(self.maxwell.free),
                "sample_volume_relative_error":float(m.sample_volume/m.sample_volume_exact-1),"cavity_volume_relative_error":float(g.volume[g.active].sum()/exact_cavity-1),
                "effective_probe_length_m":self.maxwell.effective_probe_length,"em_shape":list(self.em_grid.shape),
                "separate_channel_thermal_mesh":self.transfer is not None},
            "spectrum":[{"frequency_hz":e.frequency_hz,"dwell_weight":float(w),"s11_real":e.s11.real,"s11_imag":e.s11.imag,
                "reflected_fraction":float(abs(e.s11)**2),"sample_w":e.sample_w,"quartz_w":e.quartz_w,"metal_wall_w":e.wall_w,"aperture_w":e.aperture_w} for w,e in zip(weights,em)],
            "history":history,"assumptions":warnings,"config":copy.deepcopy(cfg)}
        summary["config"]["solver"]["power_w"]=power
        fields={"temperature_c":temperature,"e2_v2_m2":self.display_em(sum(w*e.e2 for w,e in zip(weights,em))),
            "q_sample_w_m3":self.thermal_source(em,weights,"sample",p),"q_quartz_w_m3":self.thermal_source(em,weights,"quartz",p),
            "pressure_gauge_pa":flow.pressure_pa if flow else np.zeros(g.n),"sample_fraction":m.bed,"quartz_fraction":m.quartz,
            "active":g.active,"x_m":g.axes[0],"y_m":g.axes[1],"z_m":g.axes[2]}
        if self.transfer:
            fields["channel_fraction"]=m.fluid_fraction
            channel=m.fluid_fraction>.5
            summary["temperature_c"]["channel_gas_max"]=float(temperature[channel].max())
            if flow is not None:
                tref=cfg["gas"].get("reference_temperature_c",cfg["gas"]["inlet_c"])+273.15
                summary["flow"]["diagnostics"]["ideal_gas_density_change_estimate_max"]=float(np.max(abs(tref/(temperature[channel]+273.15)-1)))
                if summary["flow"]["diagnostics"]["ideal_gas_density_change_estimate_max"]>.1:
                    warnings.append("Estimated ideal-gas density variation exceeds 10%; the constant-density channel flow is not a quantitatively justified hot-flow prediction.")
            if flow is not None and flow.outlet_mass.sum()>0:
                summary["temperature_c"]["channel_gas_outlet_mass_average"]=float(flow.outlet_mass@temperature[flow.outlet_cells]/flow.outlet_mass.sum())
        velocity=cell_velocity(g,m.fluid_fraction,temperature,cfg,flow) if flow else np.zeros((g.n,3))
        fields["superficial_velocity_m_s"]=velocity
        fields["superficial_speed_m_s"]=np.linalg.norm(velocity,axis=1)
        if flow is not None and flow.resolved_velocity is not None:
            fields["ns_velocity_m_s"]=flow.resolved_velocity
            fields["ns_speed_m_s"]=np.linalg.norm(flow.resolved_velocity,axis=1)
        # Keep complex E separately for each frequency; do not add incoherent phasors.
        for i,e in enumerate(em):
            fields[f"electric_field_{i}_real_v_m"]=self.display_em(e.e_vector.real)
            fields[f"electric_field_{i}_imag_v_m"]=self.display_em(e.e_vector.imag)
        return summary,fields
