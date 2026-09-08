"""Offline power-temperature chart built from calculated sweep rows."""
import html
from pathlib import Path


def power_chart(path, rows):
    if not rows:
        return
    xmax=max(r["incident_w"] for r in rows) or 1
    ymin=min(r["quartz_midplane_average"] for r in rows)
    ymax=max(r["sample_max"] for r in rows)
    pad=max((ymax-ymin)*.12,1);ymin-=pad;ymax+=pad
    x=lambda v:70+650*v/xmax
    y=lambda v:360-290*(v-ymin)/(ymax-ymin)
    parts=['<svg viewBox="0 0 780 420" role="img" aria-label="Calculated incident power versus temperature">']
    for i in range(6):
        px=xmax*i/5;py=ymin+(ymax-ymin)*i/5
        parts.append(f'<path d="M{x(px):.2f},70V360 M70,{y(py):.2f}H720" stroke="#dfe6e0" fill="none"/>')
        parts.append(f'<text x="{x(px):.2f}" y="385" text-anchor="middle">{px:.0f}</text><text x="60" y="{y(py)+4:.2f}" text-anchor="end">{py:.1f}</text>')
    for key,label,colour,dash in [('sample_centre_cell','Centre cell','#ad4e36',''),('sample_volume_average','Bed average','#1a8f87','5 3'),('quartz_midplane_average','Quartz midplane','#42568c','2 3')]:
        xy=' '.join(f'{x(r["incident_w"]):.2f},{y(r[key]):.2f}' for r in rows)
        parts.append(f'<polyline points="{xy}" fill="none" stroke="{colour}" stroke-width="2" stroke-dasharray="{dash}"/>')
        for r in rows:
            parts.append(f'<circle cx="{x(r["incident_w"]):.2f}" cy="{y(r[key]):.2f}" r="3" fill="{colour}"><title>{label}: {r["incident_w"]:g} W, {r[key]:.3f} °C</title></circle>')
    parts.append('<text x="400" y="412" text-anchor="middle">Incident port power / W</text><text x="16" y="200" text-anchor="middle" transform="rotate(-90,16,200)">Temperature / °C</text></svg>')
    table=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in [r['incident_w'],round(r['sample_absorbed_w'],6),round(r['sample_centre_cell'],3),round(r['quartz_midplane_average'],3),r['converged']])+'</tr>' for r in rows)
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Microwave3D power–temperature sweep</title><style>body{font-family:system-ui;max-width:900px;margin:40px auto;padding:20px;color:#173037;background:#f2f4f1}h1{font-size:28px}p{line-height:1.5}svg{width:100%;background:white}svg text{font:12px system-ui;fill:#50645d}table{width:100%;border-collapse:collapse;font-size:13px}td,th{padding:10px;border-bottom:1px solid #ccd5ce;text-align:right}.note{padding:12px;background:#fff5d9}a{color:#176e69}</style><h1>Power and temperature</h1><p class="note">Illustrative simulation, not experimental validation. These temperatures depend on the assumed probe, material table, frequency, and coarse geometry. Lines connect solved points and do not establish stable branches between them.</p><p>Centre cell: rust | Bed average: teal | Quartz midplane: blue</p>'''+''.join(parts)+'<table><thead><tr><th>Incident / W</th><th>Sample absorption / W</th><th>Centre / °C</th><th>Quartz / °C</th><th>Converged</th></tr></thead><tbody>'+table+'</tbody></table><p><a href="power-sweep.csv">Numeric sweep data</a></p></html>'
    Path(path).write_text(page,encoding='utf-8')
