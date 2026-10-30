#!/usr/bin/env python3
"""Interactive viewer (single HTML file, no external resources) of one or more
dual-engine runs: N, E, U of PPPHVA and GaiaVario and their differences, with
wheel zoom, drag pan, double-click reset and the robust sigma of the visible
window per series.

    python3 compare_html.py out.html "label1" ppphva1.txt gaiavario1.csv ["label2" ppphva2.txt gaiavario2.csv ...]
"""
import sys, json, math
import numpy as np

def unwrapper():
    st = {"prev": None, "off": 0.0}
    def f(s):
        if st["prev"] is not None and s < st["prev"] - 43200: st["off"] += 86400.0
        st["prev"] = s
        return round(s + st["off"], 1)
    return f

def load_run(ppp_path, gv_path):
    station = None; gv = {}
    uw = unwrapper()
    with open(gv_path) as f:
        for line in f:
            if line.startswith("# station ECEF"):
                station = tuple(float(v) for v in line.split()[3:6]); continue
            if line.startswith("#") or line.startswith("sod,"): continue
            p = line.rstrip("\n").split(",")
            gv[uw(float(p[0]))] = (float(p[2]), float(p[3]), float(p[4]), int(p[6]), p[15], p[16])
    x, y, z = station
    lon = math.atan2(y, x); r = math.hypot(x, y); e2 = 0.00669437999014
    lat = math.atan2(z, r * (1 - e2))
    for _ in range(5):
        Nn = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat) ** 2); h = r / math.cos(lat) - Nn
        lat = math.atan2(z, r * (1 - e2 * Nn / (Nn + h)))
    sp, cp, sl, cl = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])
    ppp = {}
    uw = unwrapper()
    with open(ppp_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.split()
            neu = R @ np.array([float(p[2]), float(p[3]), float(p[4])])
            ppp[uw(float(p[0]))] = (neu[0], neu[1], neu[2], int(p[1]))
    sods = sorted(set(ppp) | set(gv))
    t0 = sods[0]
    def mm(v): return None if v is None or not np.isfinite(v) else round(1000 * v, 2)
    rows = []
    for s in sods:
        a = ppp.get(s); b = gv.get(s)
        rows.append([round(s - t0, 1),
                     mm(a[0]) if a else None, mm(a[1]) if a else None, mm(a[2]) if a else None,
                     mm(b[0]) if b else None, mm(b[1]) if b else None, mm(b[2]) if b else None,
                     a[3] if a else None, b[3] if b else None, (b[4] or "") if b else "", (b[5] or "") if b else ""])
    return dict(sod0=t0, rows=rows)

out = sys.argv[1]
runs = []
args = sys.argv[2:]
for i in range(0, len(args), 3):
    runs.append(dict(label=args[i], **load_run(args[i + 1], args[i + 2])))

html = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>GaiaVario vs PPPHVA — interactive</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--grid:#e3e3e3;--p:#D85A30;--g:#1D9E75;--d:#534AB7}
@media (prefers-color-scheme:dark){:root{--bg:#161616;--fg:#e8e8e8;--mut:#9a9a9a;--grid:#2c2c2c}}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.4 system-ui,sans-serif}
header{padding:10px 16px;display:flex;gap:16px;align-items:center;flex-wrap:wrap;border-bottom:1px solid var(--grid)}
header b{font-weight:600}
select,button{font:inherit;padding:3px 8px}
.help{color:var(--mut);font-size:12.5px}
.panel{position:relative;margin:6px 16px}
.panel canvas{display:block;width:100%;height:150px;border:1px solid var(--grid);cursor:crosshair}
.lab{position:absolute;left:8px;top:4px;font-size:13px;background:var(--bg);padding:0 4px}
.stat{position:absolute;right:8px;top:4px;font-size:12.5px;color:var(--mut);background:var(--bg);padding:0 4px;text-align:right}
.sw{display:inline-block;width:10px;height:10px;margin-right:4px;vertical-align:middle}
#tip{position:fixed;pointer-events:none;background:var(--bg);border:1px solid var(--mut);padding:4px 6px;font-size:12px;display:none}
footer{padding:8px 16px;color:var(--mut);font-size:12.5px}
</style></head><body>
<header>
<b>GaiaVario vs PPPHVA, UPC10Z live, 10 Hz, mm per epoch</b>
<label>run <select id="run"></select></label>
<label><input type="checkbox" id="link" checked> link y-scale of N/E/U</label>
<button id="reset">reset zoom</button>
<span class="help">wheel = zoom time around the cursor · drag = pan · double-click = reset · σ = 1.4826·MAD of the visible window</span>
</header>
<div id="panels"></div>
<div id="tip"></div>
<footer id="foot"></footer>
<script>
const RUNS = __RUNS__;
const PANELS = [
 {name:"N", series:[{col:"var(--p)",idx:1,lab:"PPPHVA"},{col:"var(--g)",idx:4,lab:"GaiaVario"}]},
 {name:"E", series:[{col:"var(--p)",idx:2,lab:"PPPHVA"},{col:"var(--g)",idx:5,lab:"GaiaVario"}]},
 {name:"U", series:[{col:"var(--p)",idx:3,lab:"PPPHVA"},{col:"var(--g)",idx:6,lab:"GaiaVario"}]},
 {name:"ΔN = GaiaVario − PPPHVA", series:[{col:"var(--d)",diff:[4,1],lab:"Δ"}]},
 {name:"ΔE", series:[{col:"var(--d)",diff:[5,2],lab:"Δ"}]},
 {name:"ΔU", series:[{col:"var(--d)",diff:[6,3],lab:"Δ"}]},
];
let run = RUNS[0], rows = run.rows, x0 = 0, x1 = 1, dragging = null;
const sel = document.getElementById("run");
RUNS.forEach((r,i)=>{const o=document.createElement("option");o.value=i;o.textContent=r.label;sel.appendChild(o);});
sel.onchange = ()=>{run=RUNS[sel.value];rows=run.rows;resetX();drawAll();};
function css(v){return getComputedStyle(document.documentElement).getPropertyValue(v.slice(4,-1)).trim();}
function val(row,s){ if(s.diff){const a=row[s.diff[0]],b=row[s.diff[1]];return (a==null||b==null)?null:a-b;} return row[s.idx]; }
function resetX(){x0=rows[0][0];x1=rows[rows.length-1][0];}
resetX();
const panels=document.getElementById("panels");
const objs=PANELS.map(p=>{
  const d=document.createElement("div");d.className="panel";
  const c=document.createElement("canvas");d.appendChild(c);
  const l=document.createElement("div");l.className="lab";l.innerHTML=p.name+" &nbsp; "+p.series.map(s=>`<span class="sw" style="background:${s.col}"></span>${s.lab}`).join(" &nbsp;");d.appendChild(l);
  const st=document.createElement("div");st.className="stat";d.appendChild(st);
  panels.appendChild(d);
  c.addEventListener("wheel",e=>{e.preventDefault();const r=c.getBoundingClientRect();const f=(e.clientX-r.left)/r.width;const xc=x0+f*(x1-x0);const k=e.deltaY>0?1.25:0.8;let nx0=xc-(xc-x0)*k,nx1=xc+(x1-xc)*k;const lo=rows[0][0],hi=rows[rows.length-1][0];if(nx1-nx0<2){return;}if(nx0<lo)nx0=lo;if(nx1>hi)nx1=hi;x0=nx0;x1=nx1;drawAll();},{passive:false});
  c.addEventListener("mousedown",e=>{dragging={x:e.clientX,x0:x0,x1:x1,w:c.getBoundingClientRect().width};});
  window.addEventListener("mousemove",e=>{if(!dragging)return;const dx=(e.clientX-dragging.x)/dragging.w*(dragging.x1-dragging.x0);const lo=rows[0][0],hi=rows[rows.length-1][0];let nx0=dragging.x0-dx,nx1=dragging.x1-dx;if(nx0<lo){nx1+=lo-nx0;nx0=lo;}if(nx1>hi){nx0-=nx1-hi;nx1=hi;}x0=nx0;x1=nx1;drawAll();});
  window.addEventListener("mouseup",()=>dragging=null);
  c.addEventListener("dblclick",()=>{resetX();drawAll();});
  c.addEventListener("mousemove",e=>tip(e,c,p));
  c.addEventListener("mouseleave",()=>document.getElementById("tip").style.display="none");
  return {p,c,st};
});
function mad(v){if(v.length<5)return NaN;const s=[...v].sort((a,b)=>a-b);const m=s[s.length>>1];const d=v.map(x=>Math.abs(x-m)).sort((a,b)=>a-b);return 1.4826*d[d.length>>1];}
function visible(){let i0=0,i1=rows.length-1;while(i0<rows.length&&rows[i0][0]<x0)i0++;while(i1>0&&rows[i1][0]>x1)i1--;return [i0,i1];}
let sharedLim=null;
function drawAll(){
  const [i0,i1]=visible();
  const link=document.getElementById("link").checked;
  // common y limit for N/E/U when linked: 99.5th percentile of |values|
  sharedLim=null;
  if(link){const all=[];for(let k=0;k<3;k++)for(const s of PANELS[k].series)for(let i=i0;i<=i1;i++){const v=val(rows[i],s);if(v!=null)all.push(Math.abs(v));}all.sort((a,b)=>a-b);sharedLim=all.length?Math.max(1,all[Math.floor(all.length*0.995)]*1.3):5;}
  objs.forEach(o=>draw(o,i0,i1));
  document.getElementById("foot").textContent=`${run.label} — visible ${i1-i0+1} epochs, ${(x1-x0).toFixed(1)} s, from SOD ${(run.sod0+x0).toFixed(1)} to ${(run.sod0+x1).toFixed(1)} s. ${rows.length} epochs in the run.`;
}
function draw(o,i0,i1){
  const c=o.c,dpr=window.devicePixelRatio||1,W=c.clientWidth,H=c.clientHeight;
  if(c.width!==W*dpr||c.height!==H*dpr){c.width=W*dpr;c.height=H*dpr;}
  const g=c.getContext("2d");g.setTransform(dpr,0,0,dpr,0,0);g.clearRect(0,0,W,H);
  const ser=o.p.series.map(s=>({s,v:[]}));
  for(let i=i0;i<=i1;i++)for(const q of ser){const v=val(rows[i],q.s);if(v!=null)q.v.push(v);}
  let lim;
  if(sharedLim&&!o.p.series[0].diff)lim=sharedLim;else{const all=[];for(const q of ser)for(const v of q.v)all.push(Math.abs(v));all.sort((a,b)=>a-b);lim=all.length?Math.max(1,all[Math.floor(all.length*0.995)]*1.3):5;}
  const ym=H/2,sx=W/(x1-x0),sy=(H/2-6)/lim;
  g.strokeStyle=css("var(--grid)");g.lineWidth=1;
  for(const yy of [-lim*0.5,0,lim*0.5]){g.beginPath();g.moveTo(0,ym-yy*sy);g.lineTo(W,ym-yy*sy);g.stroke();}
  g.fillStyle=css("var(--mut)");g.font="11px system-ui";g.fillText(`+${lim.toFixed(1)} mm`,W-70,12);g.fillText(`−${lim.toFixed(1)} mm`,W-70,H-4);
  const pts=i1-i0;const asPoints=(pts/W)>3;   // dense: draw min/max bars per pixel column
  for(const q of ser){
    g.strokeStyle=css(q.s.col);g.lineWidth=1;g.beginPath();
    if(asPoints){let col=-1,mn=0,mx=0,any=false;for(let i=i0;i<=i1;i++){const v=val(rows[i],q.s);const px=Math.floor((rows[i][0]-x0)*sx);if(px!==col){if(any){g.moveTo(col+0.5,ym-mx*sy);g.lineTo(col+0.5,ym-mn*sy+0.01);}col=px;any=false;}if(v==null)continue;if(!any){mn=mx=v;any=true;}else{if(v<mn)mn=v;if(v>mx)mx=v;}}if(any){g.moveTo(col+0.5,ym-mx*sy);g.lineTo(col+0.5,ym-mn*sy+0.01);}}
    else{let pen=false;for(let i=i0;i<=i1;i++){const v=val(rows[i],q.s);if(v==null){pen=false;continue;}const px=(rows[i][0]-x0)*sx,py=ym-Math.max(-lim*1.05,Math.min(lim*1.05,v))*sy;if(!pen){g.moveTo(px,py);pen=true;}else g.lineTo(px,py);}}
    g.stroke();
    if(!asPoints&&pts<400){g.fillStyle=css(q.s.col);for(let i=i0;i<=i1;i++){const v=val(rows[i],q.s);if(v==null)continue;g.fillRect((rows[i][0]-x0)*sx-1.5,ym-v*sy-1.5,3,3);}}
  }
  o.st.innerHTML=ser.map(q=>`<span class="sw" style="background:${q.s.col}"></span>σ ${mad(q.v).toFixed(2)} mm, n ${q.v.length}`).join(" &nbsp; ");
}
function tip(e,c,p){
  const r=c.getBoundingClientRect();const t=x0+(e.clientX-r.left)/r.width*(x1-x0);
  let lo=0,hi=rows.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(rows[m][0]<t)lo=m;else hi=m;}
  const row=Math.abs(rows[lo][0]-t)<Math.abs(rows[hi][0]-t)?rows[lo]:rows[hi];
  const d=document.getElementById("tip");d.style.display="block";d.style.left=(e.clientX+12)+"px";d.style.top=(e.clientY+12)+"px";
  const f=v=>v==null?"—":v.toFixed(2);
  d.innerHTML=`SOD ${(run.sod0+row[0]).toFixed(1)} s<br>PPPHVA N ${f(row[1])} E ${f(row[2])} U ${f(row[3])} mm, sats ${row[7]??"—"}<br>GaiaVario N ${f(row[4])} E ${f(row[5])} U ${f(row[6])} mm, sats ${row[8]??"—"}`+(row[9]?`<br>rejected ${row[9]}`:"")+(row[10]?`<br>slip ${row[10]}`:"");
}
document.getElementById("reset").onclick=()=>{resetX();drawAll();};
document.getElementById("link").onchange=drawAll;
window.addEventListener("resize",drawAll);
drawAll();
</script></body></html>
"""
open(out, "w").write(html.replace("__RUNS__", json.dumps(runs, separators=(",", ":"))))
print("wrote", out, "runs:", [r["label"] for r in runs], "epochs:", [len(r["rows"]) for r in runs])
