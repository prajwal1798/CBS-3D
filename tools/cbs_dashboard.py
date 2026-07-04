#!/usr/bin/env python3
"""
CBS3D++_SI live web dashboard.

A lightweight local web monitor for a running (or finished) CBS3D_SI case.
It reads the solver's <case>_residuals.csv (written every iteration) and renders
a smooth, auto-refreshing dashboard in the browser: residual histories, CG
behaviour, time-step, and a convergence status panel.

The solver stays a headless engine; this is a separate viewer, so it works
identically for a local run or a cluster run whose output directory is mounted
or synced locally.

Usage
-----
    python tools/cbs_dashboard.py output/LidDrivenCavity3D
    python tools/cbs_dashboard.py output/LidDrivenCavity3D --tol 1e-6 --port 8000
    # then open http://127.0.0.1:8000

If a directory is given, the newest *_residuals.csv inside it is used.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import time

from flask import Flask, jsonify, render_template_string

app = Flask(__name__)

STATE = {"csv": None, "tol": 1.0e-6, "target": 0, "case": "CBS3D_SI"}


# --------------------------------------------------------------------------- #
def find_csv(path: str) -> str:
    if os.path.isfile(path):
        return path
    cands = glob.glob(os.path.join(path, "*_residuals.csv"))
    if not cands:
        cands = glob.glob(os.path.join(path, "**", "*_residuals.csv"), recursive=True)
    if not cands:
        raise SystemExit(f"No *_residuals.csv found under {path}")
    return max(cands, key=os.path.getmtime)


def read_rows(csv_path: str):
    rows = []
    try:
        with open(csv_path, "r", newline="") as f:
            for row in csv.DictReader(f):
                rec = {}
                for k, v in row.items():
                    try:
                        rec[k] = float(v)
                    except (TypeError, ValueError):
                        rec[k] = float("nan")
                rows.append(rec)
    except FileNotFoundError:
        return [], 0.0
    mtime = os.path.getmtime(csv_path) if os.path.exists(csv_path) else 0.0
    return rows, mtime


def decimate(rows, max_points=4000):
    """Keep the plot light for very long runs; always keep the last point."""
    n = len(rows)
    if n <= max_points:
        return rows
    step = n // max_points + 1
    out = rows[::step]
    if out[-1] is not rows[-1]:
        out.append(rows[-1])
    return out


# --------------------------------------------------------------------------- #
@app.route("/data")
def data():
    rows, mtime = read_rows(STATE["csv"])
    if not rows:
        return jsonify({"empty": True})

    rows_plot = decimate(rows)
    last = rows[-1]

    def col(key):
        return [r.get(key, float("nan")) for r in rows_plot]

    it = col("iteration")
    series = {k: col(k) for k in
              ["u_rel", "v_rel", "w_rel", "p_rel", "T_rel",
               "cg_relative_l2", "cg_iterations", "dt", "velocity_rel_max"]}

    vel_res = max(last.get("u_rel", 0) or 0,
                  last.get("v_rel", 0) or 0,
                  last.get("w_rel", 0) or 0)
    stale = (time.time() - mtime) > 8.0
    converged = vel_res > 0 and vel_res < STATE["tol"]

    if converged:
        status = "CONVERGED"
    elif stale:
        status = "STOPPED"
    else:
        status = "RUNNING"

    return jsonify({
        "empty": False,
        "case": STATE["case"],
        "tol": STATE["tol"],
        "target": STATE["target"],
        "status": status,
        "iteration": it,
        "series": series,
        "last": {
            "iteration": int(last.get("iteration", 0)),
            "time": last.get("time", 0.0),
            "dt": last.get("dt", 0.0),
            "u_rel": last.get("u_rel", float("nan")),
            "v_rel": last.get("v_rel", float("nan")),
            "w_rel": last.get("w_rel", float("nan")),
            "p_rel": last.get("p_rel", float("nan")),
            "T_rel": last.get("T_rel", float("nan")),
            "vel_res": vel_res,
            "cg_iterations": int(last.get("cg_iterations", 0)),
            "cg_relative_l2": last.get("cg_relative_l2", float("nan")),
        },
    })


@app.route("/")
def index():
    return render_template_string(PAGE)


# --------------------------------------------------------------------------- #
PAGE = r"""
<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CBS3D++_SI Monitor</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  :root{
    --bg:#0b0f14; --panel:#121922; --panel2:#0e151d; --line:#1e2a36;
    --ink:#e6edf3; --mut:#7d8da0; --acc:#4cc2ff; --ok:#3fb950; --warn:#d29922; --bad:#f85149;
    --u:#4cc2ff; --v:#8b7cff; --w:#36d399; --p:#ffb454; --t:#f778ba; --cg:#7d8da0;
  }
  *{box-sizing:border-box} html,body{margin:0;height:100%}
  body{background:var(--bg);color:var(--ink);
       font:14px/1.4 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
  .wrap{max-width:1280px;margin:0 auto;padding:18px}
  header{display:flex;align-items:center;justify-content:space-between;
         padding:14px 18px;background:linear-gradient(180deg,#10171f,#0d141b);
         border:1px solid var(--line);border-radius:14px}
  .brand{display:flex;align-items:center;gap:14px}
  .brand .logo{font-weight:700;letter-spacing:3px;color:var(--acc);font-size:18px}
  .brand .sub{color:var(--mut);font-size:12px}
  .pill{padding:6px 14px;border-radius:999px;font-weight:700;font-size:12px;letter-spacing:1px;
        border:1px solid var(--line)}
  .pill.RUNNING{color:#0b0f14;background:var(--acc)}
  .pill.CONVERGED{color:#0b0f14;background:var(--ok)}
  .pill.STOPPED{color:#0b0f14;background:var(--warn)}
  .pill.WAIT{color:var(--mut)}
  .grid{display:grid;gap:14px;margin-top:14px}
  .cards{grid-template-columns:repeat(6,1fr)}
  @media(max-width:980px){.cards{grid-template-columns:repeat(3,1fr)}}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:14px 16px}
  .card .k{color:var(--mut);font-size:11px;letter-spacing:1px;text-transform:uppercase}
  .card .val{font-size:22px;font-weight:700;margin-top:6px}
  .card .val small{font-size:12px;color:var(--mut);font-weight:400}
  .plots{grid-template-columns:2fr 1fr}
  @media(max-width:980px){.plots{grid-template-columns:1fr}}
  .panel{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:10px 12px}
  .panel h3{margin:6px 8px 4px;font-size:12px;letter-spacing:1px;color:var(--mut);text-transform:uppercase;font-weight:600}
  .foot{color:var(--mut);font-size:11px;text-align:center;margin:16px 0 8px}
  .dot{display:inline-block;width:7px;height:7px;border-radius:50%;background:var(--ok);margin-right:6px;
       box-shadow:0 0 0 0 rgba(63,185,80,.6);animation:p 1.6s infinite}
  @keyframes p{0%{box-shadow:0 0 0 0 rgba(63,185,80,.5)}70%{box-shadow:0 0 0 8px rgba(63,185,80,0)}100%{box-shadow:0 0 0 0 rgba(63,185,80,0)}}
</style></head>
<body><div class="wrap">
  <header>
    <div class="brand">
      <div class="logo">CBS3D++_SI</div>
      <div class="sub" id="case">— monitor</div>
    </div>
    <div><span class="dot" id="dot"></span><span class="pill WAIT" id="pill">CONNECTING</span></div>
  </header>

  <div class="grid cards">
    <div class="card"><div class="k">Iteration</div><div class="val" id="c_it">—</div></div>
    <div class="card"><div class="k">Velocity residual</div><div class="val" id="c_vel">—</div></div>
    <div class="card"><div class="k">Pressure residual</div><div class="val" id="c_p">—</div></div>
    <div class="card"><div class="k">CG iters (last)</div><div class="val" id="c_cg">—</div></div>
    <div class="card"><div class="k">Time step Δt</div><div class="val" id="c_dt">—</div></div>
    <div class="card"><div class="k">Sim time</div><div class="val" id="c_t">—</div></div>
  </div>

  <div class="grid plots">
    <div class="panel"><h3>Residual history (semilog)</h3><div id="res" style="height:430px"></div></div>
    <div class="panel"><h3>Pressure CG iterations / step</h3><div id="cg" style="height:430px"></div></div>
  </div>

  <div class="foot" id="foot">waiting for data…</div>
</div>

<script>
const COL={u_rel:'#4cc2ff',v_rel:'#8b7cff',w_rel:'#36d399',p_rel:'#ffb454',T_rel:'#f778ba',cg_relative_l2:'#7d8da0'};
const LAYOUT_BASE={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',
  font:{color:'#cfd8e3',family:'ui-monospace,Menlo,Consolas,monospace',size:11},
  margin:{l:56,r:14,t:8,b:40},showlegend:true,
  legend:{orientation:'h',y:1.12,font:{size:10}},
  xaxis:{gridcolor:'#1a242f',zeroline:false,title:'iteration'},
};
function fmt(x){if(x===null||x===undefined||isNaN(x))return '—';
  if(x!==0&&(Math.abs(x)<1e-3||Math.abs(x)>=1e4))return x.toExponential(2);
  return (Math.round(x*1e4)/1e4).toString();}

async function tick(){
  let d; try{d=await (await fetch('/data',{cache:'no-store'})).json();}catch(e){return;}
  const pill=document.getElementById('pill'), dot=document.getElementById('dot');
  if(d.empty){document.getElementById('foot').textContent='CSV found but no rows yet…';return;}

  document.getElementById('case').textContent='— '+d.case;
  pill.textContent=d.status; pill.className='pill '+d.status;
  dot.style.display = d.status==='RUNNING'?'inline-block':'none';

  const L=d.last;
  document.getElementById('c_it').innerHTML = L.iteration + (d.target? ' <small>/ '+d.target+'</small>':'');
  document.getElementById('c_vel').innerHTML = fmt(L.vel_res)+' <small>tol '+d.tol.toExponential(0)+'</small>';
  document.getElementById('c_vel').style.color = (L.vel_res<d.tol)?'#3fb950':'#e6edf3';
  document.getElementById('c_p').textContent  = fmt(L.p_rel);
  document.getElementById('c_cg').innerHTML = L.cg_iterations + ' <small>rel '+fmt(L.cg_relative_l2)+'</small>';
  document.getElementById('c_dt').textContent = fmt(L.dt);
  document.getElementById('c_t').textContent  = fmt(L.time);

  const it=d.iteration, S=d.series;
  const names={u_rel:'u',v_rel:'v',w_rel:'w',p_rel:'p',T_rel:'T',cg_relative_l2:'CG rel'};
  const traces=Object.keys(names).filter(k=>S[k]&&S[k].some(v=>v>0)).map(k=>({
    x:it,y:S[k].map(v=>v>0?v:null),name:names[k],mode:'lines',
    line:{color:COL[k],width:k==='cg_relative_l2'?1:1.8,dash:k==='cg_relative_l2'?'dot':'solid'}}));
  const tolLine={x:[it[0],it[it.length-1]],y:[d.tol,d.tol],name:'tol',mode:'lines',
    line:{color:'#f85149',width:1,dash:'dash'},hoverinfo:'skip'};
  Plotly.react('res',[...traces,tolLine],
    Object.assign({},LAYOUT_BASE,{yaxis:{type:'log',gridcolor:'#1a242f',title:'relative residual'}}),
    {displayModeBar:false,responsive:true});

  Plotly.react('cg',[{x:it,y:S.cg_iterations,mode:'lines',line:{color:'#4cc2ff',width:1.6},name:'CG iters',fill:'tozeroy',fillcolor:'rgba(76,194,255,.12)'}],
    Object.assign({},LAYOUT_BASE,{showlegend:false,yaxis:{gridcolor:'#1a242f',title:'CG iterations'}}),
    {displayModeBar:false,responsive:true});

  document.getElementById('foot').textContent =
    'updated '+new Date().toLocaleTimeString()+'  •  '+it.length+' samples  •  auto-refresh 1.5s';
}
tick(); setInterval(tick,1500);
</script>
</body></html>
"""


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="output dir (or a *_residuals.csv file)")
    ap.add_argument("--tol", type=float, default=1.0e-6, help="steady-state velocity tolerance")
    ap.add_argument("--target", type=int, default=0, help="target iteration count (optional)")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1")
    args = ap.parse_args()

    STATE["csv"] = find_csv(args.path)
    STATE["tol"] = args.tol
    STATE["target"] = args.target
    STATE["case"] = os.path.basename(STATE["csv"]).replace("_residuals.csv", "")

    print(f"[cbs_dashboard] watching {STATE['csv']}")
    print(f"[cbs_dashboard] open http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
