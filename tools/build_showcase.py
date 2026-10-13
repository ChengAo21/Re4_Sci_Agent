"""Reproduce the website gallery from archived Re4 Python solutions.

Run from the repository root: python tools/build_showcase.py
Use --case cylinder to rerun one case. All numerical algorithms, grids and final
times are those in the selected archive; only plotting and snapshot collection
are added. Generated figures are examples, not aggregate benchmark statistics.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
BASE = "src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/"
SPECS = {
    "cylinder": dict(title="Cylinder flow", subtitle="Re = 100 · staggered-grid projection", category="Flow",
        source="src/ipynb_ver/sample_results/_GPT5mini_/cylinder_for_Reproduction/sample1/solu_rev_2.dat",
        description="A fresh execution of the archived final revision. The animation contains actual velocity snapshots from the time integration; the original grid and final time are retained."),
    "sod": dict(title="Sod shock tube", subtitle="Compressible flow · t = 0.2", category="PDE",
        source=BASE+"PDEbench/SodShock/GPT_ans/sample1/solu_rev_2.dat",
        description="Density, velocity and pressure from one archived GPT final-revision solver, rerun on its original 400-cell grid."),
    "cavity": dict(title="Lid-driven cavity", subtitle="Artificial compressibility · Re = 100", category="Flow",
        source=BASE+"PDEbench/Lid_Driven/GPT_ans/sample1/solu_rev_2.dat",
        description="Velocity magnitude and streamlines from the archived artificial-compressibility implementation. The run's own stopping criterion and diagnostics are preserved."),
    "unsteady_ns": dict(title="Unsteady Navier–Stokes", subtitle="Pressure projection · t = 0.5", category="Flow",
        source=BASE+"PDEbench/Unsteady_NS/GPT_ans/sample1/solu_rev_2.dat",
        description="The two velocity components and pressure of the archived generated solver at its original final time. Displayed fields are illustrative outputs, not a claim of validation against the reference."),
    "burgers": dict(title="Viscous Burgers equation", subtitle="Nonlinear transport · t = 1", category="PDE",
        source=BASE+"PDEbench/Burgers_equ/GPT_ans/sample1/solu_rev_2.dat",
        description="Space–time evolution and the final velocity profile from an archived MUSCL-type limited solver with implicit diffusion."),
    "poisson": dict(title="Poisson / Laplace equation", subtitle="Four-hole domain · sparse linear solve", category="PDE",
        source=BASE+"PDEbench/Poisson_equ/GPT_ans/solu_rev_2.dat",
        description="The archived Poisson task is a zero-source Laplace solve on a square with four circular holes. Hole interiors are masked in the presentation."),
    "helmholtz": dict(title="Helmholtz-type equation", subtitle="Screened elliptic field · irregular domain", category="PDE",
        source=BASE+"PDEbench/Hemholtz/GPT_ans/solu_rev_2.dat",
        description="A screened elliptic solve from the archive (named Poisson–Boltzmann in the generated code), on a square with four differently sized holes."),
    "hilbert": dict(title="Hilbert linear systems", subtitle="Conditioning · direct methods · regularization", category="Linear algebra",
        source=BASE+"Hilbert/GPT_ans/sample1/solu_rev_2.dat",
        description="Forward errors and condition numbers from one generated implementation. Failed decompositions remain visible as missing values rather than being silently repaired."),
    "keyhole": dict(title="Keyhole scaling", subtitle="Dimensionless-group search · archived dataset", category="Dimensional analysis",
        source=BASE+"Keyhole/GPT_ans/solu_rev_2.dat",
        description="The archived rational-exponent search and log–log regression rerun on dataset_keyhole.csv. R² describes the fitted dataset; it is not a held-out score."),
    "rayleigh": dict(title="Rayleigh–Bénard scaling", subtitle="Candidate groups · log–log regression", category="Dimensional analysis",
        source=BASE+"Rayleigh/GPT_ans/solu_rev_2.dat",
        description="The archived dimensionless-group candidate search rerun on dataset_rayben.csv. R² describes the fitted dataset; it is not a held-out score."),
}


def read_code(path: Path) -> str:
    blocks = re.findall(r"```python\s*\n(.*?)```", path.read_text(encoding="utf-8-sig"), re.S)
    if len(blocks) != 1:
        raise ValueError(f"Expected one Python block, found {len(blocks)} in {path}")
    return blocks[0]


def theme():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":11,
        "axes.titlesize":13, "axes.labelsize":11, "axes.titleweight":"bold",
        "axes.spines.top":False, "axes.spines.right":False,
        "axes.edgecolor":"#bcc8d6", "axes.labelcolor":"#203147",
        "text.color":"#13253d", "xtick.color":"#50647a", "ytick.color":"#50647a",
        "figure.facecolor":"white", "axes.facecolor":"#f7f9fc",
        "grid.color":"#dbe3ed", "grid.alpha":.65, "savefig.dpi":160,
        "figure.dpi":100, "axes.prop_cycle":matplotlib.cycler(color=["#007f8b","#ef8c43","#6f60b6","#c34663","#3474af","#629f5f"]),
        "legend.frameon":False})


def worker(case_id: str, out: Path):
    import numpy as np
    theme()
    import matplotlib.pyplot as plt
    spec = SPECS[case_id]
    source = ROOT/spec["source"]
    code = read_code(source)
    work = ROOT/".showcase_work"/case_id
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    (out/"logs").mkdir(exist_ok=True)
    for name in ("dataset_keyhole.csv", "dataset_rayben.csv"):
        found = source.parent.parent/name
        if found.exists(): shutil.copy2(found,work/name)

    capture = {}
    originals = []
    def collect(frame, event, arg):
        if event == "return" and frame.f_code.co_filename == str(source) and frame.f_code.co_name in {"main","burgers_solver","run_experiments","solve_system","solve_laplace"}:
            capture.update(frame.f_locals.copy())
    def show(*args, **kwargs):
        for n in plt.get_fignums():
            path=work/f"shown_{len(originals):02d}.png"
            plt.figure(n).savefig(path,bbox_inches="tight")
            originals.append(path)
        plt.close("all")
    plt.show = show
    snapshots=[]
    def snapshot(u,v,xp,yp,t):
        uc=.5*(u[:-1]+u[1:]); vc=.5*(v[:,:-1]+v[:,1:])
        snapshots.append((float(t),uc.copy(),vc.copy(),xp.copy(),yp.copy()))

    ns={"__name__":"__main__", "__file__":str(source), "_showcase_snapshot":snapshot}
    if case_id=="cylinder":
        target="        if (step % progress_interval) == 0 or step == nsteps:"
        assert code.count(target)==1
        code=code.replace(target,"        if step == 1 or step == nsteps or step % max(1, nsteps // 60) == 0:\n            _showcase_snapshot(u, v, x_p, y_p, step * dt)\n"+target)
        code=code.replace("    for step in range(1, nsteps + 1):", "    _showcase_snapshot(u, v, x_p, y_p, 0.0)\n    for step in range(1, nsteps + 1):")
    old_cwd=Path.cwd(); os.chdir(work)
    log=io.StringIO(); start=time.perf_counter()
    try:
        with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            sys.setprofile(collect)
            exec(compile(code,str(source),"exec"),ns)
    finally:
        sys.setprofile(None);os.chdir(old_cwd)
        (out/"logs"/f"{case_id}.txt").write_text(log.getvalue(),encoding="utf-8")
    elapsed=time.perf_counter()-start
    values={**ns,**capture}
    plt.close("all")
    metrics={"Solver execution":f"{elapsed:.2f} s"}
    image=out/f"{case_id}.png"
    images=[]

    def finish(fig, title, subtitle):
        fig.suptitle(title,x=.06,ha="left",fontsize=19,fontweight="bold",color="#12233d")
        fig.text(.06,.9,subtitle,fontsize=10,color="#587086")
        fig.subplots_adjust(top=.81,bottom=.14,wspace=.32)
        fig.savefig(image,dpi=170,bbox_inches="tight",facecolor="white")
        plt.close(fig)

    if case_id=="cylinder":
        from matplotlib.patches import Circle
        from PIL import Image
        metrics.update({"Grid":f"{values['Nx']} × {values['Ny']}", "Final time":f"{values['T']:.0f}","Reynolds number":"100", "Snapshots":str(len(snapshots))})
        ts=[s[0] for s in snapshots]; frames=[]
        finalspeed=np.hypot(snapshots[-1][1],snapshots[-1][2])
        maxspeed=max(float(np.nanmax(np.hypot(s[1],s[2]))) for s in snapshots)
        for t,uc,vc,xp,yp in snapshots:
            fig,ax=plt.subplots(figsize=(12,6))
            im=ax.pcolormesh(xp,yp,np.hypot(uc,vc).T,cmap="viridis",vmin=0,vmax=maxspeed,shading="auto",rasterized=True)
            ax.add_patch(Circle((0,0),.5,facecolor="#ecf1f7",edgecolor="#0f2941",lw=1.1,zorder=5))
            ax.set(xlim=(-3,8),ylim=(-3,3),xlabel="x / D",ylabel="y / D",aspect="equal")
            ax.text(.015,.96,f"t = {t:05.2f}",transform=ax.transAxes,ha="left",va="top",color="#10243b",bbox=dict(boxstyle="round,pad=.35",facecolor="white",alpha=.95,edgecolor="none"),fontsize=12)
            cbar=fig.colorbar(im,ax=ax,pad=.025,fraction=.035);cbar.set_label("Speed / inlet speed")
            fig.suptitle("Flow past a cylinder",x=.075,ha="left",fontsize=21,fontweight="bold",color="#12233d")
            fig.text(.075,.91,"Re = 100  ·  220 × 120 cells  ·  actual solver snapshots",fontsize=11,color="#587086")
            fig.text(.075,.04,"Archived GPT-5 mini · final revision · original numerical settings",fontsize=9,color="#728397")
            fig.subplots_adjust(top=.84,bottom=.12,left=.075,right=.91)
            buffer=io.BytesIO();fig.savefig(buffer,format="png",dpi=100,facecolor="white");buffer.seek(0)
            frame=Image.open(buffer).convert("RGB");frames.append(frame.copy())
            if t==ts[-1]:fig.savefig(image,dpi=170,facecolor="white")
            plt.close(fig)
        frames[0].save(out/"cylinder.gif",save_all=True,append_images=frames[1:],duration=85,loop=0,optimize=True)
        np.savez_compressed(out/"cylinder_snapshots.npz",times=np.array(ts),u=np.stack([s[1] for s in snapshots]).astype(np.float32),v=np.stack([s[2] for s in snapshots]).astype(np.float32),x=snapshots[-1][3],y=snapshots[-1][4])
        metrics["Maximum final speed"]=f"{np.max(finalspeed):.4f}"
        div=values['divergence'](values['u'],values['v'],values['Nx'],values['Ny'],values['dx'],values['dy'])
        xx,yy=np.meshgrid(values['x_p'],values['y_p'],indexing='ij')
        fluid=xx**2+yy**2 > .5**2
        metrics["Max divergence, all cells"]=f"{np.max(np.abs(div)):.3e}"
        metrics["Max divergence, outside cylinder"]=f"{np.max(np.abs(div[fluid])):.3e}"
    elif case_id=="sod":
        fig,axs=plt.subplots(1,3,figsize=(12,4.4))
        for ax,key,label in zip(axs,["rho","u","p"],["Density ρ","Velocity u","Pressure p"]):
            ax.plot(values['x'],values[key],lw=2.1);ax.set(xlabel="x",ylabel=label);ax.grid(True)
        finish(fig,"Sod shock tube","400 cells  ·  t = 0.2  ·  archived final revision")
        metrics.update({"Grid":"400 cells","Final time":f"{values['t']:.3f}","Density range":f"{values['rho'].min():.3f}–{values['rho'].max():.3f}"})
    elif case_id=="cavity":
        fig,ax=plt.subplots(figsize=(7.5,6.7));u=values['u'];v=values['v'];x=values['x'];y=values['y']
        field=ax.pcolormesh(x,y,np.hypot(u,v).T,cmap="viridis",shading="auto")
        ax.streamplot(x,y,u.T,v.T,color="#d7e4ed",density=1.65,linewidth=.6,arrowsize=.65)
        ax.set(aspect="equal",xlabel="x",ylabel="y");fig.colorbar(field,ax=ax,label="Speed",pad=.03)
        finish(fig,"Lid-driven cavity","Artificial compressibility  ·  81 × 81 points  ·  Re = 100")
        metrics.update({"Grid":"81 × 81","Iterations":str(values['it']+1),"Final max update":f"{values['max_diff']:.3e}","Tolerance":f"{values['tol']:.0e}"})
    elif case_id=="unsteady_ns":
        fig,axs=plt.subplots(1,3,figsize=(13.5,4.8))
        for ax,key,label in zip(axs,["u","v","p"],["Velocity u","Velocity v","Pressure"]):
            field=ax.pcolormesh(values['x'],values['y'],values[key].T,cmap="coolwarm",shading="auto")
            ax.set(title=label,xlabel="x",ylabel="y",aspect="equal");fig.colorbar(field,ax=ax,orientation="horizontal",pad=.24,fraction=.08)
        finish(fig,"Unsteady Navier–Stokes","81 × 41 points  ·  t = 0.5  ·  archived final revision")
        metrics.update({"Grid":"81 × 41","Final time":"0.5","Time step":"0.002","Final divergence max":f"{values['div_max']:.3e}"})
    elif case_id=="burgers":
        fig,axs=plt.subplots(1,2,figsize=(12,5.2))
        im=axs[0].pcolormesh(values['x'],values['T'],values['U'],shading="auto",cmap="coolwarm");axs[0].set(xlabel="x",ylabel="t",title="Space–time evolution");fig.colorbar(im,ax=axs[0],label="u")
        axs[1].plot(values['x'],values['u0'],ls="--",color="#9baabd",label="Initial condition")
        axs[1].plot(values['x'],values['U'][-1],lw=2,label=f"t = {values['T'][-1]:.3f}");axs[1].set(xlabel="x",ylabel="u",title="Final profile");axs[1].legend();axs[1].grid(True)
        finish(fig,"Viscous Burgers equation","201 points  ·  ν = 0.01 / π  ·  archived final revision")
        metrics.update({"Grid":"201 points","Final time":f"{values['t']:.3f}","Saved states":str(len(values['T']))})
    elif case_id in {"poisson","helmholtz"}:
        fig,ax=plt.subplots(figsize=(7.5,6.6));field=values['U'] if case_id=="poisson" else values['u']
        im=ax.pcolormesh(values['X'],values['Y'],field,cmap="viridis",shading="auto",rasterized=True)
        ax.set(aspect="equal",xlabel="x",ylabel="y");fig.colorbar(im,ax=ax,label="u",pad=.03)
        finish(fig,spec['title'],"201 × 201 points  ·  masked circular holes  ·  archived final revision")
        metrics.update({"Grid":"201 × 201","Field range":f"{np.nanmin(field):.4f}–{np.nanmax(field):.4f}"})
        if 'res_norm' in values:metrics['Linear residual L2']=f"{values['res_norm']:.3e}"
    elif case_id=="hilbert":
        fig,axs=plt.subplots(1,2,figsize=(12.5,5.5));n=values['ns']
        for key,label in [('errors_lu','LU'),('errors_chol','Cholesky'),('errors_qr','Householder QR'),('errors_lu_reg','LU + regularization'),('errors_chol_reg','Cholesky + regularization'),('errors_qr_reg','QR + regularization')]:
            axs[0].semilogy(n,values[key],marker='o',ms=3,lw=1.5,label=label)
        axs[0].set(xlabel="Matrix size n",ylabel="Forward error, L∞",title="Computed solution errors");axs[0].legend(fontsize=8,loc="lower right");axs[0].grid(True)
        axs[1].semilogy(n,values['cond_numbers'],marker='o',lw=2);axs[1].set(xlabel="Matrix size n",ylabel="Estimated condition number",title="Ill-conditioning");axs[1].grid(True)
        finish(fig,"Hilbert linear systems","n = 5, 10, …, 50  ·  λ = 10⁻¹²  ·  failed factorizations remain missing")
        metrics.update({"Matrix sizes":"5–50","Regularization":"1e-12","Failed original Cholesky":str(int(np.sum(~np.isfinite(values['errors_chol']))))})
    elif case_id=="keyhole":
        fig,ax=plt.subplots(figsize=(9.5,6));xx=np.exp(values['log_Pi_final']);yy=values['e_star'];model=values['model_final'];order=np.argsort(xx)
        ax.scatter(xx,yy,s=24,alpha=.5,edgecolors="none",label="Archived dataset")
        ax.plot(xx[order],np.exp(model.predict(values['log_Pi_final'][order,None])),lw=2,color="#ef8c43",label="Log–log fit")
        ax.set(xscale="log",yscale="log",xlabel="Selected dimensionless group Π",ylabel="Keyhole aspect ratio e*");ax.grid(True,which="both");ax.legend()
        finish(fig,"Keyhole dimensionless scaling",f"{len(yy)} samples  ·  in-sample R² = {values['r2_final']:.4f}  ·  archived search")
        metrics.update({"Samples":str(len(yy)),"In-sample R²":f"{values['r2_final']:.4f}"})
    elif case_id=="rayleigh":
        if values['best_log_dimless'] is None:raise RuntimeError("Archived candidate search did not pass its own R² threshold")
        xx=np.exp(values['best_log_dimless']).ravel();yy=np.exp(values['best_log_Nu']).ravel();model=values['best_model'];order=np.argsort(xx)
        fig,ax=plt.subplots(figsize=(9.5,6));ax.scatter(xx,yy,s=23,alpha=.5,edgecolors="none",label="Archived dataset")
        ax.plot(xx[order],np.exp(model.predict(values['best_log_dimless'][order])).ravel(),color="#ef8c43",lw=2,label="Log–log fit")
        ax.set(xscale="log",yscale="log",xlabel="Selected dimensionless group Π",ylabel="Nusselt number Nu");ax.grid(True,which="both");ax.legend()
        finish(fig,"Rayleigh–Bénard dimensionless scaling",f"{len(yy)} samples  ·  in-sample R² = {values['best_r2']:.4f}  ·  archived search")
        metrics.update({"Samples":str(len(yy)),"Candidates":str(len(values['candidates'])),"In-sample R²":f"{values['best_r2']:.4f}"})

    item={"id":case_id,**spec,"image":image.relative_to(ROOT).as_posix(),"metrics":metrics,
        "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"numerical_modifications":"None; original algorithm, grid and final time retained.",
        "execution_log":(out/"logs"/f"{case_id}.txt").relative_to(ROOT).as_posix()}
    if case_id=="cylinder":item['gif']=(out/'cylinder.gif').relative_to(ROOT).as_posix()
    if case_id=="cavity" and values['max_diff']>=values['tol']:
        item['note']="The archived solver reached its iteration limit without satisfying its update tolerance. This is disclosed rather than corrected."
    if case_id in {"cylinder","unsteady_ns"}:
        item['note']="A successful execution is not an accuracy certification. Divergence diagnostics are reported from this run; no solver correction or reference-error claim is added."
    if case_id=="cylinder":
        item['note']="The unchanged archived solver has nonzero discrete divergence (see diagnostics). This is an illustrative reproduction, not a validated incompressible reference solution. Snapshot arrays are stored as float32 for display; the solver runs in its original precision."
    (out/f"{case_id}.json").write_text(json.dumps(item,ensure_ascii=False,indent=2),encoding="utf-8")
    print(f"{case_id}: {elapsed:.2f}s -> {image.relative_to(ROOT)}",flush=True)


def build(case=None):
    out=ROOT/"assets/results";out.mkdir(parents=True,exist_ok=True)
    ids=[case] if case else list(SPECS)
    failures=[]
    for case_id in ids:
        try:
            result=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--worker",case_id],cwd=ROOT,timeout=600,text=True,encoding="utf-8",errors="replace",capture_output=True)
            print(result.stdout,flush=True)
            if result.returncode:
                failures.append({"id":case_id,"error":result.stderr[-3000:]});print(result.stderr,flush=True)
        except subprocess.TimeoutExpired:
            failures.append({"id":case_id,"error":"Execution exceeded 600 s limit"})
    items=[]
    for key in SPECS:
        path=out/f"{key}.json"
        if path.exists():items.append(json.loads(path.read_text(encoding="utf-8")))
    manifest={"generated_at":datetime.now(timezone.utc).isoformat(),"provenance":"Fresh executions of selected archived final-revision Python blocks. Plots restyled; solvers and numerical parameters retained. These are individual example runs, not aggregate benchmark scores.","results":items,"failures":failures}
    (out/"manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    docs=ROOT/"docs";docs.mkdir(exist_ok=True)
    lines=["# Gallery result provenance","",manifest['provenance'],"","Rebuild with `python tools/build_showcase.py`; use `--case cylinder` (or another case ID) for one result. Dependencies: NumPy, SciPy, Matplotlib, pandas, scikit-learn, Pillow.","","All source files are preserved. The builder extracts their fenced Python blocks and executes them in `.showcase_work/<case>/`. Python plotting uses a noninteractive backend. The cylinder adds read-only snapshot collection in the original time loop; the GIF shows those computed states at fixed color limits. No physical or numerical parameter was changed.","","The gallery reports solver diagnostics and in-sample regression fits honestly. It does not certify discretization accuracy; execution success and benchmark accuracy are distinct.",""]
    for item in items:
        lines += [f"## {item['title']}","",f"Source: `{item['source']}`",f"SHA-256: `{item['source_sha256']}`","",item['description'],"",*['- '+k+': '+v for k,v in item['metrics'].items()],"",item.get('note',''),"",f"Log: `{item['execution_log']}`",""]
    if failures:
        lines += ["## Failed rebuilds","",*['- '+x['id']+': '+x['error'].replace('\n',' ') for x in failures]]
    (docs/"results_provenance.md").write_text('\n'.join(lines),encoding="utf-8")
    print(f"Built {len(items)} result cards; failures: {len(failures)}")
    return not failures


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--case",choices=list(SPECS));parser.add_argument("--worker",choices=list(SPECS))
    args=parser.parse_args()
    if args.worker:worker(args.worker,ROOT/"assets/results")
    else:sys.exit(0 if build(args.case) else 1)
