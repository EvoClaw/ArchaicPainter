"""Generate all main paper figures from v5 JSON result files."""
import json, warnings, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy import stats as sp_stats
warnings.filterwarnings("ignore")

EXEC = Path("/home/yanlin/livestock/docs/05_execution")
FIGS = Path("/home/yanlin/livestock/paper/figures")
C_AP="#E67E22"; C_POISSON="#3498DB"; C_NOMERGE="#27AE60"
C_POSONLY="#E74C3C"; C_CHB="#2980B9"; C_IBDMIX="#8E44AD"
plt.rcParams.update({
    "font.size":11,"axes.labelsize":12,"axes.titlesize":12,
    "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":9,
    "pdf.fonttype":42,"ps.fonttype":42,
    "axes.spines.top":False,"axes.spines.right":False,
})

def save_fig(fig, name):
    for ext in ("pdf","png"):
        fig.savefig(str(FIGS/f"{name}.{ext}"),bbox_inches="tight",dpi=200)
    plt.close(fig)
    print(f"  -> {name} saved")

def make_fig2():
    bench=json.load(open(EXEC/"phase4b_bench_v4.json"))
    raw=bench["raw"]; st=bench["stats"]
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,5))
    cfgs=[("archaicpainter","AP standard",C_AP),("ap_no_merge","AP no-merge",C_NOMERGE),
          ("ap_positive_only","AP pos-only",C_POSONLY),("poisson_density","Poisson",C_POISSON)]
    f1d=[[r["f1_all"] for r in raw[n] if "f1_all" in r] for n,_,_ in cfgs]
    bp=ax1.boxplot(f1d,patch_artist=True,widths=0.5,
                   medianprops=dict(color="black",lw=2),whiskerprops=dict(lw=1.2),
                   capprops=dict(lw=1.2),flierprops=dict(marker="o",ms=3,alpha=0.4))
    for patch,(_,_,c) in zip(bp["boxes"],cfgs): patch.set_facecolor(c); patch.set_alpha(0.75)
    ax1.set_xticks(range(1,5)); ax1.set_xticklabels([l for _,l,_ in cfgs],fontsize=9)
    ax1.set_ylabel("Segment-level F1"); ax1.set_ylim(-0.02,0.78)
    ax1.set_title("(a) F1 score distributions (50 replicates)")
    for i,(n,_,c) in enumerate(cfgs):
        m=st[n]["f1_mean"]
        ax1.text(i+1,m+0.01,f"{m:.3f}",ha="center",fontsize=8,color=c,fontweight="bold")
    fold=st["archaicpainter"]["f1_mean"]/st["poisson_density"]["f1_mean"]
    ax1.annotate("",xy=(4,0.64),xytext=(1,0.64),
                 arrowprops=dict(arrowstyle="<->",color="black",lw=1.2))
    ax1.text(2.5,0.655,f"{fold:.0f}x (p<10^-15)",ha="center",fontsize=8.5,fontweight="bold")
    auprcs=[st[n]["auprc_mean"] for n,_,_ in cfgs]
    astds=[st[n]["auprc_std"] for n,_,_ in cfgs]
    bars=ax2.bar(range(4),auprcs,color=[c for _,_,c in cfgs],alpha=0.8,width=0.55)
    ax2.errorbar(range(4),auprcs,yerr=astds,fmt="none",color="black",capsize=4,lw=1.2)
    ax2.set_xticks(range(4)); ax2.set_xticklabels([l for _,l,_ in cfgs],fontsize=9)
    ax2.set_ylabel("AUPRC"); ax2.set_ylim(0,0.80); ax2.set_title("(b) AUPRC comparison")
    for bar,val in zip(bars,auprcs):
        ax2.text(bar.get_x()+bar.get_width()/2.,val+0.01,f"{val:.3f}",ha="center",fontsize=9)
    fig.suptitle("Simulation benchmark (5 Mb, 20 haplotypes, 2% NEA admixture)",fontsize=11,y=1.01)
    fig.tight_layout(); save_fig(fig,"fig2_benchmark")

def make_fig3():
    real=json.load(open(EXEC/"realdata_summary_v5.json"))
    lines=[l for l in open(EXEC/"ibdmix/ibdmix_ceu20_chr21_raw.txt")
           if not l.startswith("#") and l.strip()]
    ibdmix=[(int(r[2]),int(r[3])) for r in [l.split() for l in lines[1:]]]
    ibdmix_lens=np.array([e-s for s,e in ibdmix])/1e3
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.5))
    ceu_afs=np.array(real["CEU"]["per_haplotype_ancestry"])*100
    chb_afs=np.array(real["CHB"]["per_haplotype_ancestry"])*100
    rng=np.random.default_rng(42)
    jc=rng.uniform(-0.12,0.12,len(ceu_afs)); jh=rng.uniform(-0.12,0.12,len(chb_afs))
    ax1.scatter(1+jc,ceu_afs,color=C_AP,alpha=0.55,s=20,zorder=3)
    ax1.scatter(2+jh,chb_afs,color=C_CHB,alpha=0.55,s=20,zorder=3)
    for x,afs,color in [(1,ceu_afs,C_AP),(2,chb_afs,C_CHB)]:
        m,s=np.mean(afs),np.std(afs)
        ax1.plot([x-0.18,x+0.18],[m,m],color=color,lw=2.5,zorder=5)
        ax1.errorbar(x,m,yerr=s,fmt="none",color=color,capsize=6,lw=2,zorder=5)
    ax1.set_xticks([1,2]); ax1.set_xticklabels(["CEU\n(n=40 hap)","CHB\n(n=40 hap)"])
    ax1.set_ylabel("Per-haplotype Neanderthal ancestry (%)"); ax1.set_xlim(0.5,2.5)
    top=max(max(ceu_afs),max(chb_afs)); ax1.set_ylim(0,top*1.30)
    ax1.set_title("(a) Per-haplotype NEA ancestry\nchr21:10-46.7 Mb (GRCh38)")
    for x,afs,color in [(1,ceu_afs,C_AP),(2,chb_afs,C_CHB)]:
        m,s=np.mean(afs),np.std(afs)
        ax1.text(x,top*1.19,f"{m:.2f}%+/-{s:.2f}%",ha="center",fontsize=8.5,
                 color=color,fontweight="bold")
    ax1.axhline(2.0,color="grey",ls="--",lw=1.0,alpha=0.6)
    ax1.text(2.38,2.06,"~2% genome-wide",fontsize=7.5,color="grey",ha="right")
    _,mw_p=sp_stats.mannwhitneyu(ceu_afs,chb_afs,alternative="two-sided")
    ax1.text(1.5,top*1.10,f"Mann-Whitney p={mw_p:.2f}",ha="center",fontsize=9)
    ap_lens=[]
    for l in open(EXEC/"archaicpainter_chr21_CEU_v5.bed"):
        if l.startswith("#"): continue
        p=l.strip().split("\t")
        if len(p)<3: continue
        s,e=int(p[1]),int(p[2])
        if s<19_989_084 and e>14_074_685: ap_lens.append((e-s)/1e3)
    ap_lens=np.array(ap_lens)
    bins=np.logspace(0,4,40)
    ax2.hist(ibdmix_lens,bins=bins,color=C_IBDMIX,alpha=0.6,
             label=f"IBDmix (n={len(ibdmix_lens)})",density=True)
    ax2.hist(ap_lens,bins=bins,color=C_AP,alpha=0.75,
             label=f"AP CEU (n={len(ap_lens)})",density=True)
    ax2.set_xscale("log"); ax2.set_xlabel("Segment length (kb)")
    ax2.set_ylabel("Density"); ax2.legend()
    ax2.set_title("(b) Segment-length distributions\nchr21:14-20 Mb (IBDmix LOD>=3)")
    ax2.text(0.97,0.97,
             f"IBDmix: mean={np.mean(ibdmix_lens):.1f} kb\nAP CEU: mean={np.mean(ap_lens):.1f} kb",
             transform=ax2.transAxes,fontsize=8.5,va="top",ha="right",
             bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.5))
    fig.tight_layout(); save_fig(fig,"fig3_realdata")

def make_fig4():
    sc=json.load(open(EXEC/"phase4b_scalability.json"))["results"]
    ns=np.array([r["n_haplotypes"] for r in sc],float)
    ts=np.array([r["per_hap_s"] for r in sc],float)
    fig,ax=plt.subplots(figsize=(6.5,5))
    ax.scatter(ns,ts,color=C_AP,s=60,zorder=5,label="Measured")
    slope_ph,intercept_ph,_,_,_=sp_stats.linregress(np.log(ns),np.log(ts))
    total_exp=slope_ph+1
    ns_fit=np.logspace(np.log10(ns[0]),np.log10(ns[-1]),100)
    ts_fit=np.exp(intercept_ph)*ns_fit**slope_ph
    ax.plot(ns_fit,ts_fit,"--",color="grey",lw=1.5,
            label=f"Log-log fit (total-exp={total_exp:.2f})")
    ax.scatter([1000],[ts[-1]],color="red",s=100,zorder=6,marker="*",label="n=1,000")
    ax.text(1100,ts[-1]*1.08,"0.114 s/hap",fontsize=9,color="red")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of query haplotypes (n)")
    ax.set_ylabel("Per-haplotype runtime (s)")
    ax.set_title("Near-linear scaling of ArchaicPainter\n(5 Mb region, single CPU core)")
    ax.legend(fontsize=9)
    ax.text(0.05,0.94,f"Empirical total-time exp: {total_exp:.2f}\n(theoretical O(n)=1.0)",
            transform=ax.transAxes,fontsize=9,va="top",
            bbox=dict(boxstyle="round",facecolor="lightyellow",alpha=0.8))
    fig.tight_layout(); save_fig(fig,"fig4_scalability")

def make_fig5():
    func=json.load(open(EXEC/"functional_annotation_chr21.json"))
    ai_genes=set(func["adaptive_introgression_hits"]["detected_genes"])
    ifn_genes={"IFNAR1","IFNAR2","IFNGR2","IL10RB"}
    gdet=func.get("gene_details",{})
    def load_bed(path):
        out=[]
        for l in open(path):
            if l.startswith("#"): continue
            p=l.strip().split("\t")
            if len(p)>=3: out.append((int(p[1]),int(p[2])))
        return out
    ceu_segs=load_bed(EXEC/"archaicpainter_chr21_CEU_v5.bed")
    chb_segs=load_bed(EXEC/"archaicpainter_chr21_CHB_v5.bed")
    RS=10_000_000; RE=46_709_983
    fig,axes=plt.subplots(3,1,figsize=(14,8.5),
                          gridspec_kw={"height_ratios":[3,3,2]})
    for ax,segs,color,lbl in [(axes[0],ceu_segs,C_AP,"CEU (n=40 haplotypes)"),
                               (axes[1],chb_segs,C_CHB,"CHB (n=40 haplotypes)")]:
        for s,e in segs:
            ax.barh(0,(e-s)/1e6,left=s/1e6,height=0.65,color=color,alpha=0.6,linewidth=0)
        ax.set_xlim(RS/1e6,RE/1e6); ax.set_ylim(-0.5,1.5); ax.set_yticks([])
        ax.text(0.005,0.88,lbl,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")
        ax.axvspan(33.2,33.5,color="red",alpha=0.12,zorder=0)
        ax.spines["left"].set_visible(False); ax.tick_params(axis="x",labelbottom=False)
    ax3=axes[2]; ax3.set_xlim(RS/1e6,RE/1e6); ax3.set_ylim(-0.3,2.5); ax3.set_yticks([])
    ax3.set_xlabel("Chr21 position (Mb, GRCh38)",fontsize=11)
    ax3.text(0.005,0.97,"Selected genes",transform=ax3.transAxes,
             fontsize=9,va="top",fontweight="bold")
    for gname in sorted(ai_genes|ifn_genes):
        if gname not in gdet: continue
        gs,ge=gdet[gname].get("start",0),gdet[gname].get("end",0)
        if ge<=gs: continue
        color="#C0392B" if gname in ifn_genes else "#6C3483"
        ax3.barh(0,(ge-gs)/1e6,left=gs/1e6,height=0.5,color=color,alpha=0.8)
        ax3.text((gs+ge)/2e6,0.7,gname,ha="center",va="bottom",
                 fontsize=7.5,rotation=55,color=color,fontweight="bold")
    ax3.axvspan(33.2,33.5,color="red",alpha=0.12,zorder=0)
    ax3.spines["left"].set_visible(False)
    legend_elems=[
        mpatches.Patch(color=C_AP,alpha=0.65,label="CEU introgressed"),
        mpatches.Patch(color=C_CHB,alpha=0.65,label="CHB introgressed"),
        mpatches.Patch(color="red",alpha=0.15,label="IFN receptor cluster (33.2-33.5 Mb)"),
        mpatches.Patch(color="#C0392B",label="IFN receptor genes"),
        mpatches.Patch(color="#6C3483",label="Published AI candidate genes"),
    ]
    fig.legend(handles=legend_elems,loc="upper right",ncol=1,
               fontsize=8.5,bbox_to_anchor=(0.995,0.995))
    fig.suptitle("Functional annotation: AP NEA regions, chr21 "
                 "(52 regions, 9.66 Mb; 75 named genes; 7/16 AI candidates)",fontsize=10)
    fig.subplots_adjust(hspace=0.05); save_fig(fig,"fig5_functional_annotation")

print("Generating figures from v5 results...")
make_fig2(); make_fig3(); make_fig4(); make_fig5()
print("All figures regenerated.")
