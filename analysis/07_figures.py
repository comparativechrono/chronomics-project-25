import json, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
plt.rcParams.update({"figure.dpi":130,"font.size":10,"axes.splines.top" if False else "axes.grid":False})
NAVY="#20384E"; TEAL="#2BA69A"; GOLD="#E0A42B"; GREY="#8aa0ad"; RED="#c0506a"
proms={r['ens']:r for r in json.load(open("raw/promoters.json"))}
can=json.load(open("out/eboxes_canonical.json"))
sel=json.load(open("out/selection_stats.json"))
summ=json.load(open("out/summary.json"))
import os; os.makedirs("fig",exist_ok=True)

# F1: observed vs Markov-expected canonical E-boxes (enrichment)
mk=json.load(open("out/motif_counts.json")) if os.path.exists("out/motif_counts.json") else None
if mk:
    fig,ax=plt.subplots(figsize=(4.2,3.2))
    labels=list(mk.keys()); obs=[mk[k]['obs'] for k in labels]; exp=[mk[k]['exp'] for k in labels]
    x=np.arange(len(labels)); w=0.38
    ax.bar(x-w/2,obs,w,label="observed",color=TEAL)
    ax.bar(x+w/2,exp,w,label="expected (CpG-aware)",color=GREY)
    ax.set_xticks(x); ax.set_xticklabels(labels,rotation=0); ax.set_ylabel("sites in promoters")
    ax.legend(frameon=False); ax.set_title("E-box motifs vs CpG-aware expectation")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig("fig/F1_enrichment.png"); plt.close(fig)

# F2: positional density of canonical E-boxes around TSS
rel=[x.get('rel_tss') for x in json.load(open("out/eboxes.json"))] if os.path.exists("out/eboxes.json") else []
rel=[r for r in rel if r is not None]
if rel:
    fig,ax=plt.subplots(figsize=(4.6,3.0))
    ax.hist(rel,bins=np.arange(-1000,501,100),color=NAVY)
    ax.axvline(0,color=GOLD,lw=2); ax.set_xlabel("position relative to TSS (bp)"); ax.set_ylabel("canonical E-boxes")
    ax.set_title("Where E-boxes sit in the promoter")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig("fig/F2_position.png"); plt.close(fig)

# F3: phyloP by class
fig,ax=plt.subplots(figsize=(4.4,3.2))
order=["ebox_can","cpg_ctrl","promoter_all"]; labs=["E-box (CACGTG)","other CpG","all promoter"]
data=[sel["phylop"][c]["p470"] for c in order]
means=[d["mean"] for d in data]
ax.bar(range(3),means,color=[TEAL,GOLD,GREY])
ax.set_xticks(range(3)); ax.set_xticklabels(labs,rotation=15,ha="right")
ax.set_ylabel("mean phyloP (470-way)"); ax.set_title("Conservation at E-boxes vs controls")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("fig/F3_phylop.png"); plt.close(fig)

# F4: polymorphic rate per bp by class
fig,ax=plt.subplots(figsize=(4.4,3.2))
pr=[sel["polymorphic_rate"][c]["rate"]*1000 for c in order]
ax.bar(range(3),pr,color=[TEAL,GOLD,GREY])
ax.set_xticks(range(3)); ax.set_xticklabels(labs,rotation=15,ha="right")
ax.set_ylabel("segregating SNVs per 1000 bp"); ax.set_title("Human variation at E-boxes vs controls")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("fig/F4_variation.png"); plt.close(fig)
print("figures written to fig/")
