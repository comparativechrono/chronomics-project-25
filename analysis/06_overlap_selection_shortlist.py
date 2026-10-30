import json, csv, math, numpy as np, os
from collections import defaultdict
import scipy.stats as ss

proms={r['ens']:r for r in json.load(open("raw/promoters.json"))}
seqs=json.load(open("raw/promoter_seqs.json"))
can=json.load(open("out/eboxes_canonical.json"))
pwm_sites=json.load(open("out/eboxes_pwm.json"))
M=json.load(open("out/pwm_model.json")); L=M['L']; THR=M['THR']; maxs=M['maxs']; mins=M['mins']; bg=M['bg']
pwm={}
for k,v in M['pwm'].items():
    b,i=k.split(","); pwm[(b,int(i))]=v
comp=str.maketrans("ACGT","TGCA")
def rel(sub):
    sc=0.0
    for i,b in enumerate(sub):
        if b not in "ACGT": return -1.0
        sc+=pwm[(b,i)]
    return (sc-mins)/(maxs-mins)
def best_rel(sub): return max(rel(sub), rel(sub.translate(comp)[::-1]))

phz=np.load("raw/phylop.npz",allow_pickle=True); phz={k:phz[k] for k in phz.files}
SPOPS=["AFR","AMR","EAS","EUR","SAS"]

# genes we have everything for
genes=[e for e in proms if e in seqs and (e+"|p470") in phz]
core=set(e for e in genes if proms[e]['core'])
print("genes analysed:",len(genes),"core:",len(core))

# ---- per-position class arrays, pooled selection stats ----
can_by_ens=defaultdict(list)
for x in can: can_by_ens[x['ens']].append((x['start'],x['end']))
pwm_by_ens=defaultdict(list)
for x in pwm_sites: pwm_by_ens[x['ens']].append((x['start'],x['end']))

# variant AF per position (global AF, and rare flag) from promoter_snvs.csv
var_af=defaultdict(dict)   # ens -> {pos: AF}
rows=list(csv.DictReader(open("out/promoter_snvs.csv")))
def f(x):
    try: return float(x)
    except: return float('nan')
for r in rows:
    if r['ens'] in proms:
        af=f(r['AF']); 
        # keep max AF if multiallelic duplicate pos
        var_af[r['ens']][int(r['pos'])]=max(af, var_af[r['ens']].get(int(r['pos']),0))

# accumulate position-class vectors
cls_phylop={c:{"p470":[],"p100":[]} for c in ["ebox_can","cpg_ctrl","promoter_all"]}
cls_poly={c:[0,0] for c in ["ebox_can","cpg_ctrl","promoter_all"]}  # [n_poly, n_total]
cls_af={c:[] for c in ["ebox_can","cpg_ctrl","promoter_all"]}
for e in genes:
    r=proms[e]; s=seqs[e]; ps=r['pstart']; n=len(s)
    p470=phz[e+"|p470"]; p100=phz[e+"|p100"]
    ecan=np.zeros(n,bool)
    for a,b in can_by_ens.get(e,[]):
        i0=a-ps; ecan[max(0,i0):i0+6]=True
    cpg=np.zeros(n,bool)
    for i in range(n-1):
        if s[i]=='C' and s[i+1]=='G': cpg[i]=True   # tag the C of each CpG (matches the guide's per-position code)
    vmap=var_af.get(e,{})
    for i in range(min(n,len(p470),len(p100))):
        pos=ps+i
        af=vmap.get(pos,None)
        poly=1 if (af is not None and af>0) else 0
        for c,flag in (("promoter_all",True),("ebox_can",ecan[i]),("cpg_ctrl",(cpg[i] and not ecan[i]))):
            if not flag: continue
            if not math.isnan(p470[i]): cls_phylop[c]["p470"].append(float(p470[i]))
            if not math.isnan(p100[i]): cls_phylop[c]["p100"].append(float(p100[i]))
            cls_poly[c][0]+=poly; cls_poly[c][1]+=1
            if af is not None and af>0: cls_af[c].append(af)

def summ(a): 
    a=np.array(a,float); return dict(n=len(a),mean=float(np.mean(a)),median=float(np.median(a)))
sel={"phylop":{},"polymorphic_rate":{},"maf_spectrum":{}}
for c in cls_phylop:
    sel["phylop"][c]={k:summ(v) for k,v in cls_phylop[c].items()}
    npoly,ntot=cls_poly[c]; sel["polymorphic_rate"][c]=dict(poly=npoly,total=ntot,rate=npoly/ntot if ntot else 0)
    af=np.array(cls_af[c],float)
    sel["maf_spectrum"][c]=dict(n=len(af),mean_af=float(np.mean(af)) if len(af) else 0,
                                frac_rare_lt0p5pct=float(np.mean(af<0.005)) if len(af) else 0,
                                frac_common_gt5pct=float(np.mean(af>0.05)) if len(af) else 0)
# tests: ebox_can vs cpg_ctrl
for trk in ("p470","p100"):
    u=ss.mannwhitneyu(cls_phylop["ebox_can"][trk], cls_phylop["cpg_ctrl"][trk], alternative="greater")
    sel["phylop"].setdefault("tests",{})[trk+"_ebox_gt_cpgctrl"]=dict(U=float(u.statistic),p=float(u.pvalue))
# polymorphic rate ebox vs cpg_ctrl (Fisher)
a=cls_poly["ebox_can"]; b=cls_poly["cpg_ctrl"]
odds,pf=ss.fisher_exact([[a[0],a[1]-a[0]],[b[0],b[1]-b[0]]])
sel["polymorphic_rate"]["test_ebox_vs_cpgctrl"]=dict(odds=float(odds),p=float(pf))
json.dump(sel,open("out/selection_stats.json","w"),indent=2)
print("SELECTION (phyloP470 mean): ebox=%.3f cpgctrl=%.3f all=%.3f  MWU p=%.2e"%(
    sel["phylop"]["ebox_can"]["p470"]["mean"],sel["phylop"]["cpg_ctrl"]["p470"]["mean"],
    sel["phylop"]["promoter_all"]["p470"]["mean"],sel["phylop"]["tests"]["p470_ebox_gt_cpgctrl"]["p"]))
print("polymorphic/bp: ebox=%.4f cpgctrl=%.4f all=%.4f  Fisher p=%.2e OR=%.2f"%(
    sel["polymorphic_rate"]["ebox_can"]["rate"],sel["polymorphic_rate"]["cpg_ctrl"]["rate"],
    sel["polymorphic_rate"]["promoter_all"]["rate"],pf,odds))

# ---- variant -> motif catalogue (canonical membership + PWM delta) ----
cat=[]
refmismatch=0
for r in rows:
    e=r['ens']
    if e not in seqs: continue
    pos=int(r['pos']); s=seqs[e]; ps=proms[e]['pstart']; off=pos-ps
    if off<0 or off>=len(s): continue
    if s[off]!=r['ref']: refmismatch+=1; continue
    alt=r['alt']
    # canonical membership
    in_can=any(a<=pos<=b for a,b in can_by_ens.get(e,[]))
    # PWM: best window containing off
    best=None
    for w in range(max(0,off-L+1), min(off, len(s)-L)+1):
        refw=s[w:w+L]
        altw=refw[:off-w]+alt+refw[off-w+1:]
        rr=best_rel(refw); ra=best_rel(altw)
        key=max(rr,ra)
        if best is None or key>best[0]: best=(key,rr,ra,w)
    if best is None: continue
    _,rr,ra,w=best; delta=ra-rr
    cls="none"
    if rr>=THR and ra<THR: cls="break"
    elif rr<THR and ra>=THR: cls="create"
    elif rr>=THR and ra>=THR: cls="weaken" if delta<-1e-9 else ("strengthen" if delta>1e-9 else "none")
    if in_can or cls!="none" or abs(delta)>=0.05:
        cat.append(dict(gene=r['gene'],ens=e,chrom=r['chrom'],pos=pos,vid=r['vid'],ref=r['ref'],alt=alt,
                        AF=r['AF'],**{f"AF_{p}_unrel":r[f"AF_{p}_unrel"] for p in SPOPS},
                        in_canonical_ebox=in_can,pwm_ref=round(rr,3),pwm_alt=round(ra,3),
                        pwm_delta=round(delta,3),pwm_class=cls,core=proms[e]['core']))
cat.sort(key=lambda x:(x['pwm_class'] in ("break","create"), abs(x['pwm_delta']), x['in_canonical_ebox']),reverse=True)
import csv as _csv
fields=list(cat[0].keys())
w=_csv.DictWriter(open("out/catalogue.csv","w",newline=""),fieldnames=fields); w.writeheader(); [w.writerow(x) for x in cat]
print("catalogue variants:",len(cat),"(ref-mismatch skipped:",refmismatch,")")
nb=sum(1 for x in cat if x['pwm_class']=='break'); ncr=sum(1 for x in cat if x['pwm_class']=='create')
ncanv=sum(1 for x in cat if x['in_canonical_ebox'])
print(f"  canonical-E-box variants: {ncanv};  PWM break: {nb};  PWM create: {ncr}")

# ---- shortlist: disrupting, conserved, in core/channelopathy, with frequency ----
def meanphylop(e,a,b):
    arr=phz.get(e+"|p470"); ps=proms[e]['pstart']
    if arr is None: return float('nan')
    seg=arr[a-ps:b-ps+1]; seg=seg[~np.isnan(seg)]
    return float(np.mean(seg)) if len(seg) else float('nan')
short=[]
for x in cat:
    if x['pwm_class'] not in ("break","weaken") and not x['in_canonical_ebox']: continue
    ph=meanphylop(x['ens'],x['pos'],x['pos'])
    try: af=float(x['AF'])
    except: af=float('nan')
    score=(2 if x['in_canonical_ebox'] else 0)+(1.5 if x['pwm_class']=='break' else 0)+abs(x['pwm_delta'])+(0 if math.isnan(ph) else max(0,ph))/3+(1 if x['core'] else 0)
    short.append({**x,"phylop470":round(ph,3) if not math.isnan(ph) else "","score":round(score,3)})
short.sort(key=lambda x:x['score'],reverse=True)
w=_csv.DictWriter(open("out/shortlist.csv","w",newline=""),fieldnames=list(short[0].keys())); w.writeheader(); [w.writerow(x) for x in short[:40]]
print("shortlist (top up to 40) written; total candidates:",len(short))

# ---- summary numbers for the guide ----
summary=dict(n_genes=len(genes),n_core=len(core),
             canonical_eboxes=len(can),canonical_core=sum(x['core'] for x in can),
             pwm_sites=len(pwm_sites),pwm_core=sum(x['core'] for x in pwm_sites),THR=THR,
             n_promoter_snv=len(rows),catalogue=len(cat),canon_variants=ncanv,pwm_break=nb,pwm_create=ncr,
             selection=sel)
json.dump(summary,open("out/summary.json","w"),indent=2)
print("top shortlist rows:")
for x in short[:8]:
    print("  %-8s %s:%d %s>%s %s d=%.2f AF=%s phyloP=%s core=%s"%(
        x['gene'],x['chrom'],x['pos'],x['ref'],x['alt'],x['pwm_class'],x['pwm_delta'],x['AF'],x['phylop470'],x['core']))
