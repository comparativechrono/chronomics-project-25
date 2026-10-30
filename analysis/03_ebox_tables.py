import json, math, re
seqs=json.load(open("raw/promoter_seqs.json"))
proms={r['ens']:r for r in json.load(open("raw/promoters.json"))}
comp=str.maketrans("ACGT","TGCA")
# background from promoters
from collections import Counter
pool="".join(seqs.values()); bc=Counter(c for c in pool if c in "ACGT"); tot=sum(bc.values())
bg={b:bc[b]/tot for b in "ACGT"}
# ARNTL PWM
d=json.load(open("raw/jaspar_MA0603.2.json")); pfm=d['pfm']; L=len(pfm['A']); pseudo=0.8
pwm={}
for i in range(L):
    col={b:pfm[b][i] for b in "ACGT"}; s=sum(col.values())+4*pseudo
    for b in "ACGT": pwm[(b,i)]=math.log2(((col[b]+pseudo)/s)/bg[b])
maxs=sum(max(pwm[(b,i)] for b in "ACGT") for i in range(L))
mins=sum(min(pwm[(b,i)] for b in "ACGT") for i in range(L))
def rel(sub):
    sc=0
    for i,b in enumerate(sub):
        if b not in "ACGT": return -1
        sc+=pwm[(b,i)]
    return (sc-mins)/(maxs-mins)
def best_rel(sub):  # max over both strands
    return max(rel(sub), rel(sub.translate(comp)[::-1]))
THR=0.85
can=open("out/eboxes_canonical.bed","w"); pw=open("out/eboxes_pwm.bed","w")
can_rows=[]; pwm_rows=[]
percan=Counter(); perpwm=Counter()
for e,r in proms.items():
    s=seqs.get(e)
    if not s: continue
    # canonical CACGTG (palindrome -> forward scan)
    for m in re.finditer("(?=CACGTG)",s):
        g0=r['pstart']+m.start()  # 1-based genomic start
        can.write(f"chr{r['chrom']}\t{g0-1}\t{g0+5}\t{r['sym']}|CACGTG\t0\t+\n")
        can_rows.append(dict(sym=r['sym'],ens=e,chrom=r['chrom'],start=g0,end=g0+5,strand=r['strand'],core=r['core']))
        percan[r['sym']]+=1
    # PWM sites (merge overlapping windows -> keep local-best)
    hits=[]
    for i in range(len(s)-L+1):
        sc=best_rel(s[i:i+L])
        if sc>=THR: hits.append((i,sc))
    # non-overlapping greedy by score
    hits.sort(key=lambda x:-x[1]); taken=[]
    for i,sc in hits:
        if all(abs(i-j)>=L for j,_ in taken): taken.append((i,sc))
    for i,sc in sorted(taken):
        g0=r['pstart']+i
        pw.write(f"chr{r['chrom']}\t{g0-1}\t{g0+L-1}\t{r['sym']}|pwm\t{sc:.3f}\t.\n")
        pwm_rows.append(dict(sym=r['sym'],ens=e,chrom=r['chrom'],start=g0,end=g0+L-1,relscore=round(sc,3),core=r['core']))
        perpwm[r['sym']]+=1
can.close(); pw.close()
json.dump(can_rows,open("out/eboxes_canonical.json","w"))
json.dump(pwm_rows,open("out/eboxes_pwm.json","w"))
json.dump({"THR":THR,"L":L,"maxs":maxs,"mins":mins,"pwm":{f"{b},{i}":pwm[(b,i)] for b in "ACGT" for i in range(L)},"bg":bg},open("out/pwm_model.json","w"))
print("canonical E-boxes:",len(can_rows)," PWM sites(>=%.2f):"%THR,len(pwm_rows))
print("core canonical:",sum(x['core'] for x in can_rows)," core PWM:",sum(x['core'] for x in pwm_rows))
