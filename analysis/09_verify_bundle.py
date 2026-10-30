"""Reproduce the guide's checkpoints using ONLY the shipped bundle + the guide's own code."""
import re, math, gzip, csv, io, subprocess, sys
EX="/home/claude/chronomics-project-25/data/examples"
import pysam

# ---- Part 1: promoters.fa faidx-able; count & GC ----
fa=pysam.FastaFile(f"{EX}/reference/promoters.fa")
scn1a=fa.fetch("SCN1A")
print("Part1 SCN1A length:",len(scn1a))
def read_fasta(path):
    name,seq,out=None,[],[]
    for line in open(path):
        if line.startswith(">"):
            if name: out.append((name,"".join(seq)))
            name,seq=line[1:].strip(),[]
        else: seq.append(line.strip().upper())
    if name: out.append((name,"".join(seq)))
    return out
records=read_fasta(f"{EX}/reference/promoters.fa")
print("Part1 promoters:",len(records))

# ---- Part 2: canonical CACGTG ----
eboxes=[]
for header,seq in records:
    gene=header.split()[0]; region=header.split()[1]
    chrom=region.split(":")[0]; gstart=int(region.split(":")[1].split("-")[0])
    for m in re.finditer("(?=CACGTG)",seq):
        off=m.start(); eboxes.append((chrom,gstart+off-1,gstart+off+5,gene))
print("Part2 canonical E-boxes:",len(eboxes),"(expect 205)")

# PWM
counts={}
for line in open(f"{EX}/motif/ARNTL_MA0603.2.jaspar"):
    if line.startswith(">"): continue
    counts[line[0]]=[int(x) for x in line[line.index("[")+1:line.index("]")].split()]
L=len(counts["A"]); bg={"A":0.295,"C":0.205,"G":0.205,"T":0.295}
w={}
for i in range(L):
    tot=sum(counts[b][i] for b in "ACGT")+4*0.8
    for b in "ACGT": w[(b,i)]=math.log2(((counts[b][i]+0.8)/tot)/bg[b])
best=sum(max(w[(b,i)] for b in "ACGT") for i in range(L)); worst=sum(min(w[(b,i)] for b in "ACGT") for i in range(L))
def relscore(s):
    if any(c not in "ACGT" for c in s): return -1
    return (sum(w[(c,i)] for i,c in enumerate(s))-worst)/(best-worst)
comp=str.maketrans("ACGT","TGCA")
def best_both(s): return max(relscore(s),relscore(s.translate(comp)[::-1]))
THR=0.85; pwm=0
for header,seq in records:
    i=0
    while i<=len(seq)-L:
        if best_both(seq[i:i+L])>=THR: pwm+=1; i+=L
        else: i+=1
print("Part2 PWM sites >=0.85:",pwm,"(expect 1081)")

# enrichment
from collections import Counter
pool="".join(s for _,s in records); base=Counter(c for c in pool if c in "ACGT"); N=sum(base.values())
di=Counter(pool[i:i+2] for i in range(len(pool)-1) if pool[i] in "ACGT" and pool[i+1] in "ACGT")
def P(a): return base[a]/N
def T(a,b):
    s=sum(di[a+x] for x in "ACGT"); return di[a+b]/s if s else .25
p=P("C")*T("C","A")*T("A","C")*T("C","G")*T("G","T")*T("T","G")
positions=sum(max(0,len(s)-5) for _,s in records)
print(f"Part2 enrichment: obs {len(eboxes)} exp {p*positions:.0f} ratio {len(eboxes)/(p*positions):.2f} (expect ~1.66)")

# ---- Part 3: variants VCF ----
n=int(subprocess.run(["python3","-c",
  "import pysam;print(sum(1 for _ in pysam.VariantFile('%s/variants/promoter_variants.vcf.gz')))"%EX],
  capture_output=True,text=True).stdout)
print("Part3 SNVs in VCF:",n,"(expect 22836)")
vf=pysam.VariantFile(f"{EX}/variants/promoter_variants.vcf.gz")
rec=next(iter(vf)); print("Part3 INFO has AF & AF_EUR:", "AF" in rec.info and "AF_EUR" in rec.info)
# canonical overlap
ebox_pos=set()
for chrom,s,e,g in eboxes:
    for pp in range(s+1,e+1): ebox_pos.add((chrom,pp))
vf=pysam.VariantFile(f"{EX}/variants/promoter_variants.vcf.gz")
canon=sum(1 for r in vf if (r.chrom,r.pos) in ebox_pos)
print("Part3 canonical-E-box variants:",canon,"(expect 27)")

# ---- Part 4: phylop table + selection (ebox vs control CpG) ----
import numpy as np
ph={}
with gzip.open(f"{EX}/conservation/phylop_promoters.csv.gz","rt") as f:
    rd=csv.reader(f); next(rd)
    for chrom,pos,p4,p1 in rd:
        if p4!="": ph[(chrom,int(pos))]=float(p4)
ev=[]; cv=[]
for header,seq in records:
    chrom=header.split()[1].split(":")[0]; gstart=int(header.split()[1].split(":")[1].split("-")[0])
    for i in range(len(seq)-1):
        pos=gstart+i; key=(chrom,pos)
        if key not in ph: continue
        is_ebox=key in ebox_pos; is_cpg=seq[i]=="C" and seq[i+1]=="G"
        if is_ebox: ev.append(ph[key])
        elif is_cpg: cv.append(ph[key])
print(f"Part4 phyloP mean  E-box {np.mean(ev):.2f} (expect ~2.68)  control-CpG {np.mean(cv):.2f} (expect ~0.47)")
