"""Build the student-facing data bundle in the project repo (data/examples + data/starter)."""
import json, csv, os, gzip, numpy as np, pysam, shutil
REPO="/home/claude/chronomics-project-25"
EX=f"{REPO}/data/examples"; ST=f"{REPO}/data/starter"
# clean out the workbench's stock examples/starters, keep the dirs
if os.path.exists(EX): shutil.rmtree(EX)
if os.path.exists(ST): shutil.rmtree(ST)
os.makedirs(EX+"/genes"); os.makedirs(EX+"/reference"); os.makedirs(EX+"/variants")
os.makedirs(EX+"/conservation"); os.makedirs(EX+"/motif"); os.makedirs(ST)

proms=json.load(open("raw/promoters.json")); seqs=json.load(open("raw/promoter_seqs.json"))
proms=[r for r in proms if r['ens'] in seqs]
by=lambda r:(r['chrom'],r['pstart'])

# 1) gene lists
def write_genes(path, rows):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["gene","ensembl","chrom","strand","tss","promoter_start","promoter_end","panel"])
        for r in rows: w.writerow([r['sym'],r['ens'],"chr"+r['chrom'],"+" if r['strand']==1 else "-",r['tss'],r['pstart'],r['pend'],"channelopathy" if r['core'] else "epilepsy"])
write_genes(f"{EX}/genes/brain_channelopathy.tsv",[r for r in proms if r['core']])
write_genes(f"{EX}/genes/epilepsy_all.tsv",sorted(proms,key=lambda r:r['sym']))

# 2) promoters BED (0-based)
with open(f"{EX}/genes/promoters.bed","w") as f:
    for r in sorted(proms,key=by):
        f.write(f"chr{r['chrom']}\t{r['pstart']-1}\t{r['pend']}\t{r['sym']}\t0\t{'+' if r['strand']==1 else '-'}\n")

# 3) promoter reference FASTA (forward strand) + .fai ; header carries genomic coords & strand
with open(f"{EX}/reference/promoters.fa","w") as f:
    for r in sorted(proms,key=lambda r:r['sym']):
        s=seqs[r['ens']]
        f.write(f">{r['sym']} chr{r['chrom']}:{r['pstart']}-{r['pend']} strand={'+' if r['strand']==1 else '-'} tss={r['tss']}\n")
        for i in range(0,len(s),70): f.write(s[i:i+70]+"\n")
pysam.faidx(f"{EX}/reference/promoters.fa")

# 4) variants: sites-only VCF (INFO carries AF + per-superpop AF), bgzipped + tabixed
rows=list(csv.DictReader(open("out/promoter_snvs.csv")))
SP=["AFR","AMR","EAS","EUR","SAS"]
vpath=f"{EX}/variants/promoter_variants.vcf"
with open(vpath,"w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Global allele frequency (1000G 30x, 3202 samples)">\n')
    f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Promoter gene">\n')
    for p in SP:
        f.write(f'##INFO=<ID=AF_{p},Number=A,Type=Float,Description="Allele frequency in {p} (2504 unrelated)">\n')
    # contigs
    for c in sorted(set("chr"+r['chrom'] for r in rows)): f.write(f"##contig=<ID={c}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for r in sorted(rows,key=lambda r:("chr"+r['chrom'], int(r['pos']))):
        def g(k):
            v=r.get(k,""); return v if v not in ("","None") else "."
        info=f"GENE={r['gene']};AF={g('AF')}"+ "".join(f";AF_{p}={g('AF_'+p+'_unrel')}" for p in SP)
        f.write(f"chr{r['chrom']}\t{r['pos']}\t{r['vid']}\t{r['ref']}\t{r['alt']}\t.\t.\t{info}\n")
pysam.tabix_compress(vpath, vpath+".gz", force=True); os.remove(vpath)
pysam.tabix_index(vpath+".gz", preset="vcf", force=True)

# 5) conservation: per-position phyloP for all promoters (gzipped CSV) + example-locus bedGraph
if os.path.exists("raw/phylop.npz"):
    z=np.load("raw/phylop.npz",allow_pickle=True); z={k:z[k] for k in z.files}
    with gzip.open(f"{EX}/conservation/phylop_promoters.csv.gz","wt") as f:
        f.write("chrom,pos,phyloP470,phyloP100\n")
        for r in sorted(proms,key=by):
            a470=z.get(r['ens']+"|p470"); a100=z.get(r['ens']+"|p100")
            if a470 is None: continue
            for i in range(len(a470)):
                p4=a470[i]; p1=a100[i] if a100 is not None and i<len(a100) else float('nan')
                f.write(f"chr{r['chrom']},{r['pstart']+i},{'' if p4!=p4 else round(float(p4),3)},{'' if p1!=p1 else round(float(p1),3)}\n")

# 6) PWM (JASPAR ARNTL MA0603.2) as a simple PFM text
d=json.load(open("raw/jaspar_MA0603.2.json")); pfm=d['pfm']
with open(f"{EX}/motif/ARNTL_MA0603.2.jaspar","w") as f:
    f.write(">MA0603.2\tArntl\n")
    for b in "ACGT": f.write(f"{b} [ "+" ".join(str(int(x)) for x in pfm[b])+" ]\n")

print("bundle written to",EX)
for root,_,fs in os.walk(EX):
    for fn in sorted(fs): 
        p=os.path.join(root,fn); print(f"  {os.path.relpath(p,EX)}  {os.path.getsize(p)//1024} KB")
