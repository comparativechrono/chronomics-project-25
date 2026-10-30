import pysam, json, csv, time, sys, os, gc
base="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
proms=json.load(open("raw/promoters.json"))
SPOPS=["AFR","AMR","EAS","EUR","SAS"]
fields=["gene","ens","chrom","pos","vid","ref","alt","AC","AN","AF"]+\
       [f"AF_{p}_unrel" for p in SPOPS]+[f"AC_{p}_unrel" for p in SPOPS]+[f"AN_{p}_unrel" for p in SPOPS]
outpath="out/promoter_snvs.csv"; cntpath="out/promoter_snv_counts.json"
done=set(); counts={}
if os.path.exists(outpath):
    for r in csv.DictReader(open(outpath)): done.add(r['ens'])
    if os.path.exists(cntpath): counts=json.load(open(cntpath))
todo=[r for r in proms if r['ens'] not in done]
if not todo:
    print("ALL DONE"); sys.exit(0)
bychrom={}
for r in todo: bychrom.setdefault(r['chrom'],[]).append(r)
out=open(outpath,"a" if done else "w",newline=""); w=csv.DictWriter(out,fieldnames=fields)
if not done: w.writeheader()
print(f"resume: {len(done)} done, {len(todo)} to do",flush=True)
t0=time.time(); nsnv=0
def emit(r,recs):
    global nsnv
    c=0; seen=set()
    for rec in recs:
        if rec.info.get("SVTYPE") or rec.info.get("END") is not None: continue
        for ai,alt in enumerate(rec.alts or []):
            if len(rec.ref)==1 and len(alt)==1 and alt in "ACGT" and rec.ref in "ACGT":
                key=(rec.pos,rec.ref,alt)
                if key in seen: continue
                seen.add(key)
                def gi(tag,default=""):
                    v=rec.info.get(tag); return default if v is None else (v[ai] if isinstance(v,tuple) else v)
                row=dict(gene=r['sym'],ens=r['ens'],chrom=r['chrom'],pos=rec.pos,vid=rec.id or ".",
                         ref=rec.ref,alt=alt,AC=gi("AC"),AN=rec.info.get("AN",""),AF=gi("AF"))
                for p in SPOPS:
                    row[f"AF_{p}_unrel"]=gi(f"AF_{p}_unrel"); row[f"AC_{p}_unrel"]=gi(f"AC_{p}_unrel"); row[f"AN_{p}_unrel"]=rec.info.get(f"AN_{p}_unrel","")
                w.writerow(row); nsnv+=1; c+=1
    counts[r['ens']]=c
try:
    for chrom in sorted(bychrom,key=lambda c:(len(c),c)):
        suffix=".v2.vcf.gz" if chrom in ("X","Y") else ".vcf.gz"
        f=f"{base}/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel{suffix}"
        vf=pysam.VariantFile(f)
        for r in bychrom[chrom]:
            ok=False
            for a in range(4):
                try:
                    recs=list(vf.fetch("chr"+chrom, r['pstart']-1, r['pend'])); ok=True; break
                except Exception as e:
                    try: vf.close()
                    except Exception: pass
                    gc.collect(); time.sleep(2*(a+1))
                    try: vf=pysam.VariantFile(f)
                    except Exception: time.sleep(3)
            if not ok:
                raise RuntimeError("region failed "+r['sym'])
            emit(r,recs)
        try: vf.close()
        except Exception: pass
        out.flush(); json.dump(counts,open(cntpath,"w"))
        sys.stderr.write(f"chr{chrom} ok cum_new_snv={nsnv} {time.time()-t0:.0f}s\n"); sys.stderr.flush()
except Exception as e:
    out.flush(); json.dump(counts,open(cntpath,"w"))
    sys.stderr.write(f"STOP {type(e).__name__}: {e}; saved, will resume\n")
    sys.exit(3)
out.close(); json.dump(counts,open(cntpath,"w")); print("chunk done")
