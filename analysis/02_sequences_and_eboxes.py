import requests, json, time, re
UP, DOWN = 1000, 500            # promoter = TSS-1000 .. TSS+500, strand-aware
rows=json.load(open("raw/tss.json"))
hdr={"Content-Type":"application/json","Accept":"application/json"}

# genomic promoter window (forward-strand coords), 1-based inclusive
for r in rows:
    if r['strand']==1: gs,ge=r['tss']-UP, r['tss']+DOWN
    else:              gs,ge=r['tss']-DOWN, r['tss']+UP
    r['pstart'],r['pend']=gs,ge

# fetch forward-strand genomic sequence for each promoter, batched POST (<=50 regions)
def fetch_seqs(rows):
    seqs={}
    regs=[f"{r['chrom']}:{r['pstart']}..{r['pend']}:1" for r in rows]
    for i in range(0,len(regs),50):
        chunk=regs[i:i+50]
        for attempt in range(5):
            resp=requests.post("https://rest.ensembl.org/sequence/region/human",
                               headers=hdr,data=json.dumps({"regions":chunk}),timeout=120)
            if resp.status_code==200:
                for reg,res in zip(rows[i:i+50],resp.json()):
                    seqs[reg['ens']]=res['seq'].upper()
                break
            time.sleep(2*(attempt+1))
        else:
            print("FAIL batch",i,resp.status_code,resp.text[:150])
    return seqs

seqs=fetch_seqs(rows)
print("sequences fetched:",len(seqs))

# scan forward strand for canonical E-box CACGTG (reverse-complement palindrome: forward scan suffices)
EBOX="CACGTG"
ebox_rows=[]; per_gene={}
for r in rows:
    s=seqs.get(r['ens'])
    if not s or len(s)!=(r['pend']-r['pstart']+1): 
        per_gene[r['sym']]=0; continue
    n=0
    for m in re.finditer("(?=CACGTG)",s):
        p0=m.start()                     # 0-based offset in forward seq
        g_start=r['pstart']+p0           # genomic 1-based start of the 6-mer
        g_end=g_start+5
        # position of E-box centre relative to TSS, strand-aware
        centre=g_start+2.5
        rel=(centre-r['tss']) if r['strand']==1 else (r['tss']-centre)
        ebox_rows.append(dict(sym=r['sym'],ens=r['ens'],chrom=r['chrom'],strand=r['strand'],
                              tss=r['tss'],start=g_start,end=g_end,rel_tss=round(rel,1),core=r['core']))
        n+=1
    per_gene[r['sym']]=n

json.dump(rows,open("raw/promoters.json","w"))
json.dump(seqs,open("raw/promoter_seqs.json","w"))
json.dump(ebox_rows,open("out/eboxes.json","w"))

tot=len(ebox_rows); ncore=sum(1 for e in ebox_rows if e['core'])
print(f"E-boxes (CACGTG): {tot} across {len(rows)} promoters; {ncore} in the {sum(r['core'] for r in rows)} core genes")
# expected by chance: promoter length 1500, p(CACGTG)=per-position prob under base comp
import statistics
# observed rate per promoter
rates=[per_gene[r['sym']] for r in rows]
print("mean E-boxes/promoter: %.3f ; promoters with >=1: %d (%.0f%%)"%(
      statistics.mean(rates), sum(1 for x in rates if x>0), 100*sum(1 for x in rates if x>0)/len(rates)))
# core per-gene
print("\nCore genes E-box counts:")
for r in rows:
    if r['core']: print("  %-8s %d"%(r['sym'],per_gene[r['sym']]))
