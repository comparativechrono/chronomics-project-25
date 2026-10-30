import json, time, sys, os, numpy as np, pybigtools
proms=json.load(open("raw/promoters.json"))
S3="https://genome-browser.s3.us-east-1.amazonaws.com/goldenPath/hg38"
TR={"p470":f"{S3}/phyloP470way/hg38.phyloP470way.bw","p100":f"{S3}/phyloP100way/hg38.phyloP100way.bw"}
cache="raw/phylop.npz"
done={}
if os.path.exists(cache):
    z=np.load(cache,allow_pickle=True); done={k:z[k] for k in z.files}
bw={}
def openbw(k):
    for a in range(6):
        try: return pybigtools.open(TR[k])
        except Exception: time.sleep(3*(a+1))
    raise RuntimeError("open "+k)
for k in TR: bw[k]=openbw(k)
t0=time.time(); n=0
for r in proms:
    key=r['ens']
    if key+"|p470" in done and key+"|p100" in done: continue
    for k in TR:
        for a in range(6):
            try:
                v=bw[k].values("chr"+r['chrom'], r['pstart']-1, r['pend'], fillna=float('nan'))
                done[key+"|"+k]=np.array(v,dtype=np.float32); break
            except Exception:
                time.sleep(2*(a+1)); bw[k]=openbw(k)
        else:
            print("fail",r['sym'],k)
    n+=1
    if n%80==0:
        np.savez_compressed(cache,**done); sys.stderr.write(f"{n} promoters, {time.time()-t0:.0f}s\n"); sys.stderr.flush()
np.savez_compressed(cache,**done)
print("phyloP cached for",len(proms),"promoters; arrays:",len(done))
