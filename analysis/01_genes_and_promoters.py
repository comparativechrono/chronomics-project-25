"""Phase A: gene lists (PanelApp) -> canonical TSS (Ensembl) -> promoter windows.
Writes raw/panel_90.json, raw/panel_402.json, raw/tss.json, raw/promoters.json.
Promoter = TSS-1000 .. TSS+500, strand-aware (GRCh38)."""
import requests, json, time
UP, DOWN = 1000, 500
hdr={"Content-Type":"application/json","Accept":"application/json"}

def panel_green(pid):
    d=requests.get(f"https://panelapp.genomicsengland.co.uk/api/v1/panels/{pid}/?format=json",timeout=120).json()
    out=[]
    for g in d['genes']:
        if g['confidence_level']!='3': continue
        gd=g['gene_data']; eg=gd.get('ensembl_genes',{}).get('GRch38',{})
        ens=eg[sorted(eg)[-1]]['ensembl_id'] if eg else None
        out.append((gd['gene_symbol'],ens))
    return d['name'],out

core_name,p90=panel_green(90); epi_name,p402=panel_green(402)
json.dump(p90,open("raw/panel_90.json","w")); json.dump(p402,open("raw/panel_402.json","w"))
core={e for s,e in p90 if e}; allg={}
for s,e in p402+p90:
    if e: allg[e]=s
print(f"channelopathy {len(p90)} green, epilepsy {len(p402)} green, {len(allg)} unique genes")

def post(ep,ids):
    out={}
    for i in range(0,len(ids),900):
        for a in range(5):
            r=requests.post(f"https://rest.ensembl.org{ep}",headers=hdr,data=json.dumps({"ids":ids[i:i+900]}),timeout=120)
            if r.status_code==200: out.update(r.json()); break
            time.sleep(2*(a+1))
    return out

genes=post("/lookup/id",list(allg))
canon={g:info["canonical_transcript"].split(".")[0] for g,info in genes.items() if info.get("canonical_transcript")}
tx=post("/lookup/id",list(canon.values()))
rows=[]
for g,sym in allg.items():
    ct=canon.get(g)
    if not ct or ct not in tx: continue
    t=tx[ct]; strand=t['strand']; tss=t['start'] if strand==1 else t['end']
    gs,ge=(tss-UP,tss+DOWN) if strand==1 else (tss-DOWN,tss+UP)
    rows.append(dict(sym=sym,ens=g,tx=ct,chrom=t['seq_region_name'],strand=strand,tss=tss,
                     pstart=gs,pend=ge,core=sym in core))
json.dump([dict(sym=r['sym'],ens=r['ens'],tx=r['tx'],chrom=r['chrom'],strand=r['strand'],tss=r['tss'],core=r['core']) for r in rows],open("raw/tss.json","w"))
json.dump(rows,open("raw/promoters.json","w"))
print("TSS + promoters for",len(rows),"genes;",sum(r['core'] for r in rows),"core")
