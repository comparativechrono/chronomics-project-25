import json, re
from collections import Counter
seqs=json.load(open("raw/promoter_seqs.json"))
pool="".join(seqs.values()); bc=Counter(c for c in pool if c in "ACGT"); tot=sum(bc.values())
bg={b:bc[b]/tot for b in "ACGT"}
dc=Counter(pool[i:i+2] for i in range(len(pool)-1) if pool[i] in "ACGT" and pool[i+1] in "ACGT")
trans={a+b:(dc[a+b]/sum(dc[a+c] for c in "ACGT") if sum(dc[a+c] for c in "ACGT") else .25) for a in "ACGT" for b in "ACGT"}
def comp(w): return w[::-1].translate(str.maketrans("ACGT","TGCA"))
def exp_markov(word):
    p=bg[word[0]]
    for i in range(len(word)-1): p*=trans[word[i:i+2]]
    return p
npos=sum(max(0,len(s)-5) for s in seqs.values())
out={}
for w in ["CACGTG","CACGTT","CATGTG"]:
    pal = w==comp(w)
    obs=sum(len(re.findall("(?="+w+")",s))+(0 if pal else len(re.findall("(?="+comp(w)+")",s))) for s in seqs.values())
    exp=exp_markov(w)*npos*(1 if pal else 2)
    out[w]={"obs":obs,"exp":round(exp,1),"ratio":round(obs/exp,2)}
json.dump(out,open("out/motif_counts.json","w"),indent=2)
print(out)
