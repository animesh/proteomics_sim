
from chromatography import peak

def acquire(peptides,time_axis,topn=10,dynamic_exclusion=5):
    selected=[]
    exclusion={}
    for t in time_axis:
        active=[]
        for p in peptides:
            inten=peak(p[4],p[5],t)
            if inten>1:
                active.append((p,inten))
        active.sort(key=lambda x:x[1],reverse=True)
        n=0
        for p,inten in active:
            if p[0] in exclusion and t<exclusion[p[0]]:
                continue
            exclusion[p[0]]=t+dynamic_exclusion
            selected.append((t,p,inten))
            n+=1
            if n>=topn:
                break
    return selected
