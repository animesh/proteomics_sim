from simulator import Simulator
for w in [40,20,12]:
    r=Simulator().run(100,w)
    counts=[len(v) for v in r['dia'].values()]
    print({'window':w,'median_precursors':sorted(counts)[len(counts)//2],'max_precursors':max(counts)})
