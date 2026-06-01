
import random
AA='ACDEFGHIKLMNPQRSTVWY'
MASS={a:m for a,m in zip(AA,[71,103,115,129,147,57,137,113,128,113,131,114,97,128,156,87,101,99,186,163])}

def generate_peptides(n=1000,min_len=7,max_len=25,seed=42):
    random.seed(seed)
    peps=[]
    for i in range(n):
        seq=''.join(random.choice(AA) for _ in range(random.randint(min_len,max_len)))
        mass=sum(MASS[a] for a in seq)+18
        peps.append((i,seq,mass,mass+1.007))
    return peps
