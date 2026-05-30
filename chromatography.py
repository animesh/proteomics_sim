import numpy as np
def peak(rt, abundance, t, width=8):
    sigma=width/2.355
    return abundance*np.exp(-((t-rt)**2)/(2*sigma**2))
