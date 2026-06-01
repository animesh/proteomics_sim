
import numpy as np

def peak(rt, abundance, t, width=8.0):
    sigma = width/2.355
    return abundance*np.exp(-((t-rt)**2)/(2*sigma**2))

def chromatogram_trace(rt, abundance, time_axis, width=8.0):
    return peak(rt, abundance, time_axis, width)
