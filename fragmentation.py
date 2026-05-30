
import numpy as np

def fragment_efficiency(ce, opt_ce, sigma=4):
    return np.exp(-((ce - opt_ce) ** 2) / (2 * sigma ** 2))
