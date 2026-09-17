import numpy as np
from scipy.integrate import cumtrapz
import warnings
warnings.filterwarnings("ignore")

def calc_energy(u, q):
    return q*u - cumtrapz(q, u, initial=0)
    

