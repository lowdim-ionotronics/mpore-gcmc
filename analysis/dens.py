import numpy as np
import warnings
warnings.filterwarnings("ignore")

def calc_dens(rho_factor, count_data):
    return np.sum(np.mean(count_data, axis=0))*rho_factor # np.std(count_data)*rho_factor
