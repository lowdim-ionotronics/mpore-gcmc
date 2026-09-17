import numpy as np
import warnings
warnings.filterwarnings("ignore")

def calc_cap(rho_factor, count_data, ion_charges, cap_coeff):
    correl = np.mean(count_data[:,:,np.newaxis]*count_data[:,np.newaxis,:], axis=0)
    q_correl = np.outer(ion_charges, ion_charges)
    nq_mean = np.sum(np.mean(count_data, axis=0)*ion_charges)
    nq2_mean = np.sum(correl*q_correl)
    return cap_coeff*100.0*(nq2_mean-nq_mean**2)*rho_factor