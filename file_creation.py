"""This is a collection of routines meant to be run once.

The functions here create various analytic quantities for each cosmology
in the Aemulus simulations.
"""
import numpy as np
from classy import Class

#Aemulus scale factors for test and training sims
sfs = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.])
zs = 1./sfs -1
Nz = len(zs)

#Define some constants
#k values in 1/Mpc
kmin = 1e-5
kmax = 10

def get_training_cosmos():
    """Get the training simulation cosmologies.

    Returns:
        40x8 numpy array with columns: ombh2, omch2, w0, ns, ln10As, H0, Neff, s8
    """
    return np.loadtxt("training_cosmologies.txt")

def get_testing_cosmos():
    """Get the testing simulation cosmologies.

    Returns:
        40x8 numpy array with columns: ombh2, omch2, w0, ns, ln10As, H0, Neff, s8
    """
    return np.loadtxt("testing_cosmologies.txt")

def make_linear_power_spectra():
    """Create P_lin(k) for all cosmologies.
    """
    #Define save paths

    #Define the wavenumbers
    k = np.logspace(np.log10(kmin), np.log10(kmax), num=1000) #Mpc^-1
    np.savetxt("plin/k.txt", k, header="k Mpc^-1", fmt="%e")
    
    #First do the training cosmos
    training_cos = get_training_cosmos()
    Nc = len(training_cos)
    plins = np.zeros((Nc, Nz, len(k)))
    for i in range(Nc):
        obh2, och2, w, ns, ln10As, H0, Neff, s8 = training_cos[i]
        h = H0/100.
        Omega_b = obh2/h**2
        Omega_c = och2/h**2
        Omega_m = Omega_b+Omega_c
        params = {'output': 'mPk', 'h': h, 'ln10^{10}A_s': ln10As, 'n_s': ns, 'w0_fld': w, 'wa_fld': 0.0, 'Omega_b': Omega_b, 'Omega_cdm': Omega_c, 'Omega_Lambda': 1.- Omega_m, 'N_eff': Neff, 'P_k_max_1/Mpc':10., 'z_max_pk':5. }
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()
        for j in range(len(zs)):
            plins[i,j] = np.array([cosmo.pk_lin(ki, zs[j]) for ki in k])
        print "Finished training box%d"%i
    np.save("plin/plins_training_all_mpc3", plins)

    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    Nc = len(testing_cos)
    plins = np.zeros((Nc, Nz, len(k)))
    for i in range(len(testing_cos)):
        obh2, och2, w, ns, ln10As, H0, Neff, s8 = testing_cos[i]
        h = H0/100.
        Omega_b = obh2/h**2
        Omega_c = och2/h**2
        Omega_m = Omega_b+Omega_c
        params = {'output': 'mPk', 'h': h, 'ln10^{10}A_s': ln10As, 'n_s': ns, 'w0_fld': w, 'wa_fld': 0.0, 'Omega_b': Omega_b, 'Omega_cdm': Omega_c, 'Omega_Lambda': 1.- Omega_m, 'N_eff': Neff, 'P_k_max_1/Mpc':10., 'z_max_pk':5. }
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()
        for j in range(len(zs)):
            p = np.array([cosmo.pk_lin(ki, zs[j]) for ki in k])
        print "Finished testing box%d"%i
    np.save("plin/plins_testing_all_mpc3", plins)

    #Third do the highres cosmos TODO
    return

if __name__ == "__main__":
    make_linear_power_spectra()
