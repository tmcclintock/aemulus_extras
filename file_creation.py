"""This is a collection of routines meant to be run once.

The functions here create various analytic quantities for each cosmology
in the Aemulus simulations.
"""
import numpy as np
from classy import Class
import cluster_toolkit.xi as ctxi

#Aemulus scale factors for test and training sims
sfs = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.])
zs = 1./sfs -1

#Define some constants
#k values in 1/Mpc
kmin = 1e-5
kmax = 10
k = np.logspace(np.log10(kmin), np.log10(kmax), num=1000) #Mpc^-1
#3d r values in Mpc/h
rmin = 0.1
rmax = 100.
r = np.logspace(np.log10(rmin), np.log10(rmax), num=1000) #Mpc/h

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

def calc_power_spec(cosmos, zs, lin=True):
    Nc = len(cosmos)
    Nz = len(zs)
    ps = np.zeros((Nc, Nz, len(k)))
    for i in range(Nc):
        obh2, och2, w, ns, ln10As, H0, Neff, s8 = cosmos[i]
        h = H0/100.
        Omega_b = obh2/h**2
        Omega_c = och2/h**2
        Omega_m = Omega_b+Omega_c
        params = {'output': 'mPk', 'h': h, 'ln10^{10}A_s': ln10As, 'n_s': ns, 'w0_fld': w, 'wa_fld': 0.0, 'Omega_b': Omega_b, 'Omega_cdm': Omega_c, 'Omega_Lambda': 1.- Omega_m, 'N_eff': Neff, 'P_k_max_1/Mpc':10., 'z_max_pk':5. }
        if not lin:
            params['non linear'] = 'halofit' 
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()
        for j in range(len(zs)):
            if lin:
                ps[i,j] = np.array([cosmo.pk_lin(ki, zs[j]) for ki in k])
            else:
                ps[i,j] = np.array([cosmo.pk(ki, zs[j]) for ki in k])
            continue
        print "Finished box%d"%i
    return ps

def make_linear_power_spectra():
    """Create P_lin(k) for all cosmologies. This uses CLASS.
    """
    np.savetxt("plin/k.txt", k, header="k Mpc^-1", fmt="%e")
    
    #First do the training cosmos
    training_cos = get_training_cosmos()
    plins = calc_power_spec(training_cos, zs, lin=True)
    np.save("plin/plins_training_all_mpc3", plins)

    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    plins = calc_power_spec(testing_cos, zs, lin=True)
    np.save("plin/plins_testing_all_mpc3", plins)

    #Third do the highres cosmos TODO
    print "Finished linear power spectra"
    return

def make_nonlinear_power_spectra():
    """Create P_nonlin(k) for all cosmologies. This uses halofit-takahashi as implemented in CLASS.
    """
    np.savetxt("pnl/k.txt", k, header="k Mpc^-1", fmt="%e")
    
    #First do the training cosmos
    training_cos = get_training_cosmos()
    pnls = calc_power_spec(training_cos, zs, lin=False)
    np.save("pnl/pnls_training_all_mpc3", pnls)

    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    pnls = calc_power_spec(testing_cos, zs, lin=False)
    np.save("pnl/pnls_testing_all_mpc3", pnls)

    #Third do the highres cosmos TODO
    print "Finished nonlinear power spectra"
    return

def calc_ximm(cosmos, zs, k, p):
    Nc = len(cosmos)
    Nz = len(zs)
    xis = np.zeros((Nc, Nz, len(r)))
    for i in range(Nc):
        obh2, och2, w, ns, ln10As, H0, Neff, s8 = cosmos[i]
        h = H0/100.
        kh = k/h #now h/Mpc
        for j in range(Nz):
            ph3 = p[i,j]*h**3 #now (Mpc/h)^3
            xis[i,j] = ctxi.xi_mm_at_R(r, kh, ph3)
            continue
        print "Finished ximm box %d"%i
    return xis
            
def make_linear_correlation_function():
    np.savetxt("xilin/r.txt",r,header='r [Mpc/h; comoving]')
    k = np.loadtxt("plin/k.txt") #1/Mpc

    #First do the training cosmos
    training_cos = get_training_cosmos()
    plins = np.load("plin/plins_training_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(training_cos, zs, k, plins)
    np.save("xilin/xilins_training_all", xis)
    
    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    plins = np.load("plin/plins_testing_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(testing_cos, zs, k, plins)
    np.save("xilin/xilins_testing_all", xis)

    #Third do highres cosmos TODO
    print "Finished all xi_lin calculations"
    return

def make_nonlinear_correlation_function():
    np.savetxt("xinl/r.txt",r,header='r [Mpc/h; comoving]')
    k = np.loadtxt("pnl/k.txt") #1/Mpc

    #First do the training cosmos
    training_cos = get_training_cosmos()
    pnls = np.load("pnl/pnls_training_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(training_cos, zs, k, pnls)
    np.save("xinl/xinls_training_all", xis)
    
    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    pnls = np.load("pnl/pnls_testing_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(testing_cos, zs, k, pnls)
    np.save("xinl/xinls_testing_all", xis)

    #Third do highres cosmos TODO
    print "Finished all xi_nl calculations"
    return


if __name__ == "__main__":
    #make_linear_power_spectra()
    #make_nonlinear_power_spectra()
    make_linear_correlation_function()
    make_nonlinear_correlation_function()
