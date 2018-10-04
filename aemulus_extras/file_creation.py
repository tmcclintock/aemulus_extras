"""This is a collection of routines meant to be run once.

The functions here create various analytic quantities for each cosmology
in the Aemulus simulations.
"""
import numpy as np
from classy import Class
import cluster_toolkit.xi as ctxi
import cluster_toolkit.peak_height as ctph

#Aemulus scale factors for test and training sims
sfs = np.array([0.25, 0.333333, 0.5, 0.540541, 0.588235, 0.645161, 0.714286, 0.8, 0.909091, 1.])
zs = 1./sfs -1
#Aemulus scale factors for highres sims
sfs_hr = np.array([0.165913, 0.246189, 0.323199, 0.486149, 0.513345, 0.564651,
                   0.612689, 0.683156, 0.761728, 0.80434 , 0.849337, 0.89685 ,
                   0.934222, 1.      ])
zs_hr = 1./sfs_hr - 1.

#Define some constants
#k values in 1/Mpc
kmin = 1e-5
kmax = 10
k = np.logspace(np.log10(kmin), np.log10(kmax), num=1000) #Mpc^-1
#3d r values in Mpc/h
rmin = 0.1
rmax = 200.
r = np.logspace(np.log10(rmin), np.log10(rmax), num=1000) #Mpc/h
#Mass limits for peak height, bias and mass function
Mmin = 1e11
Mmax = 1e16
M = np.logspace(np.log10(Mmin), np.log10(Mmax), num=1000) #Msun/h

def get_training_cosmos():
    """Get the training simulation cosmologies.

    Returns:
        40x8 numpy array with columns: ombh2, omch2, w0, ns, ln10As, H0, Neff, s8
    """
    return np.loadtxt("cosmologies/training_cosmologies.txt")

def get_testing_cosmos():
    """Get the testing simulation cosmologies.

    Returns:
        40x8 numpy array with columns: ombh2, omch2, w0, ns, ln10As, H0, Neff, s8
    """
    return np.loadtxt("cosmologies/testing_cosmologies.txt")

def get_highres_cosmos():
    """Get the high resolution simulation cosmologies.

    Note: these don't have sigma8 as of 10/3/2018. Zeros are appended.

    Returns:
        25x8 numpy array with columns: ombh2, omch2, w0, ns, ln10As, H0, Neff
    """
    cos = np.loadtxt("cosmologies/highres_cosmologies.txt")
    Nc = len(cos)
    return np.vstack((cos.T, np.zeros(Nc))).T

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

    #Third do the highres cosmos
    hr_cos = get_highres_cosmos()
    plins = calc_power_spec(hr_cos, zs_hr, lin=True)
    np.save("plin/plins_highres_all_mpc3", plins)
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

    #Third do the highres cosmos
    hr_cos = get_highres_cosmos()
    pnls = calc_power_spec(hr_cos, zs_hr, lin=False)
    np.save("pnl/pnls_highres_all_mpc3", pnls)
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
            xis[i,j] = ctxi.xi_mm_at_R(r, kh, ph3)#, exact=True)
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

    #Third do highres cosmos
    hr_cos = get_highres_cosmos()
    plins = np.load("plin/plins_highres_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(hr_cos, zs_hr, k, plins)
    np.save("xilin/xilins_highres_all", xis)
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

    #Third do highres cosmos 
    hr_cos = get_highres_cosmos()
    pnls = np.load("pnl/pnls_highres_all_mpc3.npy") #[Mpc]^3
    xis = calc_ximm(hr_cos, zs_hr, k, pnls)
    np.save("xinl/xinls_highres_all", xis)

    print "Finished all xi_nl calculations"
    return

def calc_peak_height(cosmos, zs, k, p):
    Nc = len(cosmos)
    Nz = len(zs)
    nus = np.zeros((Nc, Nz, len(M)))
    for i in range(Nc):
        obh2, och2, w, ns, ln10As, H0, Neff, s8 = cosmos[i]
        h = H0/100.
        Om = (obh2+och2)/h**2
        kh = k/h #now h/Mpc
        for j in range(Nz):
            ph3 = p[i,j]*h**3 #now (Mpc/h)^3
            nus[i,j] = ctph.nu_at_M(M, kh, ph3, Om)
            continue
        print "Finished peak height box %d"%i
        continue
    return nus


def make_peak_height():
    np.savetxt("peak_height/M.txt", M, header="M [Msun/h]")
    k = np.loadtxt("plin/k.txt") #1/Mpc

    #First do the training cosmos
    training_cos = get_training_cosmos()
    plins = np.load("plin/plins_training_all_mpc3.npy") #[Mpc]^3
    nus = calc_peak_height(training_cos, zs, k, plins)
    np.save("peak_height/peak_height_training_all", nus)
    
    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    plins = np.load("plin/plins_testing_all_mpc3.npy") #[Mpc]^3
    nus = calc_peak_height(testing_cos, zs, k, plins)
    np.save("peak_height/peak_height_testing_all", nus)

    #Third do the highres cosmos
    hr_cos = get_highres_cosmos()
    plins = np.load("plin/plins_highres_all_mpc3.npy") #[Mpc]^3
    nus = calc_peak_height(hr_cos, zs_hr, k, plins)
    np.save("peak_height/peak_height_highres_all", nus)
    print "Finished all peak heights"
    return

def calc_hmf(cosmos, zs):
    import aemHMF
    hmf = aemHMF.Aemulus_HMF()

    Nc = len(cosmos)
    Nz = len(zs)
    dndlMs = np.zeros((Nc, Nz, len(M)))
    for i in range(Nc):
        Ombh2, Omch2, w, ns, ln10As, H0, Neff, sig8 = cosmos[i]
        cosmo={'Obh2':Ombh2, 'Och2':Omch2, 'w0':w, 'n_s':ns, 'ln10^{10}A_s':ln10As, 'N_eff':Neff, 'H0':H0}
        hmf.set_cosmology(cosmo)
        for j in range(Nz):
            z = zs[j]
            dndlMs[i,j] = hmf.dndlM(M,z)
            continue
        print "Finished with dndlM in box%d"%i
        continue
    return dndlMs

def make_mass_function():
    #Note: this is the Aemulus mass function (McClintock+ 2018)
    np.savetxt("mass_function/M.txt", M, header="M [Msun/h]")

    #First do the training cosmos
    training_cos = get_training_cosmos()
    dndlMs = calc_hmf(training_cos, zs)
    np.save("mass_function/dndlM_training_all", dndlMs)

    #Second do the testing cosmos
    testing_cos = get_testing_cosmos()
    dndlMs = calc_hmf(testing_cos, zs)
    np.save("mass_function/dndlM_testing_all", dndlMs)

    #Third do the highres cosmos
    hr_cos = get_highres_cosmos()
    dndlMs = calc_hmf(hr_cos, zs_hr)
    np.save("mass_function/dndlM_highres_all", dndlMs)
    
    print "Finished all mass functions"
    return

def make_bias():
    #Note: this is the Aemulus bias (McClintock+ 201?)
    np.savetxt("bias/M.txt", M, header="M [Msun/h]")

    #Load the data that was computed elsewhere...
    print "Finished with bias. Computed elsewhere for now..."
    return

if __name__ == "__main__":
    #make_linear_power_spectra()
    #make_nonlinear_power_spectra()
    #make_linear_correlation_function()
    #make_nonlinear_correlation_function()
    #make_peak_height()
    #make_mass_function() #Issue in training box 4, 13, 19, 26
    #make_bias()

    print "Finished with file creation"
