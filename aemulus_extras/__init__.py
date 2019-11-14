"""This file contains the extras object, which has all of the 
analytic data for a given simulation.
"""
import inspect, os
import numpy as np
here = os.path.dirname(os.path.realpath(__file__))

class Extras(object):

    def __init__(self, index, testing=False, highres=False):
        """Constructor for the simulation extras.

        Args:
            index (int): index of the simulation of interest.

        """
        self.index = index
        if testing: name = "testing"
        elif highres: name = "highres"
        else: name = "training"

        #Load in all the data
        #First, the cosmology
        self.cosmology = np.loadtxt(here+"/cosmologies/%s_cosmologies.txt"%(name))[index]
        self.h = self.cosmology[5]/100. #Hubble constant
        #Second, the linear power spectrum
        self.k = np.loadtxt(here+"/plin/k.txt")/self.h
        self.P_lin = np.load(here+"/plin/plins_%s_all_mpc3.npy"%(name))[index]*self.h**3
        #Third, the nonlinear power spectrum
        self.P_nl = np.load(here+"/pnl/pnls_%s_all_mpc3.npy"%(name))[index]*self.h**3
        #Fourth, the linear matter correlation function
        self.r = np.loadtxt(here+"/xilin/r.txt")
        self.xi_lin = np.load(here+"/xilin/xilins_%s_all.npy"%(name))[index]
        #Fifth, the nonlinear matter correlation function
        self.xi_nl = np.load(here+"/xinl/xinls_%s_all.npy"%(name))[index]
        #Sixth, the peak height
        self.M = np.loadtxt(here+"/peak_height/M.txt")
        self.nu = np.load(here+"/peak_height/peak_height_%s_all.npy"%(name))[index]
        #Sixth, the halo mass function from aemulus
        self.dndlM = np.load(here+"/mass_function/dndlM_%s_all.npy"%(name))[index]
        #Seventh, the halo bias from aemulus
        self.bias = np.load(here+"/bias/bias_%s_all.npy"%(name))[index]
        #Eigth, the emulator predicted matter correlation function
        self.r_ximm = np.loadtxt(here+"/ximm/r.txt")
        self.xi_mm = np.load(here+"/ximm/ximms_%s_all.npy"%(name))[index]
        return

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    for i in range(40):
        e = Extras(i)#, testing=True)
        for j in range(10):
            #plt.loglog(e.k, e.P_lin[j])
            #plt.loglog(e.r, e.xi_nl[j])
            #plt.loglog(e.M, e.nu[j])
            
            #print np.argmax(e.dndlM[j][:-1] - e.dndlM[j][1:])
            #print e.M[612:616], "masses"
            #print e.dndlM[j][612:616]
            #print e.nu[j][612:616]
            plt.loglog(e.M, e.dndlM[j]) #issue with box 4-4, 6-7, 7, 
            #plt.loglog(e.M, e.bias[j])
        plt.title("box %d"%i)
        plt.ylim(1e-10,1e-1) #only for the mass function
        plt.show()
