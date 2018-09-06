"""This file contains the extras object, which has all of the 
analytic data for a given simulation.
"""
import inspect, os
import numpy as np
here = os.path.dirname(os.path.realpath(__file__))
print here, "is here"

class Extras(object):

    def __init__(self, index, testing=False):
        """Constructor for the simulation extras.

        Args:
            index (int): index of the simulation of interest.

        """
        self.index = index
        if testing: name = "testing"
        else: name = "training"

        #Load in all the data
        #First, the cosmology
        self.cosmology = np.loadtxt(here+"/%s_cosmologies.txt"%(name))[index]
        self.h = self.cosmology[5]/100. #Hubble constant
        #Second, the linear power spectrum
        self.k = np.loadtxt(here+"/plin/k.txt")/self.h
        self.P_lins = np.load(here+"/plin/plins_%s_all_mpc3.npy"%(name))[index]*self.h**3
        #Third, the nonlinear power spectrum
        self.P_nls = np.load(here+"/pnl/pnls_%s_all_mpc3.npy"%(name))[index]*self.h**3

if __name__ == "__main__":
    e = Extras(0)
    print e
