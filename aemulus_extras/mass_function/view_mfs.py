import numpy as np
import matplotlib.pyplot as plt

M = np.loadtxt("M.txt")
dn = np.load("dndlM_highres_all.npy")
for i in range(len(dn)):
    continue
    for j in range(len(dn[0])):
        plt.loglog(M, dn[i,j])
    plt.ylim(1e-12, 1e-2)
    plt.show()
    plt.clf()
    
dn = np.load("dndlM_testing_all.npy")
for i in range(len(dn)):
    continue
    for j in range(len(dn[0])):
        plt.loglog(M, dn[i,j])
    plt.ylim(1e-12, 1e-2)
    plt.show()
    plt.clf()


dn = np.load("dndlM_training_all.npy")
for i in range(len(dn)):
    for j in range(len(dn[0])):
        plt.loglog(M, dn[i,j])
    plt.ylim(1e-12, 1e-2)
    plt.show()
    plt.clf()
