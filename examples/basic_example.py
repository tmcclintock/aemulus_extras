"""
This is a basic script to make an Extras object for the aemulus sims.
This shows what all the fields are, and how to plot things.
"""

import aemulus_extras

extra = aemulus_extras.Extras(0) #, testing=True) #for a testing sim

atts = dir(extra)

print "Extra object attributes:"
for at in atts:
    if at[0:2] != "__":
        print "\t%s"%(at)

#Plot the linear correlation function at z=3 of the sim
try:
    import matplotlib.pyplot as plt
    plt.loglog(extra.k, extra.P_lin[0])
    plt.show()
except ImportError:
    print "Install matplotlib to see a plot."
    pass
