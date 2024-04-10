import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import photospline
from itertools import cycle
lines = ["-","--","-.",":"]
linecycler = cycle(lines)

for current in ["nc","em"]:

    data = np.loadtxt("M_0100MeV/dsdxdy-nu-N-%s-PDF4LHC21_mc_central.dat"%current)
    splines = {}
    for smooth in [1e-15,1e-10,1e-5]:
        splines[smooth] = photospline.SplineTable("M_0100MeV/dsdxdy-nu-N-%s-PDF4LHC21_mc_central_smooth%2.2e.fits"%(current,smooth))

    cmap = mpl.colormaps["prism"]
    energies = np.unique(data[:,0])
    energies = energies[int(len(energies)/5)::int(len(energies)/5)]
    for energy in energies:
        print(energy)
        xs = np.unique(data[:,1])
        if current=="nc":
            xs = xs[:-1:int(len(xs)/3)]
        else:
            xs = xs[-int(len(xs)/5)::int(len(xs)/15)]
        plot_anything = False
        for i,x in enumerate(xs):
            dsigdxy = data[np.logical_and(data[:,0]==energy,data[:,1]==x)]
            if sum(dsigdxy[:,3])<=0: continue
            plot_anything = True
            plt.plot(dsigdxy[:,2],dsigdxy[:,3],lw=3,color=cmap(i/len(xs)),label="x=%2.2e"%x)
            for smooth,spline in splines.items():
                plt.plot(dsigdxy[:,2],10**spline.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])]),color=cmap(i/len(xs)),ls=next(linecycler))
        if not plot_anything: continue
        plt.plot([],[],lw=3,color="black",label="True")
        for ls,smooth in zip(lines,splines.keys()):
            plt.plot([],[],color="black",ls=ls,label="Spline Fit: %2.2e smoothness"%smooth)
        plt.legend(ncol=2)
        plt.loglog()
        plt.title(r"$E_\nu =$ %2.2e GeV"%energy)
        plt.xlabel("Bjorken y")
        plt.ylabel(r"$d^2\sigma /dxdy~[{\rm cm}^{2}]$")
        plt.savefig("%s/%2.2e.pdf"%(current,energy),dpi=50)
        plt.clf()
        
        for i,x in enumerate(xs):
            dsigdxy = data[np.logical_and(data[:,0]==energy,data[:,1]==x)]
            if sum(dsigdxy[:,3])<=0: continue
            plot_anything = True
            for smooth,spline in splines.items():
                plt.plot(dsigdxy[:,2],10**spline.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])])/dsigdxy[:,3],color=cmap(i/len(xs)),ls=next(linecycler))
        if not plot_anything: continue
        for ls,smooth in zip(lines,splines.keys()):
            plt.plot([],[],color="black",ls=ls,label="Spline Fit: %2.2e smoothness"%smooth)
        plt.legend()
        plt.semilogx()
        plt.ylim(0,2)
        plt.title(r"$E_\nu =$ %2.2e GeV"%energy)
        plt.xlabel("Bjorken y")
        plt.ylabel(r"$d^2\sigma /dxdy$ (Spline/True)")
        plt.savefig("%s/ratio_%2.2e.pdf"%(current,energy),dpi=50)
        plt.clf()
