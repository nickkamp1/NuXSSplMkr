import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import photospline
from itertools import cycle
from scipy.interpolate import PchipInterpolator,LinearNDInterpolator 
lines = ["--","-.",":"]
linecycler = cycle(lines)

M = "0600"
PDFset = "GRV98lo_patched"
noHNL = True
current = "nc"
nutype = "nubar"
nutype_orig = "nutaubar" if "bar" in nutype else "nutau"
plot_orig = False

for M in ["1000"]:
    for nutype in ["nu","nubar"]:
        nutype_orig = "nutaubar" if "bar" in nutype else "nutau"
        for current in ["em"]:

            figure_dir = "M_%sMeV/figures/%s/%s/"%(M,nutype,current)
            os.makedirs(figure_dir,exist_ok=True)

            data_tot = np.loadtxt("M_%sMeV/data/sigma-%s-N-%s-%s_central.dat"%(M,nutype,current,PDFset))
            spline_tot = photospline.SplineTable("M_%sMeV/splines/sigma-%s-N-%s-%s_central.fits"%(M,nutype,current,PDFset))
            plt.plot(data_tot[:,0],data_tot[:,1],color="black",label="True")
            plt.plot(data_tot[:,0],10**spline_tot.evaluate_simple([np.log10(data_tot[:,0])]),color="red",ls="-.",label="New Spline")
            if plot_orig:
                spline_orig_tot = photospline.SplineTable("M_%sMeV/splines/Original/sigma-%s-N-nc-GRV98lo_patched_central.fits"%(M,nutype_orig))
                plt.plot(data_tot[:,0],10**spline_orig_tot.evaluate_simple([np.log10(data_tot[:,0])]),color="blue",ls="-.",label="Orig Spline")
            plt.legend()
            plt.loglog()
            plt.xlabel(r"$E_\nu$ [GeV]")
            plt.ylabel(r"$\sigma~[{\rm cm}^{2}]$")
            plt.savefig("%s/tot.pdf"%figure_dir,dpi=50)
            plt.clf()

            plt.plot(data_tot[:,0],10**spline_tot.evaluate_simple([np.log10(data_tot[:,0])])/data_tot[:,1],color="red",ls="--",label="New Spline Fit")
            if plot_orig: plt.plot(data_tot[:,0],10**spline_orig_tot.evaluate_simple([np.log10(data_tot[:,0])])/data_tot[:,1],color="blue",ls="-.",label="Original Spline Fit")
            plt.legend()
            plt.loglog()
            plt.xlabel(r"$E_\nu$ [GeV]")
            plt.ylabel(r"$\sigma$ (Spline/True)")
            plt.savefig("%s/tot_ratio.pdf"%figure_dir,dpi=50)
            plt.clf()

            data = np.loadtxt("M_%sMeV/data/dsdxdy-%s-N-%s-%s_central.dat"%(M,nutype,current,PDFset))
            data_new = data
            spline_new = photospline.SplineTable("M_%sMeV/splines/dsdxdy-%s-N-%s-%s_central.fits"%(M,nutype,current,PDFset))

            if noHNL:
                data_noHNL = np.loadtxt("M_%sMeV/data/dsdxdy-%s-N-%s-%s_central_noHNL.dat"%(M,nutype,current,PDFset))
                spline_new_noHNL = photospline.SplineTable("M_0000MeV/dsdxdy-%s-N-%s-%s_central.fits"%(nutype,current,PDFset))
                data_new = data_noHNL
                spline_new = spline_new_noHNL 
            if plot_orig: spline_orig = photospline.SplineTable("M_%sMeV/splines/Original/dsdxdy-%s-N-nc-GRV98lo_patched_central.fits"%(M,nutype_orig))

            cmap = mpl.colormaps["prism"]
            energies = np.unique(data[:,0])
            energies = energies[int(len(energies)/5)::int(len(energies)/5)]
            for energy in energies:
                print(energy)
                fig,ax = plt.subplots(2,1,figsize=(6,8),sharex=True)
                fig.subplots_adjust(hspace=0)
                xs = np.unique(data[:,1])
                if current=="nc":
                    xs = xs[:-1:int(len(xs)/5)]
                else:
                    xs = xs[-int(len(xs)/2)::int(len(xs)/50)]
                plot_anything = False
                for i,x in enumerate(xs):
                    dsigdxy = data[np.logical_and(data[:,0]==energy,data[:,1]==x)]
                    dsigdxy_new = data_new[np.logical_and(data_new[:,0]==energy,data_new[:,1]==x)]
                    if sum(dsigdxy[:,3])<=0: continue
                    if sum(dsigdxy_new[:,3])<=0: continue
                    plot_anything = True
                    true_xs,true_xs_new = dsigdxy[:,3],dsigdxy_new[:,3]
                    true_xs = np.where(true_xs<=0,1e-50,true_xs)
                    true_xs_new = np.where(true_xs_new<=0,1e-50,true_xs_new)

                    # total plots
                    ax[0].plot([],[],ls="-",color=cmap(i/len(xs)),label="x=%2.2e"%x)
                    ax[0].plot(dsigdxy[:,2],true_xs,color="black")
                    if noHNL: ax[0].plot(dsigdxy_new[:,2],true_xs_new,color="grey")
                    ax[0].plot(dsigdxy[:,2],10**spline_new.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])]),color=cmap(i/len(xs)),ls="--")
                    if plot_orig: ax[0].plot(dsigdxy[:,2],10**spline_orig.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])]),color=cmap(i/len(xs)),ls="-.")

                    # ratio plots
                    ax[1].plot(dsigdxy[:,2],true_xs/true_xs,color="black")
                    ax[1].plot(dsigdxy[:,2],10**spline_new.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])])/true_xs_new,color=cmap(i/len(xs)),ls="--")
                    if plot_orig: ax[1].plot(dsigdxy[:,2],10**spline_orig.evaluate_simple([np.log10(energy),np.log10(x),np.log10(dsigdxy[:,2])])/true_xs,color=cmap(i/len(xs)),ls="-.")

                if not plot_anything: continue
                ax[0].plot([],[],color="black",label="True (w/ HNL)")
                if noHNL: ax[0].plot([],[],color="grey",label="True (w/o HNL)")
                ax[0].plot([],[],color="black",ls="--",label="New Spline Fit")
                if plot_orig: ax[0].plot([],[],color="black",ls="-.",label="Original Spline Fit")
                l = ax[0].legend(ncol=2)
                ax[0].loglog()
                ax[1].loglog()
                l.set_title(r"$E_\nu =$ %2.2e GeV"%energy)
                ax[1].set_xlabel("Bjorken y")
                ax[0].set_ylabel(r"$d^2\sigma /dxdy~[{\rm cm}^{2}]$")
                ax[1].set_ylabel(r"$d^2\sigma /dxdy$ (Spline/True)")
                ax[0].set_ylim(1e-2*min(true_xs_new),1e2*max(true_xs_new))
                ax[1].set_ylim(1e-2,1e2)
                plt.tight_layout()
                plt.savefig("M_%sMeV/figures/%s/%s/diff_%2.2e_log.pdf"%(M,nutype,current,energy),dpi=50)
                ax[0].set_xscale("linear")
                ax[1].set_xscale("linear")
                plt.tight_layout()
                plt.savefig("%s/diff_%2.2e_linear.pdf"%(figure_dir,energy),dpi=50)
                plt.clf()
                plt.close(fig)
