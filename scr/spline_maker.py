import numpy
import operator
from functools import reduce

from glob import glob

# Fast sparse-matrix implementation
from photospline import glam_fit,ndsparse

import os

"""
Script to make structure functions spline tables.
C.A. Arguelles Delgado - aug.03.14
"""

#warnings.filterwarnings('error')

def ModLog10(x):
    if x <= 0. :
        # print(x)
        return -50
    else:
        return numpy.log10(x)

ModLog10 = numpy.vectorize(ModLog10)

def SplineFitMaker1D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 1, N = 50, outname = "", oscale = 'lin'):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the same filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x), asummes x to be the first column.
    """
    if(column < 1):
        print("Error: column < 1.")
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    of = lambda x : x
    if oscale == "log":
        of = lambda x : ModLog10(x)
    elif oscale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    #shape = tuple(numpy.unique(datas[:,i]).size for i in range(2))
    #datas = datas.reshape(shape + (datas.size/reduce(operator.mul,shape),))

    x = f(datas[:,0])
    z = of(datas[:,column])

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True)]
    order = [2]
    smooth = [1.0e-15]
    penaltyorder = [2]

    weight = numpy.ones(z.shape)
    zs,w = ndsparse.from_data(z,weight)
    result = glam_fit(zs,w,[x],knots,order,smooth,penaltyorder)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    result.write(prefix+nfilename)
    print("Done. Generated :"  + prefix+nfilename)

def SplineFitMaker2D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 2, N = 50, outname = ""):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y), asummes x y to be the first two columns.
    """
    if(column < 2):
        print("Error: column < 2.")
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    shape = tuple(numpy.unique(datas[:,i]).size for i in range(2))
    datas = datas.reshape(shape + (int(datas.size/reduce(operator.mul,shape)),))

    x = f(datas[:,0,0])
    y = f(datas[0,:,1])
    z = datas[:,:,column]

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True),numpy.linspace(y.min()-1,y.max()+1,N,endpoint = True)]
    #knots = [x,y]
    order = [2 for _ in range(2)]
    #smooth = 1.0e-5
    smooth = [1.0e-15 for _ in range(2)]
    penaltyorder = [2 for _ in range(2)]

    weight = numpy.ones(z.shape)
    #weight = 1+zz
    zs,w = ndsparse.from_data(z,weight)
    result = glam_fit(zs,w,[x,y],knots,order,smooth,penaltyorder)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    result.write(prefix+nfilename)
    print("Done. Generated :"  + prefix+nfilename)

def SplineFitMaker3D(filename, scale = 'lin', prefix = '', skip_header = 0, column = 2, N = 50, outname = "", oscale = 'lin', smooth=1e-15):
    """
    Creates a spline table from a table of numbers. Then
    saves the table using the sanem filename as given.

    Options:
    scale : linear or log. For logbicubic splines or simple bicubic splines
    prefix : prefix for outputfile
    skip_header : skipe lines in input datafile
    column : z = f(x,y,w), asummes x/y/w to be the first/second/third column.
    """
    if(column < 3):
        print("Error: column < 3.")
        exit()

    datas = numpy.loadtxt(filename, skiprows = skip_header)

    f = lambda x : x
    if scale == "log":
        f = lambda x : ModLog10(x)
    elif scale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    of = lambda x : x
    if oscale == "log":
        of = lambda x : ModLog10(x)
    elif oscale == "lin":
        pass
    else:
        print("Error: unknown scale.")
        exit()

    shape = tuple(numpy.unique(datas[:,i]).size for i in range(3))
    datas = datas.reshape(shape + (int(datas.size/reduce(operator.mul,shape)),))

    x = f(datas[:,0,0,0])
    y = f(datas[0,:,0,1])
    w = f(datas[0,0,:,2])
    z = of(datas[:,:,:,column])

    knots = [numpy.linspace(x.min()-1,x.max()+1,N,endpoint = True),
             numpy.linspace(y.min()-1,y.max()+1,N,endpoint = True),
             numpy.linspace(w.min()-1,w.max()+1,N,endpoint = True)]
    #knots = [x,y]
    order = [2 for _ in range(3)]
    #smooth = 1.0e-5
    smooth = [smooth for _ in range(3)]
    penaltyorder = [2 for _ in range(3)]

    weight = numpy.ones(z.shape)
    #weight = 1+zz
    zs,weight = ndsparse.from_data(z,weight)
    result = glam_fit(zs,weight,[x,y,w],knots,order,smooth,penaltyorder)

    # creatiing new filename and saving
    if outname == "":
        nfilename = filename.split('.').pop() + ".fits"
    else :
        nfilename = outname
    # save the result to a FITS table
    if os.path.exists(prefix+nfilename):
        os.unlink(prefix+nfilename)
    # this file can later be loaded and evaluated with the photospline C library
    result.write(prefix+nfilename)
    print("Done. Generated :"  + prefix+nfilename)

if __name__ == "__main__":
    # set inpath
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/HNL/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/HNL_test/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/HNL_modified/"
    inpath_base = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/NuXSSplMkr/bin"

    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified_2/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified_3/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified_4/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified_5/"
    # inpath_base = "/data/user/lfischer/software/NuXSSplMkr/data/tau_modified_6/"

    # # take all masses
    # inpaths = glob(os.path.join(inpath_base, 'M_*'))
    # only a specific mass
    #inpaths = [os.path.join(inpath_base, 'M_0200MeV'),
    #           os.path.join(inpath_base, 'M_0300MeV')]
    inpaths = [os.path.join(inpath_base, 'M_0100MeV')]
    # inpaths = [os.path.join(inpath_base, 'M_0200MeV')]
    # inpaths = [os.path.join(inpath_base, 'M_0000MeV')]
    # inpaths = [os.path.join(inpath_base, 'M_0400MeV')]
    # inpaths = [os.path.join(inpath_base, 'M_1777MeV')]


    print('All input directories:\n', inpaths)

    neutrino_type = ['nu']
    # neutrino_type = ['nutau','nutaubar']

    # pdf_list = ['HERAPDF15NLO_EIG_central']
    pdf_list = ['PDF4LHC21_mc_central']

    for inpath in inpaths:
        print('This input directory:\n',inpath)
        outpath = inpath + '/'
        print('Outpath: {}'.format(outpath))

        # for int_type in ["cc","nc"]:
        #     for pdf in pdf_list:
        #         for neutype in neutrino_type:
        #             filename = "sigma-"+neutype+"-N-"+int_type+"-"+pdf
        #             print("processing: "+filename)
        #             infilepath = inpath + '/' + filename + ".dat"
        #             if not os.path.isfile(infilepath):continue
        #             print('Infilepath: {}'.format(infilepath))
        #             SplineFitMaker1D(infilepath, outname = filename + ".fits",
        #                     scale = 'log',prefix = outpath, N = 65, column = 1, oscale = 'log')

        for int_type in ["em"]:
            for smooth in [1e0,1e-5,1e-10]:
                for pdf in pdf_list:
                    for neutype in neutrino_type:
                        
                        # total
                        filename = "sigma-"+neutype+"-N-"+int_type+"-"+pdf
                        print("processing: "+filename)
                        infilepath = inpath + '/' + filename + ".dat"
                        if not os.path.isfile(infilepath):continue
                        print('Infilepath: {}'.format(infilepath))
                        SplineFitMaker1D(infilepath, outname = filename + ".fits",
                                scale = 'log',prefix = outpath, N = 65, column = 1, oscale = 'log')
                        
                        # differential
                        filename = "dsdxdy-"+neutype+"-N-"+int_type+"-"+pdf
                        print("processing: "+filename)
                        infilepath = inpath + '/' + filename + ".dat"
                        if not os.path.isfile(infilepath):continue
                        print('Infilepath: {}'.format(infilepath))
                        SplineFitMaker3D(infilepath, outname = filename + "_smooth%2.2e.fits"%smooth,
                                scale = 'log',prefix = outpath, N = 65, column = 3, oscale = 'log', smooth=smooth)

#     #### Start - Original code (carguelles) ####

#     inpath = "/home/carguelles/NuXSSplMkr/data/newxs_2017/"
#     outpath = "/home/carguelles/NuXSSplMkr/fits/newxs_2017/"

#     neutrino_type = ['numu','numubar']

#     pdf_list = ['CT10nlo_central','CT10nlo_minus','CT10nlo_plus',
#                 'HERAPDF15NLO_EIG_central','HERAPDF15NLO_EIG_minus','HERAPDF15NLO_EIG_plus',
#                 'NNPDF23_nlo_as_0118_central','NNPDF23_nlo_as_0118_minus','NNPDF23_nlo_as_0118_plus']

#     #pdf_list = ['NNPDF23_nlo_as_0118_central','NNPDF23_nlo_as_0118_minus','NNPDF23_nlo_as_0118_plus']
#     pdf_list = ['HERAPDF15NLO_EIG_central']

#     #pdf = "HERAPDF15NLO_EIG_central"
#     #filename = "dsdxdy-numu-N-cc-"+pdf
#     #SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
#     #        scale = 'log',prefix = outpath, N = 50, column = 3, oscale = 'log')
# #
# #    filename = "dsdxdy-numubar-N-cc-"+pdf
# #    SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
# #            scale = 'log',prefix = outpath, N = 50, column = 3, oscale = 'log' )

# #    quit()
#     for int_type in ["cc","nc"]:
#         for pdf in pdf_list:
#             for neutype in neutrino_type:
#                 filename = "sigma-"+neutype+"-N-"+int_type+"-"+pdf
#                 print "processing: "+filename
#                 SplineFitMaker1D(inpath + filename + ".dat", outname = filename + ".fits",
#                         scale = 'log',prefix = outpath, N = 65, column = 1, oscale = 'log')

#     exit()

#     for int_type in ["cc","nc"]:
#         for pdf in pdf_list:
#             for neutype in neutrino_type:
#                 filename = "dsdxdy-"+neutype+"-N-"+int_type+"-"+pdf
#                 print "processing: "+filename
#                 SplineFitMaker3D(inpath + filename + ".dat", outname = filename + ".fits",
#                         scale = 'log',prefix = outpath, N = 65, column = 3, oscale = 'log')

#     #### End - Original code (carguelles) ####
