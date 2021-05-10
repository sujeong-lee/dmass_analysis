import numpy as np

def read_cov( filename = None ):
    
    covdat = np.genfromtxt(filename)
    cov = np.zeros((80,80))
    
    nline = covdat.shape[0]
    for k in range(nline):
        n_i, n_j, _,_,_,_,_,_, covG, covNG = covdat[k, :]
        cov[ int(n_i),int(n_j)] = covG + covNG
    return cov
    
def divide_sample_in_sysmap( cat= None, sysmap=None, cutvalue=None  ):
    hpind_high = sysmap['PIXEL'][(sysmap['SIGNAL'] > cutvalue)]
    hpind_low = sysmap['PIXEL'][~(sysmap['SIGNAL'] > cutvalue )]
    mask_high = np.in1d( cat['HPIX'], hpind_high  )
    mask_low = np.in1d( cat['HPIX'], hpind_low  )
    print 'cut value=', cutvalue
    print 'hpixel=', np.sum(sysmap['SIGNAL'] <= cutvalue), np.sum(sysmap['SIGNAL'] > cutvalue)
    print 'sample size (low, high)=', np.sum(mask_low), np.sum(mask_high)
    return cat[mask_low], cat[mask_high]

def chi2_gammat( R_amp, gammat, theory, invcov):
    dv = gammat- R_amp * theory
    chi = np.dot( dv, np.dot(invcov, dv.T))
    return chi

def fitting_gammat_amplitude( gammat=None, gammat_theory=None, cov=None, cut=0):

    import scipy
    from scipy import optimize

    gammat = gammat[cut:]
    gammat_theory = gammat_theory[cut:]
    invcov = np.linalg.inv(cov[cut:, cut:])

    xopt= scipy.optimize.fmin( chi2_gammat, 0.1, args=( gammat, gammat_theory, invcov ), ftol=1e-10, retall=True, full_output=True)
    return xopt[0]

def fitting_gammat_amplitude_jk(datajk = None, theory=None, cov = None, cut=10):
    
    Ntheta, Njk = datajk.shape
    R_amp_source = []
    #for i in range(len(filenames)):
    for i in range(Njk):
        gammat_i = datajk[:,i]
        #_,_,_,gammat_i,_,_,_,_ = np.genfromtxt(filenames[i], unpack=True)
        R_amp = fitting_gammat_amplitude( gammat=gammat_i, gammat_theory=theory, cov=cov, cut=cut)
        R_amp_source.append(R_amp)

    return R_amp_source 



