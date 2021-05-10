
import os, sys, esutil, treecorr
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
#sys.path.append('../')
#sys.path.append('/n/des/lee.5922/programs/cosmolike/DMASS-analysis/measurements/code/')
#lib_path = os.path.abspath(os.path.join('..', '..', '..','..','Dropbox','repositories','CMASS', 'code'))
#sys.path.append(lib_path)
#from cmass_modules import Cuts
sys.path.append('/n/des/lee.5922/Dropbox/repositories/CMASS/code/')
from utils import *


def measuring_boostfactor():

    lens, randoms = calling_lens_catalog()
    source = calling_source_catalog(catname='/n/des/lee.5922/data/shear/y1a1-im3shape_v5_unblind_v2_matched_v4.source_flags_select.bpz.jkindex.fits')   

    # divide source bins
    #source1 = source[ source['ZBIN_IM3'] == 0 ]
    #source2 = source[ source['ZBIN_IM3'] == 1 ]
    #source3 = source[ source['ZBIN_IM3'] == 2 ]
    #source4 = source[ source['ZBIN_IM3'] == 3 ]

    from ggl import Boostfactor_jk, construct_jk_corr_from_boostjk

    _savedir = '../measurements/ggl_debug/'

    for i in range(4):
        boostdir = _savedir+'/source'+str(i+1)+'/boostfactor/'

        src_bin = source[ source['ZBIN_IM3'] == i ]
        #Boostfactor_jk( lens=lens, source = src_bin, randoms=randoms, 
        #        nside=8, dir=boostdir, 
        #        nbins=20, minsep=2.5, maxsep=250.0, binslop=0.05)
        construct_jk_corr_from_boostjk( boostdir = boostdir )

    # plotting
    fig, ax = plt.subplots()

    for i in range(4):
        label = 'source'+str(i+1)
        boostdir = _savedir+'/'+label+'/boostfactor/'
        meanr, boost_src, err_boost_src = np.genfromtxt(boostdir+'/construct_jk/jk_mean.txt', unpack=True)
        ax.errorbar( meanr, boost_src, yerr=err_boost_src, label=label )

    ax.axhline( y=0.0, color = 'grey', ls='--')
    ax.set_xscale('log')
    ax.legend(loc='best')
    fig.savefig(_savedir+'boost.png')
    print 'fig saved to ',_savedir+'boost.png'



def measurement_ggl():
    """ Callig Catalog """

    lens, randoms = calling_lens_catalog()
    source = calling_source_catalog(catname='/n/des/lee.5922/data/shear/y1a1-im3shape_v5_unblind_v2_matched_v4.source_flags_select.bpz.jkindex.fits')

    # divide source bins
    #source1 = source[ source['ZBIN_IM3'] == 0 ]
    #source2 = source[ source['ZBIN_IM3'] == 1 ]
    #source3 = source[ source['ZBIN_IM3'] == 2 ]
    #source4 = source[ source['ZBIN_IM3'] == 3 ]

    # measuring ggl 
    _savedir = '../measurements/ggl_debug/'

    for i in range(4):
        savedir = _savedir+'/source'+str(i+1)+'/'
        src_bin = source[ source['ZBIN_IM3'] == i ]
        ggl_jk_onejk_cross( lens=lens, source = src_bin, randoms=randoms, dir=savedir)
        construct_jk_corr_from_njkcorr(dir = savedir)

    # plotting
    fig, ax = plt.subplots()
    for i in range(4):
        savedir = _savedir+'/source'+str(i+1)+'/'
        label = 'source'+str(i+1)
        meanr, gammat, err_gammat = np.genfromtxt(savedir+'/construct_jk/jk_mean.txt', unpack=True)
        ax.errorbar( meanr, gammat, yerr=err_gammat, fmt='o', capsize=5, label=label )

    #ax.axhline( y=0.0, color = 'grey', ls='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')

    ax.set_xlabel('$\\theta$ [arcmin]')
    ax.set_ylabel('$\gamma_t$')
    ax.set_yticklabels([])
    fig.savefig(_savedir+'gammat.png')
    print 'fig saved to ',_savedir+'gammat.png'



    # correct BoostFactor
    for i in range(4):
        savedir = _savedir+'/source'+str(i+1)+'/'
        construct_jk_corr_from_njkcorr(dir = savedir, boost=True)


    fig, ax = plt.subplots()
    for i in range(4):
        savedir = _savedir+'/source'+str(i+1)+'/'
        label = 'source'+str(i+1)
        #meanr, gammat, err_gammat = np.genfromtxt(savedir+'/construct_jk/jk_mean.txt', unpack=True)
        meanr, gammat_corrected, err_gammat_corrected = np.genfromtxt(savedir+'ggl_boost_corrected/construct_jk/jk_mean.txt', unpack=True)
        #ax.errorbar( meanr, gammat, yerr=err_gammat, fmt='o', capsize = 5, label=label )
        ax.errorbar( meanr, gammat_corrected, yerr=err_gammat_corrected, fmt='o', capsize = 5, label=label )

    #ax.axhline( y=0.0, color = 'grey', ls='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$\\theta$ [arcmin]')
    ax.set_ylabel('$\gamma_t$')
    ax.set_yticklabels([])
    ax.legend(loc='best')
    fig.savefig(_savedir+'gammat_boost_corrected.png')
    print 'fig saved to ',_savedir+'gammat_boost_corrected.png'


def measuring_systematics():

    lens, randoms = calling_lens_catalog()
    source = calling_source_catalog(catname='/n/des/lee.5922/data/shear/y1a1-im3shape_v5_unblind_v2_matched_v4.source_flags_select.bpz.jkindex.fits')


    savedir_sys = _savedir+'observing_condition/'
    os.system('mkdir '+savedir_sys)

    rootdir = '/n/des/lee.5922/programs/cosmolike/DMASS-analysis/sysmaps/'
    sysmap_airmass = esutil.io.read(rootdir+'Y1A1GOLD_band_r_nside4096_AIRMASS_coaddweights3_mean.fits.gz')
    sysmap_skybrite = esutil.io.read(rootdir+'Y1A1GOLD_band_r_nside4096_SKYBRITE_coaddweights3_mean.fits.gz')
    sysmap_maglim = esutil.io.read(rootdir+'Y1A1GOLD_band_r_nside4096_maglimit3__.fits.gz')
    sysmap_fwhm = esutil.io.read(rootdir+'Y1A1GOLD_band_r_nside4096_FWHM_MEAN_coaddweights3_mean.fits.gz')

    # divide source galaxies in 
    from systematics import divide_sample_in_sysmap
    source_airmass_low, source_airmass_hi \
    = divide_sample_in_sysmap( cat= source, sysmap=sysmap_airmass, cutvalue=1.26  )
    source_skybrite_low, source_skybrite_hi \
    = divide_sample_in_sysmap( cat= source, sysmap=sysmap_skybrite, cutvalue=285  )
    source_maglim_low, source_maglim_hi \
    = divide_sample_in_sysmap( cat= source, sysmap=sysmap_maglim, cutvalue=23.92  )
    source_fwhm_low, source_fwhm_hi \
    = divide_sample_in_sysmap( cat= source, sysmap=sysmap_fwhm, cutvalue=3.68  )




    bins, bs = np.linspace(0, 2, 201, retstep=True)
    binc = bins[:-1]+bs/2.
    bin_low = bins[:-1]

    labels_sys = ['airmass', 'skybrite', 'maglim', 'fwhm']
    suffix_sys = ['low', 'hi']
    labels=[ l1+'_'+l2 for l1 in labels_sys for l2 in suffix_sys ]

    srclist = [ source_airmass_low, source_airmass_hi, source_skybrite_low, source_skybrite_hi,
                source_maglim_low, source_maglim_hi, source_fwhm_low, source_fwhm_hi]

    # calculate jk ggl first
    for la, src in zip(labels, srclist):
        print 'observing condition - ggl_jk ', la
        savedir_sys_ = savedir_sys+'/'+la +'/'
        os.system('mkdir '+savedir_sys_)
        ggl_jk_onejk_cross( lens=lens, source = src, randoms=randoms, dir=savedir_sys_)
        construct_jk_corr_from_njkcorr(dir = savedir_sys_)


    # calculate nz_src
    for la, src in zip(labels, srclist):
        print 'observing condition - nz ', la
        savedir_sys_ = savedir_sys+'/'+la +'/'
        
        nz_src, _ = np.histogram(src['Z_MC'], bins=bins, density=True)
        DAT = np.column_stack(( bin_low, nz_src ))
        filename = savedir_sys_+'/'+la+'.nz'
        np.savetxt(filename, DAT, header='zlow, nz')
        

    # calculate sigma_crit_eff
    from ggl import DeltaSigma
    DS = DeltaSigma()

    z_l, nz_l = np.loadtxt('/n/des/lee.5922/programs/cosmolike/cosmosis/dmass_cat/cmass_sgc.nz_2', unpack=True)
    z_l_step = 0.5*(z_l[3]-z_l[2])
    z_l += z_l_step

    sig_crit_dics = {}
    for la, src in zip(labels, srclist):
        print 'observing condition - sigma_crit ', la
        savedir_sys_ = savedir_sys+'/'+la +'/'
        z_src_low, nz_src = np.genfromtxt(savedir_sys_+'/'+la+'.nz', unpack=True)
        z_src_step = 0.5*(z_src_low[3]-z_src_low[2])
        z_s = z_src_low + z_src_step
    
        sig_crit_inv_eff = DS.sig_crit_inv_eff(z_l, nz_l, z_s, nz_src)
        sig_crit_dics[la] = sig_crit_inv_eff
    
    import pickle
    f = open(savedir_sys+'/sig_crit_inv_eff_im3.pkl', "wb")
    pickle.dump(sig_crit_dics, f)
    f.close()



    # calculate ratio of ggl_sys
    from ggl_sys import fitting_gammat_amplitude_jk, read_cov
    cov_combined_theory=\
    read_cov( filename = '/n/des/lee.5922/programs/cosmolike/lighthouse_cov/output_source_combined/cov_im3combined_dmass_pcut_sysweight_20bins_NG_lsls_cov_Ntheta20_Ntomo1_5')
    cov_combined_theory=cov_combined_theory[40:60, 40:60]

    FITSDATA_combined = fitsio.FITS('../../simulated_data/simulated_y1_dmass_3x2pt_source_combined.fits')
    gammat_combined_theory = FITSDATA_combined['galaxy_shear_xi']['VALUE'].read()
    theta = FITSDATA_combined['galaxy_shear_xi']['ANG'].read()

    #labels_sys_ = ['airmass', 'skybrite', 'maglim', 'fwhm']
    #suffix_sys_ = ['hi', 'low']
    R_amp_source_values = {}
    for la in labels_sys:
        for suf in suffix_sys:
            fullla = la +'_'+suf
            savedir_sys_ = savedir_sys+'/'+fullla +'/'
            #savedir = rootdir+fullla+'/construct_jk/'
            
            datajk = np.genfromtxt(savedir_sys_+'/construct_jk/jk.txt')
            covjk = np.genfromtxt(savedir_sys_+'/construct_jk/jk_cov.txt')
            mean_R_amp_source=\
            fitting_gammat_amplitude_jk(datajk=datajk, theory= gammat_combined_theory , cov=cov_combined_theory*2)
            
            mean_R_amp_source = np.array(mean_R_amp_source)
            R_amp_source_values['mean_'+fullla] = np.mean(mean_R_amp_source)
            R_amp_source_values['std_'+fullla] = np.std(mean_R_amp_source) * np.sqrt((len(mean_R_amp_source) -1 ) )
            #mean_R_amp_source_values.append(mean_R_amp_source)
            #std_R_amp_source_values.append(std_R_amp_source)
            
        Rhi = R_amp_source_values['mean_'+la+'_hi']
        Rlo = R_amp_source_values['mean_'+la+'_low']
        stdhi = R_amp_source_values['std_'+la+'_hi']
        stdlo = R_amp_source_values['std_'+la+'_low']
        
        ratio = Rlo/Rhi    
        err = ratio * np.sqrt( stdlo**2/Rlo**2 + stdhi**2/Rhi**2 )
        
        R_amp_source_values['ratio_'+la] = ratio
        R_amp_source_values['err_ratio_'+la] = err
        
    import pickle
    f2 = open(savedir_sys+'/R_amp_ratio_im3.pkl', "wb")
    pickle.dump(R_amp_source_values, f2)
    f2.close()



    ## systematics plotting
    f = open(savedir_sys+'/sig_crit_inv_eff_im3.pkl', "r")
    sig_crit_dics = pickle.load(f)

    f2 = open(savedir_sys+'/R_amp_ratio_im3.pkl', "r")
    R_amp_source_values = pickle.load(f2)

    nn = 1
    fig, ax = plt.subplots()
    #for la in labels:
    for la in labels_sys:
        sig_crit_dics['ratio_'+la] = sig_crit_dics[la+'_low']/sig_crit_dics[la+'_hi']
        ax.plot(nn, sig_crit_dics['ratio_'+la], 'sk', markersize=10)
        ax.errorbar(nn, R_amp_source_values['ratio_'+la], yerr=R_amp_source_values['err_ratio_'+la], fmt='ro')
        nn += 1
    ax.set_xticks([0,1,2,3,4,5])
    ax.set_ylabel( '$\gamma_t^{\\rm low} / \gamma_t^{\\rm high}$' )
    ax.set_xticklabels(['']+labels_sys+[''])
    fig.savefig(savedir_sys+'/observing_condition_ratio.png' )
    print 'fig saved to ', savedir_sys+'/observing_condition_ratio.png'

    



################
measuring_boostfactor()
measurement_ggl()
measuring_systematics()
