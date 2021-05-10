import numpy as np
#import matplotlib.pyplot as plt
import os, sys

#filename = '/n/des/lee.5922/programs/cosmolike/MG_musigma/cov/cov_cmass_wtheta_only'
def calling_cmasscov(filename):
    cov_data = np.genfromtxt(filename)
    cov = np.zeros((500,500))
    for i in range(cov_data.shape[0]):
        ind1, ind2, _,_,_,_,_,_, co,_ = cov_data[i,:]
        cov[int(ind1)][int(ind2)] = co
        cov[int(ind2)][int(ind1)] = co
    return cov

def update_fits_cov( filename=None, template=None, cov=None):
    
    filename_update = filename.split('.fits')[0] + '_cov.fits'
    #template = '/n/des/lee.5922/programs/cosmolike/DMASS-analysis/simulated_data//simulated_y1_dmass_3x2pt_test.fits '
    print filename
    print filename_update
    
    os.system('cp '+template+' '+filename_update)
    from astropy.io.fits import getdata, update    
    _, hdr = getdata(filename_update, 1, header=True)
    updatedfits = update(filename_update, cov, header=hdr, ext=1)
    
    for ext in range(2,8):
        #data, hdr = getdata(filename, ext-1, header=True)
        data, hdr = getdata(filename, ext, header=True)
        updatedfits = update(filename_update, data, header=hdr, ext=ext)
        
        
def main():
    
    
    #covfilename = '/n/des/lee.5922/programs/cosmolike/lighthouse_cov/output_dmass_3x2pt_Y1_neff/cov_3x2pt_mcal4_dmass_pcut_sysweight_NG_Y1_neff'
    covfilename = '/n/des/lee.5922/programs/cosmolike/DMASS-analysis/data/cov/cov_3x2pt_mcal4_dmass_pcut_sysweight_NG_Y1_neff'
    cov_cosmolike=calling_cmasscov(covfilename)


    rootdir = '/n/des/lee.5922/programs/cosmolike/DMASS-analysis/simulated_data/'
    template = rootdir+'simulated_y1_dmass_3x2pt_template.fits '

    #syslist = ['baseline', 'baryons', 'iatatt', 'nolimber', 'nolimberRSD', 'magnification', 'nonlinearbiasb2']
    syslist = ['neff_baseline']

    for sys in syslist:
        filename = rootdir+'simulated_y1_dmass_3x2pt_'+sys+'.fits'
        update_fitsfile_withcov( filename=filename, template=template, cov=cov_cosmolike)
        
        