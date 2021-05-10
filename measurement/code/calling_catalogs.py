from __future__ import print_function, division, absolute_import, unicode_literals

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


def calling_source_catalog( catname = '../output_cats/source_im3.fits'):

    print ('Calling im3 shape source catalog...')
    if os.path.exists(catname):
        source_matched = esutil.io.read(catname, upper=True)

    else : 
        #columns = ['ra', 'dec', 'coadd_objects_id', 'flags', 'mask_frac', 'e1', 'e2', 'region' ]
        columns = ['ra', 'dec', 'coadd_objects_id', 'flags_select', 'e1', 'e2', 'weight' ]
        source = esutil.io.read('/n/des/lee.5922/data/shear/y1a1-im3shape_v5_unblind_v2_matched_v4.fits', 
                                #rows=1, 
                                columns=columns, 
                                upper=True)
        source = source[ (source['DEC'] < -35) & (source['FLAGS_SELECT']== 0)]
        # Add Healpix index, ring order nside=4096
        print 'Add Healpixel indeces...'
        hp_ind = hpRaDecToHEALPixel(source['RA'], source['DEC'], nside=4096, nest=False)
        source = appendColumn(cat=source, name='HPIX', value=hp_ind)
        

        # adding source bin index
        print 'Adding redshift bin index...'
        source_z_binning = esutil.io.read('/n/des/lee.5922/data/shear/y1_source_redshift_binning_v1.fits', upper=True)
        source_matched, source_z_matched = matchCatalogs(cat1=source, cat2=source_z_binning,tag='COADD_OBJECTS_ID')
        source_matched = appendColumn(cat=source_matched, name='ZBIN_IM3', value=source_z_matched['ZBIN_IM3'])

        source, source_z_matched = None, None
        #/n/des/lee.5922/data/shear/y1a1-im3shape_v5_unblind_v2_matched_v4.source_flags_select.bpz.jkindex.fits
        print 'Resulting catalog size:', source_matched.size
        print 'Saving calibrated im3shape catalog to ../cat_output/source_im3.fits'
        esutil.io.write(catname, source_matched)

    return source_matched



def calling_source_mcal_catalog(outputcatdir = '../output_cats/', combine=False, unsheared=False):

    print ('Calling mcal shape source catalog...'()
    source={}
    columns = ['ra', 'dec', 'coadd_objects_id', 'e1', 'e2', 'R11','R22', 'HPIX', 'flags_select' ]+['ZBIN_MCAL']
    src_sheared = esutil.io.read(outputcatdir+'mcal_sheared.fits', columns=columns, upper=True)
    clean_mask = (src_sheared['FLAGS_SELECT']==0) 

    src_sheared = src_sheared[clean_mask]
    source['sheared'] = {}
    source['unsheared'] = {}

    components = ['1p', '1m', '2p', '2m']
    if combine:
        mask_zbin_mcal = (src_sheared['ZBIN_MCAL'] == 0)|(src_sheared['ZBIN_MCAL'] == 1)|\
                         (src_sheared['ZBIN_MCAL'] == 2)|(src_sheared['ZBIN_MCAL'] == 3)                
        source['sheared']['cat'] = src_sheared[mask_zbin_mcal]
        if unsheared:
            for j, comp in enumerate(components):  
                columns = ['ra', 'dec', 'e1', 'e2', 'HPIX','flags_select' ] + ['ZBIN_MCAL_'+comp.upper()]
                source['unsheared'][comp] = {}
                src_unsheared = esutil.io.read(outputcatdir+'mcal_unsheared_'+comp+'.fits', columns=columns, upper=True)
                clean_mask = (src_unsheared['FLAGS_SELECT']==0) 
                mask_zbin_mcal = (src_sheared['ZBIN_MCAL'] == 0)|(src_sheared['ZBIN_MCAL'] == 1)|\
                                 (src_sheared['ZBIN_MCAL'] == 2)|(src_sheared['ZBIN_MCAL'] == 3) 
                source['unsheared'][comp] = src_unsheared[clean_mask & mask_zbin_mcal]

    else : 
        
        for i in range(4):
            source['sheared']['bin'+str(i+1)] = {}
            source['sheared']['bin'+str(i+1)]['cat'] = src_sheared[ src_sheared['ZBIN_MCAL'] == i ]
            source['unsheared']['bin'+str(i+1)] = {}

        if unsheared:
            for j, comp in enumerate(components):
                columns = ['ra', 'dec', 'coadd_objects_id', 'e1', 'e2', 'HPIX','flags_select'] + ['ZBIN_MCAL_'+comp.upper()]
                src_unsheared = esutil.io.read(outputcatdir+'mcal_unsheared_'+comp+'.fits', columns=columns, upper=True)
                clean_mask = (src_unsheared['FLAGS_SELECT']==0) 
                src_unsheared = src_unsheared[clean_mask]

                for i in range(4):
                    source['unsheared']['bin'+str(i+1)][comp]= src_unsheared[src_unsheared['ZBIN_MCAL_'+comp.upper()] == i]
                    #print i, j

        """
            if unsheared:
                for j, comp in enumerate(components):    
                    columns = ['ra', 'dec', 'coadd_objects_id', 'e1', 'e2', 'HPIX','flags_select'] + ['ZBIN_MCAL_'+comp.upper()]
                    #columns_ = columns + ['ZBIN_MCAL_'+comp.upper()]
                    #source['unsheared']['bin'+str(i+1)][comp] = {}
                    src_unsheared = esutil.io.read(outputcatdir+'mcal_unsheared_'+comp+'.fits', columns=columns, upper=True)
                    clean_mask = (src_unsheared['FLAGS_SELECT']==0) 
                    src_unsheared = src_unsheared[clean_mask]
                    source['unsheared']['bin'+str(i+1)][comp]= src_unsheared[src_unsheared['ZBIN_MCAL_'+comp.upper()] == i]
                    #print i, j
        """

    src_sheared = None
    src_unsheared = None
    
    return source





def calling_lens_catalog(catname='../cat_output/lens.fits'):

    print ('Calling DMASS catalogs and corresponding randoms')
    if os.path.exists(catname):
        dmass = esutil.io.read(catname)

    else : 
        catdir = ''.join([ c+'/' for c in catname.split('/')[:-1]])
        os.system('mkdir '+catdir)
        dmass = esutil.io.read('/n/des/lee.5922/data/dmass_cat/dmass_spt_sys_v3.fits')
        w_dmass = dmass['CMASS_PROB'] *dmass['WEIGHT0_fwhm_r']*dmass['WEIGHT1_airmass_z']
        print 'Calculatig DMASS systematic weights...'
        dmass = appendColumn(dmass, name='WEIGHT', value= w_dmass )
        dmass = dmass[ dmass['CMASS_PROB'] > 0.01 ]
        esutil.io.write(catname, dmass)
    #randoms = esutil.io.read('/n/des/lee.5922/data/dmass_cat/random_x50_dmass_spt_masked.fits')

    randoms = esutil.io.read('/n/des/lee.5922/data/dmass_cat/random_x50_dmass_spt_masked.fits')

    print ('Resulting catalog size')
    print ('DMASS=', np.sum(dmass['WEIGHT']) )
    print ('randoms=', randoms.size)
    return dmass, randoms
    


def keepGoodRegion( source, combine=True, unsheared=False, maskfile='/n/des/lee.5922/data/dmass_cat/mask/Y1LSSmask_v2_redlimcut_il22_seeil4.0_4096ring.fits'):

    print ('Masking with the Y1BAOLSS mask')
    print ('maskfile=', maskfile)
    dmass_mask_fits = esutil.io.read(maskfile)
    
    components = ['1p', '1m', '2p', '2m']
    if combine:
        mask_sheared = np.in1d( source['sheared']['cat']['HPIX'], dmass_mask_fits['PIXEL'])
        source['sheared']['cat'] = source['sheared']['cat'][mask_sheared]
        if unsheared:
            for j, comp in enumerate(components):  
                mask_unsheared = np.in1d( source['unsheared'][comp]['HPIX'], dmass_mask_fits['PIXEL'])
                source['unsheared'][comp] = source['unsheared'][comp][mask_unsheared]

    else : 
        for i in range(4):
            mask_sheared = np.in1d( source['sheared']['bin'+str(i+1)]['cat']['HPIX'], dmass_mask_fits['PIXEL'])
            source['sheared']['bin'+str(i+1)]['cat'] = source['sheared']['bin'+str(i+1)]['cat'][mask_sheared]
            print (i, mask_sheared.size, np.sum(mask_sheared) )
            if unsheared:
                for j, comp in enumerate(components): 
                    src_unsheared = source['unsheared']['bin'+str(i+1)][comp]
                    mask_unsheared = np.in1d( src_unsheared['HPIX'], dmass_mask_fits['PIXEL'])
                    source['unsheared']['bin'+str(i+1)][comp]= src_unsheared[mask_unsheared]
                    print (i, comp, mask_unsheared.size, np.sum(mask_unsheared) )

    #src_sheared = None
    src_unsheared = None
    
    return source