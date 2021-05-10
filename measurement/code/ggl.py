import numpy as np
import healpy as hp
import os, sys
sys.path.append('/n/des/lee.5922/programs/cosmolike/DMASS-analysis/measurements/code/')
from utils import hpRaDecToHEALPixel
import treecorr

def _BoostFactor( lens=None, randoms=None, filename=None, nbins=20, minsep=2.5, maxsep=250.0, binslop=0.05):
    
    # Lens and Randoms
    #w_lens_sys = lens['WEIGHT0_fwhm_r']* lens['WEIGHT1_airmass_z']
    w_lens = lens['WEIGHT'] # lens['CMASS_PROB']*w_lens_sys
    w_randoms = None #randoms['VETO']
    lens_cat = treecorr.Catalog( ra=lens['RA'], dec=lens['DEC'], w=w_lens, 
                                  ra_units='deg', dec_units='deg')

    randoms_cat = treecorr.Catalog( ra=randoms['RA'], dec=randoms['DEC'], w=w_randoms, is_rand=True,
                                  ra_units='deg', dec_units='deg')
    
    config={}
    config['nbins'] = nbins
    config['min_sep']= minsep
    config['max_sep']= maxsep
    config['bin_slop']=binslop
    config['sep_units'] = 'arcmin'
    config['bin_slop'] = binslop
    config['verbose'] = 2
    
    nn_lens = treecorr.NNCorrelation(config)
    nn_rand = treecorr.NNCorrelation(config)
    nn_lens.process(lens_cat, lens_cat)
    nn_rand.process(randoms_cat, randoms_cat)


    nn_lens = treecorr.NNCorrelation(config)
    nn_lens.process(lens_cat, source_cat)   # Compute the cross-correlation.

    nn_rand = treecorr.NNCorrelation(config)
    nn_rand.process(randoms_cat, source_cat)   # Compute the cross-correlation.

    
    N_rand = np.sum(nn_rand.weight) #
    #N_rand = randoms.size * np.sum( source['WEIGHT'] )
    N_lens = np.sum(nn_lens.weight) #
    #N_lens = np.sum(dmass_masked['CMASS_PROB']) * np.sum( source['WEIGHT'] )
    N_lens_normed = nn_lens.weight * 1./N_lens
    N_rand_normed = nn_rand.weight * 1./N_rand
    B = N_lens_normed * 1./N_rand_normed
    
    
    nn_lens.write(filename, rr=nn_rand)     # Write out to a file.
    
    
    xi, varxi = nn_lens.calculateXi(rr=nn_rand)
    #ErrB = np.sqrt( nn_lens.weight )* 1./N_lens# * 1./N_rand_normed
    #ErrB = 1./np.sqrt(N_rand_normed)
    return nn_lens.meanr, xi, varxi #B, ErrB #, nn_lens.weight, nn_rand.weight


def BoostFactor( lens=None, source=None, randoms=None, filename=None, nbins=20, minsep=2.5, maxsep=250.0, binslop=0.05):
    
    # Lens and Randoms
    #w_lens_sys = lens['WEIGHT0_fwhm_r']* lens['WEIGHT1_airmass_z']
    w_lens = lens['WEIGHT'] # lens['CMASS_PROB']*w_lens_sys
    w_randoms = None #randoms['VETO']
    lens_cat = treecorr.Catalog( ra=lens['RA'], dec=lens['DEC'], w=w_lens, 
                                  ra_units='deg', dec_units='deg')

    randoms_cat = treecorr.Catalog( ra=randoms['RA'], dec=randoms['DEC'], w=w_randoms, is_rand=True,
                                  ra_units='deg', dec_units='deg')

    source_cat = treecorr.Catalog(ra=source['RA'], dec=source['DEC'], 
                                  g1=source['E1'], g2=source['E2'], 
                                  w=source['WEIGHT'], 
                                  ra_units='deg', dec_units='deg')
    
    config={}
    config['nbins'] = nbins
    config['min_sep']= minsep
    config['max_sep']= maxsep
    config['bin_slop']=binslop
    config['sep_units'] = 'arcmin'
    config['bin_slop'] = binslop
    config['verbose'] = 1
    
    
    nn_lens = treecorr.NNCorrelation(config)
    nn_lens.process(lens_cat, source_cat)   # Compute the cross-correlation.

    nn_rand = treecorr.NNCorrelation(config)
    nn_rand.process(randoms_cat, source_cat)   # Compute the cross-correlation.
    
    nn_lens.write(filename, rr=nn_rand)     # Write out to a file.


def ggl( lens=None, source = None, randoms=None, filename=None, nbins=20, minsep=2.5, maxsep=250.0, binslop=0.05, verbose=1, cat_type='im3'):
    
    # Lens and Randoms
    #w_lens_sys = lens['WEIGHT0_fwhm_r']* lens['WEIGHT1_airmass_z']
    #w_lens = lens['CMASS_PROB']*w_lens_sys
    w_lens = lens['WEIGHT']

    w_randoms = None #randoms['VETO']
    lens_cat = treecorr.Catalog( ra=lens['RA'], dec=lens['DEC'], w=w_lens, 
                                  ra_units='deg', dec_units='deg')

    if randoms is not None : 
        randoms_cat = treecorr.Catalog( ra=randoms['RA'], dec=randoms['DEC'], is_rand=True,
                                    ra_units='deg', dec_units='deg')

    
    if cat_type == 'mcal':
        w_src = None
    elif cat_type == 'im3':
        w_src = source['WEIGHT']

    source_cat = treecorr.Catalog(ra=source['RA'], dec=source['DEC'], 
                                  g1=source['E1'], g2=source['E2'], 
                                  w = w_src, 
                                  ra_units='deg', dec_units='deg')
    config={}
    config['nbins'] = nbins
    config['min_sep']=minsep
    config['max_sep']=maxsep
    #config['bin_slop']=binslop
    config['sep_units'] = 'arcmin'
    #config['bin_slop'] = 0.05
    config['verbose'] = verbose


    ng_lens = treecorr.NGCorrelation(config)
    ng_lens.process(lens_cat, source_cat)   # Compute the cross-correlation.
    ng_lens.write(filename+'.lens')     # Write out to a file.
    #gammat_lens1 = ng_lens1.xi              # Or access the correlation function directly.

    if randoms is None : ng_rand=None
    else : 
        ng_rand = treecorr.NGCorrelation(config)
        ng_rand.process(randoms_cat, source_cat)   # Compute the cross-correlation.
        ng_rand.write(filename+'.rand') 

    ng_lens.write(filename, rg=ng_rand)     # Write out to a file.
    
    lens_cat.clear_cache()
    source_cat.clear_cache() 
    ng_lens.clear()

    if randoms is not None : 
        randoms_cat.clear_cache() 
        ng_rand.clear()
    
    """
    
    gammat, gammax, err_gammat = ng_lens.calculateXi(rg=ng_rand)
    
    meanr, meanlogr, weight, npairs= ng_lens.meanr, ng_lens.meanlogr, ng_lens.weight, ng_lens.npairs
    
    lens_cat.clear_cache()
    randoms_cat.clear_cache() 
    source_cat.clear_cache() 
    
    return 0, meanr, meanlogr, gammat, gammax, err_gammat, weight, npairs
    """

def ggl_jk_kmean( lens=None, source = None, randoms=None, 
           dir=None, 
           nbins=20, minsep=2.5, maxsep=250.0, binslop=None):
    
    njk = len( set(lens['JKindex']) )
    print 'njack=',njk
    os.system('mkdir -p '+dir+'/jk_kmean/')

    ind_jk_lens = lens['JKindex']
    ind_jk_src = source['JKindex']
    if randoms is not None : ind_jk_rand = randoms['JKindex']
    
    n = 1
    for i in range(njk):
        ngal = np.sum(ind_jk_lens != i)
        print 'njk/ntot={}/{}, ngal={}'.format(n,njk,ngal)
        filename = dir+'/jk_kmean/ggl_jk{:05}.txt'.format(i)
        
        if os.path.exists(filename): 
            print 'file exists.. ', filename
            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            #np.genfromtxt(filename, unpack=True)
            
        else : 
        
            lens_i = lens[ind_jk_lens != i]
            src_i = source[ind_jk_src != i]
            if randoms is None : rand_i = None
            else : rand_i = randoms[ind_jk_rand != i]

            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            ggl( lens=lens_i, 
                source =src_i, 
                randoms=rand_i, 
                filename=filename, nbins=nbins, 
                minsep=minsep, maxsep=maxsep, binslop=binslop)
            
            lens_i, src_i, rand_i = None, None, None
            
        #gammat_total.append( gammat )
        n += 1
    
    
    #os.system('mkdir '+dir+'../airmass_low/')
    #os.system('mkdir '+dir+'../airmass_low/jk/')
    gammat_total = []
    n = 1
    for i in range(njk):  
        filename = dir+'/jk_kmean/ggl_jk{:05}.txt'.format(i)
        print 'njk/ntot={}/{}'.format(n,njk), filename
        _,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
        np.genfromtxt(filename, unpack=True)  
        gammat_total.append( gammat )
        #os.system('cp '+filename+' '+dir+'../airmass_low/jk/.')
        n+=1

    gammat_total = np.array(gammat_total)
    gammat_avg = np.mean( gammat_total, axis = 0 )
    
    
    # covariance matrix
    xi_cov = np.zeros((gammat_avg.size, gammat_avg.size))
    for i in range(gammat_avg.size):
        for j in range(gammat_avg.size):
            xi_cov[i][j] = (njk-1.)*1./njk * np.sum( (gammat_total[:, i] - gammat_avg[i]) * (gammat_total[:,j]-gammat_avg[j] ))
            #xi_cov[i][j] = covkikj
    #inv = np.linalg.inv(xi_cov)
    np.savetxt(dir+'/jk_kmean/cov.txt', xi_cov)
    #xijkerr = np.sqrt(norm * xi_cov.diagonal())
    err_gammat = np.sqrt(xi_cov.diagonal())

    
    DAT = np.column_stack(( meanr, gammat_avg, err_gammat ))
    header = 'meanr, gammat_avg, std_gammat   \n#n_jack='+str(njk)
    np.savetxt(dir+'/jk_kmean/ggl_avg.txt', DAT, header=header)
    print 'file saved to ', dir+'/jk_kmean/ggl_avg.txt'



def ggl_jk_kmean_onejk( lens=None, source = None, randoms=None, 
           dir=None, 
           nbins=20, minsep=2.5, maxsep=250.0, binslop=None):
    
    njk = len( set(lens['JKindex']) )
    print 'njack=',njk
    os.system('mkdir -p '+dir+'/jk_kmean_onejk/')

    ind_jk_lens = lens['JKindex']
    ind_jk_src = source['JKindex']
    if randoms is not None : ind_jk_rand = randoms['JKindex']
    
    n = 1
    for i in range(njk):
        ngal = np.sum(ind_jk_lens != i)
        print 'njk/ntot={}/{}, ngal={}'.format(n,njk,ngal)
        filename = dir+'/jk_kmean_onejk/ggl_jk{:05}.txt'.format(i)
        
        if os.path.exists(filename): 
            print 'file exists.. ', filename
            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            #np.genfromtxt(filename, unpack=True)
            
        else : 
        
            lens_i = lens[ind_jk_lens == i]
            src_i = source[ind_jk_src == i]
            if randoms is None : rand_i = None
            else : rand_i = randoms[ind_jk_rand == i]

            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            ggl( lens=lens_i, 
                source =src_i, 
                randoms=rand_i, 
                filename=filename, nbins=nbins, 
                minsep=minsep, maxsep=maxsep, binslop=binslop)
            
            lens_i, src_i, rand_i = None, None, None
            
        #gammat_total.append( gammat )
        n += 1
    
    
    #os.system('mkdir '+dir+'../airmass_low/')
    #os.system('mkdir '+dir+'../airmass_low/jk/')
    gammat_total = []
    n = 1
    for i in range(njk):  
        filename = dir+'/jk_kmean_onejk/ggl_jk{:05}.txt'.format(i)
        print 'njk/ntot={}/{}'.format(n,njk), filename
        _,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
        np.genfromtxt(filename, unpack=True)  
        gammat_total.append( gammat )
        #os.system('cp '+filename+' '+dir+'../airmass_low/jk/.')
        n+=1

    gammat_total = np.array(gammat_total)
    gammat_avg = np.mean( gammat_total, axis = 0 )
    
    
    # covariance matrix

    """
    xi_cov = np.zeros((gammat_avg.size, gammat_avg.size))
    for i in range(gammat_avg.size):
        for j in range(gammat_avg.size):
            xi_cov[i][j] = (njk-1.)*1./njk * np.sum( (gammat_total[:, i] - gammat_avg[i]) * (gammat_total[:,j]-gammat_avg[j] ))
            #xi_cov[i][j] = covkikj
    #inv = np.linalg.inv(xi_cov)
    np.savetxt(dir+'/jk_kmean_onejk/cov.txt', xi_cov)
    #xijkerr = np.sqrt(norm * xi_cov.diagonal())
    err_gammat = np.sqrt(xi_cov.diagonal())

    
    DAT = np.column_stack(( meanr, gammat_avg, err_gammat ))
    header = 'meanr, gammat_avg, std_gammat   \n#n_jack='+str(njk)
    np.savetxt(dir+'/jk_kmean_onejk/ggl_avg.txt', DAT, header=header)
    print 'file saved to ', dir+'/jk_kmean_onejk/ggl_avg.txt'
    """

def ggl_jk( lens=None, source = None, randoms=None, 
           nside=8, dir=None, 
           nbins=20, minsep=2.5, maxsep=250.0, binslop=None):
    
    
    hpind_lens = hpRaDecToHEALPixel(lens['RA'], lens['DEC'], nside=nside )
    hpind_src = hpRaDecToHEALPixel( source['RA'], source['DEC'], nside=nside )
    hpind_rand = hpRaDecToHEALPixel( randoms['RA'], randoms['DEC'], nside=nside )
    
    hpix_max = np.max(np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    hpix_min = np.min( np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    
    
    njk = len( set(hpind_lens) )
    print 'njk=',njk
    os.system('mkdir '+dir+'/jk/')
    
    
    n = 1
    for i in set(hpind_lens):
        ngal = np.sum(hpind_lens != i)
        print 'njk/ntot={}/{}, ngal={}'.format(n,njk,ngal)
        filename = dir+'/jk/ggl_hpix{:05}.txt'.format(i)
        
        if os.path.exists(filename): 
            print 'file exists.. ', filename
            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            #np.genfromtxt(filename, unpack=True)
            
        else : 
        
            lens_i = lens[hpind_lens != i]
            src_i = source[hpind_src != i]
            rand_i = randoms[hpind_rand != i]

            #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
            ggl( lens=lens_i, 
                source =src_i, 
                randoms=rand_i, 
                filename=filename, nbins=nbins, 
                minsep=minsep, maxsep=maxsep, binslop=binslop)
            
            lens_i, src_i, rand_i = None, None, None
            
        #gammat_total.append( gammat )
        n += 1
    
    
    #os.system('mkdir '+dir+'../airmass_low/')
    #os.system('mkdir '+dir+'../airmass_low/jk/')
    gammat_total = []
    n = 1
    for i in set(hpind_lens):  
        filename = dir+'/jk/ggl_hpix{:05}.txt'.format(i)
        print 'njk/ntot={}/{}'.format(n,njk), filename
        _,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
        np.genfromtxt(filename, unpack=True)  
        gammat_total.append( gammat )
        #os.system('cp '+filename+' '+dir+'../airmass_low/jk/.')
        n+=1
    
    gammat_total = np.array(gammat_total)
    gammat_avg = np.mean( gammat_total, axis = 0 )
    DAT = np.column_stack(( meanr, gammat_avg ))
    np.savetxt(dir+'/jk/ggl_avg.txt', DAT)
    print 'file saved to ', dir+'/jk/ggl_avg.txt'
    
    
    # covariance matrix
    xi_cov = np.zeros((gammat_avg.size, gammat_avg.size))
    for i in range(gammat_avg.size):
        for j in range(gammat_avg.size):
            xi_cov[i][j] = (njk-1.)*1./njk * np.sum( (gammat_total[:, i] - gammat_avg[i]) * (gammat_total[:,j]-gammat_avg[j] ))
            #xi_cov[i][j] = covkikj
    #inv = np.linalg.inv(xi_cov)
    np.savetxt(dir+'/jk/cov.txt', xi_cov)
    #xijkerr = np.sqrt(norm * xi_cov.diagonal())
    

def ggl_jk_onejk( lens=None, source = None, randoms=None, 
           nside=8, dir=None, 
           nbins=20, minsep=2.5, maxsep=250.0, binslop=None):
    
    
    hpind_lens = hpRaDecToHEALPixel(lens['RA'], lens['DEC'], nside=nside )
    hpind_src = hpRaDecToHEALPixel( source['RA'], source['DEC'], nside=nside )
    hpind_rand = hpRaDecToHEALPixel( randoms['RA'], randoms['DEC'], nside=nside )
    
    hpix_max = np.max(np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    hpix_min = np.min( np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    
    
    njk = len( set(hpind_lens) )
    print 'njk=',njk
    os.system('mkdir '+dir+'/jk_onejk/')
    
    
    n = 1
    for i in set(hpind_lens):
        ngal = np.sum(hpind_lens == i)
        nsrc = np.sum(hpind_src == i)
        print 'njk/ntot={}/{}, ngal={}, nsrc={}'.format(n,njk,ngal, nsrc)
        filename = dir+'/jk_onejk/ggl_hpix{:05}.txt'.format(i)
        
        if nsrc == 0: pass
        else : 

            if os.path.exists(filename): 
                print 'file exists.. ', filename
                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                #np.genfromtxt(filename, unpack=True)
                
            else : 
            
                lens_i = lens[hpind_lens == i]
                src_i = source[hpind_src == i]
                rand_i = randoms[hpind_rand == i]

                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                ggl( lens=lens_i, 
                    source =src_i, 
                    randoms=rand_i, 
                    filename=filename, nbins=nbins, 
                    minsep=minsep, maxsep=maxsep, binslop=binslop)
                
                lens_i, src_i, rand_i = None, None, None
                
            #gammat_total.append( gammat )
        n += 1
    
    """
    #os.system('mkdir '+dir+'../airmass_low/')
    #os.system('mkdir '+dir+'../airmass_low/jk/')
    gammat_total = []
    n = 1
    for i in set(hpind_lens):  
        filename = dir+'/jk_onejk/ggl_hpix{:05}.txt'.format(i)
        print 'njk/ntot={}/{}'.format(n,njk), filename
        _,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
        np.genfromtxt(filename, unpack=True)  
        gammat_total.append( gammat )
        #os.system('cp '+filename+' '+dir+'../airmass_low/jk/.')
        n+=1
    
    gammat_total = np.array(gammat_total)
    gammat_avg = np.mean( gammat_total, axis = 0 )
    DAT = np.column_stack(( meanr, gammat_avg ))
    np.savetxt(dir+'/jk_onejk/ggl_avg.txt', DAT)
    print 'file saved to ', dir+'/jk/ggl_avg.txt'
    
    
    # covariance matrix
    xi_cov = np.zeros((gammat_avg.size, gammat_avg.size))
    for i in range(gammat_avg.size):
        for j in range(gammat_avg.size):
            xi_cov[i][j] = (njk-1.)*1./njk * np.sum( (gammat_total[:, i] - gammat_avg[i]) * (gammat_total[:,j]-gammat_avg[j] ))
            #xi_cov[i][j] = covkikj
    #inv = np.linalg.inv(xi_cov)
    np.savetxt(dir+'/jk/cov.txt', xi_cov)
    #xijkerr = np.sqrt(norm * xi_cov.diagonal())
    """

def compute_Rs(e_ix):
    mean_e1_1p = e_ix['1p'] #np.mean(source_unsheared['1p']['cat']['E1'])
    mean_e1_1m = e_ix['1m']#np.mean(source_unsheared['1m']['cat']['E1'])
    mean_e2_2p = e_ix['2p']#np.mean(source_unsheared['2p']['cat']['E2'])
    mean_e2_2m = e_ix['2m']#np.mean(source_unsheared['2m']['cat']['E2'])
    
    dgamma = 0.01*2
    Rs11_mean = ( mean_e1_1p + mean_e1_1m )/dgamma
    Rs22_mean = ( mean_e2_2p + mean_e2_2m )/dgamma
    
    Rs_mean = 0.5 * (Rs11_mean + Rs22_mean)
    return Rs_mean

def compute_Rgamma(lens, source, scalar, filename=None):
    import treecorr
    """
    Uses TreeCorr to compute the NK correlation between lens and source.
    Used to compute scale dependece responses.
    scalar: array with scalar values of some quantity of which we want to compute the
           average in angular bins.
    Returns theta and R_nk.
    """
    nk = treecorr.NKCorrelation(nbins=20, min_sep=2.5,
                                max_sep=250, sep_units='arcmin', bin_slop=0.05, 
                                num_threads=30, verbose=0)

    cat_l = treecorr.Catalog(ra=lens['RA'], dec=lens['DEC'], w=lens['WEIGHT'], ra_units='deg', dec_units='deg')
    cat_s = treecorr.Catalog(ra=source['RA'], dec=source['DEC'], k=scalar, ra_units='deg', dec_units='deg')
    nk.process(cat_l, cat_s)
    nk.write(filename)
    #theta = np.exp(nk.logr)
    #R_nk = nk.xi

    #return theta, R_nk

def _nkcorr_eix(lens, source, scalar, filename=None):
    """
    source : unsheared source
    scalar : e array
    """  
    import treecorr
    nk = treecorr.NKCorrelation(nbins=20, min_sep=2.5,
                            max_sep=250, sep_units='arcmin', bin_slop=0.05, 
                            num_threads=30, verbose=2)
    cat_l = treecorr.Catalog(ra=lens['RA'], dec=lens['DEC'], w=lens['WEIGHT'], ra_units='deg', dec_units='deg')
    cat_s = treecorr.Catalog(ra=source['RA'], dec=source['DEC'], k=scalar, ra_units='deg', dec_units='deg')
    nk.process(cat_l, cat_s)
    nk.write(filename)
    
def compute_eix(lens, source_unsheared, filename=None):
    #e_ix = {}
    components = ['1p', '1m', '2p', '2m']
    for i, comp in enumerate(components):
        source_comp = source_unsheared[comp]['cat']
        scalar = source_comp['E'+comp[0]]
        _filename = filename+'.'+comp
        _nkcorr_eix(lens, source_comp, scalar, filename=_filename)


def ggl_jk_onejk_cross( lens=None, source = None, source_unsheard=None, randoms=None, 
           nside=8, dir=None, cat_type='im3',
           nbins=20, minsep=2.5, maxsep=250.0, binslop=None):
    
    
    hpind_lens = hpRaDecToHEALPixel(lens['RA'], lens['DEC'], nside=nside )
    hpind_src = hpRaDecToHEALPixel( source['RA'], source['DEC'], nside=nside )
    hpind_rand = hpRaDecToHEALPixel( randoms['RA'], randoms['DEC'], nside=nside )

    hpix_max = np.max(np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    hpix_min = np.min( np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    
    njk = len( set(hpind_lens) )
    #njk_src = len( set(hpind_src) )
    print 'njk=',njk
    os.system('mkdir '+dir+'/jk/')
    
    n = 1
    for i in set(hpind_lens):

        neighsjk = hp.get_all_neighbours(nside, i)
        neighsjk = np.append(neighsjk, i)
        mask_src = np.in1d(hpind_src, neighsjk)
        njk_src = len(neighsjk)
        #for j in set(hpind_src):
        #for j in neighsjk:
        ngal = np.sum(hpind_lens == i)
        #nsrc = np.sum(hpind_src == j)
        nsrc = np.sum(mask_src)
        print 'njk/ntot={}/{}, nsrc={}, ngal={}, nsrc={}'.format(n,njk, nsrc, ngal, nsrc)
        #filename = dir+'/jk_onejk_cross/ggl_hpix_lens{:05}_src{:05}.txt'.format(i, j)
        filename = dir+'/jk/ggl_hpix{:05}.txt'.format(i)

        if nsrc == 0 : pass
        else : 
            if os.path.exists(filename): 
                print 'file exists.. ', filename
                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                #np.genfromtxt(filename, unpack=True)
                
            else : 
            
                lens_i = lens[hpind_lens == i]
                src_i = source[mask_src]
                rand_i = randoms[hpind_rand == i]

                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                ggl( lens=lens_i, 
                    source =src_i, 
                    randoms=rand_i, 
                    filename=filename, nbins=nbins, verbose = 0,
                    minsep=minsep, maxsep=maxsep, binslop=binslop, cat_type=cat_type)

                if cat_type=='mcal':
                    R11 = src_i['R11']
                    R22 = src_i['R22']
                    scalar = 0.5 * (R11+R22)
                    filename_R = dir+'/jk/Rgamma_hpix{:05}.txt'.format(i)
                    compute_Rgamma(lens_i, src_i, scalar, filename=filename_R)
                    filename_eix = dir+'/jk/eix_hpix{:05}.txt'.format(i)
                    compute_eix(lens, source_unsheared, filename=filename_eix)
                    print 'compute_eix'
                lens_i, src_i, rand_i = None, None, None
                
            #gammat_total.append( gammat )
            n += 1
    
    """
    #os.system('mkdir '+dir+'../airmass_low/')
    #os.system('mkdir '+dir+'../airmass_low/jk/')
    gammat_total = []
    n = 1
    for i in set(hpind_lens):  
        filename = dir+'/jk_onejk/ggl_hpix{:05}.txt'.format(i)
        print 'njk/ntot={}/{}'.format(n,njk), filename
        _,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
        np.genfromtxt(filename, unpack=True)  
        gammat_total.append( gammat )
        #os.system('cp '+filename+' '+dir+'../airmass_low/jk/.')
        n+=1
    
    gammat_total = np.array(gammat_total)
    gammat_avg = np.mean( gammat_total, axis = 0 )
    DAT = np.column_stack(( meanr, gammat_avg ))
    np.savetxt(dir+'/jk_onejk/ggl_avg.txt', DAT)
    print 'file saved to ', dir+'/jk/ggl_avg.txt'
    
    
    # covariance matrix
    xi_cov = np.zeros((gammat_avg.size, gammat_avg.size))
    for i in range(gammat_avg.size):
        for j in range(gammat_avg.size):
            xi_cov[i][j] = (njk-1.)*1./njk * np.sum( (gammat_total[:, i] - gammat_avg[i]) * (gammat_total[:,j]-gammat_avg[j] ))
            #xi_cov[i][j] = covkikj
    #inv = np.linalg.inv(xi_cov)
    np.savetxt(dir+'/jk/cov.txt', xi_cov)
    #xijkerr = np.sqrt(norm * xi_cov.diagonal())
    """


def Boostfactor_jk( lens=None, source = None, randoms=None, 
           nside=8, dir=None, 
           nbins=20, minsep=2.5, maxsep=250.0, binslop=0.05):
    
    
    hpind_lens = hpRaDecToHEALPixel(lens['RA'], lens['DEC'], nside=nside )
    hpind_src = hpRaDecToHEALPixel( source['RA'], source['DEC'], nside=nside )
    hpind_rand = hpRaDecToHEALPixel( randoms['RA'], randoms['DEC'], nside=nside )
    
    hpix_max = np.max(np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    hpix_min = np.min( np.hstack([hpind_lens, hpind_src, hpind_rand]) )
    
    
    njk = len( set(hpind_lens) )
    #njk_src = len( set(hpind_src) )
    print 'njk=',njk
    os.system('mkdir '+dir+'/jk/')
    
    
    n = 1
    for i in set(hpind_lens):

        neighsjk = hp.get_all_neighbours(nside, i)
        neighsjk = np.append(neighsjk, i)
        mask_src = np.in1d(hpind_src, neighsjk)
        njk_src = len(neighsjk)
        #for j in set(hpind_src):
        #for j in neighsjk:
        ngal = np.sum(hpind_lens == i)
        #nsrc = np.sum(hpind_src == j)
        nsrc = np.sum(mask_src)
        nrand = np.sum(hpind_rand == i)
        print 'njk/ntot={}/{}, nsrc={}, ngal={}, nrand={}'.format(n,njk, nsrc, ngal, nrand)
        #filename = dir+'/jk_onejk_cross/ggl_hpix_lens{:05}_src{:05}.txt'.format(i, j)
        filename = dir+'/jk/nn_hpix{:05}.txt'.format(i)

        if nsrc == 0 : pass
        else : 
            if os.path.exists(filename): 
                print 'file exists.. ', filename
                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                #np.genfromtxt(filename, unpack=True)
                
            else : 
            
                lens_i = lens[hpind_lens == i]
                #src_i = source[hpind_src == j]
                src_i = source[mask_src]
                rand_i = randoms[hpind_rand == i]

                #_,meanr,_,gammat,gammax,err_gammat,weight,npairs=\
                BoostFactor( lens=lens_i, 
                    source =src_i, 
                    randoms=rand_i, 
                    filename=filename, nbins=nbins, 
                    minsep=minsep, maxsep=maxsep, binslop=binslop)
                
                lens_i, src_i, rand_i = None, None, None
                
            #gammat_total.append( gammat )
            n += 1


def construct_jk_corr_from_boostjk( boostdir = None ):

    import glob
    filenamelist = glob.glob(boostdir+'/jk/nn_hpix*.txt')
    filenamelist.sort()

    DDlist = []
    RRlist = []
    for f in filenamelist:
        _,meanr,_,_,sigma,DD,RR,npairs = np.genfromtxt(f,unpack=True)
        DDlist.append( DD )
        RRlist.append( RR )

    gts =np.array(DDlist)
    grs =np.array(RRlist)

    gtnum = np.array([np.sum(np.delete(np.array(gts), i, axis = 0), axis = 0) for i in range(len(gts))])
    grnum = np.array([np.sum(np.delete(np.array(grs), i, axis = 0), axis = 0) for i in range(len(grs))])

    gt_all = gtnum/grnum -1.0
    Njk = len(gt_all)
    gt_mean = np.mean(gt_all, axis = 0)

    gt_cov = np.zeros((20,20))
    for i in range(20):
        for j in range(20):
            gt_cov[i,j] = np.sum((gt_all[:,i]-gt_mean[i])*(gt_all[:,j]-gt_mean[j])) * (Njk-1.)*1./Njk

    err_gt = np.sqrt( gt_cov.diagonal() )
    DAT = np.column_stack(( gt_all ))
    DAT2 = np.column_stack(( meanr, gt_mean, err_gt))
    header = 'Njk = {}'.format(Njk)
    header2 = 'Njk = {}'.format(Njk) + '\nmeanr   jkmean  jkerr'
    os.system('mkdir '+boostdir+'construct_jk/')
    print 'save results in ', boostdir+'construct_jk/'
    np.savetxt( boostdir+'construct_jk/jk.txt', DAT, header=header )
    np.savetxt( boostdir+'construct_jk/jk_mean.txt', DAT2, header=header2 )
    np.savetxt( boostdir+'construct_jk/jk_cov.txt', gt_cov, header=header )




def _load_jk_corr(dir, kind=None):
    
    import glob
    #filenamelist = []
    
    if kind is None : 
        filenamelist = glob.glob(dir+'ggl_hpix_*.txt')
        if len(filenamelist) == 0:filenamelist = glob.glob(dir+'ggl_jk*.txt')
    
    if kind is not None : 
        filenamelist = glob.glob(dir+'*.'+kind)
        #if len(filenamelist) == 0: filenamelist = glob.glob(dir+'*.'+kind)
        
    filenamelist.sort()
    """
        for f in filenames:
            if kind in f: filenamelist.append(f)
            else : pass
    else : 
        for f in filenames:
            if not f.endswith('.lens') and not f.endswith('.rand'):
                filenamelist.append(f)
                #filenamelist = filenames
    """
    
    gammatlist = []
    weightlist = []
    npairlist = []
    njk_lens = []
    njk_src = []
    
    for f in filenamelist:
        
        no_lens = int(f.split('/')[-1].split('lens')[1].split('_src')[0])
        no_src = int(f.split('/')[-1].split('lens')[1].split('_src')[1].split('.txt')[0])
        njk_lens.append(no_lens)
        njk_src.append(no_src)

        _,_,_,gammat,gammax,sigma,weight,npairs = np.genfromtxt(f,unpack=True)
        gammatlist.append( gammat )
        weightlist.append( weight )
        npairlist.append( npairs )

    njk_lens = len(set(np.array(njk_lens)))
    njk_src = len(set(np.array(njk_src)))
    print 'Nlens jk =', njk_lens, 
    print ' Nsrc jk =', njk_src

    gammatlist = np.array(gammatlist).reshape(njk_lens, njk_src, 20)
    weightlist = np.array(weightlist).reshape(njk_lens, njk_src, 20)
    npairlist = np.array(npairlist).reshape(njk_lens, njk_src, 20)
    
    return gammatlist, weightlist, npairlist



def load_jk_corr(dir, kind=None):
    
    import glob
    #filenamelist = []
    
    if kind is None : 
        filenamelist = glob.glob(dir+'ggl_hpix_*.txt')
        if len(filenamelist) == 0:filenamelist = glob.glob(dir+'ggl_jk*.txt')
    
    if kind is not None : 
        filenamelist = glob.glob(dir+'*.'+kind)
        #if len(filenamelist) == 0: filenamelist = glob.glob(dir+'*.'+kind)
        
    filenamelist.sort()
    """
        for f in filenames:
            if kind in f: filenamelist.append(f)
            else : pass
    else : 
        for f in filenames:
            if not f.endswith('.lens') and not f.endswith('.rand'):
                filenamelist.append(f)
                #filenamelist = filenames
    """
    
    gammatlist = []
    weightlist = []
    npairlist = []
    njk_lens = []
    njk_src = []
    
    for f in filenamelist:
        
        #no_lens = int(f.split('/')[-1].split('lens')[1].split('_src')[0])
        #no_src = int(f.split('/')[-1].split('lens')[1].split('_src')[1].split('.txt')[0])
        #njk_lens.append(no_lens)
        #njk_src.append(no_src)
        print f
        _,meanr,_,gammat,gammax,sigma,weight,npairs = np.genfromtxt(f,unpack=True)
        gammatlist.append( gammat )
        weightlist.append( weight )
        npairlist.append( npairs )

    #njk_lens = len(set(np.array(njk_lens)))
    #njk_src = len(set(np.array(njk_src)))
    #print 'Nlens jk =', njk_lens, 
    #print ' Nsrc jk =', njk_src

    gammatlist = np.array(gammatlist) #.reshape(njk_lens, njk_src, 20)
    weightlist = np.array(weightlist) #.reshape(njk_lens, njk_src, 20)
    npairlist = np.array(npairlist) #.reshape(njk_lens, njk_src, 20)
    
    return meanr, gammatlist, weightlist, npairlist

def construct_jk_corr_from_njkcorr( dir = None, boost=False ):

    _dir = dir+'/jk/'
    meanr, gts, wt, npairs = load_jk_corr( _dir , kind='lens')
    meanr, grs, wr, npairs_r = load_jk_corr( _dir , kind='rand')

    gts = gts * wt
    grs = grs * wr

    #gts = np.sum( gts * wt, axis=1)
    #grs = np.sum( grs * wr, axis=1)
    #wt = np.sum( wt, axis = 1)
    #wr = np.sum( wr, axis = 1)

    gtnum = np.array([np.sum(np.delete(np.array(gts), i, axis = 0), axis = 0) for i in range(len(gts))])
    grnum = np.array([np.sum(np.delete(np.array(grs), i, axis = 0), axis = 0) for i in range(len(grs))])
    wnum = np.array([np.sum(np.delete(np.array(wt), i, axis = 0), axis = 0) for i in range(len(gts))])
    wrnum = np.array([np.sum(np.delete(np.array(wr), i, axis = 0), axis = 0) for i in range(len(gts))])

    boost_all = 0.0
    if boost : 
        boostdir = dir+'/boostfactor/'
        boost_all = np.genfromtxt(dir+'/boostfactor/construct_jk/jk.txt' ).T
        dir = dir+'/ggl_boost_corrected/'
        os.system('mkdir '+dir)

    gt_all = (boost_all + 1.)*gtnum/wnum - grnum/wrnum
    Njk = len(gt_all)
    gt_mean = np.mean(gt_all, axis = 0)

    gt_cov = np.zeros((20,20))
    for i in range(20):
        for j in range(20):
            gt_cov[i,j] = np.sum((gt_all[:,i]-gt_mean[i])*(gt_all[:,j]-gt_mean[j])) * (Njk-1.)*1./Njk

    #gt_cov = np.array([ np.sum((gt_all[:,i]-gt_mean[i] )**2) for i in range(20) ]) * (Njk-1.)*1./Njk

    err_gt = np.sqrt( gt_cov.diagonal() )
    DAT = np.column_stack(( gt_all ))
    DAT2 = np.column_stack(( meanr, gt_mean, err_gt))
    header = 'Njk = {}'.format(Njk)
    header2 = 'Njk = {}'.format(Njk) + '\n# meanr   jkmean  jkerr'
    os.system('mkdir '+dir+'construct_jk/')
    print 'save results in ', dir+'construct_jk/'
    np.savetxt( dir+'construct_jk/jk.txt', DAT, header=header )
    np.savetxt( dir+'construct_jk/jk_mean.txt', DAT2, header=header2 )
    np.savetxt( dir+'construct_jk/jk_cov.txt', gt_cov, header=header )



                                                     


class DeltaSigma():

    def __init__(self, h0=0.678596, omm = 0.3057076, omb = 0.04845864, n_s = 0.968376):

        from astropy.cosmology import FlatLambdaCDM
        # DES cosmology
        #h0 = 0.678596
        #omm = 0.3057076
        #omb = 0.04845864
        #n_s = 0.968376

        # Andres cosmology
        #h0 = 0.6726
        #omm = 0.314
        self.cosmo = FlatLambdaCDM(H0=100.*h0, Om0=omm, Ob0=omb)

        #nz_path = '/n/des/lee.5922/programs/cosmolike/DMASS-analysis/data/'
        #z_l, nz_l = np.loadtxt(nz_path + 'cmass_full.nz', unpack=True)
        #z_l, nz_l = np.loadtxt('/n/des/lee.5922/programs/cosmolike/cosmosis/dmass_cat/cmass_sgc.nz_2', unpack=True)
        #z_l += 0.005


    def sig_crit_inv(self, zl, zs):
        import astropy.units as u
        from astropy.constants import G, c

        if zl >= zs:
            return 0 * u.kg / (u.m * u.Mpc)
        else:
            Dl = self.cosmo.angular_diameter_distance(zl)
            Ds = self.cosmo.angular_diameter_distance(zs)
            Dls = self.cosmo.angular_diameter_distance_z1z2(zl,zs)
            return ((4.*np.pi*G)/(c**2.)) * ((Dls*Dl)/Ds)

    def sig_crit(self, zl, zs):
        import astropy.units as u
        from astropy.constants import G, c

        """
        returns sigma crit in units of kg / (m Mpc)
        """
        if zl >= zs:
            return 0 * u.kg / (u.m * u.Mpc)
        else:
            Dl = self.cosmo.angular_diameter_distance(zl)
            Ds = self.cosmo.angular_diameter_distance(zs)
            Dls = self.cosmo.angular_diameter_distance_z1z2(zl,zs)
            return ((c**2.)/(4.*np.pi*G)) * (Ds/(Dls*Dl))

    def sig_crit_inv_eff(self, z_l, nz_l, z_s, nz_s):
        import astropy.units as u
        from astropy.constants import G, c
        """
        integrates over two n(z) to get an effective sigma crit (as in Prat, Sanchez et al)
        """
        #normalize nz
        l_area = np.trapz(nz_l, x=z_l)
        nz_l = nz_l/l_area
        s_area = np.trapz(nz_s, x=z_s)
        nz_s = nz_s/s_area
        
        source_integral = np.zeros(len(z_l))
        for i_l in xrange(len(z_l)):  
            print '{}/{}        \r'.format(i_l, len(z_l)),
            sig_crit_inv_s = np.array([self.sig_crit_inv(z_l[i_l], z_s[i_s]).value for i_s in xrange(len(z_s))])
            source_integral[i_l] = np.trapz(sig_crit_inv_s*nz_s, x=z_s)
        sig_crit_inv_eff = np.trapz(nz_l*source_integral, x=z_l)
        return sig_crit_inv_eff*((u.m * u.Mpc)/u.kg )







