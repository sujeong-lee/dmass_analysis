import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer


label_dict = {'cosmological_parameters--w0_fld':  r'w_{GDM}',
              'cosmological_parameters--cs2_fld': r'c_s^2',
              'cosmological_parameters--log_cs2': r'log(c_s^2)',
              'cosmological_parameters--omega_m': r'\Omega_m',
              'cosmological_parameters--omega_c': r'\Omega_c',
              'cosmological_parameters--ommh2': r'\Omega_m h^2',
              'cosmological_parameters--ombh2': r'\Omega_b h^2',
              'cosmological_parameters--omch2': r'\Omega_c h^2',
              'cosmological_parameters--h0':      r'h',
              'cosmological_parameters--omega_b': r'\Omega_b',
              'cosmological_parameters--n_s':     r'n_s',
              'cosmological_parameters--a_s':     r'A_s',
              'cosmological_parameters--w':     r'w_0',
              'cosmological_parameters--wa':     r'w_a',
              'cosmological_parameters--omnuh2':  r'\Omega_{\nu}',
              'cosmological_parameters--tau':  r'\tau',
              'cosmological_parameters--sigma_8': r'\sigma_8',
              'cosmological_parameters--s8': r'S_8',
	          'modified_gravity--sigma0':  r'\Sigma_0',
	          'modified_gravity--mu0':	  r'\mu_0',
              'modified_gravity--e_g':	  r'E_G',
	      'modified_gravity--m_sigma': r'm_{\Sigma}',
	      'modified_gravity--m_mu': r'm_{\mu}',
		'modified_gravity--m_mu2': r'm_{\mu}^2',
		'modified_gravity--m_sigma2': r'm_{\Sigma}^2',
		'modified_gravity--mu0_m_mu2': r'\mu_0 m_{\mu}^2',
		'modified_gravity--mu0_m_mu1': r'\mu_0 m_{\mu}',
		'modified_gravity--sigma0_m_sigma2': r'\Sigma_0 m_{\Sigma}^2',
		'modified_gravity--sigma0_m_sigma1': r'\Sigma_0 m_{\Sigma}',
              'intrinsic_alignment_parameters--a': r'A_{IA}',
              'intrinsic_alignment_parameters--alpha': r'\alpha_{IA}',
              'bin_bias--b1': r'b_1',
              #'bin_bias--b1': r'b_{\gamma}',
              'bin_bias--b2': 'b_2',
              'bin_bias--b3': 'b_3',
              'bin_bias--b4': 'b_4',
              'bin_bias--b5': 'b_5',
              'bias_shift--b1': r'{\Delta b}',
              'boss--b1': 'b_g',
              'boss--rcc' : r'r_{cc}',
              'shear_calibration_parameters--m1': 'm_1',
              'shear_calibration_parameters--m2': 'm_2',
              'shear_calibration_parameters--m3': 'm_3',
              'shear_calibration_parameters--m4': 'm_4',
              'shear_calibration_parameters--m5': 'm_5',
              'lens_photoz_errors--bias_1': 'z^l_1',
              'lens_photoz_errors--bias_2': 'z^l_2',
              'lens_photoz_errors--bias_3': 'z^l_3',
              'lens_photoz_errors--bias_4': 'z^l_4',
              'lens_photoz_errors--bias_5': 'z^l_5',
              'wl_photoz_errors--bias_1': 'z^s_1',
              'wl_photoz_errors--bias_2': 'z^s_2',
              'wl_photoz_errors--bias_3': 'z^s_3',
              'wl_photoz_errors--bias_4': 'z^s_4',
              'wl_photoz_errors--bias_5': 'z^s_5',
              }

    
get_label = np.vectorize(lambda label: label_dict[label] if label in label_dict.keys() else label)
    

def add_S8(data):
    data['cosmological_parameters--s8'] = data['cosmological_parameters--sigma_8']*(data['cosmological_parameters--omega_m']/0.3)**0.5
    return data

def add_rcc(data):
    if ('boss--b1' in data.keys()) & ('bin_bias--b1' in data.keys()):
        data['boss--rcc'] = data['bin_bias--b1']/data['boss--b1']
    else : pass
    return data

def add_shift( chains, shift=10 ):
    params = chains.keys()
    if 'prior' in params: params.remove('prior')
    if 'like' in params: params.remove('like')
    if 'post' in params: params.remove('post')
    if 'weight' in params: params.remove('weight')
    
    chains_copy = chains.copy()
    ppp =0.02
    for p in params:
        if p in ['modified_gravity--sigma0', 'modified_gravity--mu0']: sft = 0.05 * np.random.randint(-1*shift, shift)
        else: sft = ppp * np.random.randint(-1*shift, shift)
        values = chains[p]*(1. + sft)
        chains_copy[p] = values
    return chains_copy
    
def loading_cosmosis_chain(chainfile_name=None, burn=0, S8=True, joint=False, blind=False):
    with open(chainfile_name) as f:
        labels = f.readline()[1:-1].lower().split()

    loading_data = np.loadtxt(chainfile_name)
    if 'cosmological_parameters--s8' in labels:
        argidx = labels.index('cosmological_parameters--s8')
        #labels = labels[:argidx]
        labels.remove( 'cosmological_parameters--s8')
        nanmask = np.isnan(loading_data[0,:])
        loading_data = loading_data[:,~nanmask]

    nanmask = np.isnan(loading_data).any(axis=1)
    print np.sum(nanmask)
    data = {labels[i]: l[~nanmask] for i, l in enumerate(loading_data[burn:,:].T)}

    if S8: data = add_S8(data)
    if joint : data = add_rcc(data)
    if blind : data = add_shift(data)
    return data


def chain_concatenate(dir=None):
    import glob
    chainfiles=glob.glob(dir+'/*')
    chainfile_concatenate = dir+'chain_concatenate.txt'
    if chainfile_concatenate in chainfiles:
        chainfiles.remove(chainfile_concatenate)
    os.system('cp '+chainfiles[-1]+' '+chainfile_concatenate)
    
    f=open(chainfile_concatenate, "a")
    for i, fn in enumerate(chainfiles[1:]):
        data = np.genfromtxt(fn)
        for j, dd in enumerate(data):
            text = '\t'.join( str(d) for d in dd ) + '\n'
            f.write(text)
    f.close()
    print 'Output saved to ', chainfile_concatenate
    


def loading_IS_chain( baseline_filename = None, is_weight = None, burn = 0):
    
    base_data = loading_cosmosis_chain(baseline_filename, burn=burn)
    data = base_data.copy()
    _,_,newlike, weights = np.genfromtxt(is_weight, unpack=True)
    #weights = weights.reshape(weights.size, 1)
    data['weight'] = weights
    
    return data

def loading_cosmosis_fisher_chain(chainfile_name=None, retmean = False):
    with open(chainfile_name) as f:
        labels = np.array(f.readline()[1:-1].lower().split())
    with open(chainfile_name) as f:
        foutputs=f.readlines()
        
    means = []
    for i, fo in enumerate(foutputs):
        #if '## END_OF_PRIORS_INI' in fo: 
        if fo.startswith('#mu'):
            mean_i = float(fo.split('=')[-1])
            means.append(mean_i)
        else : pass
    means = np.array(means)
    
    # Cholesky decomposition inversion:
    Fisher = np.loadtxt(chainfile_name)
    c = np.linalg.inv(np.linalg.cholesky(Fisher))
    cov = np.dot(c.T,c)

    #cov = np.linalg.inv(np.loadtxt(chainfile_name))

    data = np.random.multivariate_normal(means, cov, size=100000)
    chain = {labels[i]: data[:,i] for i in range(data[0].size) }
    data = None
    #if S8: data = add_S8(data)
    #if joint : data = add_rcc(data)
    #if blind : data = add_shift(data)
    if retmean: 
        return chain, means
    else: return chain


def loading_cosmosis_fisher_chain_selected_params(chainfile_name=None, retmean = False, params=None, size=10000):
    with open(chainfile_name) as f:
        labels = np.array(f.readline()[1:-1].lower().split())
    with open(chainfile_name) as f:
        foutputs=f.readlines()
        
    means = []
    for i, fo in enumerate(foutputs):
        #if '## END_OF_PRIORS_INI' in fo: 
        if fo.startswith('#mu'):
            mean_i = float(fo.split('=')[-1])
            means.append(mean_i)
        else : pass
    means = np.array(means)
    
    cov = np.linalg.inv(np.loadtxt(chainfile_name))
    
    Nparams = len(params)
    argidx = [list(labels).index(params[i]) for i in range(Nparams)]
    
    cov_sub = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            cov_sub[i,j] = cov[argidx[i],argidx[j]]
    
    means_twoparams = [ means[argidx[i]] for i in range(Nparams)]
    data = np.random.multivariate_normal( means_twoparams, cov_sub, size=size)
    chain = {params[i]: data[:,i] for i in range(data[0].size) }
    data = None

    if retmean: return chain, means_twoparams
    else : return chain


def compute_FoM(chainfile_name=None, retmean = False, params=None):
    with open(chainfile_name) as f:
        labels = np.array(f.readline()[1:-1].lower().split())
    with open(chainfile_name) as f:
        foutputs=f.readlines()
        
    means = []
    for i, fo in enumerate(foutputs):
        #if '## END_OF_PRIORS_INI' in fo: 
        if fo.startswith('#mu'):
            mean_i = float(fo.split('=')[-1])
            means.append(mean_i)
        else : pass
    means = np.array(means)
    
    cov = np.linalg.inv(np.loadtxt(chainfile_name))
    
    
    argidx_w0 = list(labels).index(params[0])
    argidx_wa = list(labels).index(params[1])
    
    cov_w0wa = np.zeros((2,2))
    cov_w0wa[0,0] = cov[argidx_w0,argidx_w0]
    cov_w0wa[0,1] = cov[argidx_w0,argidx_wa]
    cov_w0wa[1,0] = cov[argidx_w0,argidx_wa]
    cov_w0wa[1,1] = cov[argidx_wa,argidx_wa]
    
    FoM1 = FoM1 = 1./np.sqrt( np.linalg.det(cov_w0wa))
    #FoM2 = np.sqrt(np.linalg.det( np.linalg.inv(cov_w0wa)) )
    print 'FoM : ', FoM1

def compute_FoM_from_chains(chains, params, kde = False ):
    
    c = ChainConsumer()
    
    chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
    #if zorder is None : zorder = [ None for i in range(len(chains)) ]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params)], name=chain_names[name_i] )
        cov = c.analysis.get_covariance(chain=name_i)
        FoM1 = 1./np.sqrt( np.linalg.det(cov[1] ))
        #FoM2 = np.sqrt(np.linalg.det(np.linalg.inv(cov[1])))
        print 'FoM (', chain_names[name_i],'):', FoM1
    

    
def _loading_cosmosis_fisher_chain(chainfile_name=None, mean = [0.3, 0.0, 0.0, 0.0, 2.1e-09, 0.0]):
    with open(chainfile_name) as f:
        labels = np.array(f.readline()[1:-1].lower().split())
    cov = np.linalg.inv(np.loadtxt(chainfile_name))
    data = np.random.multivariate_normal(mean, cov, size=1000000)
    chain = {labels[i]: data[:,i] for i in range(data[0].size) }
    data = None
    #if S8: data = add_S8(data)
    #if joint : data = add_rcc(data)
    #if blind : data = add_shift(data)    
    return chain    
    
    
def on_params(arr, params2plot):
    return np.array([arr[l] for l in params2plot]).T



def plotting_contours(chains, params2plot, truth = None, figname='test.png', 
    chain_names=None, plot_hists=True, 
    blind= None, shade_alpha = None, figsize=(7,7), extents=None,  
    shade = None, colors = None, linestyles = None, linewidths=None, 
    kde = False, 
    legend_location=None, 
    legend_color_text=True,
    legend_kwargs={'loc':'best'},
    flip=False, sigmas=[0,1,2], 
    zorder = None
    ):
    
    c = ChainConsumer()
    
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
    if zorder is None : zorder = [ None for i in range(len(chains)) ]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params2plot)], name=chain_names[name_i] )
    
    
    c.configure(
            #statistics=statistics,#[:len(chains)],
            linestyles=linestyles,#[:len(chains)], 
            label_font_size=20, tick_font_size =20,
            #summary=summary, #[: len(chains)], 
            sigmas=[0,1,2], 
            shade=shade, #[: len(chains)], \
            colors=colors, #[: len(chains)], \
            shade_alpha = shade_alpha, #[: len(chains)], 
            linewidths=linewidths, #[:len(chains)],
            #kde=[1.0,1.0,True,True,True,True][:len(chains)],
            kde = kde,
            plot_hists=plot_hists,
            flip=flip,
            legend_location=legend_location,
            legend_color_text = legend_color_text,
            legend_kwargs = legend_kwargs)
    
    fig = c.plotter.plot(filename=figname, figsize=figsize, truth=truth, 
                         blind=blind, extents=extents)

    return c

def plotting_1dhisto(chains, params2plot, truth = None, figname='test.png', 
    chain_names=None, 
    blind= None, figsize=(7,7), extents=None,  
    shade = None, colors = None, linestyles = None, linewidths=None, 
    kde = False, 
    legend_location=None, 
    legend_color_text=True,
    legend_kwargs={'loc':'best'},
    flip=False, sigmas=[0,1,2], 
    zorder = None
    ):
    
    c = ChainConsumer()
    
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
    if zorder is None : zorder = [ None for i in range(len(chains)) ]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params2plot)], name=chain_names[name_i] )
    
    
    c.configure(
            #statistics=statistics,#[:len(chains)],
            linestyles=linestyles,#[:len(chains)], 
            label_font_size=20, tick_font_size =20,
            summary=True, #[: len(chains)], 
            #shade=shade, #[: len(chains)], \
            colors=colors, #[: len(chains)], \
            #shade_alpha = shade_alpha, #[: len(chains)], 
            linewidths=linewidths, #[:len(chains)],
            kde = kde,
            legend_location=legend_location,
            legend_color_text = legend_color_text,
            legend_kwargs = legend_kwargs)
    
    fig = c.plotter.plot(filename=figname, figsize=figsize, truth=truth, 
                         blind=blind, extents=extents)

    return c    

def plotting_errorbars(chains, params2plot, truth=None, include_truth_chain=False, figname='test.png',
    extra_parameter_spacing=1.3, vertical_spacing_ratio = 1.5, blind=None, colors = None, extents=None,
    chain_names=None, figsize=1.5, kde=None, errorbar=True, shadecolor="#FB8C00" ):
    """
    keep : list. choose in the second chain
    params_fid : fiducial value of chain2
    """
    c = ChainConsumer()
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
        
    if colors is None : colors = ['black']
    if kde is None:  kde = []
    if truth is None : truth=chain_names[0]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params2plot)], name=chain_names[name_i])
        if colors is None : colors.append('black')
        if kde is None: kde.append(1)
         
        
    c.configure(label_font_size=15, tick_font_size =30,
                colors=colors, kde = kde)

    #c.configure_truth(ls=":", color="#FB8C00")
    c.configure_truth(ls=":", color=shadecolor)
    c.plotter.plot_summary(filename = figname, errorbar=errorbar, blind=blind,
                       truth=truth, 
                       include_truth_chain=include_truth_chain,
                       figsize=figsize, extents=extents, 
                       extra_parameter_spacing=extra_parameter_spacing, 
                          vertical_spacing_ratio = vertical_spacing_ratio)
    print "plot save to ", figname
    return c

def _plotting_errorbars(chains, params2plot, truth=None, figname='test.png',
    kde = None, 
    extra_parameter_spacing=1.3,
    chain_names=None, figsize=1.5,
           
                      
                      ):
    """
    keep : list. choose in the second chain
    params_fid : fiducial value of chain2
    """
    c = ChainConsumer()
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
        
    colors = ['red']
    if kde is None: kde = []
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params2plot)], name=chain_names[name_i])
        colors.append('black')
        if kde is None: kde.append(1)
         
    c.configure(label_font_size=20, tick_font_size =20,
                colors=colors, kde = kde)

    c.configure_truth(ls=":", color="#FB8C00")
    c.plotter.plot_summary(filename = figname, errorbar=True,
                       truth=chain_names[0], include_truth_chain=True,figsize=figsize, 
                       extra_parameter_spacing=extra_parameter_spacing)
    print "plot save to ", figname
    return c



def plotting_chainwalkers(chains, params2plot, truth = None, 
                          figname=None, chain_names=None, figsize=(7,7)):
    
    c = ChainConsumer()
    
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=['$'+l+'$' for l in get_label(params2plot)], name=chain_names[name_i] )
    
    fig = c.plotter.plot_walks(filename=figname, truth=truth, convolve=100)
    return c


def calculate_summary_statistics(chains, params2plot, truth=None, kde=None, blind=False,
    chain_names=None):

    c = ChainConsumer()
    if chain_names is None : chain_names = ['chain'+str(i+1) for i in range(len(chains)) ]
        
    colors = ['red']
    if kde is None: kde = []
    parameters = ['$'+l+'$' for l in get_label(params2plot)]
    for name_i, data in enumerate(chains) :
        c.add_chain(on_params(data, params2plot), weights=data['weight'] if 'weight' in data.keys() else None,
            parameters=parameters, name=chain_names[name_i])
        colors.append('black')
        if kde is None: kde.append(1)
        
    c.configure(bar_shade=False, kde=kde)
    #c.plotter.plot_summary()
    
    table = c.analysis.get_summary( )
    for tind, tab in enumerate(table):
        print chain_names[tind], ':'
        #for kind in keep:
        for p in parameters:
            values = tab[p]
            #print values
            mean = values[1]
            low = mean-values[0]
            upp = values[2]-mean
            #print low, mean, upp
	    if blind:
		mean=0.0
		low=0.0
		upp=0.0
		print 'blinded  -- ', p , ':', '{:0.3f}'.format(mean)+'^{+'+'{:0.3f}'.format(upp)+'}_{-'+'{:0.3f}'.format(low)+'}'
	    else : 
            	print '  -- ', p , ':', '{:0.2f}'.format(mean)+'^{+'+'{:0.2f}'.format(upp)+'}_{-'+'{:0.2f}'.format(low)+'}'
            
    latex_table = c.analysis.get_latex_table(caption="Results for the tested models", label="tab:example")
    
    return latex_table


def plotting_measurement( rootdir, blind=True, type='gammat', ylabel='', yscale='log', yaxhline=None):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    for i in range(4):
        savedir = rootdir+'/source'+str(i+1)+'/'
        meanr, gt, err_gt=\
        np.genfromtxt(savedir+'/finalize/meanjk_'+type+'.txt', unpack=True)
        ax.errorbar(meanr, gt, yerr=err_gt, fmt='o', label=labels[i], capsize=5 )

    if yaxhline is not None:
        ax.axhline(y=yaxhline, ls='--', color='grey')

    ax.axvspan(1, 27.0, alpha=0.2, color='grey')
    ax.set_xscale('log')
    ax.set_yscale(yscale)
    ax.set_xlabel('$\\theta$ [arcmin]')
    ax.set_ylabel(ylabel)
    ax.set_xlim(2.3, 300)
    if blind: 
        ax.set_yticklabels([])
        ax.tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
        ax.tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
    ax.legend()
    os.system('mkdir '+rootdir+'/fig/')
    plt.tight_layout()
    #fig.savefig(rootdir+'/fig/'+type+'.pdf')
    print 'fig saved to', rootdir+'/fig/'+type+'.pdf'


def plotting_measurement_split( rootdir, blind=True, types='gammat', input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False,
yextent = None, data_labels = None, mask_ang=[27,27,27,27], prefactor=1, figname=None):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,4, figsize=(16,3.5) )
    ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    #mask_mg = [28, 28, 35,35]
    #mask_lcdm = [27, 27, 27,27]
    
    #mask_ang = mask_mg
    #mask_ang = mask_lcdm
    fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = ['tab:blue', 'tomato', None]
        
    chi2_null_4_tot = np.zeros(4)
    chi2_null_12_tot = np.zeros(4)	
    for i in range(4):

        for ti, type in enumerate(types):
    
            savedir = rootdir+'/source'+str(i+1)+'/'
            meanr, gt, err_gt=\
            np.genfromtxt(savedir+'/finalize/meanjk_'+type+'.txt', unpack=True)

	    covgt = np.genfromtxt(savedir+'/finalize/covjk_'+type+'.txt')
	    chi2_null_4_tot[i] = calculate_chi2_null( meanr, gt, covgt, angcut=9, blind=False )
	    chi2_null_12_tot[i] = calculate_chi2_null(meanr, gt, covgt, angcut=27, blind=False )
	    #chi2_null_4 = calculate_chi2_null( meanr, gt, None, err=err_gt, angcut=9, blind=False )
            #chi2_null_12 = calculate_chi2_null(meanr, gt, None, err=err_gt, angcut=27, blind=False )
	    
            if thetagamma : 
                gt = gt*meanr
                err_gt = err_gt * meanr

            if input_err is not None :
                err_gt = input_err[meanr.size*i: meanr.size*(i+1)]

            dlabel=None
            if data_labels != None:
                #print data_labels[ti]
                if i == 0 : dlabel= data_labels[ti]
            #if ti != 0 : label=None
            ax[i].errorbar(meanr*(1 + 0.05*ti), prefactor*gt, yerr=prefactor*err_gt, fmt=fmts[ti], markerfacecolor=mfc[ti],  color=lc[ti], capsize=5, label=dlabel)

	    

        if yaxhline is not None:
            ax[i].axhline(y=yaxhline, ls='--', color='k', lw=1)

        ax[i].text(0.95, 0.95, labels[i], transform=ax[i].transAxes,fontsize=17, horizontalalignment='right', verticalalignment='top', backgroundcolor='w' )
        ax[i].axvspan(1, mask_ang[i], alpha=0.2, color='grey')
        ax[i].axvline(x=25, alpha=0.3, color='k', ls='--', lw=3)
        #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
        ax[i].set_xscale('log')
        ax[i].set_yscale(yscale, fontsize=17)
        ax[i].set_xlabel(xlabel,fontsize=17)
        ax[i].set_xlim(2.3, 300)
        ax[i].tick_params(labelsize=15)

        if yextent is not None : ax[i].set_ylim(yextent)
        #ax[i].set_ylabel(ylabel)

        if i != 0: ax[i].set_yticklabels([])
        if blind: 
            ax[i].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
            ax[i].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
            ax[i].set_yticklabels([])
        

        #ax[i].legend(loc=4)
    ax[0].set_ylabel(ylabel,fontsize=17)
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0, wspace=0)
    
     # labels along the bottom edge are off

    chi2_null_4_tot = np.sum(chi2_null_4_tot)
    chi2_null_12_tot = np.sum(chi2_null_12_tot)
    print '4 Mpc/h: chi2={:0.2f}/{}'.format(chi2_null_4_tot, 14*4)
    print '12Mpc/h: chi2={:0.2f}/{}'.format(chi2_null_12_tot, 10*4)

    os.system('mkdir '+rootdir+'/fig/')
    if figname is None : 
        fig.savefig(rootdir+'/fig/'+type+'_split.pdf')
        print 'fig saved to', rootdir+'/fig/'+type+'_split.pdf'
    else : 
        fig.savefig(figname)
        print 'fig saved to ', figname
    



def plotting_measurement_from_fits( fitsname=None, figname=None, blind=True, types=['gammat'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    
    fits = fitsio.FITS(fitsname)
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    #data = fits[4]['VALUE'].read()
    #ang = fits[4]['ANG'].read()
    #bin_no = fits[4]['BIN2'].read()
    
    """
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    mask2d = mask[d1] * mask[d2]

    cov_masked = cov[mask2d].reshape(Nma, Nma)
    data_masked = data[mask]
    """
    #fig, ax = plt.subplots()
    #ax.imshow( np.log10(cov_masked)  )
    
    #covinv = np.linalg.inv(cov_masked)
    #chi = np.sqrt(  np.dot( np.dot( data_masked, covinv), data_masked.T)  )
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory['galaxy_shear_xi']['VALUE'].read()
        #theory_masked = theory[mask]
        ang_theory = fits_theory['galaxy_shear_xi']['ANG'].read()
        bin_no_theory = fits_theory['galaxy_shear_xi']['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,4, figsize=(18,4) )
    ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            #mask_i = (ang > angcut) & (bin_no == bin_i+1)
            mask_i = (bin_no == bin_i+1)
            Nma = np.sum(mask_i)

            d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
            mask2d_i = mask_i[d1] * mask_i[d2]
            cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
            err_masked_i = np.sqrt(cov_masked_i.diagonal())
            data_masked_i = data[mask_i]
            ang_masked_i = ang[mask_i]

            #covinv = np.linalg.inv(cov_masked_i)
            #chi_i = np.sqrt(  np.dot( np.dot( data_masked_i, covinv), data_masked_i.T)  )

            if thetagamma :
                data_masked_i = data_masked_i * ang_masked_i
                err_masked_i = err_masked_i * ang_masked_i

            if input_err is not None :
                err_masked_i = input_err[mask_i] 

            
            if theory is not None:
                mask_theory_i = (bin_no_theory == bin_i+1)
                theory_masked_i = theory[mask_theory_i]
                ang_theory_masked_i = ang_theory[mask_theory_i]
                #diff_dv = theory_masked_i - data_masked_i
                #diff_chi_i = np.sqrt(  np.dot( np.dot( diff_dv, covinv), diff_dv.T)  )
                #print diff_chi_i
               

            dlabel=None
            if data_labels != None:
                if bin_i == 0 : dlabel= data_labels[bin_i]

                    
            if theory is not None:
                ax[bin_i].plot(ang_theory_masked_i, theory_masked_i, ls='-', lw=1, color='black')
                
            ax[bin_i].errorbar(ang_masked_i*(1 + 0.05*ti), data_masked_i, yerr=err_masked_i, fmt=fmts[ti], markerfacecolor=mfc[ti],  color=lc[ti], capsize=5, label=dlabel)

                
            if yaxhline is not None:
                ax[bin_i].axhline(y=yaxhline, ls='--', color='grey')

            
            ax[bin_i].text(0.95, 0.95, labels[bin_i], transform=ax[bin_i].transAxes,fontsize=17, 
                           horizontalalignment='right', verticalalignment='top', backgroundcolor='w' )
            ax[bin_i].axvspan(1, angcut, alpha=0.2, color='grey')
            #ax[bin_i].axvline(x=angcut, alpha=0.3, color='k', ls='-', lw=3)
            ax[bin_i].axvline(x=25, alpha=0.3, color='k', ls='--', lw=3)
            #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
            ax[bin_i].set_xscale('log')
            ax[bin_i].set_yscale(yscale, fontsize=17)
            ax[bin_i].set_xlabel(xlabel,fontsize=17)
            ax[bin_i].set_xlim(2.3, 290)
            ax[bin_i].tick_params(labelsize=15)

            if yextent is not None : ax[bin_i].set_ylim(yextent)

            #if bin_i != 0: ax[bin_i].set_yticklabels([])
            if blind: 
                ax[bin_i].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
                ax[bin_i].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
                ax[bin_i].set_yticklabels([])


            ax[bin_i].legend(loc=4)
            ax[0].set_ylabel(ylabel,fontsize=17)
            
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0, wspace=0.15)

    if figname is not None: fig.savefig(figname)
        


def plotting_measurement_from_fits_in_onepanel( fitsname=None, figname=None, blind=True, types=['gammat'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    fits = fitsio.FITS(fitsname)
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    #data = fits[4]['VALUE'].read()
    #ang = fits[4]['ANG'].read()
    #bin_no = fits[4]['BIN2'].read()
    
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    mask2d = mask[d1] * mask[d2]

    cov_masked = cov[mask2d].reshape(Nma, Nma)
    data_masked = data[mask]
    
    #fig, ax = plt.subplots()
    #ax.imshow( np.log10(cov_masked)  )
    
    #covinv = np.linalg.inv(cov_masked)
    #chi = np.sqrt(  np.dot( np.dot( data_masked, covinv), data_masked.T)  )
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory['galaxy_shear_xi']['VALUE'].read()
        #theory_masked = theory[mask]
        ang_theory = fits_theory['galaxy_shear_xi']['ANG'].read()
        bin_no_theory = fits_theory['galaxy_shear_xi']['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    #fig, ax = plt.subplots(1,4, figsize=(18,4) )
    fig, ax = plt.subplots(figsize=(6,5) )

    #ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'o', 'd', 'd']
    mfc = [None, 'w', None, 'w']
    lc = ['steelblue', 'orangered', 'orange', 'black', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            #mask_i = (ang > angcut) & (bin_no == bin_i+1)
            mask_i = (bin_no == bin_i+1)
            Nma = np.sum(mask_i)

            d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
            mask2d_i = mask_i[d1] * mask_i[d2]
            cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
            err_masked_i = np.sqrt(cov_masked_i.diagonal())
            data_masked_i = data[mask_i]
            ang_masked_i = ang[mask_i]

            #covinv = np.linalg.inv(cov_masked_i)
            #chi_i = np.sqrt(  np.dot( np.dot( data_masked_i, covinv), data_masked_i.T)  )

            if thetagamma :
                data_masked_i = data_masked_i * ang_masked_i
                err_masked_i = err_masked_i * ang_masked_i

            if input_err is not None :
                err_masked_i = input_err[mask_i] 

            
            if theory is not None:
                mask_theory_i = (bin_no_theory == bin_i+1)
                theory_masked_i = theory[mask_theory_i]
                ang_theory_masked_i = ang_theory[mask_theory_i]
                #diff_dv = theory_masked_i - data_masked_i
                #diff_chi_i = np.sqrt(  np.dot( np.dot( diff_dv, covinv), diff_dv.T)  )
                #print diff_chi_i
               

            dlabel=None
            if data_labels != None:
                if bin_i == 0 : dlabel= data_labels[bin_i]

                    
            if theory is not None:
                ax.plot(ang_theory_masked_i, theory_masked_i, ls='-', lw=1, color=lc[bin_i])
                
            ax.errorbar(ang_masked_i*(0.99 + 0.02*(-1)**bin_i), data_masked_i, yerr=err_masked_i, fmt=fmts[bin_i], 
            markerfacecolor=mfc[bin_i],  color=lc[bin_i], capsize=5, label=labels[bin_i])

                
        if yaxhline is not None:
            ax.axhline(y=yaxhline, ls='--', color='grey')

        
        #ax.text(0.95, 0.95, labels[bin_i], transform=ax.transAxes,fontsize=17, 
        #               horizontalalignment='right', verticalalignment='top' )
        #ax.axvspan(1, angcut, alpha=0.2, color='grey')
        #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
        ax.set_xscale('log')
        ax.set_yscale(yscale, fontsize=17)
        ax.set_xlabel(xlabel,fontsize=17)
        ax.set_xlim(2.3, 290)
        ax.set_ylim(None, 8e-03)
        ax.tick_params(labelsize=15)

        if yextent is not None : ax.set_ylim(yextent)

        #if bin_i != 0: ax[bin_i].set_yticklabels([])
        if blind: 
            ax.tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
            ax.tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
            ax.set_yticklabels([])
    
        ax.axvspan(1, angcut, alpha=0.2, color='grey')
        #ax.axvspan(1, 27, alpha=0.15, color='grey')
        #ax.axvline(x=angcut, color='k', ls='-', lw=3, alpha=0.3  )
        ax.axvline(x=25, color='k', ls='--', lw=3, alpha=0.3 )
        #ax.legend(loc=4)
        ax.legend(loc=1, fontsize=15, frameon=False)
        ax.set_ylabel(ylabel,fontsize=17)
            
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0, wspace=0.15)

    if figname is not None: fig.savefig(figname)

def chi2_from_fits( fitsname=None, figname=None, blind=True, types=['gammat'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    fits = fitsio.FITS(fitsname)
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    #data = fits[4]['VALUE'].read()
    #ang = fits[4]['ANG'].read()
    #bin_no = fits[4]['BIN2'].read()
    
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    mask2d = mask[d1] * mask[d2]

    cov_masked = cov[mask2d].reshape(Nma, Nma)
    data_masked = data[mask]
    
    #fig, ax = plt.subplots()
    #ax.imshow( np.log10(cov_masked)  )
    
    covinv = np.linalg.inv(cov_masked)
    chi = np.sqrt(  np.dot( np.dot( data_masked, covinv), data_masked.T)  )
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory['galaxy_shear_xi']['VALUE'].read()
        theory_masked = theory[mask]
        ang_theory = fits_theory['galaxy_shear_xi']['ANG'].read()
        bin_no_theory = fits_theory['galaxy_shear_xi']['BIN2'].read()
    
    diff_dv = theory_masked - data_masked
    #diff_chi = np.sqrt(  np.dot( np.dot( diff_dv, covinv), diff_dv.T)  )
    diff_chi = np.sum(diff_dv**2/cov_masked.diagonal())
    print 'total:', '{:0.1f}/{:0.0f}, {:0.1f}'.format(diff_chi, diff_dv.size, diff_chi*1./diff_dv.size)
    
    import matplotlib.pyplot as plt
    #fig, ax = plt.subplots(1,4, figsize=(18,4) )
    fig, ax = plt.subplots(figsize=(6,5) )

    #ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'o', 'd', 'd']
    mfc = [None, 'w', None, 'w']
    lc = ['steelblue', 'orangered', 'orange', 'black', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            mask_i = (ang > angcut) & (bin_no == bin_i+1)
            #mask_i = (bin_no == bin_i+1)
            Nma = np.sum(mask_i)

            d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
            mask2d_i = mask_i[d1] * mask_i[d2]
            cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
            err_masked_i = np.sqrt(cov_masked_i.diagonal())
            data_masked_i = data[mask_i]
            ang_masked_i = ang[mask_i]

            covinv = np.linalg.inv(cov_masked_i)
            #chi_i = np.sqrt(  np.dot( np.dot( data_masked_i, covinv), data_masked_i.T)  )

            if thetagamma :
                data_masked_i = data_masked_i * ang_masked_i
                err_masked_i = err_masked_i * ang_masked_i

            if input_err is not None :
                err_masked_i = input_err[mask_i] 

            
            if theory is not None:
                mask_theory_i = (ang > angcut) & (bin_no_theory == bin_i+1)
                theory_masked_i = theory[mask_theory_i]
                ang_theory_masked_i = ang_theory[mask_theory_i]
                diff_dv = theory_masked_i - data_masked_i
                #diff_chi_i = np.sqrt(  np.dot( np.dot( diff_dv, covinv), diff_dv.T)  )
                diff_chi_i = np.sum(diff_dv**2/cov_masked_i.diagonal())
                print labels[bin_i], '{:0.1f}/{:0.0f} {:0.2f}'.format(diff_chi_i, diff_dv.size, diff_chi_i*1./diff_dv.size)
               

            dlabel=None
            if data_labels != None:
                if bin_i == 0 : dlabel= data_labels[bin_i]

                    
            if theory is not None:
                ax.plot(ang_theory_masked_i, theory_masked_i, ls='-', lw=1, color=lc[bin_i])
                
            ax.errorbar(ang_masked_i*(0.99 + 0.02*(-1)**bin_i), data_masked_i, yerr=err_masked_i, fmt=fmts[bin_i], 
            markerfacecolor=mfc[bin_i],  color=lc[bin_i], capsize=5, label=labels[bin_i])

                
        if yaxhline is not None:
            ax.axhline(y=yaxhline, ls='--', color='grey')

        
        #ax.text(0.95, 0.95, labels[bin_i], transform=ax.transAxes,fontsize=17, 
        #               horizontalalignment='right', verticalalignment='top' )
        #ax.axvspan(1, angcut, alpha=0.2, color='grey')
        #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
        ax.set_xscale('log')
        ax.set_yscale(yscale, fontsize=17)
        ax.set_xlabel(xlabel,fontsize=17)
        ax.set_xlim(2.3, 290)
        ax.set_ylim(None, 8e-03)
        ax.tick_params(labelsize=15)

        if yextent is not None : ax.set_ylim(yextent)

        #if bin_i != 0: ax[bin_i].set_yticklabels([])
        if blind: 
            ax.tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
            ax.tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
            ax.set_yticklabels([])

        ax.axvspan(1, angcut, alpha=0.2, color='grey')
        #ax.axvspan(1, 27, alpha=0.15, color='grey')
        #ax.axvline(x=angcut, color='k', ls='-', lw=3, alpha=0.3  )
        ax.axvline(x=25, color='k', ls='--', lw=3, alpha=0.3 )
        #ax.legend(loc=4)
        ax.legend(loc=1, fontsize=15, frameon=False)
        ax.set_ylabel(ylabel,fontsize=17)
            
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0, wspace=0.15)

    if figname is not None: fig.savefig(figname)  

def plotting_measurement_split_diff( rootdir, blind=True, types='gammat', input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False,
yextent = None, data_labels = None, mask_ang=[27,27,27,27], figname=None):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,4, figsize=(16,3.5) )
    ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'

    fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        

    
    for i in range(4):

        #for ti, type in enumerate(types):
        ti = 0

        savedir = rootdir+'/source'+str(i+1)+'/'
        meanr, gt, err_gt=\
        np.genfromtxt(savedir+'/finalize/meanjk_'+types[0]+'.txt', unpack=True)
        meanr2, gt2, err_gt2=\
        np.genfromtxt(savedir+'/finalize/meanjk_'+types[1]+'.txt', unpack=True)
	    
        """
        if thetagamma : 
            gt = gt*meanr
            err_gt = err_gt * meanr

        if input_err is not None :
            err_gt = input_err[meanr.size*i: meanr.size*(i+1)]
        
        dlabel=None
        if data_labels != None:
            #print data_labels[ti]
            if i == 0 : dlabel= data_labels[ti]
        #if ti != 0 : label=None
        """
        
        #fracdiff = 1. - gt/gt2
        #fracerr = np.sqrt(2)*err_gt/gt
        diff = gt-gt2
        differr = err_gt
        ax[i].errorbar(meanr, diff, yerr=differr, fmt=fmts[ti], markerfacecolor=mfc[ti],  color=lc[ti], capsize=5)

        if yaxhline is not None:
            ax[i].axhline(y=yaxhline, ls='--', color='grey')

        ax[i].text(0.95, 0.95, labels[i], transform=ax[i].transAxes,fontsize=17, horizontalalignment='right', verticalalignment='top' )
        ax[i].axvspan(1, mask_ang[i], alpha=0.2, color='grey')
        #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
        ax[i].set_xscale('log')
        ax[i].set_yscale(yscale, fontsize=17)
        ax[i].set_xlabel(xlabel,fontsize=17)
        ax[i].set_xlim(2.3, 300)
        ax[i].tick_params(labelsize=15)

        if yextent is not None : ax[i].set_ylim(yextent)
        #ax[i].set_ylabel(ylabel)

        if i != 0: ax[i].set_yticklabels([])
        if blind: 
            ax[i].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
            ax[i].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
            ax[i].set_yticklabels([])
        

        ax[i].legend(loc=4)
    ax[0].set_ylabel(ylabel,fontsize=17)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0)
    
     # labels along the bottom edge are off

    os.system('mkdir '+rootdir+'/fig/')
    if figname is None : 
        fig.savefig(rootdir+'/fig/'+type+'_split.pdf')
        print 'fig saved to', rootdir+'/fig/'+type+'_split.pdf'
    else : 
        fig.savefig(figname)
        print 'fig saved to ', figname


def plotting_gammat_measurement_from_fits( fitsname=None, figname=None, blind=True, types=['gammat'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    fits = fitsio.FITS(fitsname)
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory['galaxy_shear_xi']['VALUE'].read()
        #theory_masked = theory[mask]
        ang_theory = fits_theory['galaxy_shear_xi']['ANG'].read()
        bin_no_theory = fits_theory['galaxy_shear_xi']['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,4, figsize=(10,2.8) )
    ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            #mask_i = (ang > angcut) & (bin_no == bin_i+1)
            mask_i = (bin_no == bin_i+1)
            Nma = np.sum(mask_i)

            d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
            mask2d_i = mask_i[d1] * mask_i[d2]
            cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
            err_masked_i = np.sqrt(cov_masked_i.diagonal())
            data_masked_i = data[mask_i]
            ang_masked_i = ang[mask_i]


            if thetagamma :
                data_masked_i = data_masked_i * ang_masked_i
                err_masked_i = err_masked_i * ang_masked_i

            if input_err is not None :
                err_masked_i = input_err[mask_i] 

            
            if theory is not None:
                mask_theory_i = (bin_no_theory == bin_i+1)
                theory_masked_i = theory[mask_theory_i]
                ang_theory_masked_i = ang_theory[mask_theory_i]
               

            dlabel=None
            if data_labels != None:
                if bin_i == 0 : dlabel= data_labels[bin_i]

                    
            if theory is not None:
                ax[bin_i].plot(ang_theory_masked_i, theory_masked_i, ls='--', color='black')
                
            ax[bin_i].errorbar(ang_masked_i*(1 + 0.05*ti), data_masked_i, yerr=err_masked_i, fmt=fmts[ti], markerfacecolor=mfc[ti],  color=lc[ti], capsize=5, label=dlabel)

                
            if yaxhline is not None:
                ax[bin_i].axhline(y=yaxhline, ls='--', color='grey')

            
            ax[bin_i].text(0.95, 0.95, labels[bin_i], transform=ax[bin_i].transAxes,fontsize=15, 
                           horizontalalignment='right', verticalalignment='top' )
            ax[bin_i].axvspan(1, angcut[bin_i], alpha=0.3, color='grey')
            #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
            ax[bin_i].set_xscale('log')
            ax[bin_i].set_yscale(yscale, fontsize=17)
            ax[bin_i].set_xlabel(xlabel,fontsize=17)
            ax[bin_i].set_xlim(2.3, 290)
            ax[bin_i].set_ylim(2e-06, 5e-03)
            ax[bin_i].tick_params(labelsize=15)

            if yextent is not None : ax[bin_i].set_ylim(yextent)

            
            #if blind: 
            ax[bin_i].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
            ax[bin_i].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
            if bin_i != 0: ax[bin_i].set_yticklabels([])


            ax[bin_i].legend(loc=4)
            ax[0].set_ylabel(ylabel,fontsize=17)
            
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)

    if figname is not None: fig.savefig(figname)


def plotting_gammat_measurement_from_fits( fitsname=None, figname=None, blind=True, types=['gammat'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    fits = fitsio.FITS(fitsname)
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory['galaxy_shear_xi']['VALUE'].read()
        #theory_masked = theory[mask]
        ang_theory = fits_theory['galaxy_shear_xi']['ANG'].read()
        bin_no_theory = fits_theory['galaxy_shear_xi']['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,4, figsize=(10,2.8) )
    ax = ax.ravel()

    labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            #mask_i = (ang > angcut) & (bin_no == bin_i+1)
            mask_i = (bin_no == bin_i+1)
            Nma = np.sum(mask_i)

            d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
            mask2d_i = mask_i[d1] * mask_i[d2]
            cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
            err_masked_i = np.sqrt(cov_masked_i.diagonal())
            data_masked_i = data[mask_i]
            ang_masked_i = ang[mask_i]


            if thetagamma :
                data_masked_i = data_masked_i * ang_masked_i
                err_masked_i = err_masked_i * ang_masked_i

            if input_err is not None :
                err_masked_i = input_err[mask_i] 

            
            if theory is not None:
                mask_theory_i = (bin_no_theory == bin_i+1)
                theory_masked_i = theory[mask_theory_i]
                ang_theory_masked_i = ang_theory[mask_theory_i]
               

            dlabel=None
            if data_labels != None:
                if bin_i == 0 : dlabel= data_labels[bin_i]

                    
            if theory is not None:
                ax[bin_i].plot(ang_theory_masked_i, theory_masked_i, ls='--', lw = 1, color='black')
                
            ax[bin_i].errorbar(ang_masked_i, data_masked_i, yerr=err_masked_i, fmt=fmts[ti], markerfacecolor=mfc[ti],  color=lc[ti], capsize=5, label=dlabel)

                
            if yaxhline is not None:
                ax[bin_i].axhline(y=yaxhline, ls='--', color='grey')

            
            ax[bin_i].text(0.95, 0.95, labels[bin_i], transform=ax[bin_i].transAxes,fontsize=15, 
                           horizontalalignment='right', verticalalignment='top' )
            ax[bin_i].axvspan(1, angcut[bin_i], alpha=0.3, color='grey')
            #ax[i].axvline(x=meanr, color='green', lw=27, alpha=0.5)
            ax[bin_i].set_xscale('log')
            ax[bin_i].set_yscale(yscale, fontsize=17)
            ax[bin_i].set_xlabel(xlabel,fontsize=17)
            ax[bin_i].set_xlim(2.3, 290)
            ax[bin_i].set_ylim(1e-06, 5e-03)
            ax[bin_i].tick_params(labelsize=15, which='both', direction='in')

            if yextent is not None : ax[bin_i].set_ylim(yextent)

            
            if blind: 
                ax[bin_i].tick_params( which='both', axis='y', length = 0, width = 0, color = 'blue')

            if bin_i != 0: ax[bin_i].set_yticklabels([])


            ax[bin_i].legend(loc=4)
            ax[0].set_ylabel(ylabel,fontsize=17)
            
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)

    if figname is not None: fig.savefig(figname)



def plotting_xipm_measurement_from_fits( fitsname=None, figname=None, blind=True, types=['shear_xi_plus'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    section = 'shear_xi_plus'
    fits = fitsio.FITS(fitsname)
    data = fits[section]['VALUE'].read()
    ang = fits[section]['ANG'].read()
    bin_no1 = fits[section]['BIN1'].read()
    bin_no2 = fits[section]['BIN2'].read()
    cov = fits['COVMAT'].read()[:200, :200]
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory[section]['VALUE'].read()
        ang_theory = fits_theory[section]['ANG'].read()
        bin_no1_theory = fits_theory[section]['BIN1'].read()
        bin_no2_theory = fits_theory[section]['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(6,4, figsize=(12,9) )
    #ax = ax.ravel()

    #labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            for bin_j in np.arange(bin_i,4):
                
                #mask_i = (ang > angcut) & (bin_no == bin_i+1)
                mask_ij = (bin_no1 == bin_i+1) * (bin_no2 == bin_j+1)
                #mask_j = (bin_no == bin_j+1)
                Nma = np.sum(mask_ij)

                d1, d2 = np.mgrid[0:mask_ij.size, 0:mask_ij.size]
                mask2d_ij = mask_ij[d1] * mask_ij[d2]
                cov_masked_ij = cov[mask2d_ij].reshape(Nma, Nma)
                err_masked_ij = np.sqrt(cov_masked_ij.diagonal())
                data_masked_ij = data[mask_ij]
                ang_masked_ij = ang[mask_ij]

            
                if theory is not None:
                    mask_theory_ij = (bin_no1_theory == bin_i+1) * (bin_no2_theory == bin_j+1)
                    theory_masked_ij = theory[mask_theory_ij]
                    ang_theory_masked_ij = ang_theory[mask_theory_ij]


                #dlabel=None
                #if data_labels != None:
                #    if bin_i == 0 : dlabel= data_labels[bin_i]

                if theory is not None:
                    ax[5-bin_i][3-bin_j].plot(ang_theory_masked_ij, theory_masked_ij, ls='--', lw=1, color='black')

                ax[5-bin_i][3-bin_j].errorbar(ang_masked_ij, data_masked_ij, 
                                          yerr=err_masked_ij, fmt=fmts[ti], markerfacecolor=mfc[ti],  
                                          color=lc[ti], capsize=5, label=None)


                if yaxhline is not None:
                    ax[5-bin_i][3-bin_j].axhline(y=yaxhline, ls='--', color='grey')


                ax[5-bin_i][3-bin_j].text(0.95, 0.95, '$({},{})$'.format(bin_i+1, bin_j+1), 
                                      transform=ax[5-bin_i][3-bin_j].transAxes,fontsize=17, 
                                      horizontalalignment='right', verticalalignment='top' )
                ax[5-bin_i][3-bin_j].axvspan(1, angcut[5-bin_i][3-bin_j], alpha=0.3, color='grey')
                ax[5-bin_i][3-bin_j].set_xscale('log')
                ax[5-bin_i][3-bin_j].set_yscale(yscale, fontsize=17)
                ax[5-bin_i][3-bin_j].set_xlabel(xlabel,fontsize=17)
                ax[5-bin_i][3-bin_j].set_xlim(2.3, 290)
                ax[5-bin_i][3-bin_j].set_ylim(1e-07, 9e-05)
                ax[5-bin_i][3-bin_j].tick_params(labelsize=15, direction='in', which='both' )
                if bin_j != 3: ax[5-bin_i][3-bin_j].set_yticklabels([])
                ax[5-bin_i][3-bin_j].legend(loc=4)
                ax[5-bin_i][0].set_ylabel(ylabel,fontsize=17)

                #if yextent is not None : ax[bin_i+2][bin_j].set_ylim(yextent)


                if blind: 
                    ax[5-bin_i][3-bin_j].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
                    ax[5-bin_i][3-bin_j].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
                
    
    # --------------------------------------------
    # xi minus starts
    
    section = 'shear_xi_minus'
    ylabel='$\\xi_{-} (\\theta)$'

    #fits = fitsio.FITS(fitsname)
    data = fits[section]['VALUE'].read()
    ang = fits[section]['ANG'].read()
    bin_no1 = fits[section]['BIN1'].read()
    bin_no2 = fits[section]['BIN2'].read()
    cov = fits['COVMAT'].read()[200:400, 200:400]
    
    
    if theory is not None:
        #fits_theory = fitsio.FITS(theory)
        theory = fits_theory[section]['VALUE'].read()
        ang_theory = fits_theory[section]['ANG'].read()
        bin_no1_theory = fits_theory[section]['BIN1'].read()
        bin_no2_theory = fits_theory[section]['BIN2'].read()
    

    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(5):
            for bin_j in np.arange(bin_i, 4):
  
                mask_ij = (bin_no1 == bin_i+1) * (bin_no2 == bin_j+1)
                Nma = np.sum(mask_ij)

                d1, d2 = np.mgrid[0:mask_ij.size, 0:mask_ij.size]
                mask2d_ij = mask_ij[d1] * mask_ij[d2]
                cov_masked_ij = cov[mask2d_ij].reshape(Nma, Nma)
                err_masked_ij = np.sqrt(cov_masked_ij.diagonal())
                data_masked_ij = data[mask_ij]
                ang_masked_ij = ang[mask_ij]

            
                if theory is not None:
                    mask_theory_ij = (bin_no1_theory == bin_i+1) * (bin_no2_theory == bin_j+1)
                    theory_masked_ij = theory[mask_theory_ij]
                    ang_theory_masked_ij = ang_theory[mask_theory_ij]

                #dlabel=None
                #if data_labels != None:
                #    if bin_i == 0 : dlabel= data_labels[bin_i]

                #if theory is not None:
                ax[bin_i][bin_j].plot(ang_theory_masked_ij, theory_masked_ij, ls='--', lw=1, color='black')
                ax[bin_i][bin_j].errorbar(ang_masked_ij, data_masked_ij, 
                                          yerr=err_masked_ij, fmt=fmts[ti], markerfacecolor='w' ,  
                                          color=None, capsize=5, label=None)


                if yaxhline is not None:
                    ax[bin_i][bin_j].axhline(y=yaxhline, ls='--', color='grey')


                ax[bin_i][bin_j].text(0.95, 0.95, '$({},{})$'.format(bin_i+1, bin_j+1), 
                                      transform=ax[bin_i][bin_j].transAxes,fontsize=17, 
                                      horizontalalignment='right', verticalalignment='top' )
                ax[bin_i][bin_j].axvspan(1, angcut[bin_i][bin_j], alpha=0.3, color='grey')
                ax[bin_i][bin_j].set_xscale('log')

                
                ax[bin_i][bin_j].set_yscale(yscale, fontsize=17)
                #ax[bin_i][bin_j].set_xlabel(xlabel,fontsize=17)
                ax[bin_i][bin_j].set_xlim(2.3, 290)
                ax[bin_i][bin_j].set_ylim(1e-07, 8e-05)
                ax[bin_i][bin_j].tick_params(labelsize=15, which='both', direction='in', labelright=True, right=True )
                
                #if yextent is not None : ax[bin_i][bin_j].set_ylim(yextent)

                #if blind: 
                #ax[bin_i][bin_j].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
                #ax[bin_i][bin_j].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
                ax[bin_i][bin_j].set_xticklabels([])
                if bin_j != 3: ax[bin_i][bin_j].set_yticklabels([])


                ax[bin_i][bin_j].legend(loc=4)
                ax[bin_i][3].set_ylabel(ylabel,fontsize=17)
                ax[bin_i][bin_j].yaxis.set_ticks_position('right') 
                ax[bin_i][bin_j].yaxis.set_label_position('right') 

                if bin_j == (bin_i + 1) :
                    ax[bin_j][bin_i].remove()
                
    ax[4][3].remove()       
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)

    if figname is not None: fig.savefig(figname)




def plotting_xipm_measurement_from_fits2( fitsname=None, figname=None, blind=True, types=['shear_xi_plus'], input_err=None,
ylabel='', yscale='log', yaxhline=None, scalecut=True, thetagamma=False, theory=None, fmts = None, 
yextent = None, data_labels = None, angcut=1):
    
    
    import fitsio
    
    section = 'shear_xi_plus'
    fits = fitsio.FITS(fitsname)
    data = fits[section]['VALUE'].read()
    ang = fits[section]['ANG'].read()
    bin_no1 = fits[section]['BIN1'].read()
    bin_no2 = fits[section]['BIN2'].read()
    cov = fits['COVMAT'].read()[:200, :200]
    
    
    if theory is not None:
        fits_theory = fitsio.FITS(theory)
        theory = fits_theory[section]['VALUE'].read()
        ang_theory = fits_theory[section]['ANG'].read()
        bin_no1_theory = fits_theory[section]['BIN1'].read()
        bin_no2_theory = fits_theory[section]['BIN2'].read()
    
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(5,4, figsize=(12,8) )
    #ax = ax.ravel()

    #labels = [ '$0.20 < z_s < 0.43$' ,  '$0.43 < z_s < 0.63$' ,  '$0.63 < z_s < 0.90$' ,  '$0.90 < z_s < 1.30$'  ]
    xlabel = '$\\theta$ ${\\rm [arcmin]}$'
    if fmts is None: fmts = ['o', 'd', 's']
    mfc = [None, 'w', 'w']
    lc = [None, 'tomato', None]
        
    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            for bin_j in np.arange(bin_i,4):
                
                #mask_i = (ang > angcut) & (bin_no == bin_i+1)
                mask_ij = (bin_no1 == bin_i+1) * (bin_no2 == bin_j+1)
                #mask_j = (bin_no == bin_j+1)
                Nma = np.sum(mask_ij)

                d1, d2 = np.mgrid[0:mask_ij.size, 0:mask_ij.size]
                mask2d_ij = mask_ij[d1] * mask_ij[d2]
                cov_masked_ij = cov[mask2d_ij].reshape(Nma, Nma)
                err_masked_ij = np.sqrt(cov_masked_ij.diagonal())
                data_masked_ij = data[mask_ij]
                ang_masked_ij = ang[mask_ij]

            
                if theory is not None:
                    mask_theory_ij = (bin_no1_theory == bin_i+1) * (bin_no2_theory == bin_j+1)
                    theory_masked_ij = theory[mask_theory_ij]
                    ang_theory_masked_ij = ang_theory[mask_theory_ij]


                #dlabel=None
                #if data_labels != None:
                #    if bin_i == 0 : dlabel= data_labels[bin_i]

                #if theory is not None:
                ax[4-bin_i][3-bin_j].plot(ang_theory_masked_ij, theory_masked_ij, ls='--', lw=1, color='black')
                ax[4-bin_i][3-bin_j].errorbar(ang_masked_ij, data_masked_ij, 
                                          yerr=err_masked_ij, fmt=fmts[ti], markerfacecolor=mfc[ti],  
                                          color=lc[ti], capsize=5, label=None)


                if yaxhline is not None:
                    ax[4-bin_i][3-bin_j].axhline(y=yaxhline, ls='--', color='grey')


                ax[4-bin_i][3-bin_j].text(0.95, 0.95, '$({},{})$'.format(bin_i+1, bin_j+1), 
                                      transform=ax[4-bin_i][3-bin_j].transAxes,fontsize=17, 
                                      horizontalalignment='right', verticalalignment='top' )
                ax[4-bin_i][3-bin_j].axvspan(1, angcut[4-bin_i][3-bin_j], alpha=0.2, color='grey')
                ax[4-bin_i][3-bin_j].set_xscale('log')
                ax[4-bin_i][3-bin_j].set_yscale(yscale, fontsize=17)
                ax[4-bin_i][3-bin_j].set_xlabel(xlabel,fontsize=17)
                ax[4-bin_i][3-bin_j].set_xlim(2.3, 290)
                ax[4-bin_i][3-bin_j].set_ylim(1e-07, 9e-05)
                ax[4-bin_i][3-bin_j].tick_params(labelsize=15, direction='in', which='both' )
                if bin_j != 3: ax[4-bin_i][3-bin_j].set_yticklabels([])
                ax[4-bin_i][3-bin_j].legend(loc=4)
                ax[4-bin_i][0].set_ylabel(ylabel,fontsize=17)

                

                #if yextent is not None : ax[bin_i+2][bin_j].set_ylim(yextent)


                if blind: 
                    ax[4-bin_i][3-bin_j].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
                    ax[4-bin_i][3-bin_j].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
                

    # --------------------------------------------
    # xi minus starts
    
    section = 'shear_xi_minus'
    ylabel='$\\xi_{-} (\\theta)$'

    #fits = fitsio.FITS(fitsname)
    data = fits[section]['VALUE'].read()
    ang = fits[section]['ANG'].read()
    bin_no1 = fits[section]['BIN1'].read()
    bin_no2 = fits[section]['BIN2'].read()
    cov = fits['COVMAT'].read()[200:400, 200:400]
    
    
    if theory is not None:
        #fits_theory = fitsio.FITS(theory)
        theory = fits_theory[section]['VALUE'].read()
        ang_theory = fits_theory[section]['ANG'].read()
        bin_no1_theory = fits_theory[section]['BIN1'].read()
        bin_no2_theory = fits_theory[section]['BIN2'].read()
    

    for ti, type in enumerate(types):
        # For each bin
        for bin_i in range(4):
            for bin_j in np.arange(bin_i, 4):
  
                mask_ij = (bin_no1 == bin_i+1) * (bin_no2 == bin_j+1)
                Nma = np.sum(mask_ij)

                d1, d2 = np.mgrid[0:mask_ij.size, 0:mask_ij.size]
                mask2d_ij = mask_ij[d1] * mask_ij[d2]
                cov_masked_ij = cov[mask2d_ij].reshape(Nma, Nma)
                err_masked_ij = np.sqrt(cov_masked_ij.diagonal())
                data_masked_ij = data[mask_ij]
                ang_masked_ij = ang[mask_ij]

            
                if theory is not None:
                    mask_theory_ij = (bin_no1_theory == bin_i+1) * (bin_no2_theory == bin_j+1)
                    theory_masked_ij = theory[mask_theory_ij]
                    ang_theory_masked_ij = ang_theory[mask_theory_ij]

                #dlabel=None
                #if data_labels != None:
                #    if bin_i == 0 : dlabel= data_labels[bin_i]

                #if theory is not None:
                ax[bin_i][bin_j].plot(ang_theory_masked_ij, theory_masked_ij, ls='--', lw=1, color='black')
                ax[bin_i][bin_j].errorbar(ang_masked_ij, data_masked_ij, 
                                          yerr=err_masked_ij, fmt=fmts[ti], markerfacecolor='w' ,  
                                          color=None, capsize=5, label=None)


                if yaxhline is not None:
                    ax[bin_i][bin_j].axhline(y=yaxhline, ls='--', color='grey')


                ax[bin_i][bin_j].text(0.95, 0.95, '$({},{})$'.format(bin_i+1, bin_j+1), 
                                      transform=ax[bin_i][bin_j].transAxes,fontsize=17, 
                                      horizontalalignment='right', verticalalignment='top' )
                ax[bin_i][bin_j].axvspan(1, angcut[bin_i][bin_j], alpha=0.2, color='grey')
                ax[bin_i][bin_j].set_xscale('log')

                
                ax[bin_i][bin_j].set_yscale(yscale, fontsize=17)
                #ax[bin_i][bin_j].set_xlabel(xlabel,fontsize=17)
                ax[bin_i][bin_j].set_xlim(2.3, 290)
                ax[bin_i][bin_j].set_ylim(1e-07, 8e-05)
                ax[bin_i][bin_j].tick_params(labelsize=15, which='both', direction='in', labelright=True, right=True )
                
                #if yextent is not None : ax[bin_i][bin_j].set_ylim(yextent)

                #if blind: 
                #ax[bin_i][bin_j].tick_params( which='major', axis='y', length = 0, width = 0, color = 'blue')
                #ax[bin_i][bin_j].tick_params( which='minor', axis='y', length = 0, width = 0, color = 'blue')
                ax[bin_i][bin_j].set_xticklabels([])
                if bin_j != 3: ax[bin_i][bin_j].set_yticklabels([])


                ax[bin_i][bin_j].legend(loc=4)
                ax[bin_i][3].set_ylabel(ylabel,fontsize=17)
                ax[bin_i][bin_j].yaxis.set_ticks_position('right') 
                ax[bin_i][bin_j].yaxis.set_label_position('right') 

                #if bin_j == (bin_i + 1) :
                #    ax[bin_j][bin_i].spines['top'].set_linewidth(2)
                #    ax[bin_j][bin_i].spines['right'].set_linewidth(2)
                #    #ax[bin_j][bin_i].remove()


    for bin_i in range(4):
        for bin_j in np.arange(bin_i, 4):
            if bin_j == bin_i :
                ax[bin_i+1][bin_j].spines['top'].set_linewidth(3)
                ax[bin_i][bin_j].spines['left'].set_linewidth(3)

    ax[0][0].spines['left'].set_linewidth(1)          
    #ax[4][3].remove()       
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)

    if figname is not None: fig.savefig(figname)




def calculate_gammat_SNR( fitsname=None, angcut=9, blind=True ):
    
    import fitsio
    fits = fitsio.FITS(fitsname)
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    data = fits['galaxy_shear_xi']['VALUE'].read()
    ang = fits['galaxy_shear_xi']['ANG'].read()
    bin_no = fits['galaxy_shear_xi']['BIN2'].read()
    
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    mask2d = mask[d1] * mask[d2]

    cov_masked = cov[mask2d].reshape(Nma, Nma)
    data_masked = data[mask]
    
    #fig, ax = plt.subplots()
    #ax.imshow( np.log10(cov_masked)  )
    
    covinv = np.linalg.inv(cov_masked)
    chi = np.sqrt(  np.dot( np.dot( data_masked, covinv), data_masked.T)  )
    #chi = np.sum(data_masked**2/cov_masked.diagonal())
    
    print 'total:  ',
    if blind: print 'SNR=--, dof={}'.format(Nma) 
    else : print 'SNR={:0.2f}, dof={}'.format(chi,Nma)
        
    
    # For each bin
    for bin_i in range(1,5):
        mask_i = (ang > angcut) & (bin_no == bin_i)
        Nma = np.sum(mask_i)
    
        d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
        mask2d_i = mask_i[d1] * mask_i[d2]
        cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
        data_masked_i = data[mask_i]
        
        covinv = np.linalg.inv(cov_masked_i)
        chi_i = np.sqrt(  np.dot( np.dot( data_masked_i, covinv), data_masked_i.T)  )
        #chi_i = np.sum(data_masked_i**2/cov_masked_i.diagonal())
        print 'source{}:'.format(bin_i),
        if blind: print 'SNR=--, dof={}'.format(Nma) 
        else : print 'SNR={:0.2f}, dof={}'.format(chi_i,Nma)


def calculate_chi2_diff( fitsname=None, theory=None, section='galaxy_shear_xi', angcut=9, blind=True ):
    
    import fitsio
    fits = fitsio.FITS(fitsname)
    cov = fits['COVMAT'].read()[400:-20, 400:-20]
    if section == 'galaxy_xi': cov = fits['COVMAT'].read()[-20:, -20:]
    data = fits[section]['VALUE'].read()
    ang = fits[section]['ANG'].read()
    bin_no = fits[section]['BIN2'].read()
    
    
    fits_theory = fitsio.FITS(theory)
    data_theory = fits_theory[section]['VALUE'].read()    
    
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    #mask2d = mask[d1] * mask[d2]
    cov_masked = cov[:,mask][mask,:]
    #cov_masked = cov[mask2d].reshape(Nma, Nma)
    data_masked = data[mask]
    theory_masked = data_theory[mask]
    #fig, ax = plt.subplots()
    #ax.imshow( np.log10(cov_masked)  )
        
    covinv = np.linalg.inv( cov_masked ) 
    diff_dv = (data_masked - theory_masked)
    chi =  np.dot( np.dot( diff_dv, covinv), diff_dv.T)  
    #chi = np.sum(diff_dv**2 / cov_masked.diagonal())
    
    print 'total:  ',
    if blind: print 'chi2=--, dof={}'.format(Nma) 
    else : print 'chi2={:3.2f}/{}'.format(chi,Nma)
        
    
    # For each bin
    for bin_i in range(1,5):
        mask_i = (ang > angcut) & (bin_no == bin_i)
        Nma = np.sum(mask_i)
    
        #d1, d2 = np.mgrid[0:mask_i.size, 0:mask_i.size]
        #mask2d_i = mask_i[d1] * mask_i[d2]
        #cov_masked_i = cov[mask2d_i].reshape(Nma, Nma)
	cov_masked_i = cov[:,mask_i][mask_i,:]
        data_masked_i = data[mask_i]
        theory_masked_i = data_theory[mask_i]
        
        covinv = np.linalg.inv(cov_masked_i)
        diff_dv_i = (data_masked_i - theory_masked_i)
        chi_i =  np.dot( np.dot( diff_dv_i, covinv), diff_dv_i.T)  
        #chi_i = np.sum(diff_dv_i**2 / cov_masked_i.diagonal())
        
        print 'source{}:'.format(bin_i),
        if blind: print 'chi2=--, dof={}'.format(Nma) 
        else : print 'chi2={:3.2f}/{}'.format(chi_i,Nma)
            
            
def calculate_chi2_null( ang, data, cov, err=None, angcut=5, blind=True):
    
    mask = ang > angcut
    Nma = np.sum(mask)
    
    d1, d2 = np.mgrid[0:mask.size, 0:mask.size]
    mask2d =  mask[d1] * mask[d2]

    data_masked = data[mask]
    diff_datav = data_masked 
    #-1.
    
    if err is not None: 
        err_masked = err[mask]
        chi2 = np.sum(data_masked**2 * 1./err_masked**2)
    
    else : 
        cov_masked = cov[mask2d].reshape(Nma, Nma)
        covinv = np.linalg.inv(cov_masked)
        chi2 = np.dot( np.dot( data_masked, covinv), data_masked.T) 
        #chi2 = np.sum( diff_datav**2/cov_masked.diagonal() )
    
    if blind: print 'chi2=--, dof={}'.format(Nma) 
    else : print 'chi2={:0.2f}/{}'.format(chi2,Nma)

    return chi2
