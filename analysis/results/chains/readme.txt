# description: chains of DMASS modified gravity analysis
#  
# if you use the chain files in this directory, cite these papers below:
# http://arxiv.org/abs/1906.01136 
# http://arxiv.org/abs/2104.11319 
# http://arxiv.org/abs/2104.14515
# 
# 
des+ext/

dwgsb_run2_chain.txt            : DMASS 3x2pt + BOSS
dwgsbp_run2_chain.txt           : DMASS 3x2pt + BOSS + Planck
dwgsbp_src34_run2_chain.txt     : DMASS 3x2pt (no 1st&2nd gammat signals) + BOSS + Planck 
dwgsbp_redmagic_run2_chain.txt  : redMaGiC 3x2pt + BOSS + Planck 


ext-only/

boss_run1_chain.txt             : BOSS only
boss_planck_chain.txt           : BOSS + Planck


systematics/
*This directory contains chains for robustness tests shown in Figure 5 in arxiv:2104.14515)

* DMASS 3x2pt + BOSS 
wgsb_base.txt                                     : baseline chain
is_weights_wgsb_run5_magnification_alpha3.24.txt  : Importance Sampling weights for magnification
is_weights_wgsb_run5_nolimberRSD.txt              : Importance Sampling weights for nolimberRSD
is_weights_wgsb_run5_nolimber.txt                 : Importance Sampling weights for nolimber
wgsb_iatatt_chain.txt                             : multinest chain for iatatt

* DMASS 3x2pt + BOSS + Planck
wgsbp_base.txt                                     : baseline chain 
is_weights_wgsbp_run2_magnification_alpha3.24.txt  : Importance Sampling weights for magnification
is_weights_wgsbp_run2_nolimberRSD.txt              : Importance Sampling weights for nolimberRSD  
is_weights_wgsbp_run2_nolimber.txt                 : Importance Sampling weights for nolimber
wgsbp_iatatt_chain.txt                             : multinest chain for iatatt 