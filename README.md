# DMASS analysis


## Relevant papers
Please cite relevant papers if you use the DMASS galaxy sample or any other products from the sample. 
* catalog: https://arxiv.org/abs/1906.01136
* galaxy-galaxy lensing measurement: https://arxiv.org/abs/2104.11319 
* modified gravity analysis: https://arxiv.org/abs/2104.14515

## DMASS and random catalogs

The DMASS selection algorithm can be found in 
* https://github.com/sujeong-lee/DMASS.git

The DMASS catalog and corresponding randoms are available at
* DMASS: [download](https://drive.google.com/uc?export=download&id=1XABi761R4OLsWxmQZGW03t40jOWvgC_q) 
* randoms: [download](https://drive.google.com/uc?export=download&id=1mek4JB6PiKK0S0rpuUGnPlEuSnN2GZPD)

The DMASS catalog includes `CMASS_PROB`, the selection weights for the DES Y1 GOLD galaxies that weights GOLD galaxies to produce a statistical match to the BOSS CMASS sample. See Section 3.4 of [the catalog paper](https://arxiv.org/abs/1906.01136) for details. Systematics weight `WEIGHT_SYS` has been characterized in Section 4 of the paper. These weights should be applied for any two-point function-related computations.

To use the catalogs, follow the steps below:

1. Exclude contaminated regions and noisy sources as follows: 
   > for dmass  : `(VETO != 0) && (CMASS_PROB > 0.01)` <br>
   > for random : `(VETO != 0)` <br>
   
   `VETO != 0` keeps galaxies in the regions where the impact of survey properties (observing condition, depth) is negligible. The cut `CMASS_PROB > 0.01` removes low probability galaxies that are less likely to be a CMASS. The cut was determined to yield the same number density as CMASS. More details can be found in Section 3.5 of [the catalog paper](https://arxiv.org/abs/1906.01136). We recommend applying the same cut as stated here but leave it as the user's choice.    

2. Apply weights to the DMASS galaxies as follows: 
    > `weight = CMASS_PROB * WEIGHT_SYS`


## Measurements and chains
#### galaxy-galaxy lensing and cross-correlation coefficient: 
- All products of [the galaxy-galaxy lensing measurement paper](https://arxiv.org/abs/2104.14515) can be found in  `measurement` folder
- measurement of 3x2pt datavector (ggl + clustering + shear) : `measurement/results/fits/measurement_y1_dmass_3x2pt.fits` 
- 2pt-related outputs : `measurement/results/2pt/` 
- chains : `measurement/results/chains/`
- Please see [jupyter notebook](https://github.com/sujeong-lee/DMASS-analysis-publish/blob/master/notebook/DMASS-GGL%20results.ipynb) for relevant figures

#### modified gravity results: 
- All products of [the modified gravity paper](https://arxiv.org/abs/2104.14515) can be found in  `analysis` folder
- chains: `analysis/chains/`
- Please see [jupyter notebook](https://github.com/sujeong-lee/DMASS-analysis-publish/blob/master/notebook/DMASS-MG%20Results.ipynb) for relevant figures




