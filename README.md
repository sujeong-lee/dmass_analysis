# DMASS analysis


## Relevant papers
Please cite relevant papers if you use the DMASS galaxy sample or any other products from the sample. 
* catalog: https://arxiv.org/abs/1906.01136
* galaxy-galaxy lensing and galaxy clustering measurement: https://arxiv.org/abs/2104.11319 
* modified gravity analysis: https://arxiv.org/abs/2104.14515

## DMASS and random catalogs

The DMASS selection algorithm can be found in 
* https://github.com/sujeong-lee/DMASS.git

The DMASS catalog and corresponding randoms are available at
* DMASS: [download](https://drive.google.com/uc?export=download&id=1OtOby6uA-hdR5CsMBBGvgLcHNmYAqh1R) 
* randoms: [download](https://drive.google.com/uc?export=download&id=1X2DNTdQyQLPlWlTkFCqYiLt3NRm5Eiu-)

1. Probability and systematics weights: The `CMASS_PROB` column contains a CMASS membership probability for each source. This column should be applied to sources as a weight. Systematics weight `WEIGHT` has been characterized in Section 4 in [the catalog paper](https://arxiv.org/abs/1906.01136). These two columns should be used as `weight = CMASS_PROB * WEIGHT` for any two-point computations.   
2. Excluding low probability galaxies: In the catalog paper, we remove galaxies lower than a probability threshold `CMASS_PROB < 0.01` to reduce noise and yield the same number density as CMASS. More details can be found in Section 3.5 in [the catalog paper](https://arxiv.org/abs/1906.01136).  We recommend users to apply the same cut but leave it as user's choice.      


## Measurements and chains
#### galaxy-galaxy lensing and cross-correlation coefficient: 

- datavector : `measurement/results/fits/measurement_y1_dmass_3x2pt.fits` 
- chains : `measurement/results/chains/`

#### modified gravity chains: 
- `analysis/chains/`




