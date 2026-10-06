from magiconfig import MagiConfig
from model_fcdc_base import config as common
from svjHelper import masses_snowmass, fcdc_configs_NcNf1
from copy import deepcopy

# generate new configs for each input
config = MagiConfig()
for mpi in [3.1, 5.0, 10.0, 20.0]:
    for ratio in [0.1, 0.5, 1.0, 1.6, 2.0]:
        common = masses_snowmass(config=common, scale=mpi/ratio, mpi_over_scale=ratio)
        common2 = deepcopy(common)

        cfg = fcdc_configs_NcNf1(common=common)
        name = f'comp_mpi{mpi:.1f}_ratio{ratio:.1f}'.replace('.','p')
        setattr(config,name,cfg)

        common2.spectrum = 'fcdcSimp'
        cfg2 = fcdc_configs_NcNf1(common=common2)
        name2 = name.replace('comp','simp')
        setattr(config,name2,cfg2)
