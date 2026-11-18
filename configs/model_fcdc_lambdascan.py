import copy
from magiconfig import MagiConfig
from model_fcdc_base import config as base
from svjHelper import masses_snowmass, fcdc_config, fcdc_config_simp, fcdc_rinv, gchi_lhcdm

mpi = 6
ratios = [0.3, 0.6, 1.0, 1.4, 1.6, 1.7, 2.0, 2.5]
Nc, Nf, Ns = 3, 3, 1

config = MagiConfig()
for ratio in ratios:
    tag = 'r{:g}'.format(ratio).replace('.', 'p')
    common = masses_snowmass(config=copy.deepcopy(base), scale=mpi/ratio, mpi_over_scale=ratio)
    _, full = fcdc_config(common=common, Nc=Nc, Nf=Nf, Ns=Ns)
    setattr(config, f'{tag}_full', full)
    simp_common = copy.deepcopy(common)
    simp_common.spectrum = 'fcdcSimp'
    simp_common.Nc, simp_common.Nf = Nc, Nf
    simp_common.gchi = gchi_lhcdm(gDM=1.0, Nc=Nc, Nf=Nf)
    _, simp = fcdc_config_simp(common=simp_common, rinv=fcdc_rinv(Nf=Nf, Ns=Ns))
    setattr(config, f'{tag}_simp', simp)
