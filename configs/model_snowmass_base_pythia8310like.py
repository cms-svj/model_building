from magiconfig import MagiConfig
from svjHelper import gchi_lhcdm

config = MagiConfig()
config.channel = 's'
config.mmed = 1000
config.Nc = 3
config.Nf = 3
# StringFlav:mesonUDvector = 0.5 in Monash tune (and CP5 tune)
# StringFlav:mesonUDvector sets V/PS, while probVector sets V/(V+PS)
config.pvector = 0.333
# this rinv is only applied to diagonal pions; overall rinv value is (6+k)/9
k = 1
config.rinv = k/3
config.spectrum = 'snowmass'
config.gq = 0.25
config.gchi = gchi_lhcdm(gDM=1.0, Nc=config.Nc, Nf=config.Nf)
#
config.hv_setabsigma = 2
config.hv_alund = 0.3
config.hv_blund = 0.0078130
config.hv_sigmalund = 5.0595
config.hv_rfact = [1.,1.,1.,1.,1.,1.,1.,1.]
config.recoil_strategy_rf = 1
config.hv_stop_mass = 9.0          # 1.5 * m(pi_d)
config.hv_stop_new_flav = 2.0
config.hv_stop_smear = 0.2
config.hv_mjoin = 0.3
config.hv_mjoin_junction = 1.0
config.hv_mstring_min = 1.0
config.hv_constituent_mass = 0.325
config.hv_closed_m2max = 25.0
config.hv_ktolerance = 1.0

