from magiconfig import MagiConfig
from svjHelper import mqconst_snowmass, mrho_snowmass

config = MagiConfig()
config.channel = 's'
config.mmed = 1000
config.Nc = 3
config.Nf = 3
config.scale = 10
config.mpi = 0.6*config.scale
config.mq = mqconst_snowmass(mpi=config.mpi, scale=config.scale)
config.mrho = mrho_snowmass(mpi=config.mpi, scale=config.scale)
config.pvector = 0.5
# this rinv is only applied to diagonal pions; overall rinv value is (6+k)/9
k = 1
config.rinv = k/3
config.spectrum = 'snowmass'
