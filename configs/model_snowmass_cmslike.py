from magiconfig import MagiConfig
from svjHelper import mqconst_snowmass, mrho_snowmass

config = MagiConfig()
config.channel = 's'
config.mmed = 1000
config.Nc = 3
config.Nf = 3
config.mpi = 20
config.scale = config.mpi/0.6
config.mq = mqconst_snowmass(mpi=config.mpi, scale=config.scale)
config.mrho = mrho_snowmass(mpi=config.mpi, scale=config.scale)
config.pvector = 0.5
config.rinv = 0.3
config.spectrum = 'snowmass_cmslike'
