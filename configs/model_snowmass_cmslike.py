from model_snowmass_base import config
from svjHelper import masses_snowmass

mpi = 20
mpi_over_scale = 0.6
config = masses_snowmass(config=config,scale=mpi/mpi_over_scale,mpi_over_scale=mpi_over_scale)

config.rinv = 0.3
config.spectrum = 'snowmass_cmslike'
