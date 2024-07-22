from model_snowmass_base import config
from svjHelper import masses_snowmass

config = masses_snowmass(config=config,scale=50,mpi_over_scale=0.6)
