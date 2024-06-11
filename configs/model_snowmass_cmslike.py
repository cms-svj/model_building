from magiconfig import MagiConfig
import math

config = MagiConfig()
config.channel = 's'
config.mmed = 1000
config.Nc = 3
config.Nf = 3
config.mpi = 20
config.scale = config.mpi/0.6
config.mq = (config.mpi/5.5)**2/config.scale + config.scale
config.mrho = config.scale*math.sqrt(5.76+1.5*config.mpi**2/config.scale**2)
config.pvector = 0.5
config.rinv = 0.3
config.spectrum = 'snowmass'