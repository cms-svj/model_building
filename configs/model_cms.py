from magiconfig import MagiConfig
import math

config = MagiConfig()
config.channel = 's'
config.mmed = 1000
config.Nc = 2
config.Nf = 2
config.scale = 10
config.mpi = 20
config.mrho = config.mpi
config.scale = 3.2*math.pow(config.mpi,0.8)
config.mq = config.mpi/2.
config.pvector = 0.75
config.rinv = 0.3
config.spectrum = 'cms'
