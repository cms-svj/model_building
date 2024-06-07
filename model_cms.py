from magiconfig import MagiConfig
import math

model = MagiConfig()
model.mmed = 1000
model.Nc = 2
model.Nf = 2
model.scale = 10
model.mpi = 20
model.scale = 3.2*math.pow(model.mpi,0.8)
model.mq = model.mpi/2.
model.pvector = 0.75
model.rinv = 0.3
model.spectrum = 'cms'
