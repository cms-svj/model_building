from magiconfig import MagiConfig
import math

model = MagiConfig()
model.channel = 's'
model.mmed = 1000
model.Nc = 3
model.Nf = 3
model.scale = 10
model.mpi = 0.6*model.scale
model.mq = (model.mpi/5.5)**2/model.scale + model.scale
model.mrho = model.scale*math.sqrt(5.76+1.5*model.mpi**2/model.scale**2)
model.pvector = 0.5
k = 1
model.rinv = (6+k)/9
model.spectrum = 'snowmass'
