from magiconfig import MagiConfig

config = MagiConfig()
config.Nc = 3
config.Nf = 3
config.channel = 's'
config.delphes = 'cards/delphes_card_CMS.tcl'
config.dir = 'models'
config.events = 1000
config.mmed = 1000.0
config.mpi = 20.0
config.mq = 33.73002754820937
config.mrho = 83.66600265340756
config.pvector = 0.5
config.pythia = ['cards/CMS_Common.txt', 'cards/CMS_Tune_CP5.txt']
config.quiet = False
config.rinv = 0.3
config.scale = 33.333333333333336
config.spectrum = 'snowmass'
config.steps = ['all']
config.verbose = True