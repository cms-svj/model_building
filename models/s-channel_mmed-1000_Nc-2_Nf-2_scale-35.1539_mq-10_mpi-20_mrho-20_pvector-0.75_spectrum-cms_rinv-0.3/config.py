from magiconfig import MagiConfig

config = MagiConfig()
config.Nc = 2
config.Nf = 2
config.channel = 's'
config.delphes = 'cards/delphes_card_CMS.tcl'
config.dir = 'models'
config.events = 1000
config.mmed = 1000.0
config.mpi = 20.0
config.mq = 10.0
config.mrho = 20.0
config.pvector = 0.75
config.pythia = ['cards/CMS_Common.txt', 'cards/CMS_Tune_CP5.txt']
config.quiet = False
config.rinv = 0.3
config.scale = 35.15393738579578
config.spectrum = 'cms'
config.steps = ['all']
config.verbose = True