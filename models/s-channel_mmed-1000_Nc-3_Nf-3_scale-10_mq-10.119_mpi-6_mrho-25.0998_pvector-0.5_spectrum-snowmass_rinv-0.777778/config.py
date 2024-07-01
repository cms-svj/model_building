from magiconfig import MagiConfig

config = MagiConfig()
config.Nc = 3
config.Nf = 3
config.channel = 's'
config.delphes = 'cards/delphes_card_CMS.tcl'
config.dir = 'models'
config.events = 1000
config.mmed = 1000.0
config.mpi = 6.0
config.mq = 10.11900826446281
config.mrho = 25.099800796022265
config.pvector = 0.5
config.pythia = ['cards/CMS_Common.txt', 'cards/CMS_Tune_CP5.txt']
config.quiet = False
config.rinv = 0.7777777777777778
config.scale = 10.0
config.spectrum = 'snowmass'
config.steps = ['all']
config.verbose = True