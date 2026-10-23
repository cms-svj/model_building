import os
from copy import deepcopy

from model_cms import config as _base_config
from svjHelper import scale_cms

# Start from the standard CMS config
config = deepcopy(_base_config)

# Allow overriding of variables
config.channel = os.getenv("CHANNEL", getattr(config, "channel", "s"))

config.mmed = int(os.getenv("MMED", getattr(config, "mmed", 1000)))
config.Nc   = int(os.getenv("NC",   getattr(config, "Nc", 2)))
config.Nf   = int(os.getenv("NF",   getattr(config, "Nf", 2)))

config.mpi  = float(os.getenv("MPI",  getattr(config, "mpi", 20.0)))
config.mrho = float(os.getenv("MRHO", getattr(config, "mrho", config.mpi)))

config.scale = float(os.getenv("SCALE", getattr(config, "scale", scale_cms(mpi=config.mpi))))
config.mq    = float(os.getenv("MQ",    getattr(config, "mq", config.mpi / 2.0)))

config.pvector = float(os.getenv("PVECTOR", getattr(config, "pvector", 0.75)))
config.rinv    = float(os.getenv("RINV",    getattr(config, "rinv", 0.30)))

config.spectrum = os.getenv("SPECTRUM", getattr(config, "spectrum", "cms"))


print()
print("I am using my CMS-based master model")
print()

print("Configuration:")
print(f"  channel  = {config.channel}")
print(f"  mmed     = {config.mmed}")
print(f"  Nc       = {config.Nc}")
print(f"  Nf       = {config.Nf}")
print(f"  mpi      = {config.mpi}")
print(f"  mrho     = {config.mrho}")
print(f"  scale    = {config.scale}")
print(f"  mq       = {config.mq}")
print(f"  pvector  = {config.pvector}")
print(f"  rinv     = {config.rinv}")
print(f"  spectrum = {config.spectrum}")
print()
