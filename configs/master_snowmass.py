import os
from model_snowmass_base import config as _base_config
from svjHelper import masses_snowmass, scale_cms

# Start from a known-good config -> It wasn't filling at first, so I had to do this
config = _base_config

# Rebuild the standard snowmass mass/scale wiring so required fields exist
mpi = float(os.getenv("MPI", "20"))
mpi_over_scale = float(os.getenv("MPI_OVER_SCALE", "0.6"))
config = masses_snowmass(
    config=config,
    scale=mpi / mpi_over_scale,
    mpi_over_scale=mpi_over_scale
)

# Allow overriding of variables
config.channel  = os.getenv("CHANNEL", getattr(config, "channel", "s"))
config.mmed     = int(os.getenv("MMED", getattr(config, "mmed", 1000)))
config.Nc       = int(os.getenv("NC", getattr(config, "Nc", 3)))
config.Nf       = int(os.getenv("NF", getattr(config, "Nf", 3)))

config.mpi      = float(os.getenv("MPI", getattr(config, "mpi", 20.0)))
config.mrho     = float(os.getenv("MRHO", getattr(config, "mrho", config.mpi)))

config.scale    = float(os.getenv("SCALE", getattr(config, "scale", scale_cms(mpi=config.mpi))))
config.mq       = float(os.getenv("MQ", getattr(config, "mq", config.mpi / 2.0)))

config.pvector  = float(os.getenv("PVECTOR", getattr(config, "pvector", 0.5)))
config.rinv     = float(os.getenv("RINV", getattr(config, "rinv", 0.30)))

config.spectrum = os.getenv("SPECTRUM", getattr(config, "spectrum", "independent"))




print()
print("I am using the snowmass model")
print()