# PYTHIA configuration file.
# Generated on Mon Jul 1 02:49:17 PM CDT 2024 with the user supplied options:
# --with-hepmc2=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/HepMC/2.06.11/x86_64-el9-gcc12-opt
# --with-lhapdf6=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/lhapdf/6.5.3/x86_64-el9-gcc12-opt
# --with-python
# --with-gzip

# Install directory prefixes.
PREFIX_BIN=/uscms/home/jakiac/nobackup/PURSUE/CMSSW_14_1_0_pre3/src/model_building/install/pythia8/bin
PREFIX_INCLUDE=/uscms/home/jakiac/nobackup/PURSUE/CMSSW_14_1_0_pre3/src/model_building/install/pythia8/include
PREFIX_LIB=/uscms/home/jakiac/nobackup/PURSUE/CMSSW_14_1_0_pre3/src/model_building/install/pythia8/lib
PREFIX_SHARE=/uscms/home/jakiac/nobackup/PURSUE/CMSSW_14_1_0_pre3/src/model_building/install/pythia8/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/12.1.0-57c96/x86_64-centos9/bin/g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -DGZIP
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.11-d5a39/x86_64-el9-gcc12-opt/
HEPMC2_INCLUDE=-I/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.11-d5a39/x86_64-el9-gcc12-opt/include
HEPMC2_LIB=-L/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.11-d5a39/x86_64-el9-gcc12-opt/lib -Wl,-rpath,/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.11-d5a39/x86_64-el9-gcc12-opt/lib -lHepMC

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-el9-gcc12-opt/bin/
LHAPDF6_INCLUDE=-I/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-el9-gcc12-opt/include
LHAPDF6_LIB=-L/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-el9-gcc12-opt/lib -Wl,-rpath,/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-el9-gcc12-opt/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=true
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=-lz

PYTHON_USE=true
PYTHON_CONFIG=python-config
PYTHON_BIN=
PYTHON_INCLUDE=-I/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/include/python3.9 -I/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/include/python3.9
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
