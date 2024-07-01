export PYTHIA8=${MODEL_BUILDING}/install/pythia8
export PATH=${PYTHIA8}/bin/:${PATH}
export PYTHIA8DATA=$(pythia8-config --xmldoc)
export PYTHIA8RUNNER=${PYTHIA8}/examples/main42
