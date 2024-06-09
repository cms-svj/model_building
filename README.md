# Dark QCD Model Building

## Installation

```bash
git clone git@github.com:cms-svj/model_building
cd model_building
./install.sh
```

## Environment

```bash
source init.sh
```

## Predefined model configurations

A few predefined model configuration files are provided in the [configs](./configs) folder:
1. The CMS model used in [EXO-19-020](https://arxiv.org/abs/2112.11125) (largely based on [arXiv:1503.00009](https://www.arxiv.org/abs/1503.00009) and [arXiv:1707.05326](https://arxiv.org/abs/1707.05326)).
2. One of the Snowmass benchmark models from [arXiv:2203.09503](https://arxiv.org/abs/2203.09503).

## Running

An example run command is:
```bash
./run_model -C configs/model_cms.py --run 10 --verbose
```

This command will, in order:
1. load the CMS model configuration
2. generate the Pythia and Delphes cards for this configuration
3. run Pythia and Delphes

If the `--run` command is omitted, step 3 will be skipped.

The input model configuration can be modified using command-line arguments, and the resulting configuration file will be generated along with the Pythia and Delphes cards.

## Analysis

The notebook [example_MT.ipynb](./example_MT.ipynb) shows an example of analyzing the generated events using [coffea](https://github.com/CoffeaTeam/coffea/):
* applying selections
* computing new quantities
* filling histograms
* making plots

