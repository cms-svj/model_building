import os, sys, imp
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter, _rgetattr
from multiprocessing import Pool
from svjHelper import svjHelper
from common import config_path
config_path()

def model_done(config_name, obj, outdir):
    # get config name from helper
    helper = svjHelper.build(config_name, obj_name=obj)
    matches = [d for d in os.listdir(outdir) if d.endswith(helper.name())]
    return len(matches)==1 and os.path.exists(os.path.join(outdir, matches[0], "Hists.pkl"))

def run_objs(config_name, top_obj, args, suffix, outdir, nproc, dryrun, force):
    def fullname(obj):
        return f"{top_obj}.{obj}"

    os.makedirs(outdir, exist_ok=True)
    configs_fcdc = imp.load_source("configs_fcdc", config_name)
    objs = []
    for obj in vars(_rgetattr(configs_fcdc, top_obj)):
        if not force and model_done(config_name, fullname(obj), outdir):
            print(f'skipping {obj}: Hists.pkl already exists (use --force to rerun)')
        else:
            objs.append(obj)

    cmds = []
    for obj in objs:
        logName = f'model_{obj}'
        if suffix: logName += f'_{suffix}'
        cmd = f'./run_model helper -C {config_name} -O {fullname(obj)} --dir {outdir} {args} > {outdir}/{logName}.log 2>&1'
        cmds.append(cmd)

    for (obj,cmd) in zip(objs,cmds):
        print(obj)
        print(cmd)
        # serial processing
        if not dryrun and nproc<=1:
            os.system(cmd)

    # parallel processing
    if not dryrun and nproc>1:
        with Pool(nproc) as pool:
            exit_codes = pool.map(os.system, cmds)
        print('\n'.join(f'{i}: {j}' for i,j in zip(objs,exit_codes)))

if __name__=="__main__":
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsRawHelpFormatter
    )
    parser.add_argument("--config-name", type=str, default="configs/model_fcdc_10.py", help="config file containing multiple models")
    parser.add_argument("--obj", type=str, default="config", help="name of top-level config object")
    parser.add_argument("--args", type=str, default="--steps all --events 1000 --verbose", help="arguments to pass to run_model")
    parser.add_argument("--suffix", type=str, default="", help="suffix for log file names")
    parser.add_argument("--outdir", type=str, default="models/fcdc", help="outdir for logfile and events")
    parser.add_argument("--nproc", type=int, default=1, help="run in parallel (if > 1)")
    parser.add_argument("--dryrun", default=False, action="store_true", help="print commands but don't run")
    parser.add_argument("-f", "--force", default=False, action="store_true", help="rerun models even if output (Hists.pkl) already exists")
    args = parser.parse_args()
    run_objs(args.config_name, args.obj, args.args, args.suffix, args.outdir, args.nproc, args.dryrun, args.force)
