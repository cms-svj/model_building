import os, sys, imp, re
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter
from multiprocessing import Pool

config_dir = os.path.join(os.getcwd(), "configs")
if config_dir not in sys.path:
    sys.path.append(config_dir)

def model_done(obj, outdir):
    # match config name (e.g. Nc7Nf7Ns1, rinv0p250) to model dir name parts (e.g. ..._Nc-7_Nf-7_..._Ns-1, ..._rinv-0.25)
    # simp config names round rinv to 3 decimals, so compare values at that precision
    parts = {key: float(val.replace('p','.')) for key,val in re.findall(r'([A-Za-z]+)([\dp]+)', obj)}
    def matches_dir(d):
        dparts = dict(p.split('-',1) for p in d.split('_') if '-' in p)
        return all(key in dparts and round(float(dparts[key]),3)==val for key,val in parts.items())
    matches = [d for d in os.listdir(outdir) if matches_dir(d)]
    return len(matches)==1 and os.path.exists(os.path.join(outdir, matches[0], "Hists.pkl"))

def run_objs(config_name, args, suffix, outdir, nproc, dryrun, force):
    os.makedirs(outdir, exist_ok=True)
    configs_fcdc = imp.load_source("configs_fcdc", config_name)
    objs = []
    for obj in vars(configs_fcdc.config):
        if not force and model_done(obj, outdir):
            print(f'skipping {obj}: Hists.pkl already exists (use --force to rerun)')
        else:
            objs.append(obj)

    cmds = []
    for obj in objs:
        logName = f'model_{obj}'
        if suffix: logName += f'_{suffix}'
        cmd = f'./run_model helper -C {config_name} -O config.{obj} --dir {outdir} {args} > {outdir}/{logName}.log 2>&1'
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
    parser.add_argument("--args", type=str, default="--steps all --events 1000 --verbose", help="arguments to pass to run_model")
    parser.add_argument("--suffix", type=str, default="", help="suffix for log file names")
    parser.add_argument("--outdir", type=str, default="models/fcdc/", help="outdir for logfile and events")
    parser.add_argument("--nproc", type=int, default=1, help="run in parallel (if > 1)")
    parser.add_argument("--dryrun", default=False, action="store_true", help="print commands but don't run")
    parser.add_argument("-f", "--force", default=False, action="store_true", help="rerun models even if output (Hists.pkl) already exists")
    args = parser.parse_args()
    run_objs(args.config_name, args.args, args.suffix, args.outdir, args.nproc, args.dryrun, args.force)
