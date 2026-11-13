import os
import hist
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mplhep as hep
import pickle
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter
from glob import glob
import itertools
from common import set_plot_style, resolve_models, accumulate_data

def all_samples(study):
    samples = []
    study_version = int(study[0])
    if study_version==1:
        samples = [
            {"name": "FCDC", "models": glob("models/fcdc/s-channel_mmed-1000_Nc-*_Nf-*_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-*_Ns-*")},
            {"name": "simp", "models": glob("models/fcdc/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdcSimp_gq-0.25_gchi-0.333333_rinv-*")},
            {"name": "FCDC (3-body)", "models": glob("models/fcdc/s-channel_mmed-1000_Nc-*_Nf-*_scale-3.52941_mq-3.8666_mpi-6_mrho-11.2139_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-*_Ns-*")},
            {"name": "simp (3-body)", "models": glob("models/fcdc/s-channel_mmed-1000_Nc-3_Nf-3_scale-3.52941_mq-3.8666_mpi-6_mrho-11.2139_pvector-0.333_spectrum-fcdcSimp_gq-0.25_gchi-0.333333_rinv-*")},
        ]
    elif study_version==2:
        model = study[1]
        samples = [
            {"name": "FCDC", "models": glob(f"models/fcdc_var/comp_{model}/*/")},
            {"name": "simp", "models": glob(f"models/fcdc_var/simp_{model}/*/")},
        ]
    elif study_version==3:
        # pythia version comparison: each version is compared to 8.310 ("ref") in shift/chain/scatter plots
        # color by version, marker by base version (open = default, filled = 8.310-like settings), linestyle by model
        base = "/store/user/easmith/svj/models/fcdc-pythia83{v}/s-channel_mmed-1000_Nc-*_Nf-*_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-{spec}_gq-0.25_gchi-*"
        versions = ["13", "13-like8310", "15", "15-like8310", "18", "18-like8310"]
        markers = {"13": "^", "15": "v", "18": "D"}
        for label, spec, linestyle in [("FCDC", "fcdc", "solid"), ("Simple", "fcdcSimp", "dotted")]:
            ref = f"{label} 8.310"
            samples.append({"name": ref, "models": resolve_models(base.format(v="10", spec=spec)), "style": dict(color="black", marker="o", linestyle=linestyle)})
            samples.extend([
                {"name": f"{label} 8.3{v}", "models": resolve_models(base.format(v=v, spec=spec)), "ref": ref,
                 "style": dict(color=style["color"], marker=markers[v[:2]], linestyle=linestyle, fillstyle="full" if "like" in v else "none")}
                for v, style in zip(versions, custom_cycler)
            ])
    else:
        raise ValueError(f"Unknown study version {study_version}")

    # account for two complete models with rinv=0.75
    for sample in samples:
        if "simp" in sample["name"]:
            idx = next((i for i, val in enumerate(sample["models"]) if 'rinv-0.75' in val), None)
            if idx is not None:
                sample["models"].insert(idx+1, sample["models"][idx])

    return samples

custom_cycler = set_plot_style()

def get_xvals(x, meta):
    # account for difference between complete and simplified models
    if "rinvpred" in x and x not in meta:
        return meta[x.replace("rinvpred","rinv")]
    return meta[x]

def hist_stats(h):
    # for quantities without analysis values (e.g. jet substructure), estimate them from the histogram
    vals, centers = h.values(), h.axes[0].centers
    mean = vals @ centers / vals.sum()
    stdev = np.sqrt(vals @ (centers-mean)**2 / vals.sum())
    return dict(mean=mean, stdev=stdev, stderr=stdev/np.sqrt(vals.sum()))

def process_data(data, samples, x, qname, forcex):
    processed = []
    for sample, models in data.items():
        xvals = np.array([get_xvals(x,model['meta']) for model in models])
        names = np.array([os.path.basename(os.path.dirname(model['file'])) for model in models])
        stats = [model['meta'].get(qname) or hist_stats(model['hist'][qname]) for model in models]
        means, stdevs, stderrs, hists = zip(*[
            (s['mean'], s['stdev'], s['stderr'], model['hist'].get(qname,None)) for s,model in zip(stats,models)
        ])
        # order by rinv_pred (then by model name, to match models across samples)
        order = np.lexsort((names, xvals))
        processed_dict = dict(
            sample = sample,
            ref = samples[sample].get('ref'),
            style = samples[sample].get('style'),
            names = names[order],
            xvals = xvals[order],
            means = np.array(means)[order],
            stdevs = np.array(stdevs)[order],
            stderrs = np.array(stderrs)[order],
            hists = [hists[i] for i in order],
        )
        hists_missing = [h is None for h in processed_dict['hists']]
        if all(hists_missing):
            processed_dict['hists'] = []
        elif any(hists_missing):
            raise RuntimeError(f"Missing histogram {qname} in:\n"+','.join([model[file] for i,model in enumerate(models) if hists_missing[i]]))
        processed.append(processed_dict)
    if forcex:
        if forcex not in data:
            raise ValueError(f"Unknown forcex sample {forcex}")
        forcexvals = next((pd['xvals'] for pd in processed if pd['sample']==forcex))
        for processed_dict in processed:
            processed_dict['xvals'] = forcexvals
    # compare to the same models in the reference sample
    by_sample = {pd['sample']: pd for pd in processed}
    for processed_dict in processed:
        if processed_dict['ref'] is None:
            continue
        ref = by_sample[processed_dict['ref']]
        if not np.array_equal(ref['names'], processed_dict['names']):
            raise RuntimeError(f"Models in {processed_dict['sample']} don't match reference {ref['sample']}")
        processed_dict['ref_means'] = ref['means']
        processed_dict['shifts'] = 100*(processed_dict['means']/ref['means'] - 1)
    return processed

# types comparing each sample to its reference sample (e.g. study 3: newer pythia versions vs 8.310)
ref_types = ['shift', 'chain', 'scatter']
# dark jet collections from dark partons to reco, for the chain plot
chain_collections = [("DP", "DP"), ("DH", "DH"), ("DH stable", "DHIV"), ("DH visible", "DHV"), ("Reco", "")]

# helper to make a plot
def make_plot(type, data, x, xlabel, qname, outdir, offset, plabel, tag=""):
    fig, ax = plt.subplots(figsize=(8,6))
    # iterator for manual control
    props = iter(custom_cycler)
    # advance iterator if requested
    next(itertools.islice(props, offset, offset), None)
    # use the sample's own style if it has one
    get_style = lambda entry: {'fillstyle': 'none', **(entry['style'] or next(props))}
    ylim = None
    ylabel = None
    if type=='stat':
        for entry in data:
            # extract label
            if ylabel is None and len(entry['hists'])>0:
                ylabel = entry['hists'][0].axes[0].label
            style = get_style(entry)
            # means
            line, = ax.plot(entry['xvals'], entry['means'], label=entry['sample'], **style)
            # stderr as errorbar
            ax.errorbar(entry['xvals'], entry['means'], yerr=entry['stderrs'], fmt='none', ecolor=style['color'], capsize=3, **style)
            # stdev as filled
            ax.fill_between(entry['xvals'], entry['means']-entry['stdevs'], entry['means']+entry['stdevs'], color=style['color'], alpha=0.15)
    elif type=='violin':
        # get global max hist height for each xval
        all_hists = []
        samples = []
        labels = []
        bin_edges = []
        heights = []
        centers = []
        max_per_label = []
        for entry in data:
            if len(entry['hists'])==0:
                continue

            # these are assumed to be shared across all samples
            if len(labels)==0:
                labels = entry['xvals']
                bin_edges = entry['hists'][0].axes[0].edges
                ylabel = entry['hists'][0].axes[0].label
                ylim = (bin_edges[0], bin_edges[-1])
                heights = np.diff(bin_edges)
                centers = bin_edges[:-1] + heights/2
            # extract values from histograms
            all_hists.append({x:h.values() for x,h in zip(entry['xvals'],entry['hists'])})
            samples.append(entry)

        # get global max hist height for each xval
        for label in labels:
            max_per_label.append(max([max(h[label]) for h in all_hists]))
        widths = np.array(max_per_label)
        x_locs = np.cumsum(widths) - 0.5 * widths

        for hists,entry in zip(all_hists, samples):
            style = get_style(entry)
            # plot all hists from this sample
            style.update({
                'fc': 'none',
                'ec': style['color'],
                'label': entry['sample'],
            })
            style.pop('marker')
            style.pop('fillstyle')
            for x_loc, label in zip(x_locs, labels):
                hist = hists[label]
                lefts = x_loc - 0.5 * hist
                ax.barh(centers, hist, height=heights, left=lefts, **style)
                # only apply label once
                style.pop('label', None)
            ax.set_xticks(x_locs, [f'{label:.3f}' for label in labels], rotation=45)
    elif type in ['shift', 'chain']:
        # median shift over models sharing an x value (chain: jet collection), 16-84% of models as filled
        for entry in data:
            style = get_style(entry)
            xs = np.unique(entry['xvals'])
            groups = [entry['shifts'][entry['xvals']==xv] for xv in xs]
            ax.plot(xs, [np.nanmedian(g) for g in groups], label=entry['sample'], **style)
            ax.fill_between(xs, [np.nanpercentile(g,16) for g in groups], [np.nanpercentile(g,84) for g in groups], color=style['color'], alpha=0.15)
        ax.axhline(0, color='gray', linewidth=1)
        ylabel = "Mean shift vs reference [%]"
        if type=='chain':
            ax.set_xticks(range(len(chain_collections)), [label for label,_ in chain_collections], rotation=30)
            xlabel = "Jet collection"
    elif type=='scatter':
        # per-model mean vs the same model in the reference sample
        for entry in data:
            ax.plot(entry['ref_means'], entry['means'], label=entry['sample'], **(get_style(entry) | dict(linestyle='none')))
        # y=x over the data range
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, color='black', linestyle=':')
        ax.set_xlim(lims)
        xlabel = f"Mean ({data[0]['ref']})"
        ylabel = "Mean"
    if type in ref_types:
        ax.set_title(next((entry['hists'][0].axes[0].label for entry in data if entry['hists']), qname))
    if type=='stat':
        ax.axline((0,0), slope=1, color='black', linestyle=':')
    if xlabel is None: xlabel = x
    ax.set_xlabel(xlabel)
    if ylim: ax.set_ylim(*ylim)
    if ylabel is None: ylabel = qname
    ax.set_ylabel(ylabel)
    handles, labels = ax.get_legend_handles_labels()
    if plabel is not None:
        handles.append(Rectangle((0,0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0, label=plabel))
    ax.legend(handles=handles, framealpha=0.5)
    plt.savefig(f'{outdir}/{type}{tag}_{qname}.pdf',bbox_inches='tight')
    plt.close(fig)

def make_all_plots(outdir, types, sample_list, study, x, y, xlabel, forcex, offset, plabel):
    samples = all_samples(study)

    unknown_samples = [s for s in sample_list if not any([s==sm['name'] for sm in samples])]
    if unknown_samples:
        raise ValueError("Unknown sample(s) requested:",','.join(unknown_samples))

    selected_samples = samples
    if sample_list:
        selected_samples = [sample for sample in samples if sample["name"] in sample_list]
    data = accumulate_data(selected_samples)

    samples_by_name = {sample["name"]: sample for sample in selected_samples}
    refs = list(dict.fromkeys(sample["ref"] for sample in selected_samples if sample.get("ref")))
    if not refs and any(t in ref_types for t in types):
        raise ValueError(f"{','.join(ref_types)} need samples with a reference (e.g. study 3)")

    os.makedirs(outdir, exist_ok=True)
    for qname in y:
        processed = process_data(data, samples_by_name, x, qname, forcex)
        for plot_type in types:
            if plot_type not in ref_types:
                make_plot(plot_type, processed, x, xlabel, qname, outdir, offset, plabel)
                continue
            entries = [entry for entry in processed if entry['ref']]
            if plot_type=='chain':
                if not qname.startswith("Jet12_"):
                    continue
                # combine shifts from each jet collection, using the collection index as x
                chain = [process_data(data, samples_by_name, x, f"{pre}{qname}", forcex) for _,pre in chain_collections]
                entries = [dict(entry, xvals=np.repeat(np.arange(len(chain)), len(entry['names'])), shifts=np.concatenate([c[i]['shifts'] for c in chain]))
                           for i,entry in enumerate(processed) if entry['ref']]
            # one plot per reference sample
            for ref in refs:
                make_plot(plot_type, [entry for entry in entries if entry['ref']==ref], x, xlabel, qname, outdir, offset, plabel, tag=f"_vs_{ref.replace(' ','_')}")

if __name__=="__main__":
    allowed_types = ['stat', 'violin'] + ref_types
    qtys_default = [
        'stable_invisible_fraction',
        'alpha_3body',
        'kappa',
        'DHIVJet12_rinv_proj',
        'AvgDHIVJet_rinv_proj',
        'DiDHIVJet_rinv_proj',
        'DHIVJet12_rinv_proj_global',
        'DHIVJet12_rinv_shape',
        'AvgDHIVJet_rinv_shape',
        'DiDHIVJet_rinv_shape',
        'DHIVJet12_rinv_shape_global',
    ]

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsRawHelpFormatter
    )
    parser.add_argument("--dir", type=str, default="All_metaplots", help="output directory")
    parser.add_argument("--types", type=str, default=['stat', 'violin'], nargs='*', choices=allowed_types, help="plot types (shift, chain, scatter need a reference sample, e.g. study 3)")
    parser.add_argument("--samples", type=str, default=[], nargs='*', help="list of samples to plot")
    parser.add_argument("--study", type=str, default=[], nargs='+', help="study version + misc arguments (to select samples)")
    parser.add_argument("-x", type=str, default='rinvpred', help="x variable")
    parser.add_argument("--xlabel", type=str, default=None, help="x axis label")
    parser.add_argument("--forcex", type=str, default=None, help="force use of x values from specified sample")
    parser.add_argument("-y", type=str, default=qtys_default, nargs='*', help="y variable(s)")
    parser.add_argument("--offset", type=int, default=0, help="offset for color/style cycler")
    parser.add_argument("--plabel", type=str, default=None, help="plot label (e.g. model info)")
    args = parser.parse_args()

    make_all_plots(args.dir, args.types, args.samples, args.study, args.x, args.y, args.xlabel, args.forcex, args.offset, args.plabel)
