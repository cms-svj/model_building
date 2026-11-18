import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import pickle
from magiconfig import ArgumentParser, ArgumentDefaultsRawHelpFormatter
from common import resolve_models

# --summary: samples with the same group (model configuration) are overlaid in the same panel; group: panel title
groups = {
    "Nc8Nf8Ns1": r"$N_c=N_f=8,\ N_s=1$ ($r_{\rm inv}=0.22$)",
    "Nc3Nf3Ns1": r"$N_c=N_f=3,\ N_s=1$ ($r_{\rm inv}=0.50$)",
    "Nc8Nf8Ns6": r"$N_c=N_f=8,\ N_s=6$ ($r_{\rm inv}=0.86$)",
}
samples = [
    {"name": r"Pythia 8.310", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.313", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8313/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.313-like8310", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8313-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.315", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8315/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.315-like8310", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8315-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.318", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8318/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.318-like8310", "group": "Nc8Nf8Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8318-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-1"},
    {"name": r"Pythia 8.310", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8310/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.313", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8313/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.313-like8310", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8313-like8310/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.315", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8315/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.315-like8310", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8315-like8310/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.318", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8318/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.318-like8310", "group": "Nc3Nf3Ns1", "model": "/store/user/easmith/svj/models/fcdc-pythia8318-like8310/s-channel_mmed-1000_Nc-3_Nf-3_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.333333_Ns-1"},
    {"name": r"Pythia 8.310", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.313", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8313/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.313-like8310", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8313-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.315", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8315/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.315-like8310", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8315-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.318", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8318/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
    {"name": r"Pythia 8.318-like8310", "group": "Nc8Nf8Ns6", "model": "/store/user/easmith/svj/models/fcdc-pythia8318-like8310/s-channel_mmed-1000_Nc-8_Nf-8_scale-10_mq-10.119_mpi-6_mrho-25.0998_pvector-0.333_spectrum-fcdc_gq-0.25_gchi-0.125_Ns-6"},
]

# stylistic options
mpl.rcParams.update({
    "axes.labelsize" : 18,
    "legend.fontsize" : 16,
    "xtick.labelsize" : 14,
    "ytick.labelsize" : 14,
    "font.size" : 18,
    "legend.frameon": True,
})
# based on https://github.com/mpetroff/accessible-color-cycles (10-color set)
# blue, orange, red, gray, purple, brown, cyan, tan, slate, dark orange (cyan moved up: brown vs dark orange too similar)
colors = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#92dadd", "#b9ac70", "#717581", "#e76300"]
# follows samples order within each group; color by version: 8.310/8.313/8.315/8.318 blue/orange/red/gray, 8.313/8.315/8.318-like8310 purple/brown/cyan
version_colors = [colors[i] for i in [0,1,4,2,5,3,6]]
custom_cycler = mpl.cycler(color=version_colors)

hists = {}      # Contains the lists of histos for all models, keyed by (group, name)

for sample in samples:
    with open(f'{resolve_models(sample["model"])[0]}/Hists.pkl', "rb") as inp:
        hists[(sample["group"], sample["name"])] = pickle.load(inp)['hist']     # Dict Contains all the histos for 1 model

# quick run (--summary): only these histograms
summary = ["dark_pion_initial_pt", "dark_rho_initial_pt", "kappa", "n_pion", "n_rho", "DHJet12_nconst", "Jet12_girth", "Jet12_mass"]

# helper to make a plot: all samples overlaid in one panel, or (grid) one panel per group with a ratio to the group's first sample below
def make_plot(hname,outdir,liny=False,grid=False):
    if any(hname not in h for h in hists.values()): return
    if grid: panels = [(groups[group], [(l,h[hname]) for (g,l),h in hists.items() if g==group]) for group in groups]
    else: panels = [("", [(l,h[hname]) for (g,l),h in hists.items()])]
    nrows = 2 if grid else 1
    fig, axs = plt.subplots(nrows, len(panels), figsize=(8*len(panels),8 if grid else 6), squeeze=False, sharex="col", gridspec_kw=dict(height_ratios=[3,1][:nrows], hspace=0.05))
    for ax, (title, items) in zip(axs[0], panels):
        ax.set_prop_cycle(custom_cycler)
        for l,h in items:
            hep.histplot(h,density=True,ax=ax,label=l,flow="none",yerr=0)
        ax.set_xlim(h.axes[0].edges[0],h.axes[0].edges[-1])
        if not liny: ax.set_yscale("log")
        ax.set_ylabel("Arbitrary units")
        ax.set_title(title, fontsize=16)
        ax.legend(framealpha=0.5, ncol=2, fontsize=12, loc="upper right")
    if grid:
        for ax, rax, (title, items) in zip(axs[0], axs[1], panels):
            rax.set_prop_cycle(custom_cycler)
            ref = items[0][1].values()/items[0][1].values().sum()
            for l,h in items:
                with np.errstate(divide="ignore", invalid="ignore"):
                    rax.stairs(np.where(ref > 0, h.values()/h.values().sum()/ref, np.nan), h.axes[0].edges)
            rax.set_xlabel(ax.get_xlabel()); ax.set_xlabel("")
            rax.set_ylim(0.5, 1.5); rax.axhline(1, color="black", lw=1)
            rax.set_ylabel("/ " + items[0][0].replace("Pythia ", ""), fontsize=14)
    # raise the y range so the data max sits just below the legend
    fig.canvas.draw()
    for ax in axs[0]:
        leg_bottom = ax.get_legend().get_window_extent().transformed(ax.transAxes.inverted()).y0 - 0.02
        ymin, ymax = ax.get_ylim()[0], ax.dataLim.y1
        if liny: ax.set_ylim(ymin, ymin + (ymax-ymin)/leg_bottom)
        else: ax.set_ylim(ymin, ymin*(ymax/ymin)**(1/leg_bottom))
    plt.savefig('{}/{}.pdf'.format(outdir,hname),bbox_inches='tight')
    plt.close(fig)

def make_all_plots(outdir, liny, hnames, grid=False):
    os.makedirs(outdir,exist_ok=True)
    for hname in hnames:
        make_plot(hname,outdir,liny=any([x in hname for x in liny]),grid=grid)

if __name__=="__main__":
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsRawHelpFormatter
    )
    parser.add_argument("--dir", type=str, default="All_plots", help="output directory")
    parser.add_argument("--liny", type=str, default=[], nargs='*', help="plot histograms with matching names in linear scale")
    parser.add_argument("--summary", default=False, action="store_true", help="only plot the histograms in the summary list")
    args = parser.parse_args()

    # summary: one panel per group with ratio, linear in y; otherwise: all samples overlaid in one plot
    if args.summary: make_all_plots(args.dir, summary, summary, grid=True)
    else: make_all_plots(args.dir, args.liny, next(iter(hists.values())))
