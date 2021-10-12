import fanc
import fanc.plotting as fancplot
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
plt.rcParams['pdf.fonttype'] = 42

# from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(level=logging.INFO)

sample_names = ["wt_G1_none_1", "CTCFdegron_G1_auxin0_1", "SCC1degron_G1_auxin0_1",
                "SCC1degron_prometaphase_noRNAi_1", "CTCFdegron_G1_auxin120_1", "SCC1degron_G1_auxin180_1", ]

resolution = "50kb"
logging.info("Using Hi-C files at %s resolution", str(resolution))

hic_dict = {}
for name in sample_names:
    hic_dict[name] = fanc.load(os.path.join("data", name, "hic", "downsampled",
                               name + "_" + resolution + ".hic"), mode="r")


region = "chr2:120,000,000-128,000,000"
region = fanc.GenomicRegion.from_string(region)


fig = plt.figure(figsize=(3, 2), dpi=300)
gs = GridSpec(ncols=4, nrows=2,
              width_ratios=[10, 10, 10, 1.5],
              hspace=0.3, wspace=0.1)


# 1. plot Hi-C

for i, sample in enumerate(sample_names):
    if i < 3:
        row = 0
        col = i
    else:
        row = 1
        col = i - 3

    ax_hic = plt.subplot(gs[row, col])
    p_hic = fancplot.SquareMatrixPlot(hic_dict[sample],
                           norm="lin", title="", vmin=0, vmax=0.02,
                           show_colorbar=False,
                           draw_tick_legend=False, draw_minor_ticks=False,
                           draw_tick_labels=i==0)

    p_hic.plot(region)
    ax_hic.set_xticks([region.start, region.end])
    ax_hic.set_yticks([region.start, region.end])
    ax_hic.set_yticklabels(["", ""])

cax_hic = plt.subplot(gs[:, 3])
p_hic.add_colorbar(ax=cax_hic, aspect=80)
p_hic.colorbar.ax.minorticks_off()
p_hic.colorbar.set_label("Normalised contact probability")

fig.savefig("figures/final/Wutz_examples.pdf")
