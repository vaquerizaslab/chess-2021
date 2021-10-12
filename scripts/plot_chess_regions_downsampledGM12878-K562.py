import fanc
import fanc.plotting as fancplot
import logging
import re
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import pybedtools as pbt
import sys

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]
output_file = sys.argv[2]


logging.info("Working on input regions file %s", str(input_file))
logging.info("Will write output to %s", str(output_file))


sample_names = ["GM12878_3", "GM12878_5", "K562_1", "K562_2"]

resolution = re.match(".*_([0-9]*kb)", os.path.basename(input_file)).groups()[0]
logging.info("Using Hi-C files at %s resolution", str(resolution))

hic_dict = {}
for name in sample_names:
    hic_dict[name] = fanc.load(os.path.join("data", name, "hic", "downsampled",
                               name + "_" + resolution + ".hic"), mode="r")


def plot_region(region):
    fig = plt.figure(figsize=(9, 6), dpi=300)
    gs = GridSpec(ncols=5, nrows=3,
                  width_ratios=[10, 10, 10, 10, 2],
                  hspace=0.3, wspace=0.1)
 

    # 1. plot Hi-C

    for i, sample in enumerate(sample_names):
        ax_hic = plt.subplot(gs[0, i])
        p_hic = fancplot.SquareMatrixPlot(hic_dict[sample],
                               norm="lin", title=sample, vmin=0, vmax=0.03,
                               show_colorbar=False,
                               draw_tick_legend=False, draw_minor_ticks=False,
                               draw_tick_labels=i==0)

        p_hic.plot(region)
        ax_hic.set_xticks([region.start, region.end])
        ax_hic.set_yticks([region.start, region.end])
        ax_hic.set_yticklabels(["", ""])

    cax_hic = plt.subplot(gs[0, 4])
    p_hic.add_colorbar(ax=cax_hic, aspect=40)
    p_hic.colorbar.ax.minorticks_off()
    p_hic.colorbar.set_label("Normalised contact\nprobability")


   #  2. plot Hi-C log2 obs/exp

    for i, sample in enumerate(sample_names):
        ax_oe = plt.subplot(gs[1, i])
        p_oe = fancplot.SquareMatrixPlot(hic_dict[sample], colormap="coolwarm",
                              oe=True, log=True,
                              vmin=-1.5, vmax=1.5,
                              show_colorbar=False,
                              draw_tick_legend=False, draw_minor_ticks=False,
                              draw_tick_labels=False)

        p_oe.plot(region)
        ax_oe.set_xticks([region.start, region.end])
        ax_oe.set_yticks([region.start, region.end])
        ax_oe.set_yticklabels(["", ""])

    cax_oe = plt.subplot(gs[1, 4])
    p_oe.add_colorbar(ax=cax_oe, aspect=40)
    p_oe.colorbar.ax.minorticks_off()
    p_oe.colorbar.set_label("Log2(observed/expected)")

    # 4. plot uncorrected Hi-C

    for i, sample in enumerate(sample_names):
        ax_hic = plt.subplot(gs[2, i])
        p_hic = fancplot.SquareMatrixPlot(hic_dict[sample],
                               norm="lin", title=sample, matrix_norm=False, vmax=40,
                               show_colorbar=False,
                               draw_tick_legend=False, draw_minor_ticks=False,
                               draw_tick_labels=i==0)

        p_hic.plot(region)
        ax_hic.set_xticks([region.start, region.end])
        ax_hic.set_yticks([region.start, region.end])
        ax_hic.set_yticklabels(["", ""])

    cax_hic = plt.subplot(gs[2, 4])
    p_hic.add_colorbar(ax=cax_hic, aspect=40)
    p_hic.colorbar.ax.minorticks_off()
    p_hic.colorbar.set_label("Read pairs")


bed = pbt.BedTool(input_file)

with PdfPages(output_file) as pdf:
    for region in bed:
        region = "{}:{}-{}".format(region.chrom, region.start, region.end)
        logging.info(region)
        region = fanc.GenomicRegion.from_string(region)
        plot_region(region)
        pdf.savefig()
