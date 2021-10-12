import fanc
import fanc.plotting as fancplot
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 10

logging.basicConfig(level=logging.INFO)

sample_names = ["control", "DLBCL"]

resolution = "25kb"
logging.info("Using Hi-C files at %s resolution", str(resolution))

hic_dict = {}
for name in sample_names:
    hic_dict[name] = fanc.load(os.path.join("data", name, "hic",
                               name + "_" + resolution + ".hic"), mode="r")


fig = plt.figure(figsize=(7, 2), dpi=300)
gs = GridSpec(ncols=7, nrows=2,
              width_ratios=[10, 10, 10, 10, 10, 10, 2],
              hspace=0.3, wspace=0.1)

regions = ["chr2:31500000-33500000", "chr2:146000000-148000000", "chr2:182500000-184500000"]

for index, region in enumerate(regions):
    index = index * 2
    region = fanc.GenomicRegion.from_string(region)
    
    # 1. plot Hi-C

    for i, sample in enumerate(sample_names):
        ax_hic = plt.subplot(gs[0, index + i])
        p_hic = fancplot.SquareMatrixPlot(hic_dict[sample],
                               norm="lin", vmin=0, vmax=0.03,
                               title=sample,
                               show_colorbar=False,
                               draw_tick_legend=False, draw_minor_ticks=False,
                               draw_tick_labels=False)

        p_hic.plot(region)
        ax_hic.set_xticks([region.start, region.end])
        ax_hic.set_yticks([region.start, region.end])
        ax_hic.set_yticklabels(["", ""])


   #  2. plot Hi-C log2 obs/exp

    for i, sample in enumerate(sample_names):
        ax_oe = plt.subplot(gs[1, index + i])
        p_oe = fancplot.SquareMatrixPlot(hic_dict[sample], colormap="coolwarm",
                              oe=True, log=True,
                              vmin=-1.5, vmax=1.5,
                              show_colorbar=False,
                              draw_tick_legend=False, draw_minor_ticks=False,
                              draw_tick_labels=i==0)

        p_oe.plot(region)
        ax_oe.set_xticks([region.start, region.end])
        ax_oe.set_yticks([region.start, region.end])
        ax_oe.set_yticklabels(["", ""])

cax_hic = plt.subplot(gs[0, 6])
p_hic.add_colorbar(ax=cax_hic, aspect=40)
p_hic.colorbar.ax.minorticks_off()
p_hic.colorbar.set_label("Normalised contact\nprobability")

cax_oe = plt.subplot(gs[1, 6])
p_oe.add_colorbar(ax=cax_oe, aspect=40)
p_oe.colorbar.ax.minorticks_off()
p_oe.colorbar.set_label("Log2(obs/exp)")

fig.savefig("figures/final/DLBCL_examples.pdf")

