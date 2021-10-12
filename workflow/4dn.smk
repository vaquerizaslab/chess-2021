import pandas as pd
from urllib import request
from fanc.tools.general import str_to_int
import csv
import fanc
import os
import sys


GM12878_metadata = pd.read_csv("4DNES3JX38V5_processed_files_2021-09-18-15h-40m.tsv", sep="\t", comment='#'
    ).assign(sample_name=lambda row: (row['Biosource'] + "_" + row['Bio Rep No'].astype(str) + "_" + row['Tech Rep No'].astype(str))
    ).set_index("sample_name", drop=False)
K562_metadata = pd.read_csv("4DNESI7DEJTM_processed_files_2021-09-18-15h-41m.tsv", sep="\t", comment='#'
    ).assign(sample_name=lambda row: (row['Biosource'] + "_" + row['Bio Rep No'].astype(str) + "_" + row['Tech Rep No'].astype(str))
    ).set_index("sample_name", drop=False)

fourDN_metadata = pd.concat([GM12878_metadata, K562_metadata]).query('`File Format` == "pairs"')

fourDN_merged = list(set(fourDN_metadata['Biosource'] + "_" + fourDN_metadata['Bio Rep No'].astype(str)))

Wutz_metadata = pd.read_csv("metadata_2021-09-27-15h-22m.tsv", sep="\t", comment='#'
    ).assign(sample_name=lambda row: (row['Condition Name'] + "_" + row['Rep'].astype(str))
    ).set_index("sample_name", drop=False)

fourDN_GM12878_K562_mcools_metadata = pd.read_csv("4DN_GM12878_K562_mcools.tsv", sep="\t", comment='#'
    ).assign(sample_name=lambda row: (row['Biosource'] + "_all")
    ).set_index("sample_name", drop=False)

fourDN_mcools_metadata = pd.concat([Wutz_metadata, fourDN_GM12878_K562_mcools_metadata], axis=0)

rule all:
    input:
        expand("data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}.pairs", zip,
            cell_line=fourDN_metadata['Biosource'],
            bio_rep=fourDN_metadata['Bio Rep No'],
            tech_rep=fourDN_metadata['Tech Rep No'],
            id=fourDN_metadata['File Accession']),
        expand("data/{sample_name}/hic/{sample_name}.hic", sample_name=fourDN_merged),
        expand("data/{sample_name}/hic/{sample_name}_{res}.hic",
            sample_name=fourDN_merged, res=["5kb", "10kb", "25kb"]),
        expand("data/chess/{sample1}_vs_{sample2}/{chr}_window{window}_step{step}_{res}.txt",
            sample1=fourDN_merged, sample2=fourDN_merged,
            window="2Mb", step="500kb", res=["25kb"], chr=["chr16", "chr2"]),
        expand("data/chess/{sample1}_vs_{sample2}/{chr}_window{window}_step{step}_{res}.txt",
            sample1="GM12878_1", sample2="GM12878_merged",
            window="2Mb", step="500kb", res=["25kb"], chr=["chr16", "chr2"]),
        expand("data/boundaries/{sample}_{res}_{w}.bw",
            sample=fourDN_merged + ["GM12878_merged"], res="25kb", w=[4, 6, 8, 10, 15]),
        expand("data/{sample}/hic/{sample}_{res}_marginals.bed",
            sample=fourDN_merged + ["GM12878_merged"], res="25kb"),
        # merged .mcools
        expand("data/chess/{sample1}_vs_{sample2}/{chr}_window{window}_step{step}_{res}.txt",
            sample1="GM12878_all", sample2="K562_all",
            window="2Mb", step="500kb", res=["25kb"], chr=["chr16", "chr2"]),
        expand("data/boundaries/{sample}_{res}_{w}.bw",
            sample=fourDN_GM12878_K562_mcools_metadata.index, res="25kb", w=[4, 6, 8, 10, 15]),
        expand("data/{sample}/hic/{sample}_{res}_marginals.bed",
            sample=fourDN_GM12878_K562_mcools_metadata.index, res="25kb"),
        ## downsampled
        expand("data/chess/{sample1}_vs_{sample2}/downsampled/{chr}_window{window}_step{step}_{res}.txt",
            sample1=["GM12878_3", "GM12878_5", "K562_1", "K562_2"], 
            sample2=["GM12878_3", "GM12878_5", "K562_1", "K562_2"],
            window="2Mb", step="500kb", res=["25kb"], chr=["chr16", "chr2"]),
        expand("data/boundaries/downsampled/{sample}_{res}_{w}.bw",
            sample=["GM12878_3", "GM12878_5", "K562_1", "K562_2"], res="25kb", w=[4, 6, 8, 10, 15]),
        expand("data/{sample}/hic/downsampled/{sample}_{res}_marginals.bed",
            sample=["GM12878_3", "GM12878_5", "K562_1", "K562_2"], res="25kb"),
        expand("figures/chess/{comparison}_window2Mb_step500kb_25kb_ssim10_SN90.pdf",
            comparison=["K562_1_vs_GM12878_3", "K562_1_vs_GM12878_5",
            "K562_2_vs_GM12878_3", "K562_2_vs_GM12878_5"]),
        "figures/chess/GM12878_vs_K562_window2Mb_step500kb_25kb_negative_set.pdf",

        # Wutz et al. data
        expand("data/{condition}_{bio_rep}/hic/{condition}_{bio_rep}_{res}.hic", zip,
            condition=Wutz_metadata["Condition Name"],
            bio_rep=Wutz_metadata["Rep"],
            res=["25kb", "50kb"]),
        expand("data/boundaries/{sample}_{res}_{w}.bw",
            sample=Wutz_metadata.index, res="50kb", w=[4, 6, 8, 10, 15]),
        expand("data/{sample}/hic/{sample}_{res}_marginals.bed",
            sample=Wutz_metadata.index, res="50kb"),
        expand("data/chess/{sample1}_vs_{sample2}/{chr}_window{window}_step{step}_{res}.txt",
            sample1=Wutz_metadata.index, sample2=Wutz_metadata.index,
            window="4Mb", step="1Mb", res=["50kb"], chr=["chr16", "chr2"]),
        expand("data/{sample}/hic/{sample}_{res}_hic_stats.txt",
            sample=Wutz_metadata.index, res=["25kb", "50kb"]),
        ## downsampled         
        expand("data/chess/{sample1}_vs_{sample2}/downsampled/{chr}_window{window}_step{step}_{res}.txt",
            sample1=Wutz_metadata.index,
            sample2=Wutz_metadata.index,
            window="4Mb", step="1Mb", res=["50kb"], chr=["chr16", "chr2"]),
        expand("data/boundaries/downsampled/{sample}_{res}_{w}.bw",
            sample=Wutz_metadata.index, res="50kb", w=[4, 6, 8, 10, 15]),
        expand("data/{sample}/hic/downsampled/{sample}_{res}_marginals.bed",
            sample=Wutz_metadata.index, res="50kb"),

rule get_4dn_data:
    output:
        "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}.pairs.gz"
    run:
        sample_name = "_".join([wildcards.cell_line, wildcards.bio_rep, wildcards.tech_rep])
        url = fourDN_metadata.loc[sample_name, "Open Data URL"]
        file_id = fourDN_metadata.loc[sample_name, "File Accession"]
        request.urlretrieve(url, output[0])

rule convert_pairs:
    input:
        pairs = "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}.pairs.gz",
        genome = "data/hs_MboI_fragments.bed"
    output:
        pairs = "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}.pairs",
        redist_plot = "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}_filtered_re_dist.png",
        ligation_plot = "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}_filtered_ligation_error.png",
        stats_plot = "data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}_filtered_stats.png",
    shell:
        "fanc pairs -tmp "
        "-g {input.genome} "
        "--filter-inward 5000 --filter-outward 5000 "
        "--filter-self-ligations --filter-pcr-duplicates 1 "
        "--statistics-plot {output.stats_plot} "
        "--re-dist-plot {output.redist_plot} "
        "--ligation-error-plot {output.ligation_plot} "
        "{input.pairs} {output.pairs} "


def get_4dn_pairs_input(wildcards):
    ids = list(fourDN_metadata.loc[(fourDN_metadata['Biosource'] == wildcards.cell_line) &
                                         (fourDN_metadata['Bio Rep No'] == int(wildcards.bio_rep)), 'File Accession'])

    tech_reps = list(fourDN_metadata.loc[(fourDN_metadata['Biosource'] == wildcards.cell_line) &
                                         (fourDN_metadata['Bio Rep No'] == int(wildcards.bio_rep)), 'Tech Rep No'])
    pairs = expand("data/4dn/{cell_line}_{bio_rep}_{tech_rep}/pairs/{id}.pairs", zip,
        cell_line=[wildcards.cell_line for i in ids],
        bio_rep=[wildcards.bio_rep for i in ids],
        tech_rep=tech_reps,
        id=ids)
    return(pairs)


rule merge_4dn_pairs:
    input:
        get_4dn_pairs_input
    output:
        "data/{cell_line}_{bio_rep}/pairs/{cell_line}_{bio_rep}.pairs"
    wildcard_constraints:
        cell_line="|".join(set(fourDN_metadata['Biosource']))
    script:
        "../scripts/merge_4dn_pairs.py"


rule get_pairs_stats:
    input:
        pairs = "data/{sample_name}/pairs/{sample_name}.pairs"
    output:
        "data/{sample_name}/{sample_name}_hic_stats.txt"
    wildcard_constraints:
        sample_name="|".join(fourDN_merged)
    shell: "python scripts/get_stats.py data/{wildcards.sample_name}/pairs/ {output}"

rule make_merged_hic:
    input:
        pairs = "data/{sample_name}/pairs/{sample_name}.pairs",
        stats = "data/{sample_name}/{sample_name}_hic_stats.txt",
    output:
        "data/{sample_name}/hic/{sample_name}.hic"
    wildcard_constraints:
        sample_name="|".join(fourDN_merged)
    shell:
        "fanc hic -tmp {input.pairs} {output}"

rule bin_hic:
    input:
        "data/{sample_name}/hic/{sample_name}.hic"
    output:
        "data/{sample_name}/hic/{sample_name}_{res}.hic"
    threads: 8
    wildcard_constraints:
        res="|".join(["25kb"]),
        sample_name="|".join(fourDN_merged)
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --normalise  "
        "{input} {output}"

rule merge_replicates:
    input:
        expand("data/{sample_name}/hic/{sample_name}.hic",
            sample_name=["GM12878_" + str(i) for i in range(2, 10)])
    output:
        "data/GM12878_merged/hic/GM12878_merged.hic"
    shell:
        "fanc hic {input} {output}"

rule subset_hic:
    input:
        "data/{sample_name}/hic/{sample_name}_{res}.hic"
    output:
        "data/{sample_name}/hic/{sample_name}_{res}_{chr}.hic"
    shell:
        "fanc subset {input} {output} {wildcards.chr}"

rule chess_make_pairs:
    output:
        "data/chess/{genome}_pairs_window{window}_step{step}.bedpe"
    params:
        window_bp = lambda wildcards: str_to_int(wildcards.window),
        step_bp = lambda wildcards: str_to_int(wildcards.step)
    shell:
        "chess pairs {wildcards.genome} {params.window_bp} {params.step_bp} {output} && "
        "sed -i -e 's/chr//g' {output}" # remove chr prefixes!

rule chess_sim:
    input:
        matrix1 = "data/{sample1}/hic/{sample1}_{res}_{chr}.hic",
        matrix2 = "data/{sample2}/hic/{sample2}_{res}_{chr}.hic",
        pairs = "data/chess/hg38_pairs_window{window}_step{step}.bedpe"
    output:
        "data/chess/{sample1}_vs_{sample2}/{chr}_window{window}_step{step}_{res}.txt"
    threads: 8
    shell:
        "chess sim "
        "{input.matrix1} {input.matrix2}  "
        "{input.pairs} {output} "
        "-p {threads}"

# insulation score calulations

rule insulation_score:
    input: 
        "data/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/boundaries/{sample}_{res}.ii"
    params:
        window_sizes = lambda wildcards: [int(re.findall('([0-9]+)kb', wildcards.res)[0]) * 1000 * w for w in [4, 6, 8, 10, 15]]
    wildcard_constraints:
        sample = "|".join(Wutz_metadata.index),
    shell: 
        "fanc insulation {input} {output} -w {params.window_sizes}"

rule insulation_bw:
    input:
        "data/boundaries/{sample}_{res}.ii"
    output:
        "data/boundaries/{sample}_{res}_{w}.bw"
    params:
        window_size = lambda wildcards: int(re.findall('([0-9]+)kb', wildcards.res)[0]) * 1000 * int(wildcards.w)
    wildcard_constraints:
        sample = "|".join(Wutz_metadata.index),
    run:
        import fanc
        ii = fanc.load(input[0])
        ii.to_bigwig(output[0], params.window_size)

rule export_marginals:
    input:
        "data/{sample}/hic/{sample}_{res}.hic"
    output:
        "data/{sample}/hic/{sample}_{res}_marginals.bed"
    wildcard_constraints:
        sample = "|".join(Wutz_metadata.index),
    shell:
        "python scripts/export_marginals.py {input} {output}"

rule get_4dn_mcools:
    output:
        "data/{cell_line}_all/hic/{id}.mcool"
    run:
        sample_name = wildcards.cell_line + "_all"
        url = fourDN_mcools_metadata.loc[sample_name, "Open Data URL"]
        file_id = fourDN_mcools_metadata.loc[sample_name, "File Accession"]
        request.urlretrieve(url, output[0])

# rule get_wutz_4dn_data:
#     output:
#         "data/{condition}_{bio_rep}/hic/{id}.mcool"
#     run:
#         sample_name = "_".join([wildcards.condition, wildcards.bio_rep])
#         url = Wutz_metadata.loc[sample_name, "Open Data URL"]
#         file_id = Wutz_metadata.loc[sample_name, "File Accession"]
#         request.urlretrieve(url, output[0])


def get_mcool(wildcards):
    sample_name = wildcards.sample_name
    id = fourDN_mcools_metadata.loc[sample_name, "File Accession"]
    return(f"data/{sample_name}/hic/{id}.mcool")


rule convert_mcool:
    input:
        get_mcool
    output:
        "data/{sample_name}/hic/{sample_name}_{res}.hic"
    wildcard_constraints:
        sample_name = "|".join(fourDN_mcools_metadata.index),
        res = "25kb|50kb"
    shell:
        "fanc from-cooler {input}@{wildcards.res} {output} -tmp"


def write_stats(input_file, output_file):
    """Write stats to file."""
    hic = fanc.load(input_file, mode='r')
    total = sum(e.weight for e in hic.edges(lazy=True, norm=False))

    with open(output_file, "a+") as out:
        out = csv.writer(out)
        out.writerow(["total", total])


rule get_hic_stats:
    input:
        "data/{sample_name}/hic/{sample_name}_{res}.hic"
    output:
        "data/{sample_name}/hic/{sample_name}_{res}_hic_stats.txt"
    run:
        write_stats(str(input), str(output))


rule downsample_hic:
    input:
         "data/{sample_name}/hic/{sample_name}_50kb.hic"
    output:
         "data/{sample_name}/hic/downsampled/{sample_name}_50kb.hic"
    params:
        ref = 147993483
    wildcard_constraints:
        sample_name = "|".join(Wutz_metadata.index),
        res = "50kb"
    shell:
        "fanc hic -tmp --downsample {params.ref} "
        "--normalise --norm-method KR {input} {output} "

rule downsample_hic2:
    input:
         "data/{sample_name}/hic/{sample_name}_25kb.hic"
    output:
         "data/{sample_name}/hic/downsampled/{sample_name}_25kb.hic"
    params:
        ref = 219937477
    wildcard_constraints:
        sample_name = "|".join(["GM12878_3", "GM12878_5", "K562_1", "K562_2"]),
    shell:
        "fanc hic -tmp --downsample {params.ref} "
        "--normalise --norm-method KR {input} {output} "

rule subset_downsampled_hic:
    input:
        "data/{sample_name}/hic/downsampled/{sample_name}_{res}.hic"
    output:
        "data/{sample_name}/hic/downsampled/{sample_name}_{res}_{chr}.hic"
    shell:
        "fanc subset {input} {output} {wildcards.chr}"

rule chess_sim_downsampled:
    input:
        matrix1 = "data/{sample1}/hic/downsampled/{sample1}_{res}_{chr}.hic",
        matrix2 = "data/{sample2}/hic/downsampled/{sample2}_{res}_{chr}.hic",
        pairs = "data/chess/hg38_pairs_window{window}_step{step}.bedpe"
    output:
        "data/chess/{sample1}_vs_{sample2}/downsampled/{chr}_window{window}_step{step}_{res}.txt"
    threads: 8
    shell:
        "chess sim "
        "{input.matrix1} {input.matrix2}  "
        "{input.pairs} {output} "
        "-p {threads}"

rule insulation_score_downsampled:
    input:
        "data/{sample}/hic/downsampled/{sample}_{res}.hic"
    output:
        "data/boundaries/downsampled/{sample}_{res}.ii"
    params:
        window_sizes = lambda wildcards: [int(re.findall('([0-9]+)kb', wildcards.res)[0]) * 1000 * w for w in [4, 6, 8, 10, 15]]
    wildcard_constraints:
        sample = "|".join(list(Wutz_metadata.index) + ["GM12878_3", "GM12878_5", "K562_1", "K562_2"]),
    shell:
        "fanc insulation {input} {output} -w {params.window_sizes}"

rule insulation_bw_downsampled:
    input:
        "data/boundaries/downsampled/{sample}_{res}.ii"
    output:
        "data/boundaries/downsampled/{sample}_{res}_{w}.bw"
    params:
        window_size = lambda wildcards: int(re.findall('([0-9]+)kb', wildcards.res)[0]) * 1000 * int(wildcards.w)
    wildcard_constraints:
        sample = "|".join(list(Wutz_metadata.index) + ["GM12878_3", "GM12878_5", "K562_1", "K562_2"]),
    run:
        import fanc
        ii = fanc.load(input[0])
        ii.to_bigwig(output[0], params.window_size)

rule export_marginals_downsampled:
    input:
        "data/{sample}/hic/downsampled/{sample}_{res}.hic"
    output:
        "data/{sample}/hic/downsampled/{sample}_{res}_marginals.bed"
    wildcard_constraints:
        sample = "|".join(list(Wutz_metadata.index) + ["GM12878_3", "GM12878_5", "K562_1", "K562_2"]),
    shell:
        "python scripts/export_marginals.py {input} {output}"


rule plot_chess_regions:
    input:
        bed = "data/chess/{region_set}.bed"
    output:
        "figures/chess/{region_set}.pdf",
    shell:
        "python scripts/plot_chess_regions_downsampledGM12878-K562.py {input.bed} {output}"
