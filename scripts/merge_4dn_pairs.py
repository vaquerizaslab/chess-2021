import fanc
import os
import logging
logging.basicConfig(level=logging.INFO)

pairs_list = [f for f in snakemake.input]
for f in pairs_list:
    logging.info(f)

output_file = str(snakemake.output)
directory = os.path.dirname(output_file)
if not os.path.exists(directory):
    os.makedirs(directory)

pairs = [fanc.load(file_name) for file_name in pairs_list]
logging.info(f"Merging and saving to {output_file}")
merged = fanc.ReadPairs.merge(pairs, file_name=output_file)