import os
import subprocess

treads = subprocess.run( ["grep '^>' " + snakemake.input.spplited + " | wc -l"],stdout=subprocess.PIPE, shell=True)
totalReads =  treads.stdout.decode('utf-8').strip()
tr = int(totalReads)
     #"split_libraries_fastq.py -m {input.mapFile} -i {input.rFile} "
     #"-o {params.outDirRC} -b {input.bcFile} -q {config[split][q]} "
     #"--barcode_type {config[split][barcode_type]} {config[split][extra_params]} --rev_comp_mapping_barcodes --rev_comp"
#if we run split and split rc with retain_unassigned_reads we can duplicate a lot of sequences
extra_params = snakemake.config["split"]["extra_params"]
#extra_params = extra_params.replace("--retain_unassigned_reads", "")
subprocess.run( [snakemake.config["qiime"]["path"]+"split_libraries_fastq.py -m " + snakemake.input.mapFile + " -i "+ snakemake.input.rFile + " -o " + snakemake.params.outDirRC +
 " -b " + snakemake.input.bcFile + " -q " + snakemake.config["split"]["q"] + " -r " + snakemake.config["split"]["r"] + " --barcode_type " +
 snakemake.config["split"]["barcode_type"] +" " + extra_params + " --retain_unassigned_reads -s " + str(tr) ],  stdout=subprocess.PIPE, shell=True)
#before changes
#subprocess.run( ["split_libraries_fastq.py -m " + snakemake.input.mapFile + " -i "+ snakemake.input.rFile + " -o " + snakemake.params.outDirRC +
# " -b " + snakemake.input.bcFile + " -q " + snakemake.config["split"]["q"] + " --barcode_type " + snakemake.config["split"]["barcode_type"] +
# " " + extra_params + " --rev_comp_mapping_barcodes --rev_comp -s " + str(tr) ],  stdout=subprocess.PIPE, shell=True)
