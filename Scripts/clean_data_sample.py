import glob, os, shutil

#shutil.rmtree('testdir')

with open(snakemake.output[0], "w") as log:
  if snakemake.config["KEEP_TMP"]  == "T":
    log.write("No intermediate files were removed\n") 
    log.write("Config file value: " + snakemake.config["KEEP_TMP"])
  else:   
  #barcodes
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq_corrected"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/reads.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/barcodes.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/barcodes.fastq_corrected"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/reads.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  #demultiplexed data
    dir=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/primer_removed/*.fastq.gz"
    for f in glob.glob(dir):
      os.remove(f)
      log.write("removed file: "+f+"\n")
    dir=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/reads_discarded_primer"
    #if os.path.exists(dir): os.rmdir(dir); log.write("removed directory: "+dir+"\n")
    if os.path.exists(dir): shutil.rmtree(dir); log.write("removed directory: "+dir+"\n")
  #dada2 filtered fastq files(only for ASV)
    dir=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/filtered"
    #if os.path.exists(dir): os.rmdir(dir); log.write("removed directory: "+dir)
    if os.path.exists(dir): shutil.rmtree(dir); log.write("removed directory: "+dir+"\n")
  
  #peared
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")      
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.discarded.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.unassembled.forward.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.unassembled.reverse.fastq"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  
  #Split reads
    dir=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/split*/*.fna"
    for f in glob.glob(dir):
      os.remove(f)
      log.write("removed file: "+f+"\n")
  
  #intermediate files
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n") 
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.fna"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")                                                                            
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_removed.fna"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")  
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_statistics.txt"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")  
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs.unassigned.fna"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  log.close()
print("done!\nLog file:"+snakemake.output[0])

