import glob, os

with open(snakemake.output[0], "w") as log:
  if snakemake.config["KEEP_TMP"]  == "T":
    log.write("No intermediate files were removed\n") 
    log.write("Config file value: " + snakemake.config["KEEP_TMP"])
  else:   
  #derep
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/seqs_fw_rev_combined_derep.fasta"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/seqs_fw_rev_combined_derep.uc"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  #OTU                        
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/seqs_fw_rev_combined_derep_clusters.uc"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/seqs_fw_rev_combined_derep_otus.txt"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  
  #OTU taxonomy
    file=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable.tmp.biom"
    if os.path.exists(file): os.remove(file); log.write("removed file: "+file+"\n")
  log.close()
print("done!\nLog file:"+snakemake.output[0])

