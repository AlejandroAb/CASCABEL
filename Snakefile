"""
CASCABEL
Version: 7.0.0
Authors:
 - Alejandro Abdala
 - Julia Engelmann
 - Wietse Reitsma

Last update: 25/11/2025
"""
run=config["RUN"]

def selectInput():
    if (config["align_vs_reference"]["align"] == "T" and config["primers"]["remove"] != "F"):
        return ["{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.oneline.fna"]
    elif (config["align_vs_reference"]["align"] == "T" and config["primers"]["remove"] == "F"):
        return ["{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.oneline.fna"]
    elif (config["align_vs_reference"]["align"] != "T" and config["primers"]["remove"] != "F" and config["ANALYSIS_TYPE"]=="OTU"):
        return ["{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"]
    else:
        return ["{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"]

rule all:
    input:
        expand("{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".pdf" if config["pdfReport"] == "T" else 
               "{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        if config["ANALYSIS_TYPE"]=="OTU" else
         expand("{PROJECT}/runs/{run}/report_{sample}_dada2.pdf" if config["pdfReport"] == "T" else
               "{PROJECT}/runs/{run}/report_{sample}_dada2.html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run),
        expand("{PROJECT}/runs/{run}/report_"+config["assignTaxonomy"]["tool"].lower()+".zip",PROJECT=config["PROJECT"], run=run,sample=config["LIBRARY"])
        if config["ANALYSIS_TYPE"]=="OTU" else
        expand("{PROJECT}/runs/{run}/report_dada2.zip",PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)

if len(config["LIBRARY"])==1 and config["demultiplexing"]["demultiplex"] == "T" and len(config["input_files"])<2 and config["LIBRARY_LAYOUT"] != "SE":
    rule init_structure:
        input:
            fw = config["fw_reads"],
            rv = config["rv_reads"],
            metadata = config["metadata"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        shell:
            "Scripts/init_sample.sh "+config["PROJECT"]+" "+config["LIBRARY"][0]+" {input.metadata} {input.fw} {input.rv}"
elif len(config["LIBRARY"])==1 and config["demultiplexing"]["demultiplex"] == "T" and len(config["input_files"])<2 and config["LIBRARY_LAYOUT"] == "SE":
    rule init_structure:
        input:
            fw = config["fw_reads"],
            metadata = config["metadata"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        shell:
            "Scripts/init_sample_SE.sh "+config["PROJECT"]+" "+config["LIBRARY"][0]+" {input.metadata} {input.fw}"
elif len(config["LIBRARY"])==1 and config["demultiplexing"]["demultiplex"] == "F" and len(config["input_files"])<2 and config["LIBRARY_LAYOUT"] != "SE":
    rule init_structure:
        input:
            fw = config["fw_reads"],
            rv = config["rv_reads"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        shell:
            "Scripts/init_sample_dmx.sh "+config["PROJECT"]+" "+config["LIBRARY"][0]+"  {input.fw} {input.rv}"

elif len(config["LIBRARY"])==1 and config["demultiplexing"]["demultiplex"] == "F" and len(config["input_files"])<2 and config["LIBRARY_LAYOUT"] == "SE":
    rule init_structure:
        input:
            fw = config["fw_reads"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        shell:
            "Scripts/init_sample_dmx_SE.sh "+config["PROJECT"]+" "+config["LIBRARY"][0]+"  {input.fw}"

elif config["demultiplexing"]["demultiplex"] == "T" and config["LIBRARY_LAYOUT"] != "SE" and len(config["input_files"])>2:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        script:
            "Scripts/init_sample.py"
elif config["demultiplexing"]["demultiplex"] == "T" and config["LIBRARY_LAYOUT"] == "SE" and len(config["input_files"])>2:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        script:
            "Scripts/init_sample_SE.py"

elif config["LIBRARY_LAYOUT"] == "SE" and len(config["input_files"])>2:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        script:
            "Scripts/init_sample_SE.py"
else: 
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/init_structure.benchmark"
        script:
            "Scripts/init_sample.py"


rule fast_qc:
    """
    Runs QC on raw reads
    """
    input:
        r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" and config["LIBRARY_LAYOUT"] != "SE" else 
        "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["LIBRARY_LAYOUT"] != "SE"  else []
    output:
        o1="{PROJECT}/samples/{sample}/qc/fastqc/fw_fastqc.html",
        o2="{PROJECT}/samples/{sample}/qc/fastqc/rv_fastqc.html" if config["LIBRARY_LAYOUT"] != "SE" else [],
        s1="{PROJECT}/samples/{sample}/qc/fastqc/fw_fastqc/summary.txt",
        s2="{PROJECT}/samples/{sample}/qc/fastqc/rv_fastqc/summary.txt" if config["LIBRARY_LAYOUT"] != "SE" else [],
        r=report(directory("{PROJECT}/samples/{sample}/qc/fastqc/"), patterns=["{name}_fastqc.html"], category="2. Quality", subcategory="Raw reads",labels={"sample":"{sample}","Tool":"FastQC", "Read":"{name}"},caption="resources/caption/fast_qc.rst",)
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/fast_qc.benchmark"
    shell:
        #"{config[fastQC][command]} {input.r1} {input.r2} --extract {config[fastQC][extra_params]} -t {config[fastQC][threads]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/"
        "{config[fastQC][command]} {input} --extract {config[fastQC][extra_params]} -t {config[fastQC][threads]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/fastqc/"

rule sequali:
    """
    Runs QC with sequali on raw reads
    """
    input:
        r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" and config["LIBRARY_LAYOUT"] != "SE" else 
        "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["LIBRARY_LAYOUT"] != "SE" else []
        if config["LIBRARY_LAYOUT"] != "SE" else []
    output:
        o1="{PROJECT}/samples/{sample}/qc/sequali/sequali.html",
        s2="{PROJECT}/samples/{sample}/qc/sequali/sequali.json",
        r=report(directory("{PROJECT}/samples/{sample}/qc/sequali/"), htmlindex="sequali.html", category="2. Quality", subcategory="Raw reads",labels={"sample":"{sample}","Tool":"Sequali", "Read":"fw & rv"},)
    params:
        outdir="{PROJECT}/samples/{sample}/qc/sequali/"
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/sequali.benchmark"
    shell:
        "sequali --outdir {params.outdir} --html  sequali.html --json sequali.json -t {config[sequali][threads]}  {config[sequali][extra_params]}  {input}"

rule validateQC:
    """
    Interpret FastQC output and stops on interactive mode, if too many
    errors.
    """
    input:
        fws="{PROJECT}/samples/{sample}/qc/fastqc/fw_fastqc/summary.txt",
        rvs="{PROJECT}/samples/{sample}/qc/fastqc/rv_fastqc/summary.txt" if config["LIBRARY_LAYOUT"] != "SE" else [],
        fhtml="{PROJECT}/samples/{sample}/qc/fastqc/fw_fastqc.html",
        rhtml="{PROJECT}/samples/{sample}/qc/fastqc/rv_fastqc.html" if config["LIBRARY_LAYOUT"] != "SE" else [],
        rfw="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        rrv="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F"  and config["LIBRARY_LAYOUT"] == "PE" 
        else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["LIBRARY_LAYOUT"] == "PE" else []
    output:
        "{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt"
    params:
        base_dir="{PROJECT}/samples/{sample}/",
        threshold= str(config["fastQC"]["qcLimit"])
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/validateQC.benchmark"
    shell:
        "python Scripts/validateFastQC.py -f {input.fws} -r {input.rvs} -fr {input.rfw} -rr {input.rrv} -o {output} -t {params.threshold}"
        if config["LIBRARY_LAYOUT"] != "SE" else
        "python Scripts/validateFastQC.py -f {input.fws}  -fr {input.rfw}  -o {output} -t {params.threshold}" #validate indexes

rule bc_mapping_validation:
    """
    Check if the mapping file is ok: 
      - Header line needs to start with "#"
      - Looks for ...
      - Creates 1 or 2 fasta files with barcodes
    """
    input:
        mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/bc_mapping_validation.benchmark"
    params:
        dmx="--demultiplex-check" if config["demultiplexing"]["demultiplex"].lower() == "t" else "--no-demultiplex-check",
        prm="--primer-check" if config["primers"]["remove"] == "CFG" else "--no-primer-check"
    output:
        log="{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log", # antes .log pero creo carpeta .log??
        fw="{PROJECT}/metadata/{sample}/fw_barcodes.fasta", #applies for PE and SE
        rv="{PROJECT}/metadata/{sample}/rv_barcodes.fasta" if config["demultiplexing"]["demultiplex"].lower() == "t" and config["LIBRARY_LAYOUT"] == "PE" 
        else []
    shell:
        "python Scripts/validateMappingFile.py -l " + config["LIBRARY_LAYOUT"] +" -o {output.log} {params.dmx} {params.prm} -f {output.fw} -r {output.rv} {input.mapp}"
        if config["demultiplexing"]["demultiplex"].lower() == "t" and config["LIBRARY_LAYOUT"] == "PE" else 
        "python Scripts/validateMappingFile.py -l " + config["LIBRARY_LAYOUT"] +" -o {output.log} {params.dmx} {params.prm} -f {output.fw}  {input.mapp}"
        if config["demultiplexing"]["demultiplex"].lower() == "t" and config["LIBRARY_LAYOUT"] == "SE" else 
        "python Scripts/validateMappingFile.py -l " + config["LIBRARY_LAYOUT"] +" -o {output.log}  {params.dmx} {params.prm} {input.mapp}"

rule create_temp_mapping_file:
    input:
        mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
    output:
        tsv_mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.tsv"
    shell:
        "ln -s $PWD\"/\"{input} {output}"

rule create_yaml_metadata_tbl:
    input:
        "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.tsv"
    output:
        "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.yaml"
    shell:
        "datavzrd suggest --files {input} --separators $'\\t' |  awk -v s=\"{wildcards.sample}\"  "+
        "'NR==1{{print \"__use_yte__: true\"}}"+
        " NR==2{{print \"name: \" s  }} "+
        " NR==3{{print $0  }} "+
        #" NR==4{{print \"  Metadata \" s}} "+
        " NR==4{{print $0}} "+
        " NR==5{{print \"    path: ?input.table\" }}"+
        " NR>5{{print $0}}' |  perl -0777 " +
        " -pe 's/.*plot:\\n.*heatmap.*\\n.*scale: ordinal.*\\n.*clamp.*\\n.*color-scheme.*\\n.*range.*\\n.*aux-domain.*\\n//g' >  {output}"
        #" NR>5{{if($0 ~ \"path:\"){{}}else{{print $0}}}}'> {output}"



rule datavzrd_metadata:
    input:
        config="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.yaml",
        table="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.tsv"
    output:
        report(
            directory("{PROJECT}/metadata/{sample}_report/"),
            htmlindex="index.html",
            category="1. Project information: {PROJECT}",
            subcategory="Run: " + run,
            labels={"sample": "{sample}"}, 
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"  


#rule project_information:
#    output:
#        tsv_data="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.tsv"
#    shell:
#        "ln -s $PWD\"/\"{input} {output}"

rule cutadapt_demux:
    """
    This rule runs cutadapt on raw reads to perform demultiplexing.
    Before running this quality on raw reads must be checked, and mapping
    file should be validated.
    """
    input:
        r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" and config["LIBRARY_LAYOUT"] == "PE" else  
        "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["gzip_input"] == "T" and config["LIBRARY_LAYOUT"] == "PE" else [],
        tmp1="{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
        tmp3="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else [],
        map_log="{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log", # antes .log pero creo carpeta .log??
        rep_bc="{PROJECT}/metadata/{sample}_report/",
        fw_bc="{PROJECT}/metadata/{sample}/fw_barcodes.fasta", #applies for PE and SE
        rv_bc="{PROJECT}/metadata/{sample}/rv_barcodes.fasta" if config["LIBRARY_LAYOUT"] == "PE" else []
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/cutadapt_demux.benchmark"
    params:
        odir="{PROJECT}/samples/{sample}/cutadapt/",
        fix_opt="  --pair-adapters",
        cpus=" -j " + str(config["demultiplexing"]["cpus"])        
    output:
        log="{PROJECT}/samples/{sample}/cutadapt/cutadapt_s1.log",
        ufw="{PROJECT}/samples/{sample}/cutadapt/u_1_fw.fastq.gz",
        urv="{PROJECT}/samples/{sample}/cutadapt/u_1_rv.fastq.gz" if config["LIBRARY_LAYOUT"] == "PE" else []
    shell:
        "cutadapt -e " + str(config["demultiplexing"]["bc_mismatch"]) + "{params.fix_opt} {params.cpus}" +
        " -g ^file:{input.fw_bc} -G ^file:{input.rv_bc} " +
        " -o  {params.odir}{{name}}_1.fastq.gz -p {params.odir}{{name}}_2.fastq.gz " +
        " --untrimmed-output {output.ufw} --untrimmed-paired-output {output.urv} "  + str(config["demultiplexing"]["cutadapt_extra_params"]) + 
        " {input.r1} {input.r2} > {output.log}"  if config["LIBRARY_LAYOUT"] == "PE" else
         "cutadapt -e " + str(config["demultiplexing"]["bc_mismatch"]) + " {params.cpus}" +
        " -g ^file:{input.fw_bc} " +
        " -o  {params.odir}{{name}}_1.fastq.gz " +
        " --untrimmed-output {output.ufw} " + str(config["demultiplexing"]["cutadapt_extra_params"]) +
        " {input.r1} > {output.log}"
                        
rule demux_unassigned:
    """
    This rule take all the unassigned reads, and performs a second round of cutadapt
    but now using the reverse file as the forward and the other way around.
    """
    input:
        fw_bc="{PROJECT}/metadata/{sample}/fw_barcodes.fasta", 
        rv_bc="{PROJECT}/metadata/{sample}/rv_barcodes.fasta"  if config["LIBRARY_LAYOUT"] == "PE" else [],
        ufw="{PROJECT}/samples/{sample}/cutadapt/u_1_fw.fastq.gz",
        urv="{PROJECT}/samples/{sample}/cutadapt/u_1_rv.fastq.gz" if config["LIBRARY_LAYOUT"] == "PE" else []
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/demux_unassigned.benchmark"
    params:
        odir="{PROJECT}/samples/{sample}/cutadapt/",
        fix_opt="  --pair-adapters",
        cpus=" -j " + str(config["demultiplexing"]["cpus"])
    output:
        log="{PROJECT}/samples/{sample}/cutadapt/cutadapt_s2.log",
        ufw="{PROJECT}/samples/{sample}/demultiplexed/unassigned_R1.fastq.gz",
        urv="{PROJECT}/samples/{sample}/demultiplexed/unassigned_R2.fastq.gz"  if config["LIBRARY_LAYOUT"] == "PE" else []
    shell:
        "cutadapt -e " + str(config["demultiplexing"]["bc_mismatch"]) + "{params.fix_opt} {params.cpus}" +
        " -g ^file:{input.fw_bc} -G ^file:{input.rv_bc} " +
        " -o  {params.odir}{{name}}_fw.fastq.gz -p {params.odir}{{name}}_rv.fastq.gz " +
        " --untrimmed-output {output.ufw} --untrimmed-paired-output {output.urv} " + str(config["demultiplexing"]["cutadapt_extra_params"])+
        " {input.urv} {input.ufw} > {output.log}"  if config["LIBRARY_LAYOUT"] == "PE" else 
        "cutadapt -e " + str(config["demultiplexing"]["bc_mismatch"]) + " {params.cpus}" +
        " -g ^file:{input.fw_bc} " +
        " -o  {params.odir}{{name}}_fw.fastq.gz " +
        " --untrimmed-output {output.ufw}  --revcomp " + str(config["demultiplexing"]["cutadapt_extra_params"]) +
        " {input.ufw}> {output.log}"

rule concat_demux:
    """
    This rule takes the fastq files from both cutadapt runs and concatenates them into one single file.
    At the end this rule produces: "{PROJECT}/samples/{sample}/demultiplexed/{samp}.fastq.gz 
    """
    input:
        "{PROJECT}/samples/{sample}/cutadapt/cutadapt_s2.log"
    params:
        indir="{PROJECT}/samples/{sample}/cutadapt/",
        outdir="{PROJECT}/samples/{sample}/demultiplexed/"
    output:
        sum="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
        #wf="{PROJECT}/samples/{sample}/demultiplexed/summary_wf.txt"  
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/concat_demux.benchmark"
    shell:
        "Scripts/concat_demux.sh {params.indir} {params.outdir} > {output.sum}" if config["LIBRARY_LAYOUT"] == "PE" else 
        "EVALUATE HERE- Rev comp 2nd results and then concat: TBD"


rule tag_jump_analysis:
    input:
        mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
        fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
        rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" and config["LIBRARY_LAYOUT"] == "PE" 
        else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["gzip_input"].upper() == "T" and config["LIBRARY_LAYOUT"] == "PE"
        else []
    output:
        "{PROJECT}/samples/{sample}/barcodes/sample_matrix.tsv"
    shell:
        "{config[java][command]} -jar Scripts/BarcodeCorrector.jar  -b {input.mapp} -rf <(zcat -f {input.fw}) -rr <(zcat -f {input.rv}) -m "+ str(config["demultiplexing"]["bc_mismatch"]) +
        " -x {output} -l "+ str(config["demultiplexing"]["barcode_length"])+ " "  + str(config["demultiplexing"]["tagjump_extra_params"])
        if config["gzip_input"] == "T" else   
        "{config[java][command]} -jar Scripts/BarcodeCorrector.jar  -b {input.mapp} -rf  {input.fw} -rr {input.rv} -m "+ str(config["demultiplexing"]["bc_mismatch"])+
        " -x {output} -l "+ str(config["demultiplexing"]["barcode_length"])+ " " + str(config["demultiplexing"]["tagjump_extra_params"])

rule create_yaml_sample_matrix_tbl:
        input:
            "{PROJECT}/samples/{sample}/barcodes/sample_matrix.tsv"
        output:
            "{PROJECT}/samples/{sample}/barcodes/matrix.yaml"
        shell:
            "cat resources/datavzrd/matrix_base.yaml {input} | awk -F'\\t' 'NR<13{{print $0}} "+
            "NR==13{{for(i=2;i<=NF;i++){{print \"        \"$i\":\\n"+
            "          plot:\\n            heatmap:\\n              scale: linear\\n"+
            "              color-scheme: greens\"}}}}' >{output}"

rule datavzrd_sample_matrix:
        input:
            config="{PROJECT}/samples/{sample}/barcodes/matrix.yaml",
            table="{PROJECT}/samples/{sample}/barcodes/sample_matrix.tsv"
        params:
            num="3. " if config["ANALYSIS_TYPE"] == "ASV" else "4. "
        output:
            report(
                directory("{PROJECT}/samples/{sample}/barcodes/tagjump/"),
                htmlindex="index.html",
                category="3. Demultiplex",
                subcategory="Tag jump analysis",
                labels={"sample": "{sample}","Type":"matrix"}, 
            ) if config["ANALYSIS_TYPE"] == "ASV" else 
            report(
                directory("{PROJECT}/samples/{sample}/barcodes/tagjump/"),
                htmlindex="index.html",
                category="4. Demultiplex",
                subcategory="Tag jump analysis",
                labels={"sample": "{sample}","Type":"matrix"}, 
            ),
        wrapper:
            "v4.7.2/utils/datavzrd" 
rule create_heat_map:
    '''
    TBD
    This rule creates a heatmap, but first we need to create the sample matix, for this we first
    need to run cutadapt in a different way, when we create the barcode files, we need to give a different name
    to the reverse barcode and then it will try to create all the combinations, then we will need to process
    the results to make correct pairs and summarize results. 
    '''
    input:
        matrix="{PROJECT}/samples/{sample}/barcodes/sample_matrix.tsv",
        tag="{PROJECT}/samples/{sample}/barcodes/tagjump/" 
    output:
        report("{PROJECT}/samples/{sample}/barcodes/heat_map.png", category="3. Demultiplex", subcategory="Tag jump analysis",labels={"sample": "{sample}","Type":"plot"},)
        if config["ANALYSIS_TYPE"] == "ASV" else
        report("{PROJECT}/samples/{sample}/barcodes/heat_map.png", category="4. Demultiplex", subcategory="Tag jump analysis",labels={"sample": "{sample}","Type":"plot"},)

    params:
        ghmap="{PROJECT}/samples/{sample}/barcodes/heat_map_golay.png"
        #num="3. " if config["ANALYSIS_TYPE"] == "ASV" else "4. "
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/create_heat_map.benchmark"
    shell:
        "Rscript Scripts/heatMapDemux.R $PWD {input.matrix} {output} {params.ghmap}"


rule validate_demux:
    """
    This rule validates the demultiplexing procss 
    """
    input:
        raw="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz", 
        demux_log="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt",
        heatmap="{PROJECT}/samples/{sample}/barcodes/heat_map.png" if config["demultiplexing"]["create_tag_jump_analysis"].upper() == "T" else []
    params:
        outdir="{PROJECT}/samples/{sample}/demultiplexed/",
        interactive="--interactive " if config["interactive"].upper() == "T" else " ",
        tagjump="--tagjump " if config["demultiplexing"]["create_tag_jump_analysis"].upper() == "T" else " "
    output:
        summary="{PROJECT}/samples/{sample}/demultiplexed/summary_demux.tsv", 
        log="{PROJECT}/samples/{sample}/demultiplexed/demux_validation.log"              
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/validate_demux.benchmark"
    shell:
        "python Scripts/validateDemultiplex.py -r {input.raw} -o {output.summary} -l {output.log} {params.interactive} {params.tagjump} {params.outdir} "

rule create_yaml_demux_reads_tbl:
    output:
        "{PROJECT}/samples/{sample}/demultiplexed/demultiplex.yaml"
    shell:
        "cp resources/datavzrd/demultiplex.yaml {output}"

rule datavzrd_demux_reads:
    input:
        config="{PROJECT}/samples/{sample}/demultiplexed/demultiplex.yaml",
        table="{PROJECT}/samples/{sample}/demultiplexed/summary_demux.tsv"
    
    output:
        report(
            directory("{PROJECT}/samples/{sample}/demultiplexed/results/"),
            htmlindex="index.html",
            category="3. Demultiplex",
            subcategory="Sample distribution",
            labels={"sample": "{sample}", "Description":"Demultiplex summary table"},
        ) if config["ANALYSIS_TYPE"] == "ASV" else
        report(
            directory("{PROJECT}/samples/{sample}/demultiplexed/results/"),
            htmlindex="index.html",
            category="4. Demultiplex",
            subcategory="Sample distribution",
            labels={"sample": "{sample}", "Description":"Demultiplex summary table"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"

if config["demultiplexing"]["demultiplex"].lower() != "t":
    rule skip_demux_reads:
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" and config["LIBRARY_LAYOUT"] == "PE" else  
            "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz" if config["gzip_input"] == "T" and config["LIBRARY_LAYOUT"] == "PE" else [],
            tmp1="{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
            #tmp1="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
            #tmp2="{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" and config["LIBRARY_LAYOUT"] == "PE" else [],
            tmp3="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else []
        params:
            dir="{PROJECT}/samples/{sample}/demultiplexed/",
            ext=".fastq" if config["gzip_input"].lower() == "f" else ".fastq.gz"
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/summary_demux.tsv"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/skip_demux_reads.benchmark"
        shell:
            "echo 'No read demultiplexing' > {output} "
            "ln -s $PWD/{input.r1} {params.dir}{wildcards.sample}_1{ext}"
            "ln -s $PWD/{input.r2} {params.dir}{wildcards.sample}_2{ext}"
            if config["LIBRARY_LAYOUT"] == "PE" else
            "echo 'No read demultiplexing' > {output} "
            "ln -s $PWD/{input.r1} {params.dir}{wildcards.sample}_1{ext}" 
            
rule pear_reads:
    """
    Pair/extend paired-end data.
    """
    input:
        "{PROJECT}/samples/{sample}/demultiplexed/summary_demux.tsv",
        "{PROJECT}/samples/{sample}/demultiplexed/results/"            
    output:
        "{PROJECT}/samples/{sample}/peared/pear.log"
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/pear_reads.benchmark"
    params:
        #"{PROJECT}/samples/{sample}/peared/seqs",
        in_dir="{PROJECT}/samples/{sample}/demultiplexed/",
        out_dir="{PROJECT}/samples/{sample}/peared/"
    shell:
        "Scripts/assembly_files.sh {params.in_dir} {config[pear][t]} {config[pear][v]}  {config[pear][p]} {config[pear][j]} {config[pear][q]} {params.out_dir}  {config[pear][extra_params]} > {output} 2>&1"

rule validate_read_assembly:
    input:
        "{PROJECT}/samples/{sample}/peared/pear.log"
    output:
        summary="{PROJECT}/samples/{sample}/peared/pear_validation.tsv", #Per sample log file
        log="{PROJECT}/samples/{sample}/peared/pear_validation.out"  #All totals
    params:
        indir="{PROJECT}/samples/{sample}/peared/",
        interactive="--interactive" if config["interactive"].upper() == "T" else " "            
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/validate_read_assembly.benchmark"
    shell:
        "python Scripts/validatePearLog.py {params.interactive}  -m {config[pear][prcpear]}  -o {output.summary} -l {output.log} {params.indir} "

rule create_yaml_peared_reads_tbl:
    output:
        "{PROJECT}/samples/{sample}/peared/peared.yaml"
    shell:
        "cp resources/datavzrd/pear.yaml {output}"

rule datavzrd_peared_reads:
    input:
        config="{PROJECT}/samples/{sample}/peared/peared.yaml",
        table="{PROJECT}/samples/{sample}/peared/pear_validation.tsv"
    output:
        report(
            directory("{PROJECT}/samples/{sample}/peared/results/"),
            htmlindex="index.html",
            category="3. Assembled reads",
            labels={"Pear results": "{sample}"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"
    
rule quality_assembled_reads:
    input: 
        "{PROJECT}/samples/{sample}/peared/pear_validation.tsv",
        "{PROJECT}/samples/{sample}/peared/results/"
    output:
        out="{PROJECT}/samples/{sample}/peared/pear_quality.out",
        log="{PROJECT}/samples/{sample}/peared/pear_quality.log",
        rfq=report(directory("{PROJECT}/samples/{sample}/peared/fastqc/"), 
            patterns=["{name}.assembled_fastqc.html"], 
            category="2. Quality", subcategory="Assembled reads",
            labels={"sample":"{sample}","Tool":"FastQC","Reads":"{name}"},
            caption="resources/caption/fastqc.pear.rst") if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
        rsq=report(directory("{PROJECT}/samples/{sample}/peared/sequali/"), 
            patterns=["sequali.{name}.html"], 
            category="2. Quality", subcategory="Assembled reads",
            labels={"sample":"{sample}","Tool":"Sequali","Reads":"{name}"},
            caption="resources/caption/sequali.pear.rst") if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else []
        #qc= "" "{PROJECT}/samples/{sample}/qc/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else []
    params:
        fastqs="{PROJECT}/samples/{sample}/peared/",
        fo="{PROJECT}/samples/{sample}/peared/fastqc",
        so="{PROJECT}/samples/{sample}/peared/sequali"
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/quality_assembled_reads.benchmark"
    shell:
        "python Scripts/FastQC_assembled_reads.py -fo {params.fo} -so {params.so} -cf  {config[fastQC][threads]} -cs  {config[sequali][threads]} -o {output.out} {params.fastqs} > {output.log} 2>&1"

rule create_assigned_reads_unpaired:
    """
    This rule gzip all the demultiplexed unpaired files place them at:
    demultiplexed/unpaired/{id}_1|2.fastq.gz
    """
    input:
        "{PROJECT}/samples/{sample}/peared/pear_validation.out" if config["fastQCPear"].lower() != "t" else
        "{PROJECT}/samples/{sample}/peared/pear_quality.out"   
    output: 
        log="{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
        #think about using this OUTPUT (el de abajo)
        #splitLog="{PROJECT}/samples/{sample}/split_library_log.txt"
    params:
        in_dir="{PROJECT}/samples/{sample}/peared/",
        out_dir="{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
    benchmark:
        "{PROJECT}/samples/{sample}/benchmark/create_assigned_reads_unpaired.benchmark"
    shell:
        "Scripts/gzip_unpaired.sh {params.in_dir} {params.out_dir} {output.log}"


 ######OTU#################
 #Here we continue with mainly OTU workflow rules
  ##################################
if config["LIBRARY_LAYOUT"].lower() == "pe":
    rule fastqToFasta:
        """
        This rule take all the demultiplexed paired files and create one single fasta
        file and renames the headers according to the ID.
        """
        input:
            "{PROJECT}/samples/{sample}/peared/pear_validation.out" if config["fastQCPear"].lower() != "t" else
            "{PROJECT}/samples/{sample}/peared/pear_quality.out" 
            #expand("{PROJECT}/samples/{sample}/peared/pear_validation.out", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run) if config["fastQCPear"].lower() != "t" else
            #expand("{PROJECT}/samples/{sample}/peared/pear_quality.out", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        output:
            seqs="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna", 
            #seqs="{PROJECT}/samples/{sample}/seqs.fna", #antes "{PROJECT}/samples/{sample}/splitLibs/seqs.fna"
            seqs_un="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna",
            #seqs_un="{PROJECT}/samples/{sample}/seqs.unassembled.fna",
            #splitLog="{PROJECT}/samples/{sample}/split_library_log.txt" #this was commented for the error:
            #Not all output files of rule fastqToFasta contain the same wildcards.
            splitLog="{PROJECT}/runs/{run}/{sample}_data/split_library_log.txt"
        params:
            #"{PROJECT}/samples/{sample}/peared/seqs",
            in_dir="{PROJECT}/samples/{sample}/peared/"
            #in_dir="dfsdf"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/fastqToFasta.{sample}.benchmark"
        shell:
            "Scripts/convert_2_fasta.sh {params.in_dir} {output.splitLog} {output.seqs} {output.seqs_un} {config[CHAR_TO_PAIR]} {input}" 
        
       
else: #SE   
#HERE WE NEED TO ADD SOME QUALITY FILTERING!    
    rule skip_demultiplexing_fq2fasta_OTU_SE:
            input:
                raw="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                tmp="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else [],
                tmp2="{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else []
            output:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
            params:
                command="zcat " if config["gzip_input"] == "T" else "cat "
            benchmark:
                 "{PROJECT}/runs/{run}/benchmark/skip_demultiplexing_fq2fasta_OTU_SE.{sample}.benchmark"
            shell:
                #to change vsearch --fastq_filter N_debres_r1.txt -fastaout test.txt
                #"fq2fa {input.extended_reads} {output}"
                #"sed -n '1~4s/^@/>/p;2~4p' {input.extended_reads} >  {output}"
                "{params.command} {input.extended_reads} | sed -n '1~4s/^@/>/p;2~4p'  | "
                "awk  '{{if($0 ~ \"^>\"){{seq=seq+1;print \">{wildcards.sample}_\"seq\" \"substr($1,2)}}else{{print $0}}}}' > {output}"
                #if config["gzip_input"] == "T" and config["LIBRARY_LAYOUT"] == "SE" else
                #"sed -n '1~4s/^@/>/p;2~4p' {input.extended_reads} | "
                #"awk  '{{if($0 ~ \"^>\"){{seq=seq+1;print \">{wildcards.sample}_\"seq\" \"substr($1,2)}}else{{print $0}}}}' > {output}"
     
#####################################################
#                   ASV Rules                       #
#  Here we continue with only ASV workflow rules    # 
##################################################### 
if config["ANALYSIS_TYPE"] == "ASV":
    if config["demultiplexing"]["demultiplex"].lower() == "f" and config["LIBRARY_LAYOUT"] != "SE" and config["UNPAIRED_DATA_PIPELINE"] != "T":
        rule skip_dmx_file_creation_ASV_PE:
            input:
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
                tmp="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else [],
                tmp2="{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else []

            params:
                "{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq"
                if config["gzip_input"] == "F" else
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq.gz "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq.gz"

 
    elif config["demultiplexing"]["demultiplex"].lower() == "f" and config["LIBRARY_LAYOUT"] == "SE":
        rule skip_dmx_file_creation_ASV_SE:
            input:
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                #tmp2="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
                #tmp3="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else [],
                tmp1="{PROJECT}/samples/{sample}/qc/fq_internal_validation.txt" if config["QC"].lower() == "fastqc" or config["QC"].lower() == "both" else [],
                tmp3="{PROJECT}/samples/{sample}/qc/sequali/sequali.html" if config["QC"].lower() == "sequali" or config["QC"].lower() == "both" else []
                
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
                #if config["primers"]["remove"].lower() != "f" else
                #"{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                if config["gzip_input"].lower() == "f" else
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq.gz"

    elif config["demultiplexing"]["demultiplex"].lower() == "f" and config["LIBRARY_LAYOUT"] != "SE" and config["UNPAIRED_DATA_PIPELINE"] == "T":
         rule skip_dmx_file_creation_ASV_Unpaired:
            input:
                r1="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.fastq",
                r2="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.fastq" 
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
                #if config["primers"]["remove"].lower() != "f" else
                #"{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq "
# -------
# - Primer removal section.
# When we demultiplex, we search for barcode combination in both sides, and with that we 
# place the reads in correct orientation in the case they aren't, so when we don't demultiplex
# we still need to do two passes of cutadapt (rev comp), in this sence there are two blocks of
# rules, remove primers for demultiplexed files by Cascabel and the second block for reads
# demux by other tools. 
# ------- 
    if config["demultiplexing"]["demultiplex"].lower() == "t":
        if config["primers"]["remove"].lower() == "cfg":
            rule remove_primers_dmx_files_cfg_asv:
                input:
                #fw=expand("{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run),
                #rv=expand("{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run) #if int(config["demultiplexing"]["bc_mismatch"]) < 1
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
                params:
                    "{PROJECT}/samples/{sample}/demultiplexed" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired" ,
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    config["LIBRARY_LAYOUT"].upper(),
                    "{PROJECT}/runs/"+ run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/remove_primers_dmx_files_cfg_asv.benchmark" 
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                script:
                    "Scripts/removePrimersDemultiplex_cfg.py"
                    
        elif config["primers"]["remove"].lower() == "metadata":
            rule remove_primers_dmx_files_metadata_asv:
                input:
                    config="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                    map_log="{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log", #if we choose to remove primers and no demux there was no mapping file validation, now there is 
                    rep_bc="{PROJECT}/metadata/{sample}_report/"
                    #rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    #else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
                params:
                    "{PROJECT}/samples/{sample}/demultiplexed"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired",
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    config["LIBRARY_LAYOUT"].upper(),
                    "{PROJECT}/runs/"+ run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/remove_primers_dmx_files_metadata_asv.benchmark"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                script:
                    "Scripts/removePrimersDemultiplex.py"
        else:
            rule skip_remove_primers_dmx_files_asv:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/skip_remove_primers_dmx_files.benchmark"
                #here we only need to symlink a wk output file the next step for the ASV WF (asvFilter.R) has an option to see if we 
                #removed primers or not and then take the reads either from the demx or the demx/primer_removed directory.
                shell:
                    "ln -S {input} {output}"
  
    else: #we start from raw fastq files from user.
        if config["primers"]["remove"].lower() == "cfg":
            rule remove_primers_fq_files_cfg_no_demux_asv:
                '''
                The script on this rule reverse complement the reads with no primers on the first round.
                '''
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
                params:                
                    "{PROJECT}/samples/{sample}/demultiplexed"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired",
                    config["primers"]["extra_params"],
                    #"fastq"  if config["gzip_input"] == "F" or (config["ANALYSIS_TYPE"] == "ASV" and config["UNPAIRED_DATA_PIPELINE"] == "T") else "fastq.gz",
                    "fastq"  if config["gzip_input"] == "F" else "fastq.gz",
                    config["LIBRARY_LAYOUT"],
                    "{PROJECT}/runs/"+run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/remove_primers_fq_files_cfg_no_demux_asv.benchmark"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                script:
                    "Scripts/removePrimersNoCascabelDemultiplex_cfg.py"
        elif config["primers"]["remove"].lower() == "metadata":
            rule remove_primers_fq_files_metadata_asv:
                input:
                    config="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt" if config["demultiplexing"]["demultiplex"] == "T"
                    else config["metadata"],
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                    map_log="{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log", # inforce mapping file validation and report gen.
                    rep_bc="{PROJECT}/metadata/{sample}_report/"

                params:
                    "{PROJECT}/samples/{sample}/demultiplexed"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired",
                    config["primers"]["extra_params"],
                    #"fastq"  if config["gzip_input"] == "F"  or (config["ANALYSIS_TYPE"] == "ASV" and config["UNPAIRED_DATA_PIPELINE"] == "T")  else "fastq.gz",
                    "fastq"  if config["gzip_input"] == "F"  else "fastq.gz",
                    config["LIBRARY_LAYOUT"],
                    "{PROJECT}/runs/"+run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/remove_primers_fq_files_metadata_asv.benchmark"
         
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                script:
                    "Scripts/removePrimersDemultiplex.py"
        else:#este es nuevo y hay q ver!! nadie produce los inputs de arriba
            rule skip_remove_primers_dmx_files_asv:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
                benchmark:
                    "{PROJECT}/samples/{sample}/benchmark/skip_remove_primers_dmx_files.benchmark"
                #here we only need to symlink a wk output file the next step for the ASV WF (asvFilter.R) has an option to see if we 
                #removed primers or not and then take the reads either from the demx or the demx/primer_removed directory.
                shell:
                    "ln -S {input} {output}"

    rule create_yaml_primer_reads_tbl:
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/primer.yaml"
        shell:
            "cp resources/datavzrd/primer.yaml {output}"

    rule datavzrd_primer_reads:
        input:
            config="{PROJECT}/samples/{sample}/demultiplexed/primer.yaml",
            table="{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
        output:
            report(
                directory("{PROJECT}/samples/{sample}/demultiplexed/primer_results/"),
                htmlindex="index.html",
                category="4. Remove primers",
                labels={"Cutadat results": "Primer removed: {sample}"},
            ),
        wrapper:
            "v4.7.2/utils/datavzrd" 
#------ End Primer removal section --------------#                   
    rule dada2Filter:
        input:
            fq_files=expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
            p_report=expand("{PROJECT}/samples/{sample}/demultiplexed/primer_results/" if config["primers"]["remove"].lower() == "metadata" or config["primers"]["remove"].lower() == "cfg"
            else [] , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
            d_report=expand("{PROJECT}/samples/{sample}/demultiplexed/results/" if config["demultiplexing"]["demultiplex"].lower() == "t"
            else [], PROJECT=config["PROJECT"],sample=config["LIBRARY"])            
        output:
            summary="{PROJECT}/runs/{run}/asv/filter_summary.out",
            r=report("{PROJECT}/runs/{run}/asv/QA_plots.pdf", category="2. Quality",subcategory="Dada2 plots",labels={"Plot":"QA Plots", "Strand":"FW & RV"},) if config["dada2_filter"]["generateQAplots"].lower() == "t"
            else []
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/dada2Filter.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/asvFilter.R $PWD " + str(config["dada2_filter"]["generateQAplots"]) + " " + str(config["dada2_filter"]["truncFW"]) + " " + str(config["dada2_filter"]["truncRV"]) + " "+str(config["dada2_filter"]["maxEE_FW"]) + " "+str(config["dada2_filter"]["maxEE_RV"]) + " " +str(config["dada2_filter"]["cpus"]) + " \"" +str(config["dada2_filter"]["extra_params"]) + "\" " +str(config["interactive"])+ " {output.summary} " +config["primers"]["remove"] +" {input} " 

    rule validate_dada2Filter:
        input:
            "{PROJECT}/runs/{run}/asv/filter_summary.out"
        output:
            "{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/validate_dada2Filter.benchmark"
        script:
            "Scripts/validateFilterASV.py"
    rule subset_for_errors_fw:
        input:
            expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
             "{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        output:
            temp("{PROJECT}/runs/{run}/asv/subset.fastq") 
        params:
            dir="{PROJECT}/runs/{run}/asv/"
            #ssize=config["dada2_asv"]["subset_size"]
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/subset_for_errors_fw.benchmark"
        shell:
            "cat test_out_demux/samples/NIOZ102/demultiplexed/filtered/*1.fastq.gz "
            "| gzip -d |  vsearch --fastx_subsample - --sample_size 200000  --fastaout test_out_demux/samples/NIOZ102/demultiplexed/filtered/subsample.fastq.gz"
             
    rule learn_errors:
        input:
            reads=expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
            flow="{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        output:
            "{PROJECT}/runs/{run}/asv/errors.RData",
            report("{PROJECT}/runs/{run}/asv/fw_err.pdf", category="2. Quality",subcategory="Dada2 plots",labels={"Plot":"Error plots", "Strand":"FW"}) if config["dada2_asv"]["generateErrPlots"].lower() == "t"
            else [],
            report("{PROJECT}/runs/{run}/asv/rv_err.pdf", category="2. Quality",subcategory="Dada2 plots",labels={"Plot":"Error plots", "Strand":"RV"}) if config["dada2_asv"]["generateErrPlots"].lower() == "t"
            else []
        params:
            dir="{PROJECT}/runs/{run}/asv/"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/learn_errors.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/errDada2.R $PWD " +str(config["dada2_asv"]["cpus"]) + " "
            " "+ str(config["dada2_asv"]["generateErrPlots"])  + " {params.dir} "  + " "
            " "+ str(config["dada2_asv"]["nbases"]) +" " + str(config["binned_q_scores"])  + " "
            + "{input.reads}"
    

    rule run_dada2:
        input:
            reads=expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
            wf="{PROJECT}/runs/{run}/asv/filter_summary.validation.txt",
            err="{PROJECT}/runs/{run}/asv/errors.RData"
        output:
            "{PROJECT}/runs/{run}/asv/stats_dada2.tsv",
            "{PROJECT}/runs/{run}/asv/representative_seq_set.fasta",
            "{PROJECT}/runs/{run}/asv/representative_bimeras.fasta",
            temp("{PROJECT}/runs/{run}/asv/dada2_asv_table.txt"),
            temp("{PROJECT}/runs/{run}/asv/dada2_asv_bimera_table.txt") 
        params:
            dir="{PROJECT}/runs/{run}/asv/",
            script="Scripts/asvDada2_bd.R"  if config["big_data_wf"].lower()  == "t"
            else "Scripts/asvDada2.R"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/run_dada2.benchmark"
        shell:
            "{config[Rscript][command]} {params.script} $PWD " #1 
            " "+ str(config["dada2_asv"]["pool"]) + " " #2
            " "+ str(config["dada2_asv"]["cpus"])  + " " #3
            " \""+ str(config["dada2_asv"]["extra_params"]) + "\" " #5 - #4
            " "+ str(config["dada2_asv"]["chimeras"]) + " " #10 - #5
            " {params.dir} " #6 - #6 (move chimera - 10 to here)  
            " "+ str(config["rm_reads"]["shorts"]) + " "#7
            " "+ str(config["rm_reads"]["longs"]) + " " #8
            " "+ str(config["rm_reads"]["offset"])  + " "#9
            # " "+ str(config["dada2_asv"]["chimeras"]) + " " #10 - #5
            #" "+ str(config["dada2_taxonomy"]["db"]) + " " #11
            #" "+ str(config["dada2_taxonomy"]["add_sps"]["db_sps"])  + " "#12
            #" "+ str(config["dada2_taxonomy"]["add_sps"]["add"]) + " "  #13
            #" \""+ str(config["dada2_taxonomy"]["extra_params"]) + "\" " #14
            " "+ str(config["dada2_merge"]["minOverlap"]) +" " #15 - 10
            " "+ str(config["dada2_merge"]["maxMismatch"]) +" "#16 -11
            " "+ str(config["UNPAIRED_DATA_PIPELINE"]) +" " #17 - 12
            #" \""+ str(config["dada2_taxonomy"]["add_sps"]["extra_params"]) + "\" " #18
            " "+ str(config["interactive"]) + " " #19 - 13
            " "+ str(config["rm_reads"]["non_interactive_behaviour"]) + " " #20 - 14
            " "+ str(config["dada2_asv"]["chimeras_method"]) + " " #15  
            "{input.reads}" #21 - 16

    rule create_yaml_dada2_tbl:
        output:
            "{PROJECT}/runs/{run}/asv/dada2_stats.yaml"
        shell:
            "cp resources/datavzrd/dada2_stats.yaml {output}"

    rule datavzrd_dada2_stats:
        input:
            config="{PROJECT}/runs/{run}/asv/dada2_stats.yaml",
            table="{PROJECT}/runs/{run}/asv/stats_dada2.tsv"
        output:
            report(
                directory("{PROJECT}/runs/{run}/asv/dada2_stats/"),
                htmlindex="index.html",
                category="5. ASV results",
                subcategory="Dada2 stats",
                labels={"table": "Dada2 statistics"},
            ),
        wrapper:
            "v4.7.2/utils/datavzrd" 
    rule assignTaxonomyDada2:
        input:
            rep_seqs="{PROJECT}/runs/{run}/asv/representative_seq_set.fasta",
            dd2_stats="{PROJECT}/runs/{run}/asv/dada2_stats/"
        output:
            report("{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt", category="6. Taxonomy", subcategory="1. Taxonomy assignation results",labels={"Representative sequences":"ASV rep. seqs"},),
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.bootstrap.txt"
        params:
            dir="{PROJECT}/runs/{run}/asv/"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/assignTaxonomyDada2.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/dada2AssignTaxo.R $PWD "
            "{input.rep_seqs} " +" "
            " \"" + str(config["dada2_taxonomy"]["extra_params"]) + "\" "
            " \"" + str(config["dada2_taxonomy"]["add_sps"]["extra_params"]) + "\" "
            " " + str(config["dada2_taxonomy"]["seed"])+" "
            " " + str(config["dada2_taxonomy"]["db"]) + " "
            " " + str(config["dada2_taxonomy"]["add_sps"]["db_sps"])  + " "
            " " + str(config["dada2_taxonomy"]["add_sps"]["add"]) + " "
            "{output[0]} "
            "{output[1]} "
            " "+ str(config["dada2_asv"]["cpus"]) + " "
            " asv. "
            
    if config["dada2_asv"]["chimeras"].lower() == "t" and config["dada2_asv"]["chimeras_taxonomy"].lower() == "t":
        rule assignTaxonomyBimeras:
            input:
                "{PROJECT}/runs/{run}/asv/representative_bimeras.fasta"
            output:
                "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeras_tax_assignments.txt",
                "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeras_tax_assignments.bootstrap.txt"
            params:
                dir="{PROJECT}/runs/{run}/asv/"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/assignTaxonomyBimeras.benchmark"
            shell:
                "{config[Rscript][command]} Scripts/dada2AssignTaxo.R $PWD "
                "{input} " +" "
                " \"" + str(config["dada2_taxonomy"]["extra_params"]) + "\" "
                " \"" + str(config["dada2_taxonomy"]["add_sps"]["extra_params"]) + "\" "
                " " + str(config["dada2_taxonomy"]["seed"])+" "
                " " + str(config["dada2_taxonomy"]["db"]) + " "
                " " + str(config["dada2_taxonomy"]["add_sps"]["db_sps"])  + " "
                " " + str(config["dada2_taxonomy"]["add_sps"]["add"]) + " "
                "{output[0]} "
                "{output[1]} "
                " "+ str(config["dada2_asv"]["cpus"]) + " "
                " bimera. " 
    else: 
        rule skipAssignTaxonomyBimeras:
            output:
                "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeras_tax_assignments.txt",
                "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeras_tax_assignments.bootstrap.txt"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/skipAssignTaxonomyBimeras.benchmark"
            shell:
                "touch {output}" 

    rule asv_table:
        input:
            "{PROJECT}/runs/{run}/asv/dada2_asv_table.txt" 
        output:
            "{PROJECT}/runs/{run}/asv/asv_table.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/asv_table.benchmark"
        shell:
            "cat {input} | awk '{{if(NR==1){{header=\"#OTUID\";for(i=1;i<=NF;i++){{header=header\"\\t\"$i}};print header}}else{{print $0}}}}'|   awk '{{ for (i=1; i<=NF; i++){{ a[NR,i] = $i }} }} NF>p {{ p = NF }} END {{ for(j=1; j<=p; j++) {{ str=a[1,j]; for(i=2; i<=NR; i++){{ str=str\"\\t\"a[i,j]; }} print str }} }}' > {output}"

    rule asv_tax_table:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
            "{PROJECT}/runs/{run}/asv/asv_table.txt" 
        output:
            report("{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt", category="5. ASV results", subcategory="ASV length filtered", labels={"File": "ASV table"},)

        benchmark:
            "{PROJECT}/runs/{run}/benchmark/asv_tax_table.benchmark"
        shell:
            "cat {input[0]} | awk  -F \"\\t\" 'NR==FNR{{if(NR>1){{tax=$2;for(i=3;i<=NF;i++){{tax=tax\";\"$i}};h[$1]=tax;}}next;}} {{if(FNR==1){{print $0\"\\ttaxonomy\"}}else{{print $0\"\\t\"h[$1]}}}}' -  {input[1]} > {output}" 
    rule bimera_table:
        input:
            "{PROJECT}/runs/{run}/asv/dada2_asv_bimera_table.txt"
        output:
            "{PROJECT}/runs/{run}/asv/bimera_table.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/bimera_table.benchmark"
        shell:
            "cat {input} | awk '{{if(NR==1){{header=\"#OTUID\";for(i=1;i<=NF;i++){{header=header\"\\t\"$i}};print header}}else{{print $0}}}}'|   awk '{{ for (i=1; i<=NF; i++){{ a[NR,i] = $i }} }} NF>p {{ p = NF }} END {{ for(j=1; j<=p; j++) {{ str=a[1,j]; for(i=2; i<=NR; i++){{ str=str\"\\t\"a[i,j]; }} print str }} }}' > {output}"

    rule bimera_tax_table:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeras_tax_assignments.txt",
            "{PROJECT}/runs/{run}/asv/bimera_table.txt"
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeraTable.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/bimera_tax_table.benchmark"
        shell:
            "cat {input[0]} | awk  -F \"\\t\" 'NR==FNR{{if(NR>1){{tax=$2;for(i=3;i<=NF;i++){{tax=tax\";\"$i}};h[$1]=tax;}}next;}} {{if(FNR==1){{print $0\"\\ttaxonomy\"}}else{{print $0\"\\t\"h[$1]}}}}' -  {input[1]} > {output}"

    rule asv_to_biom:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt" 
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/asv_to_biom.benchmark"
        shell:
            "{config[biom][command]} convert -i {input[0]} -o {output} --table-type \"OTU table\" --to-hdf5 --process-obs-metadata taxonomy "

    rule filter_asv:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt"
        output:
            report("{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt", category="5. ASV results", subcategory="ASV no singletons", labels={"File": "ASV table no singletons"},)
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/filter_asv.benchmark"
        params: 
            min_abundance=config["filterOtu"]["n"],
            min_sample=config["filterOtu"]["s"]
        shell:
            "cat {input} | " + 
            "awk  -v n={params.min_abundance} -v s={params.min_sample} -F\"\\t\" 'NR==1{{print $0;tax=-1; "+
            "for(i=1;i<=NF;i++){{if($i==\"taxonomy\"){{tax=i}}}}}} "+
            "NR>1{{samples=0;total=0; "+
            "for(i=2;i<tax;i++){{ " +
            "  total+=$i;" +
            "  if($i>0){{samples++}}}}; "+
            "if(samples>s && total>n){{print $0}}}}' > {output}"
    
    #Convert to/from the BIOM table format.
    rule convert_filter_asv:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.biom"
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/convert_filter_asv.benchmark"
        shell:
            "{config[biom][command]} convert -i {input} -o {output} {config[biom][tableType]} "
            "{config[biom][headerKey]} {config[biom][outFormat]} {config[biom][extra_params]}"
             

else: #OTU Workflow
    if config["align_vs_reference"]["align"] == "T":
        rule align_vs_reference:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna"
            output:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align"
            params:
                "{PROJECT}/runs/{run}/{sample}_data/"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/align_vs_reference.benchmark"
            shell:
               mapp= "cd {params} && "
                "{config[align_vs_reference][mothur_cmd]} '#align.seqs(fasta=seqs_fw_rev_accepted.fna, reference={config[align_vs_reference][dbAligned]}, processors={config[align_vs_reference][cpus]})'"
    
    if config["primers"]["remove"] != "F":
        """
        If we are going to remove primers/adapters and they do not came from our own
        demultiplexing, it is likely to have the sequences in both direction, FW and RV
        thus, we need to rev com the sequences, concatenate them and then when we run cutadapt
        always in this step, discard-untrimmed so this way we end up with the sequences in the correct orientation
        """
        if config["demultiplexing"]["demultiplex"] != "T" :
            rule rev_com_seqs:
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna" 
                output:
                    temp("{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.fna")
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/align_vs_reference.{sample}.benchmark"
                shell:
                    "vsearch --fastx_revcomp {input} --fastaout {output}"
            rule concat_seqs:
                input:
                    sq="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna",
                    rc="{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.fna"
                output:
                    temp("{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.concat.fna")
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/concat_seqs.{sample}.benchmark"
                shell:
                    "cat {input.sq} {input.rc} > {output}"
    
            rule remove_primers_cutadapt_OTU_wf:
                '''
                This rule remove primers for the OTU workflow from fasta files NO demultiplexed by Cascabel.
                For the moment it only works this way for metadata. Need to test how to implement it with 
                config file information.
                TODO: else on the 3rd input need to be fixed.                
                '''
                input:
                    #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                    #else 
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna",
                    config["metadata"]  if config["primers"]["remove"].lower() == "metadata"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.concat.fna"
                output:
                    out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna",
                    summ="{PROJECT}/runs/{run}/{sample}_data/cutadapt.{sample}.summary.tsv"
                    #log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/remove_primers_cutadapt_OTU_wf.{sample}.benchmark"
                params:
                    config["primers"]["extra_params"],
                    "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt",
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                    "{PROJECT}/runs/{run}/report_files/{sample}.summary.tsv_",
                    "{sample}",
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                script:
                    "Scripts/remove_adapters_by_sample.py"
    
        else:   #We demux
            if config["primers"]["remove"].lower() == "metadata":
                rule remove_primers_metadata_OTU:
                    input:
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                        else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                        else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna",
                        "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt" if config["demultiplexing"]["demultiplex"] == "T" 
                        else config["metadata"],
                        "{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log" if config["demultiplexing"]["demultiplex"] == "T"
                        else [],
                        "{PROJECT}/metadata/{sample}_report/" if config["demultiplexing"]["demultiplex"] == "T"
                        else [],
                    output:
                        out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna",
                        summ="{PROJECT}/runs/{run}/{sample}_data/cutadapt.{sample}.summary.tsv"
                       # log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                    benchmark:
                        "{PROJECT}/runs/{run}/benchmark/remove_primers_metadata_OTU.{sample}.benchmark"
                    params:
                        config["primers"]["extra_params"],
                        "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt",
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                        "{PROJECT}/runs/{run}/report_files/{sample}.summary.tsv_",
                        "{PROJECT}/runs/{run}/{sample}_data/cutadapt_tmp/",
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                    script:
                         "Scripts/remove_adapters_v2.py" # && ln -s ../../report_files/cutadapt.{wildcards.sample}.summary.tsv {params[4]}   "
            elif config["primers"]["remove"].lower() == "cfg":
                rule remove_primers_cfg_OTU:
                    input:
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                        else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["UNPAIRED_DATA_PIPELINE"].lower() == "f"
                        else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.unassembled.fna"
                        
                    output:
                        out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna",
                        summ="{PROJECT}/runs/{run}/{sample}_data/cutadapt.{sample}.summary.tsv"
                        #log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
    
                    benchmark:
                        "{PROJECT}/runs/{run}/benchmark/remove_primers_cfg_OTU.{sample}.benchmark"
                    params:
                        config["primers"]["extra_params"],
                        "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt", 
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                        "{PROJECT}/runs/{run}/report_files/{sample}.summary.tsv_",
                        "{PROJECT}/runs/{run}/{sample}_data/cutadapt_tmp/",
                        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                    script:
                        "Scripts/remove_adapters_v2.py"          
        rule create_yaml_primer_reads_tbl:
            output:
                "{PROJECT}/runs/{run}/{sample}_data/primer.yaml"
            shell:
                "cp resources/datavzrd/primer.yaml {output}"

        rule datavzrd_primer_reads:
            input:
                config="{PROJECT}/runs/{run}/{sample}_data/primer.yaml",
                table="{PROJECT}/runs/{run}/{sample}_data/cutadapt.{sample}.summary.tsv"
            output:
                report(
                    directory("{PROJECT}/runs/{run}/{sample}_data/primer_results/"),
                    htmlindex="index.html",
                    category="5. Remove primers",
                    subcategory="Primer removal stats",
                    labels={"Cutadapt results": "Primer removed: {sample}"},
                ),
            wrapper:
                "v4.7.2/utils/datavzrd" 
    
    if config["align_vs_reference"]["align"] == "T":
        rule degap_alignment:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["primers"]["remove"] == "F"
                else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"
            output:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna" if config["primers"]["remove"] == "F"
                else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/degap_alignment.{sample}.benchmark"
            shell:
                "degapseq {input} {output}"
    
        rule single_line_fasta:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna" if config["primers"]["remove"] == "F"
                else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna"
            output:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.oneline.fna" if config["primers"]["remove"] == "F"
                else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.oneline.fna"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/single_line_fasta.{sample}.benchmark"
            shell:
                "{config[java][command]} -cp Scripts FastaOneLine -f {input} -m 1 --write-discarded -o {output}"
    
    
    #Creates file with sequence length ditribution
    rule histogram:
        input:
            #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_nc.fna",
            #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" #for tmp
            *selectInput(),
            "{PROJECT}/runs/{run}/{sample}_data/primer_results/" if config["primers"]["remove"].lower() != "f" else []
        output:
            temp("{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_lengths.txt"),
            temp("{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_hist.txt")
        benchmark:
                "{PROJECT}/runs/{run}/benchmark/histogram.{sample}.benchmark"
        shell:
            "cat {input[0]} | grep -v '^>' | awk '{{print length}}' > {output[0]} "
            "&&  sort -g {output[0]} | uniq -c > {output[1]}"
    
    rule histogram_chart:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_hist.txt",
            "{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_lengths.txt"
        output:
            report("{PROJECT}/runs/{run}/report_files/seqs_dist_hist.{sample}.png", category="5. Remove primers", subcategory="Read length", labels={"Read length histogram":"Sample: {sample}"},),
            #"{PROJECT}/runs/{run}/report_files/seqs_dist_hist.{sample}.png",
            "{PROJECT}/runs/{run}/{sample}_data/seqs_statistics.txt"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/",
            "{PROJECT}/runs/{run}/report_files/",
            "{sample}"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/histogram_chart.{sample}.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/histogram.R $PWD {input[0]} {input[1]} {params[0]} {output[0]}"
    
    #remove to long and to short sequences
    rule remove_short_long_reads:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_statistics.txt",
            "{PROJECT}/runs/{run}/report_files/seqs_dist_hist.{sample}.png",
            "{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_hist.txt",
            *selectInput()
            #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_nc.fna"
        #params:
        #    "{PROJECT}/runs/{run}/{sample}_data/splitLibs"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",
            "{PROJECT}/runs/{run}/{sample}_data/filter.log"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/remove_short_long_reads.{sample}.benchmark"
        script:
            "Scripts/rmShortLong.py"
    if config["chimera"]["search"] == "T":
            rule search_chimera_vsearch:
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
                output:
                    "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.summary.txt"
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/search_chimera_vsearch.{sample}.benchmark"
                shell:
                    "vsearch --{config[chimera][method]} {input} --threads {config[chimera][threads]} {config[chimera][extra_params]} --uchimeout {output}"
            rule filter_chimera_vsearch:
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.summary.txt"
                output:
                    "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.txt"
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/filter_chimera_vsearch.{sample}.benchmark"
                shell:
                    "cat {input} | awk '$18==\"Y\"{{print $2\"\\t\"$1}}' > {output}"
            rule remove_chimera:           
                """
                This script counts the sequence in input[0], the chimeras in input[2] if it is non interactive it
                directly removes the chimeras, otherwise ask the user, if the selection is NO, it will rename
                input[0] as output[0]
                """
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",
                #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna" if config["align_vs_reference"]["align"] == "T" and config["cutAdapters"] == "T"
                #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna" elif config["align_vs_reference"]["align"] == "T" and config["cutAdapters"] != "T"
                #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna" elif config["align_vs_reference"]["align"] != "T" and config["cutAdapters"] == "T"
                #else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna",
                #*selectInput(),
                    "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.txt"
                output:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta",
                    "{PROJECT}/runs/{run}/{sample}_data/chimera/chimera.log"
                benchmark:
                    "{PROJECT}/runs/{run}/benchmark/remove_chimera.{sample}.benchmark"
                script:
                    "Scripts/remove_chimera.py"          
    rule combine_filtered_samples:
        input:
            allFiltered = expand("{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T" 
            else  "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        output:
            "{PROJECT}/runs/{run}/seqs_fw_rev_combined.fasta",
            "{PROJECT}/runs/{run}/samples.log"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/combine_filtered_samples.benchmark"
        script:
            "Scripts/combineAllReads.py"
    
    if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch") or config["pickOTU"]["m"] == "swarm" :
    #here make two steps
        rule dereplicate:
            input:
                "{PROJECT}/runs/{run}/seqs_fw_rev_combined.fasta"
            output:
                #"{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_otus.txt", <- picking OTUs method
                "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.fasta",
                "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.uc"
            params:
                "{PROJECT}/runs/{run}/derep/"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/dereplicate.benchmark"
            shell:
                "{config[derep][vsearch_cmd]} --derep_fulllength {input} --output {output[0]} --uc {output[1]} --strand {config[derep][strand]} "
                "--fasta_width 0 --minuniquesize {config[derep][min_abundance]} --sizeout" 
                
    if config["pickOTU"]["m"] == "swarm" :
        rule cluster_OTUs_swarm:
            input:
                "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.fasta" 
            output:
                swarms="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt",
               # rep_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set_swarm.fasta", (in case we want to generate reps -w {output.rep_seqs})
                uc="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.uc",
                rep_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set_tmp.fasta" 
            params:
                otuDir="{PROJECT}/runs/{run}/otu/"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/cluster_OTUs_swarm.benchmark"
            #-i  {params.otuDir}swarm.struct
            shell:
                "swarm  -s {params.otuDir}swarm.stats -d {config[pickOTU][s]} -z "
                "-o {output.swarms}  -u {output.uc}  -t {config[pickOTU][cpus]} "
                "-w {output.rep_seqs} "
                "{config[pickOTU][extra_params]} < {input} "

    if config["pickOTU"]["m"] == "vsearch":
        rule cluster_OTUs_vsearch:
            input:
                "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.fasta" 
            output:
                otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt",
                uc="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.uc",
                rep_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set_tmp.fasta" 
            params:
                otuDir="{PROJECT}/runs/{run}/otu/",
                clusters="{PROJECT}/runs/{run}/otu/cluster"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/cluster_OTUs_vsearch.benchmark"
            shell:
                "vsearch --cluster_size {input} "
                "--threads {config[pickOTU][cpus]} --id  {config[pickOTU][s]} "
                "--uc {output.uc} "
                "--strand plus --sizein --sizeout --xsize "
                "--fasta_width 0  "
                #"--relabel OTU "
                "{config[pickOTU][extra_params]} "
                "--centroids {output.rep_seqs} "
                "--otutabout {output.otus} "
                #"--clusters {params.clusters}"
    rule remap_clusters:
        '''
        During the dereplicaition we use "sample" ids so it is important NOT
        to relabel clusters with OTU ids, then the remapping fails or do weird staff.
        So we need to keep original ids.
        '''
        input:
            uc_otu="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.uc",
            uc_derep="{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.uc"
        output:
            map="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt",
            log="{PROJECT}/runs/{run}/otu/remap.log",
            otu="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt.biom.txt"
        benchmark:
           "{PROJECT}/runs/{run}/benchmark/cluster_OTUs_vsearch.benchmark"
        shell:
           "{config[java][command]} -cp Scripts/ClusterMapper/build/classes clustermapper.ClusterMapper uc2uc "
           "-uc {input.uc_derep} -uc2 {input.uc_otu} -o {output.map} --otu-tbl  --full-uc --relabel -l OTU -lidx 1 > {output.log}"
    rule rename_sequences:
        '''
        This rule rename the representative sequences fasta file. 
        As we dont ask for relabel we still have original ids on both (vsearch and swarm)
        The remapping process takes the original id and assign OTU number according to 
        the remapped OTU table:otu/seqs_fw_rev_combined_remapped_otus.txt
        '''
        input:
            rep_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set_tmp.fasta",
            map="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt"
        output:
            renamed_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set.fasta"
        benchmark:
           "{PROJECT}/runs/{run}/benchmark/rename_sequences.benchmark"
        shell:
           "awk -F\"\\t\" 'NR==FNR{{h[$2]=$1;next}} "
           "{{if($0 ~ \"^>\"){{gsub(\">\",\"\",$1); split($1,id,\";\"); "
           " if(h[id[1]]){{print(\">\"h[id[1]]\" \"$1)}} "
           " else{{print \"ERROR! No id found for: \"id[1]}}}}else{{print $0}}}}' "
           "{input.map} {input.rep_seqs} > {output.renamed_seqs}" 

    if  config["assignTaxonomy"]["tool"].lower() == "blast":
        rule assign_taxonomy:
            input:
                "{PROJECT}/runs/{run}/otu/representative_seq_set.fasta"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.out")
            params:
                reference="-db " +config["assignTaxonomy"]["blast"]["blast_db"] if len(str(config["assignTaxonomy"]["blast"]["blast_db"])) > 1
                else "-subject " +config["assignTaxonomy"]["blast"]["fasta_db"]
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/assign_taxonomy.benchmark"
            shell:
                "{config[assignTaxonomy][blast][command]} {params.reference} -query {input} -evalue {config[assignTaxonomy][blast][evalue]} "
                "-outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -num_threads {config[assignTaxonomy][blast][jobs]} "
                "-max_target_seqs {config[assignTaxonomy][blast][max_target_seqs]} -perc_identity {config[assignTaxonomy][blast][identity]} "
                "{config[assignTaxonomy][blast][extra_params]} -out {output[0]} "
        rule complete_blast_out:
            """
             Add missing ids to blast out
            """
            input:
                blastout="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.out",
                #otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt" if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch")
                map="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt.biom.txt"
                #or config["pickOTU"]["m"] == "swarm" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
            output:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.complete.out"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/complete_blast_out.benchmark"
            shell:
                "cat {input.blastout} | cut -f1 | sort | uniq | grep -v -w -F -f - {input.map} "
                "| awk '$0 !~ \"^#\"{{print $1\"\\tUnassigned\\t-\\t-\\t-\\t-\"}}' | cat {input.blastout} - > {output}"
        rule filter_best_hits_blast:
            """
             Before LCA only keep results with id and coverage eq. to the best hit.
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.complete.out"
            output:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.filtered.out"
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/filter_best_hits_blast.benchmark"
            shell:
                "cat {input}  "
                "| awk -v curr_s=\"\" ' {{if($1 != curr_s){{ id=$3;  cov=$4;  curr_s=$1;  print $0 }}else{{  if($3>=id && $4 >= cov){{   print $0  }} }}  }}' "
                " > {output}"
        
        rule prepare_blast_for_stampa:
            """
             Take completed blast output file and make some reformat in order to
             fulfill stampa format.
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.filtered.out"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.blast.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/prepare_blast_for_stampa.benchmark"
            shell:
                "cat {input}  | cut -f2 | sort | uniq | grep -F -w -f -  {config[assignTaxonomy][blast][mapFile]} | "
                "awk 'NR==FNR {{h[$1] = $2; next}} {{print $1\"\\t\"$3\"\\t\"$2\" \"h[$2]}}' FS=\"\\t\" - FS=\"\\t\" {input} "
                " > {output}"
        rule skip_stampa:
            """
             skip stampa by selecting the best hit. It keeps the accesions from the other hits
             and add them as metadata with the identity. If we use stampa the identities are the same
             for all the best hits, for blast the identity can be different and we only report the best one.
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.blast.out"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.blast.no_stampa.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/skip_stampa.benchmark"
            shell:
                #"cat {input}  | sed 's/ /\\t/1' | awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}}{{if(!h[$1]){{if($3 == \"*\"){{print $1,1,$2,\"Unassigned\",$3;h[$1]=$1}}else{{print $1,1,$2,$4,$3;h[$1]=$1}}}}}}' > {output}"
                "cat {input}  | sed 's/ /\t/1' | awk -F'\\t' -v current='' 'BEGIN{{OFS=\"\\t\"}}{{if(length(current)>0 && current != $1 && length(line)> 0){{print line\"\\t\"h[current];line=\"\";}};if(!h[$1]){{if($3 == \"*\"){{print $1,1,$2,\"Unassigned\",$3;h[$1]=$3;}}else{{line=$1\"\\t1\\t\"$2\"\\t\"$4;h[$1]=$3;current=$1}}}}else{{h[$1]=h[$1]\";\"$3}}}}END{{if(length(line)>0){{print line\"\\t\"h[current]}}}}' > {output}"
    
        rule run_stampa:
            """
             compute lca using stampa merge script
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.blast.out"
            params:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.blast.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/run_stampa.benchmark"
            shell:
                "Scripts/stampa_merge.py {params} {config[assignTaxonomy][blast][taxo_separator]}"
        rule normalize_taxo_out:
            """
             Normalize the output in terms of its format and names in order to be able
             to continue with the pipeline
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.blast.out" if config["assignTaxonomy"]["map_lca"].lower() == "t"
                else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.blast.no_stampa.out" 
            output:
                report("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt", category="7. Taxonomy", subcategory="1. Taxonomy assignation results",labels={"Representative sequences":"OTU rep. seqs"},)
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/normalize_taxo_out.benchmark"
            shell:
                "cat {input} |  awk -F\"\\t\" '{{print $1\"\\t\"$4\"\\t\"$3\"\\t\"$5}}' | sed 's/N;o;_;h;i;t/Unassigned/' > {output}"
    
    elif  config["assignTaxonomy"]["tool"].lower() == "vsearch":
        rule assign_taxonomy:
            input:
                "{PROJECT}/runs/{run}/otu/representative_seq_set.fasta"
            output:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_vsearch.out"
            params:
                reference="-db " +config["assignTaxonomy"]["blast"]["blast_db"] if len(str(config["assignTaxonomy"]["blast"]["blast_db"])) > 1
                else "-subject " +config["assignTaxonomy"]["blast"]["fasta_db"]
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/assign_taxonomy.benchmark"
            shell:
                "{config[assignTaxonomy][vsearch][command]}  --usearch_global {input} --db {config[assignTaxonomy][vsearch][db_file]} "
                "--dbmask none --qmask none --rowlen 0 --id {config[assignTaxonomy][vsearch][identity]} "
                "--iddef {config[assignTaxonomy][vsearch][identity_definition]}  --userfields query+id{config[assignTaxonomy][vsearch][identity_definition]}+target "
                "--threads {config[assignTaxonomy][vsearch][jobs]} {config[assignTaxonomy][vsearch][extra_params]} "
                " --maxaccepts {config[assignTaxonomy][vsearch][max_target_seqs]} --output_no_hits --userout {output[0]} "
    
        rule mapp_vsearch_taxo:
            """
             Take vsearch output and make some reformat in order to
             fulfill stampa format.
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_vsearch.out"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/taxons.txt")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/mapp_vsearch_taxo.benchmark"
            shell:
                "cat {input}  | cut -f3 | sort | uniq | grep -F -w -f -  {config[assignTaxonomy][vsearch][mapFile]} > {output} "
    
        rule prepare_vsearch_for_stampa:
            """
             Take completed blast output file and make some reformat in order to
             fulfill stampa format.
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/taxons.txt",
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_vsearch.out"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.vsearch.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/prepare_vsearch_for_stampa.benchmark"
            shell:
                "echo '*\\tUnassigned' | cat {input[0]} - | awk 'NR==FNR {{h[$1] = $2; next}} {{print $1\"\\t\"$2\"\\t\"$3\" \"h[$3]}}' FS=\"\\t\" - FS=\"\\t\" {input[1]} "
                " > {output}"
    
        rule skip_stampa:
            """
             skip stampa by selecting the best hit. It keeps the accesions from the other hits
             and add them as metadata with the identity. If we use stampa the identities are the same 
             for all the best hits, for blast the identity can be different and we only report the best one. 
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.vsearch.out"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.no_stampa.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/skip_stampa.benchmark"
            shell:
                #"cat {input}  | sed 's/ /\\t/1' | awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}}{{if(!h[$1]){{if($3 == \"*\"){{print $1,1,$2,\"Unassigned\",$3;h[$1]=$1}}else{{print $1,1,$2,$4,$3;h[$1]=$1}}}}}}' > {output}"
                "cat {input}  | sed 's/ /\t/1' | awk -F'\\t' -v current='' 'BEGIN{{OFS=\"\\t\"}}{{if(length(current)>0 && current != $1 && length(line)> 0){{print line\"\\t\"h[current];line=\"\";}};if(!h[$1]){{if($3 == \"*\"){{print $1,1,$2,\"Unassigned\",$3;h[$1]=$3;}}else{{line=$1\"\\t1\\t\"$2\"\\t\"$4;h[$1]=$3;current=$1}}}}else{{h[$1]=h[$1]\";\"$3}}}}END{{if(length(line)>0){{print line\"\\t\"h[current]}}}}' > {output}"
    
        rule run_stampa:
            """
             compute lca using stampa merge script
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.vsearch.out"
            params:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/"
            output:
                temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.out")
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/run_stampa.benchmark"
            shell:
                "Scripts/stampa_merge.py {params} {config[assignTaxonomy][vsearch][taxo_separator]}"
        rule normalize_taxo_out:
            """
             Normalize the output in terms of format and names in order to be able
             to continue with the pipeline
            """
            input:
                "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.out" if config["assignTaxonomy"]["map_lca"].lower() == "t"
                else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.no_stampa.out"  
            output:
                report("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt", category="7. Taxonomy", subcategory="1. Taxonomy assignation results",labels={"Representative sequences":"OTU rep. seqs"},)
            benchmark:
                "{PROJECT}/runs/{run}/benchmark/normalize_taxo_out.benchmark"
            shell:
                "cat {input} |  awk -F\"\\t\" '{{print $1\"\\t\"$4\"\\t\"$3\"\\t\"$5}}' | sed 's/N;o;_;h;i;t/Unassigned/' > {output}"


    rule map_taxonomy_otu_table:
        input:
            tax="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt",
            otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt.biom.txt" # if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch") 
           # or config["pickOTU"]["m"] == "swarm" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
        output:
            report("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt", category="6. OTU results", labels={"Type": "OTU table"},)
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/map_taxonomy_otu_table.benchmark"
        shell:
            "cat {input.tax} |  awk -F'\\t' 'NR==FNR{{split($1,record,\";\"); id=record[1];gsub(\" \",\"\",$2); h[id]=$2;next}} {{if(FNR==1){{print $0\"\\ttaxonomy\"}}else{{print $0\"\\t\"h[$1]}}}}' - {input.otus} > {output}"
    
    rule make_otu_table:
        input:
            #tax="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt",
            #otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt.biom.txt" # if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch") 
           # or config["pickOTU"]["m"] == "swarm" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt" 
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.tmp.biom" 
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/make_otu_table.benchmark"
        shell:
            #"{config[qiime][path]}make_otu_table.py -i {input.otus} -t {input.tax} -o {output} {config[makeOtu][extra_params]}"
            "{config[biom][command]} convert -i {input[0]} -o {output} --table-type \"OTU table\" --to-hdf5 --process-obs-metadata taxonomy "
    
    rule prapare_obs_metadata:
        '''
        Prepare observation metadata 
        '''
        input:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt"
        output:
            temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otu.metadata.tsv")
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/prapare_obs_metadata.benchmark"
        shell:
            "cat {input} | awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\";print \"#OTUID\\tIdentity\\tACCs\"}}{{print $1,$3,$4}}' > {output}" 
    
    rule add_obs_metadata:
        input:
            metadata="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otu.metadata.tsv",
            otu="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.tmp.biom"
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/add_obs_metadata.benchmark"
        shell:
            "{config[biom][command]} add-metadata -i {input.otu} -o {output} --observation-metadata-fp {input.metadata} --float-fields Identity"
            
    rule filter_otu:
        input:
            #"{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom"
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt"     
        output:
            report("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt", category="6. OTU results", labels={"Type": "OTU table no singletons"},)
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/filter_otu.benchmark"
        params: 
            min_abundance=config["filterOtu"]["n"],
            min_sample=config["filterOtu"]["s"]

        shell:
            #"{config[qiime][path]}filter_otus_from_otu_table.py -i {input} -o {output} -n {config[filterOtu][n]} {config[filterOtu][extra_params]}"
            "cat {input} | "
            "awk  -v n={params.min_abundance} -v s={params.min_sample} -F\"\\t\" 'NR==1{{print $0;tax=-1; "
            "for(i=1;i<=NF;i++){{if($i==\"taxonomy\"){{tax=i}}}}}} "
            "NR>1{{samples=0;total=0; "
            "for(i=2;i<tax;i++){{ "
            "  total+=$i;"
            "  if($i>0){{samples++}}}}; "
            "if(samples>=s && total>=n){{print $0}}}}' > {output}"

    
    #Convert to/from the BIOM table format.
    rule convert_filter_otu:
        input:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt"
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.biom"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/convert_filter_otu.benchmark"
        shell:
            "{config[biom][command]} convert -i {input} -o {output} --table-type \"OTU table\" --to-hdf5 --process-obs-metadata taxonomy "   
                 
#------------------------------------------------
#-- Taxonomy summarization and OTU tabe filtering for both ASV and OTU wfs 
#----------------------------------------------------

#TODO change summarize taxa from qiime
rule summarize_taxa:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/otuTable_L2.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/asvTable_L2.txt"
    params:
        prefix="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/otuTable"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/asvTable",
        row_names="--print-row-names" if config["summTaxa"]["row_names"].lower() == "t"
        else ""
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/summarize_taxa.benchmark"
    shell:
        "python Scripts/summarize_taxonomy.py -l {config[summTaxa][levels]} -o {params.prefix} {params.row_names} {input}"    
rule summarize_taxa_no_singletons:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/otuTable_noSingletons_L2.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/asvTable_noSingletons_L2.txt"
    params:
        prefix="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/otuTable_noSingletons"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/asvTable_noSingletons",
        row_names="--print-row-names" if config["summTaxa"]["row_names"].lower() == "t"
        else ""
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/summarize_taxa_no_singletons.benchmark"
    shell:
        "python Scripts/summarize_taxonomy.py -l {config[summTaxa][levels]} -o {params.prefix} -t {config[summTaxa][taxo_separator]} {params.row_names} {input}"

if  config["krona"]["report"].casefold() == "t" or config["krona"]["report"].casefold() == "true":
    rule krona:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt" if config["ANALYSIS_TYPE"] == "ASV" and config["krona"]["otu_table"].casefold() != "singletons"  
            else "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt" 
            if config["ANALYSIS_TYPE"] == "ASV"
            else
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt" 
            if config["krona"]["otu_table"].casefold() != "singletons"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt"
        output:
            #"{PROJECT}/runs/{run}/report_files/krona_report.dada2.html"
            #if config["ANALYSIS_TYPE"] == "ASV" else
            #"{PROJECT}/runs/{run}/report_files/krona_report."+config["assignTaxonomy"]["tool"].lower()+".html",
            #report(directory("{PROJECT}/runs/{run}/report_files/"),htmlindex="krona_report.dada2.html") if config["ANALYSIS_TYPE"] == "ASV" else
            #report(directory("{PROJECT}/runs/{run}/report_files/"),htmlindex="krona_report."+config["assignTaxonomy"]["tool"].lower()+".html")
            report("{PROJECT}/runs/{run}/report_files/krona_report.dada2.html" ,htmlindex="krona_report.dada2.html", category="6. Taxonomy", subcategory="2. Krona: Taxonomy plot", labels={"Taxonomy assignation":"Krona plot"},)
            if config["ANALYSIS_TYPE"] == "ASV" else
            report("{PROJECT}/runs/{run}/report_files/krona_report."+config["assignTaxonomy"]["tool"].lower()+".html", htmlindex="krona_report."+config["assignTaxonomy"]["tool"].lower()+".html", category="7. Taxonomy", subcategory="2. Krona: Taxonomy plot", labels={"Taxonomy assignation":"Kona plot"},)

        params:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/" 
            if config["ANALYSIS_TYPE"] == "ASV"  else 
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/" 
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/krona_report.benchmark"
        script:
            "Scripts/otu2krona.py"
    
    #rule krona_report:
    #    input:
    #        "{PROJECT}/runs/{run}/report_files/krona_report.dada2.html"
    #        if config["ANALYSIS_TYPE"] == "ASV" else
    #        "{PROJECT}/runs/{run}/report_files/krona_report."+config["assignTaxonomy"]["tool"].lower()+".html"
    #    output:
    #        report(directory("{PROJECT}/runs/{run}/krona_report/"),htmlindex="krona_report.dada2.html") if config["ANALYSIS_TYPE"] == "ASV" else
    #        report(directory("{PROJECT}/runs/{run}/krona_report/"),htmlindex="krona_report."+config["assignTaxonomy"]["tool"].lower()+".html")
    #    shell:
    #        "cp "

else:
    rule skip_krona:
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/krona_report.benchmark"
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/skip.krona" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/skip.krona"
        shell:
            "touch {output}"
#OTU map-based filtering: Keep all sequences that show up in an OTU map.
rule filter_rep_seqs:
    input:
        fastaRep="{PROJECT}/runs/{run}/asv/representative_seq_set.fasta"
        if config["ANALYSIS_TYPE"] == "ASV" else
        "{PROJECT}/runs/{run}/otu/representative_seq_set.fasta",
        otuNoSingleton="{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt"
        if config["ANALYSIS_TYPE"] == "ASV" else
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt"
    output:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta"
         if config["ANALYSIS_TYPE"] == "ASV" else
         "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_noSingletons.fasta"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/filter_rep_seqs.benchmark"
    shell:
        "seqkit grep -i -f <(cat {input.otuNoSingleton} | grep -v \"^#\" | cut -f1 ) --line-width 0 -o {output} {input.fastaRep}"

if config["alignRep"]["align"] == "T":
#Align representative sequences
    rule align_rep_seqs:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta" if config["ANALYSIS_TYPE"] == "ASV" 
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_noSingletons.fasta"
        output:
            aligned="{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/representative_seq_set_noSingletons_aligned.fasta",
            log="{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_log.txt" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/representative_seq_set_noSingletons_log.txt"
        params:
            outdir="{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/" if config["ANALYSIS_TYPE"] == "ASV" 
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/align_rep_seqs.benchmark"
        shell:
            "{config[qiime][path]}align_seqs.py -m {config[alignRep][m]} -i {input} -o {params.outdir} {config[alignRep][extra_params]}"

#This step should be applied to generate a useful tree when aligning against a template alignment (e.g., with PyNAST)
    rule filter_alignment:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/representative_seq_set_noSingletons_aligned.fasta"
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.fasta" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.fasta"
        params:
            outdir="{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/filter_alignment.benchmark"
        shell:
            "{config[qiime][path]}filter_alignment.py -i {input} -o {params.outdir} {config[filterAlignment][extra_params]}"

#Many downstream analyses require that the phylogenetic tree relating the OTUs in a study be present.
#The script make_phylogeny.py produces this tree from a multiple sequence alignment
    rule make_tree:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.fasta" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.fasta"
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.tre" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.tre"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/make_tree.benchmark"
        shell:
            "{config[qiime][path]}make_phylogeny.py -i {input} -o {output} -t {config[makeTree][method]} {config[makeTree][extra_params]}"

#----------Align and filtering section ----#
#####---Both workflows statistics   
if config["demultiplexing"]["demultiplex"] == "T" and  config["ANALYSIS_TYPE"] != "ASV":
    rule count_samples_final:
        input:
            fasta="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/count_samples_final.{sample}.benchmark"
        shell:
            "cat {input.fasta} | grep '^>' |  cut -d'_' -f1 | sed 's/>//g' "
            "| sort | uniq -c | sort -nr | awk '{{print $1\"\\t\"$2}}' "
            "| awk 'NR==FNR{{h[$2]=$1; next}} {{print $1\"\\t\"h[$1]}}' - {input.metadata} | grep -v \"#\" > {output}"

elif config["ANALYSIS_TYPE"] == "ASV" and config["demultiplexing"]["demultiplex"] == "T":

    rule count_samples_final:
        input:
            "{PROJECT}/runs/{run}/asv/filter_summary.out",
            "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/count_samples_final.{sample}.benchmark"
        shell:
            "cat {input[0]} | awk 'NR==FNR{{if(NR>1){{h[$1]=$2;}}next}}{{if(FNR>1){{print $1\"\t\"h[$1]}}}}' - {input[1]} > {output}"

elif config["ANALYSIS_TYPE"] == "ASV" and config["demultiplexing"]["demultiplex"] != "T":

    rule count_samples_final:
        input:
            "{PROJECT}/runs/{run}/asv/filter_summary.out",
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/count_samples_final.{sample}.benchmark"
        shell:
            "cat {input} | awk 'NR>1{{print $1\"\\t\"$2}}' > {output}"


else:
   rule count_samples_final:
        input:
            fasta="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/count_samples_final.{sample}.benchmark"
        shell:
            "echo {wildcards.sample}\\t 100 > {output}"


rule distribution_chart:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
    output:
        #"{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered2.{sample}.dist.png"
        report("{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered.{sample}.dist.png", category="3. Demultiplex", subcategory="Sample distribution", labels={"sample": "{sample}","Description": "Plot: number of reads per sample"})
        if config["ANALYSIS_TYPE"] == "ASV" else
        report("{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered.{sample}.dist.png", category="4. Demultiplex", subcategory="Sample distribution", labels={"sample": "{sample}","Description": "Plot: number of reads per sample"})

    params:
        d="{PROJECT}/runs/{run}/report_files/",
        num="3. " if config["ANALYSIS_TYPE"] == "ASV" else "4. "
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/distribution_chart.{sample}.benchmark"
    shell:
        #"Scripts/sampleDist.py"
        "python  Scripts/sampleDist.manual.py {input} {output}"

if config["ANALYSIS_TYPE"]=="ASV":
    rule create_sample_log:
        input:
            ff=expand("{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt",PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run) 
        output:
            "{PROJECT}/runs/{run}/samples.log"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/create_sample_log.benchmark"
        script:
            "Scripts/combineAllReads_asv.py"

# -------
# - Both workflows: create md5sum files   
# -------
if config["LIBRARY_LAYOUT"] != "SE":
    rule fastqCheckSum_PE:
        input:
            #fw="{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            #else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
            #rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            #else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
            fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
        params:
            "{PROJECT}/samples/{sample}/demultiplexed/" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/fastqCheckSum_PE.benchmark"
        shell:
            "for file in {params}*.fastq.gz; do"
            "  fq=$(echo $file | awk -F\"/\" '{{print $NF}}');" 
            "  check=$(md5sum $file | cut -f1 -d \" \");"
            "  echo -e $check\"\\t\"$fq; "
            "done > {output}"   
else:
    rule fastqCheckSum_SE:
        input:
            fw="{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv"
        params:
            "{PROJECT}/samples/{sample}/demultiplexed/" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/benchmark/fastqCheckSum_SE.benchmark"
        shell:
            "for file in {params}*.fastq.gz; do"
            "  fq=$(echo $file | awk -F\"/\" '{{print $NF}}');" 
            "  check=$(md5sum $file | cut -f1 -d \" \");"
            "  echo -e $check\"\\t\"$fq; "
            "done > {output}"   

rule create_version_table:
    output:
        "{PROJECT}/runs/{run}/run.versions.tsv"
    script:
        "Scripts/all_versions_asv.py" if config["ANALYSIS_TYPE"].lower() == "asv"
        else "Scripts/all_versions_otu.py"

rule create_yaml_versions_tbl:
    output:
        "{PROJECT}/runs/{run}/versions.yaml"
    shell:
        "cp resources/datavzrd/versions.yaml {output}"

rule datavzrd_versions:
    input:
        config="{PROJECT}/runs/{run}/versions.yaml",
        table="{PROJECT}/runs/{run}/run.versions.tsv"
    output:
        report(
            directory("{PROJECT}/runs/{run}/versions/"),
            htmlindex="index.html",
            #category="1. Project Information",
            category="1. Project information: {PROJECT}",
            subcategory="Tools and versions",
            labels={"Run": "{run}", "Description":"Tools used for this run"},
        ),
    wrapper:
        "v4.7.2/utils/datavzrd"


if config["ANALYSIS_TYPE"] != "ASV": 
    rule report_all:
       input:
            a="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt",
#            if config["ANALYSIS_TYPE"] == "OTU" else
#            "{PROJECT}/runs/{run}/asv/asv_table.biom",
            b="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/otuTable_L2.txt",
#            if config["ANALYSIS_TYPE"] == "OTU" else
#            "{PROJECT}/runs/{run}/asv/summary/asv_table_L6.txt",
            c="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt",
            d="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/otuTable_noSingletons_L2.txt",
#            if config["ANALYSIS_TYPE"] == "OTU" else 
#            "{PROJECT}/runs/{run}/asv/stats_dada2.txt",
            e="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.tre"
            if config["alignRep"]["align"] == "T"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_noSingletons.fasta",
            f="{PROJECT}/runs/{run}/report_files/krona_report."+config["assignTaxonomy"]["tool"].lower()+".html"
            if config["krona"]["report"] == "T" and config["ANALYSIS_TYPE"] != "ASV" else 
            "{PROJECT}/runs/{run}/report_files/krona_report.dada2.html"
            if config["krona"]["report"] == "T" and config["ANALYSIS_TYPE"] == "ASV" else
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/skip.krona",
            g="{PROJECT}/runs/{run}/samples.log",
            h="{PROJECT}/runs/{run}/versions/"
       output:
            temp("{PROJECT}/runs/{run}/reporttmp_all.html")
       benchmark:
            "{PROJECT}/runs/{run}/benchmark/report_all.benchmark"
       shell:
            #"Scripts/report_all_v2.py"
            "touch {output[0]}"
else:
    rule report_all_asv:
       input:
            a="{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom",
            b="{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/asvTable_L2.txt",
            d="{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/asvTable_noSingletons_L2.txt",
            c="{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.tre"
            if config["alignRep"]["align"] == "T"
            else "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta",
 #           if config["alignRep"]["align"] == "T"
 #           else "{PROJECT}/runs/{run}/asv/dada2_representatives.fasta"
 #           e="{PROJECT}/runs/{run}/asv/dada2_representatives.fasta",
 #           if config["ANALYSIS_TYPE"] == "ASV"
 #           else
 #           "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"]+"/representative_seq_set_noSingletons.fasta",
            f="{PROJECT}/runs/{run}/report_files/krona_report.dada2.html"
            if config["krona"]["report"] == "T" else
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/skip.krona",
            g="{PROJECT}/runs/{run}/samples.log",
            h="{PROJECT}/runs/{run}/asv/taxonomy_dada2/bimeraTable.txt"
            if config["dada2_asv"]["chimeras"].lower() == "t" and config["dada2_asv"]["chimeras_taxonomy"].lower() == "t" else
            "{PROJECT}/runs/{run}/asv/representative_bimeras.fasta",
            i="{PROJECT}/runs/{run}/versions/"
       output:
            temp("{PROJECT}/runs/{run}/reporttmp_all.html")
       benchmark:
            "{PROJECT}/runs/{run}/benchmark/report_all_asv.benchmark"
       shell:
            #"Scripts/report_all_asv.py"
            "touch {output[0]}"


rule tune_report_all:
    input:
        "{PROJECT}/runs/{run}/reporttmp_all.html"
    output:
        "{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv_report_dada2.html"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/tune_report_all.benchmark"
    shell:
        #"Scripts/tuneReport.py"
        "touch {output[0]}"
if config["pdfReport"].casefold() == "t":
    rule translate_to_pdf:
        input:
            "{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
            if config["ANALYSIS_TYPE"] == "OTU" else
            "{PROJECT}/runs/{run}/asv_report_dada2.html"
        output:
            "{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".pdf"
            if config["ANALYSIS_TYPE"] == "OTU" else
            "{PROJECT}/runs/{run}/asv_report_dada2.pdf"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/translate_to_pdf.benchmark"
        shell:
            "{config[wkhtmltopdf_command]}  {input} {output}"
            #"touch {output[0]}"

rule report:
    input:
        report_all="{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv_report_dada2.html",
        dist_chart="{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered.{sample}.dist.png",
        #dmxFiles="{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.tsv" if config["UNPAIRED_DATA_PIPELINE"] != "T"
        #else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.tsv",
        md="{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
        else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/report.html"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/report.{sample}.benchmark"
    shell:
        #"Scripts/report_v2.py" if config["ANALYSIS_TYPE"] == "OTU" and config["LIBRARY_LAYOUT"] != "SE" 
        #else "Scripts/report_v2_SE.py" if config["ANALYSIS_TYPE"] == "OTU" and config["LIBRARY_LAYOUT"] == "SE"  
        #else "Scripts/report_asv.py"
        "touch {output}"
rule tune_report:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/report.html"
    output:
        "{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/report_{sample}_dada2.html"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/tune_report.{sample}.benchmark"
    script:
        "Scripts/tuneReport.py"
if config["pdfReport"].casefold() == "t":
    rule translate_pdf_final_report:
        input:
            toTranslate="{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html"
            if config["ANALYSIS_TYPE"] == "OTU" else
            "{PROJECT}/runs/{run}/report_{sample}_dada2.html",
            tmp="{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".pdf"
            if config["ANALYSIS_TYPE"] == "OTU" else
            "{PROJECT}/runs/{run}/asv_report_dada2.pdf"
        output:
            "{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".pdf"
            if config["ANALYSIS_TYPE"] == "OTU" else
            "{PROJECT}/runs/{run}/report_{sample}_dada2.pdf"
        benchmark:
            "{PROJECT}/runs/{run}/benchmark/translate_pdf_final_report.{sample}.benchmark"
        shell:
            "{config[wkhtmltopdf_command]} {input.toTranslate} {output}"

rule clean_data_sample:
#"""
#In this case, we use only dummy input files to maintain the
#workflow. The script clean all the intermedita files.
#the cleaning steps are divid into two rules du to wildcard restriction
#with the samples. 
#"""
    input:
        "{PROJECT}/runs/{run}/{sample}_data/report.html"
    output:
        "{PROJECT}/runs/{run}/clean_files_{sample}.log"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/clean_data_sample.{sample}.benchmark"
    script:
        "Scripts/clean_data_sample.py"

rule clean_data_project:
    input:
        expand("{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        if config["ANALYSIS_TYPE"] == "OTU" else
        expand("{PROJECT}/runs/{run}/report_{sample}_dada2.html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run),
        expand("{PROJECT}/runs/{run}/clean_files_{sample}.log", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
    output:
        "{PROJECT}/runs/{run}/clean_files.log"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/clean_data_project.benchmark"
    script:
        "Scripts/clean_data_otu.py"



rule create_portable_report:
    input:
        expand("{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        if config["ANALYSIS_TYPE"] == "OTU" else
        expand("{PROJECT}/runs/{run}/report_{sample}_dada2.html", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run),
        "{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv_report_dada2.html",
        "{PROJECT}/runs/{run}/clean_files.log"
    output:
        "{PROJECT}/runs/{run}/report_"+config["assignTaxonomy"]["tool"].lower()+".zip"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/report_dada2.zip"
    params:
        "{PROJECT}/runs/{run}/report_files"
    benchmark:
        "{PROJECT}/runs/{run}/benchmark/create_portable_report.benchmark"
    shell:
        "zip -r {output} {params} {input[0]} {input[1]}"

          
