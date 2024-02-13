"""
CASCABEL
Version: 6.0.2
Authors:
 - Julia Engelmann
 - Alejandro Abdala
 - Wietse Reitsma

Last update: 09/02/2024
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
        shell:
            "Scripts/init_sample_dmx.sh "+config["PROJECT"]+" "+config["LIBRARY"][0]+"  {input.fw} {input.rv}"

elif len(config["LIBRARY"])==1 and config["demultiplexing"]["demultiplex"] == "F" and len(config["input_files"])<2 and config["LIBRARY_LAYOUT"] == "SE":
    rule init_structure:
        input:
            fw = config["fw_reads"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
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
        script:
            "Scripts/init_sample.py"
elif config["demultiplexing"]["demultiplex"] == "T" and config["LIBRARY_LAYOUT"] == "SE" and len(config["input_files"])>2:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        script:
            "Scripts/init_sample_SE.py"

elif config["LIBRARY_LAYOUT"] == "SE" and len(config["input_files"])>2:
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
        script:
            "Scripts/init_sample_SE.py"
else: 
    rule init_structure:
        input:
            file_list = config["input_files"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"  if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        script:
            "Scripts/init_sample.py"

if config["LIBRARY_LAYOUT"] != "SE":
    rule fast_qc:
        """
        Runs QC on raw reads
        """
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        output:
            o1="{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            o2="{PROJECT}/samples/{sample}/qc/rv_fastqc.html",
            s1="{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt",
            s2="{PROJECT}/samples/{sample}/qc/rv_fastqc/summary.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/qc/fq.benchmark"
        shell:
            "{config[fastQC][command]} {input.r1} {input.r2} --extract {config[fastQC][extra_params]} -t 10 -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/"

    rule validateQC:
        """
        Interpret FastQC output and stops on interactive mode, if too many
        errors.
        """
        input:
            "{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt",
            "{PROJECT}/samples/{sample}/qc/rv_fastqc/summary.txt",
            "{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            "{PROJECT}/samples/{sample}/qc/rv_fastqc.html",
            "{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            "{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
        output:
            "{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
            "{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
        script:
            "Scripts/validateQC.py"

    rule pear:
        """
        Paire/extend paired-end data
        """
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
            tmp1="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
            tmp2="{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
        output:
            "{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq",
            "{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.fastq",
            "{PROJECT}/samples/{sample}/peared/pear.log",
            "{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.fastq"

        benchmark:
            "{PROJECT}/samples/{sample}/peared/pear.benchmark"
        params:
            "{PROJECT}/samples/{sample}/peared/seqs"
        shell:
            "{config[pear][command]} -f {input.r1} -r {input.r2} -o {params[0]} "
            "-t {config[pear][t]} -v {config[pear][v]} -j {config[pear][j]} -p {config[pear][p]} {config[pear][extra_params]} > "
            "{output[2]}"

    if config["UNPAIRED_DATA_PIPELINE"] == "T" or config["demultiplexing"]["add_unpair"] == "T":
        rule identify_unpaired_fw:
            """
            Identify reads that they were not paired.
            """
            input:
                fw="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.fastq"
            output:
                fwo="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.out"
            shell:
                "cat {input.fw} | awk '{{if((NR-1)%4==0){{header=$1}}else if((NR-2)%4==0){{seq=$0}}else if(NR%4==0){{print header\"\\t\"seq\"\\t\"$0}}}}' > {output.fwo}"
        rule identify_unpaired_rv:
            """
            Is like the rule above, but this RC the sequence and reverse the quality values!!
            """
            input:
                rv="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.fastq"
            output:
                rvo="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.out"
            shell:
                "cat {input.rv} | awk '{{if((NR-1)%4==0){{header=$1}}else if((NR-2)%4==0){{seq=$0}}else if(NR%4==0){{print header\"\\t\"seq\"\\t\"$0}}}}' > {output.rvo}"
             #"cat {input.rv} | awk 'BEGIN{{a[\"T\"]=\"A\";a[\"A\"]=\"T\";a[\"C\"]=\"G\";a[\"G\"]=\"C\";a[\"N\"]=\"N\"}}"
             #"{{if((NR-1)%4==0){{header=$1}}else if((NR-2)%4==0){{seq=$0}}else if(NR%4==0){{rc=\"\";rq=\"\";for(i=length(seq);i>0;i--){{k=substr(seq,i,1);rc=rc a[k];rq=rq substr($0,i,1)}} "
             #"printf \"%s\\t%s\\t%s\\n\",header,rc,rq}}}}' >   {output.rvo}"

        rule pair_unpaired:
            input:
                fw="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.out",
                rv="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.out"
            output:
                "{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq"
            shell:
                #"cat {input.fw} | awk -v char={config[CHAR_TO_PAIR]} 'NR==FNR{h[$1]=$2;next}{print \">\"$1\"\\n\"$2 char h[$1]}' - {input.rv} > {output}
                "cat {input.fw} | awk -v ch=\"{config[CHAR_TO_PAIR]}\" 'BEGIN{{for(i=1;i<=length(ch);i++){{chq=chq \"{config[QUALITY_CHAR]}\"}}}} "
                "NR==FNR{{seqs[$1]=$2;qa[$1]=$3;next}}{{print $1\"\\n\"seqs[$1] ch $2\"\\n+\\n\"qa[$1] chq $3}}' - {input.rv} > {output}"

    rule validatePear:
        """
        Validate the percentage of paired reads.
        """
        input:
            "{PROJECT}/samples/{sample}/peared/pear.log"
        output:
            "{PROJECT}/samples/{sample}/peared/pear.log.validation"
        script:
            "Scripts/validatePear.py"

    if config["fastQCPear"] == "T" and config["UNPAIRED_DATA_PIPELINE"] != "T":
        rule fastQCPear:
            input:
                r1="{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq"
            output:
                o1="{PROJECT}/samples/{sample}/peared/qc/seqs.assembled_fastqc.html",
                s2="{PROJECT}/samples/{sample}/peared/qc/seqs.assembled_fastqc/summary.txt"
            params:
                "{PROJECT}/samples/{sample}/peared/qc/"
            benchmark:
                "{PROJECT}/samples/{sample}/peared/qc/fq.benchmark"
            shell:
                "{config[fastQC][command]} {input.r1} -t 10  --extract -o {params}"

        rule validateFastQCPear:
            input:
                "{PROJECT}/samples/{sample}/peared/qc/seqs.assembled_fastqc/summary.txt",
                "{PROJECT}/samples/{sample}/peared/qc/seqs.assembled_fastqc.html",
                "{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq"
            output:
                "{PROJECT}/samples/{sample}/peared/qc/fq_fw_internal_validation_pear.txt"
            script:
                "Scripts/validatePearedQC.py"
    else:
        rule skipFastQCPear:
            input:
                "{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" else
       	        "{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq"
            output:
                "{PROJECT}/samples/{sample}/peared/qc/fq_fw_internal_validation_pear.txt"
            shell:
                "touch {output}"           

else: #SE Workflow only runs QC and validate QC all the above rules 

    rule fast_qc:
        """
        Runs QC on raw reads
        """
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
        output:
            o1="{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            s1="{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/qc/fq.benchmark"
        shell:
            "{config[fastQC][command]} {input.r1} --extract {config[fastQC][extra_params]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/"

    rule validateQC_SE:
        """
        Interpret FastQC output and stops on interactive mode, if too many
        errors.
        """
        input:
            "{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt",
            "{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            "{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
        output:
            "{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt"
        script:
            "Scripts/validateQC_SE.py"

if config["demultiplexing"]["demultiplex"] == "T":

    rule bc_mapping_validation:
        """
        Check if the mapping file is ok. Header line needs to start with #
        """
        input:
            mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        output:
            "{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log" # antes .log pero creo carpeta .log??
        benchmark:
            "{PROJECT}/metadata/bc_validation/{sample}/validation.benchmark"
        params:
            "{PROJECT}/metadata/bc_validation/{sample}/"
        shell:
            "{config[qiime][path]}validate_mapping_file.py -o {params} -m {input.mapp}"
#validate bc validation log file stop WF in failure
    rule validateBCV:
        input:
            "{PROJECT}/metadata/bc_validation/{sample}/sampleList_mergedBarcodes_{sample}.log",
            "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        output:
            "{PROJECT}/metadata/bc_validation/{sample}/validation.log"
        params:
            "{PROJECT}/metadata/bc_validation/{sample}/"
        script:
            "Scripts/validateBCV.py"

    if config["demultiplexing"]["add_unpair"] == "T" and  config["LIBRARY_LAYOUT"] != "SE":
        rule concat_assembly:
            input: 
                pair="{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq",
                unpair="{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq"
            output:
                temp("{PROJECT}/samples/{sample}/peared/seqs.assembled.PAIRED_UNPAIRED.fastq")
            shell:
                "cat {input.pair} {input.unpair} > {output}" 

    if  config["LIBRARY_LAYOUT"] != "SE":
        rule extract_barcodes_pe:
            """
            Extracts the barcodes from the fastq reads. This rules applies to the PAIRED END mode.
            For choosing the correct "assembly" it uses the peared/seqs.assembled.PAIRED_UNPAIRED.fastq (Paired + Unpaired) if the WF is
            different to the UNPAIRED WF and if the add_unpair option is set to T.
            It selects only the paired sequences if the workflow is different to the unpaired WF and ther option add_unpair is set to F.
            Finally, it uses only the unpaired sequences if the Workflow is the UNPAIRED no matter what the user choose for the add_unpair option.
            Now
            """
            input:
                tmp3="{PROJECT}/metadata/bc_validation/{sample}/validation.log", 
                tmp2="{PROJECT}/samples/{sample}/peared/pear.log.validation",  #see how we can change this
                assembly="{PROJECT}/samples/{sample}/peared/seqs.assembled.PAIRED_UNPAIRED.fastq" if config["demultiplexing"]["add_unpair"] == "T"
                else "{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq",
                #if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"     
                #else "{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" 
                #else "{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq", 
                tmpinput="{PROJECT}/samples/{sample}/peared/qc/fq_fw_internal_validation_pear.txt" 
            output:
                "{PROJECT}/samples/{sample}/barcodes/barcodes.fastq",
                "{PROJECT}/samples/{sample}/barcodes/reads.fastq"
            params:
                "{PROJECT}/samples/{sample}/barcodes/"
            benchmark:
                "{PROJECT}/samples/{sample}/barcodes/barcodes.benchmark"
            shell:
                "{config[qiime][path]}extract_barcodes.py -f {input.assembly} -c {config[ext_bc][c]} "
                "{config[ext_bc][bc_length]} {config[ext_bc][extra_params]} -o {params}"
        rule extract_barcodes_unpaired:
            """
            Now, we always run a separated rule to (for the moment) always create the unpaired reads 
            """
            input:
                tmp3="{PROJECT}/metadata/bc_validation/{sample}/validation.log", 
                tmp2="{PROJECT}/samples/{sample}/peared/pear.log.validation",  #see how we can change this
                assembly="{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq", 
                tmpinput="{PROJECT}/samples/{sample}/peared/qc/fq_fw_internal_validation_pear.txt" 
            output:
                "{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq",
                "{PROJECT}/samples/{sample}/barcodes_unpaired/reads.fastq"
            params:
                "{PROJECT}/samples/{sample}/barcodes_unpaired/"
            benchmark:
                "{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.benchmark"
            shell:
                "{config[qiime][path]}extract_barcodes.py -f {input.assembly} -c {config[ext_bc][c]} "
                "{config[ext_bc][bc_length]} {config[ext_bc][extra_params]} -o {params}"
    else: #The workflow is SE
        rule extract_barcodes_se:
            input:
                assembly="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                tmp="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
                tmp3="{PROJECT}/metadata/bc_validation/{sample}/validation.log"
            output:
                #"{PROJECT}/runs/{run}/{sample}_data/barcodes/barcodes.fastq",
                #"{PROJECT}/runs/{run}/{sample}_data/barcodes/reads.fastq"
                "{PROJECT}/samples/{sample}/barcodes/barcodes.fastq",
                "{PROJECT}/samples/{sample}/barcodes/reads.fastq"
            params:
                "{PROJECT}/samples/{sample}/barcodes/"
            benchmark:
                "{PROJECT}/samples/{sample}/barcodes/barcodes.benchmark"
            shell:
                "{config[qiime][path]}extract_barcodes.py -f {input.assembly} -c {config[ext_bc][c]} "
                "{config[ext_bc][bc_length]} {config[ext_bc][extra_params]} -o {params}"

#If allow mismatch correct bar code
    #if config["demultiplexing"]["bc_mismatch"] > 0:
    rule correct_barcodes:
        input:
            #bc="{PROJECT}/runs/{run}/{sample}_data/barcodes/barcodes.fastq",
            bc="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq",
            mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
            reads="{PROJECT}/samples/{sample}/barcodes/reads.fastq"
        output:
            ob="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq_corrected",
            ore="{PROJECT}/samples/{sample}/barcodes/reads.fastq_corrected",
            l="{PROJECT}/samples/{sample}/barcodes/demux.log",
            mx="{PROJECT}/samples/{sample}/barcodes/sample_matrix.txt",
            rc="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq_corrected_toRC" 
        params:
            hmap="{PROJECT}/samples/{sample}/barcodes/heat_map.png",
            ghmap="{PROJECT}/samples/{sample}/barcodes/heat_map_golay.png"
        benchmark:
            "{PROJECT}/samples/{sample}/barcodes/barcodes_corrected.benchmark"
        shell:
            "java -jar Scripts/BarcodeCorrector.jar  -fb {input.bc} -fr {input.reads} -b  {input.mapp} -m  "  + str(config["demultiplexing"]["bc_mismatch"]) + ""
            " -o {output.ob} -or {output.ore} -rc -x {output.mx} " +  str(config["demultiplexing"]["bcc_params"]) + " > {output.l} "
    rule correct_barcodes_unpaired:
        input:
            #bc="{PROJECT}/runs/{run}/{sample}_data/barcodes/barcodes.fastq",
            bc="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq",
            mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
            reads="{PROJECT}/samples/{sample}/barcodes_unpaired/reads.fastq"
        output:
            ob="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq_corrected",
            ore="{PROJECT}/samples/{sample}/barcodes_unpaired/reads.fastq_corrected",
            l="{PROJECT}/samples/{sample}/barcodes_unpaired/demux.log",
            mx="{PROJECT}/samples/{sample}/barcodes_unpaired/sample_matrix.txt",
            rc="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq_corrected_toRC" 
        params:
            hmap="{PROJECT}/samples/{sample}/barcodes_unapired/heat_map.png",
            ghmap="{PROJECT}/samples/{sample}/barcodes_unapired/heat_map_golay.png"
        benchmark:
            "{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes_corrected.benchmark"
        shell:
            "java -jar Scripts/BarcodeCorrector.jar  -fb {input.bc} -fr {input.reads} -b  {input.mapp} -m  "  + str(config["demultiplexing"]["bc_mismatch"]) + ""
            " -o {output.ob} -or {output.ore} -rc -x {output.mx} " +  str(config["demultiplexing"]["bcc_params"]) + " > {output.l} "
    if config["demultiplexing"]["create_tag_pairs_heatmap"] == "T":
        rule create_heat_map:
            input:
                "{PROJECT}/samples/{sample}/barcodes/sample_matrix.txt"
            output:
                "{PROJECT}/samples/{sample}/barcodes/heat_map.png"
            params:
                ghmap="{PROJECT}/samples/{sample}/barcodes/heat_map_golay.png"
            shell:
                "Rscript Scripts/heatMapDemux.R $PWD {input} {output} {params.ghmap}"
 
#split libraries - demultiplex
    rule split_libraries:
        input:
            rFile="{PROJECT}/samples/{sample}/barcodes/reads.fastq_corrected", # if config["demultiplexing"]["bc_mismatch"] else "{PROJECT}/runs/{run}/{sample}_data/barcodes/reads.fastq",
            mapFile="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
            bcFile="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq_corrected", # if config["demultiplexing"]["bc_mismatch"] else "{PROJECT}/runs/{run}/{sample}_data/barcodes/barcodes.fastq",
            #tmp3="{PROJECT}/metadata/bc_validation/{sample}/validation.log"
        output:
            seqs="{PROJECT}/samples/{sample}/splitLibs/seqs.fna",
            spliLog="{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt",
        params:
            outDir="{PROJECT}/samples/{sample}/splitLibs",
        benchmark:
            "{PROJECT}/samples/{sample}/splitLibs/splitLibs.benchmark"
        shell:
            "{config[qiime][path]}split_libraries_fastq.py -m {input.mapFile} -i {input.rFile} "
            "-o {params.outDir} -b {input.bcFile} -q {config[split][q]} -r {config[split][r]} "
            "--retain_unassigned_reads --barcode_type {config[split][barcode_type]} {config[split][extra_params]}"
    rule split_libraries_unpaired:
        input:
            rFile="{PROJECT}/samples/{sample}/barcodes_unpaired/reads.fastq_corrected", # if config["demultiplexing"]["bc_mismatch"] else "{PROJECT}/runs/{run}/{sample}_data/barcodes/reads.fastq",
            mapFile="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
            bcFile="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq_corrected", # if config["demultiplexing"]["bc_mismatch"] else "{PROJECT}/runs/{run}/{sample}_data/barcodes/barcodes.fastq",
            #tmp3="{PROJECT}/metadata/bc_validation/{sample}/validation.log"
        output:
            seqs="{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.fna",
            spliLog="{PROJECT}/samples/{sample}/splitLibs_unpaired/split_library_log.txt",
        params:
            outDir="{PROJECT}/samples/{sample}/splitLibs_unpaired",
        benchmark:
            "{PROJECT}/samples/{sample}/splitLibs_unpaired/splitLibs.benchmark"
        shell:
            "{config[qiime][path]}split_libraries_fastq.py -m {input.mapFile} -i {input.rFile} "
            "-o {params.outDir} -b {input.bcFile} -q {config[split][q]} -r {config[split][r]} "
            "--retain_unassigned_reads --barcode_type {config[split][barcode_type]} {config[split][extra_params]}"
#This condition does not go any more, because NOW we always run the barcode correction tool to orientate the reads.
 
#    if config["LIBRARY_LAYOUT"] != "SE" and config["demultiplexing"]["bc_mismatch"] == 0:
#        rule get_unassigned:
#            '''
#            This IF was running always for a PE workflow, but now, the bc correction tool
#            test the reverse complement of the unassigneds and ouput the reads in the expected order
#            so now, we only get the unassigneds and continue with the old "methodology" when there is
#            no barcode correction performed.
#            '''
#            input:
#                split="{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.fna",
#                assembly="{PROJECT}/samples/{sample}/peared/seqs.assembled.PAIRED_UNPAIRED.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"     
#                else "{PROJECT}/runs/{run}/{sample}_data/peared/seqs.assembled.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" 
#                else "{PROJECT}/runs/{run}/{sample}_data/peared/seqs.assembled.UNPAIRED.fastq" 
#            output:
#                temp("{PROJECT}/runs/{run}/{sample}_data/splitLibs/unassigned.fastq")
#            shell:
#                "cat {input.split} | grep \"^>Unassigned\" |  sed 's/>Unassigned_[0-9]* /@/g' | "
#                "sed 's/ .*//' | grep -F --no-group-separator -w -A3  -f - {input.assembly}  > {output}"
#    elif config["LIBRARY_LAYOUT"] == "SE" and config["demultiplexing"]["bc_mismatch"] == 0:
#        rule uncompress_singleEND:
#            input:
#                "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
#            output:
#                temp("{PROJECT}/samples/{sample}/rawdata/fw.tmp.fastq")
#            shell:
#                "gzip -cd {input} > {output}"
#        rule get_unassigned:
#            input:
#                split="{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.fna",
#                assembly="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.tmp.fastq"
#            output:
#                temp("{PROJECT}/runs/{run}/{sample}_data/splitLibs/unassigned.fastq")
#            shell:
#                "cat {input.split} | grep \"^>Unassigned\" |  sed 's/>Unassigned_[0-9]* /@/g' | "
#                "sed 's/ .*//' | grep -F -w -A3  -f - {input.assembly} |  sed '/^--$/d' > {output}"
#
#    rule rc_unassigned:
#        input:
#            "{PROJECT}/runs/{run}/{sample}_data/splitLibs/unassigned.fastq"
#        output:
#            temp("{PROJECT}/runs/{run}/{sample}_data/splitLibs/unassigned.reversed.fastq")
#        shell:
#            "vsearch --fastx_revcomp {input} --fastqout {output}"
#    rule extract_barcodes_unassigned:
#        input:
#            assembly="{PROJECT}/runs/{run}/{sample}_data/splitLibs/unassigned.reversed.fastq"
#        output:
#            "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.fastq",
#            "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/reads.fastq"
#        params:
#            "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/"
#        benchmark:
#            "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.benchmark"
#        shell:
#            "{config[qiime][path]}extract_barcodes.py -f {input.assembly} -c {config[ext_bc][c]} "
#            "{config[ext_bc][bc_length]} {config[ext_bc][extra_params]} -o {params}"
#If allow mismatch correct bar code
#    if config["demultiplexing"]["bc_mismatch"]:
#        rule correct_barcodes_unassigned:
#            input:
#                bc="{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.fastq",
#                mapp="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
#            output:
#                "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.fastq_corrected"
#            benchmark:
#                "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes_corrected.benchmark"
#            shell:
#                "java -jar Scripts/BarcodeCorrector.jar barcodecorrector.BarcodeCorrector -fb {input.bc} -b  {input.mapp} -m  "  + str(config["demultiplexing"]["bc_mismatch"])
#                #"{config[Rscript][command]} Scripts/errorCorrectBarcodes.R $PWD {input.mapp} {input.bc} "  + str(config["demultiplex"]["bc_mismatch"])
#
#    rule split_libraries_rc:
#        """
#        This rule will call a script in order to execute the librarie splitting for the
#        RC sequences. It still need the seqs.fna file bz this file contains the exact
#        number of sequences assigned during the first split, then the script takes that
#        number and start to assign reads for the new splitting starting at the previous
#        number...
#        """
#
#        input:
#            spplited="{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.fna",
#            rFile="{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/reads.fastq",
#            mapFile="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt",
#            bcFile="{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.fastq_corrected" if config["demultiplexing"]["bc_mismatch"]
#            else "{PROJECT}/runs/{run}/{sample}_data/barcodes_unassigned/barcodes.fastq"
#        output:
#            seqsRC=temp("{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.fna"), #marc as tmp
#            spliLog="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/split_library_log.txt"
#        params:
#            outDirRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC"
#        benchmark:
#            "{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/splitLibs.benchmark"
#        script:
#            "Scripts/splitRC.py"

    rule remove_unassigned:
        input:
            "{PROJECT}/samples/{sample}/splitLibs/seqs.fna"
        output:
            "{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna"
        shell:
            "cat {input} | grep -P -A1 --no-group-separator '(?!>Unass)^>'  > {output}"
    rule remove_unassigned_unpaired:
        input:
            "{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.fna"
        output:
            "{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.no_unassigneds.fna"
        shell:
            "cat {input} | grep -P -A1 --no-group-separator '(?!>Unass)^>'  > {output}"            

    rule create_assigned_reads:
        input:
            rc="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq_corrected_toRC",
            assigned="{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
            else "{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna"
        output:
            temp("{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.ori.txt")
        shell:
            "cat {input.rc} | awk '{{print substr($1,2)}}'| grep -F -w -v -f - {input.assigned} | grep \"^>\" > {output} || true"
    rule create_assigned_reads_unpaired:
        input:
            rc="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq_corrected_toRC",
            assigned="{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.no_unassigneds.fna" 
        output:
            temp("{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.assigned.ori.txt")
        shell:
            "cat {input.rc} | awk '{{print substr($1,2)}}'| grep -F -w -v -f - {input.assigned} | grep \"^>\" > {output} || true"


    rule create_assigned_reads_rc:
        input:
            rc="{PROJECT}/samples/{sample}/barcodes/barcodes.fastq_corrected_toRC",
            assigned="{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
            else "{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna"
        output:
            temp("{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.rc.txt")
        shell:
            "cat {input.rc} | awk '{{print substr($1,2)}}'| grep -F -w -f - {input.assigned} > {output} || true"
    rule create_assigned_reads_rc_unpaired:
        input:
            rc="{PROJECT}/samples/{sample}/barcodes_unpaired/barcodes.fastq_corrected_toRC",
            assigned="{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.no_unassigneds.fna"
        output:
            temp("{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.assigned.rc.txt")
        shell:
            "cat {input.rc} | awk '{{print substr($1,2)}}'| grep -F -w -f - {input.assigned} > {output} || true"



#    rule remove_unassigned_rv:
#        input:
#            splitRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.fna"
#        output:
#            "{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.no_unassigneds.fna"
#        shell:
#            "cat {input} | grep -P -A1 '(?!>Unass)^>' | sed '/^--$/d' > {output}"

    rule create_unassigned_file:
        """
        With the new demux tool, we alwyas run the bc correction regardless of the bc_mismatch, so reads
        always end up in the correct orientation and thus we dont use any more the results from splitLibRC
        """
        input:
            #splitRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.fna" if config["demultiplexing"]["bc_mismatch"] == 0
            #else "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.fna"
            "{PROJECT}/samples/{sample}/splitLibs/seqs.fna"
        output:
            "{PROJECT}/samples/{sample}/splitLibs/seqs.unassigned.fna"
        shell:
            "cat {input} | grep -A1 --no-group-separator \"^>Unassigned\"  > {output}"

    if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T" and config["LIBRARY_LAYOUT"] != "SE":
        rule find_unassembled_fw_ids:
             """
             This rule creates a file with all the new ids assigned to unassembled reads.
             The file seqs.unassembled.forward.out is equivalent to the seqs.unassembled.reverse.out in this context
             since we only need the fadstq ids which should be exactly the same in both files.
             
             From whatever we had assigned, we remove the unassembled reads. 
             """
             input:
                 fasta="{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna",
                 mapa="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.out"
             output:
                 temp("{PROJECT}/samples/{sample}/splitLibs/seqs.unassembled_ids.txt")
             shell: 
                 "cat {input.mapa} | cut -f1 | sed 's/@//' | grep -F -w -f - {input.fasta} | cut -f1 -d\" \" | sed 's/>//' > {output} || true" 
        rule find_unassembled_rv_ids:
             """
             Similar rule as find_unassembled_fw_ids. Here we can 
             use either the reverse or forward, they should be the same      
             """
             input:
                 fasta="{PROJECT}/samples/{sample}/splitLibsRC/seqs.no_unassigneds.fna",
                 mapa="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.out"
             output:
                 temp("{PROJECT}/samples/{sample}/splitLibsRC/seqs.unassembled_ids.txt")
             shell: 
                 "cat {input.mapa} | cut -f1 | sed 's/@//' | grep -F -w -f - {input.fasta} | cut -f1 -d\" \" | sed 's/>//' > {output} || true"
        rule remove_unassembled:
             '''
             In the output I have added the if else, so if we perform the demultiplex, no need to do any
             reverse complement search for barcodes.
             '''
             input:
                 fasta="{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna",
                 ids="{PROJECT}/samples/{sample}/splitLibs/seqs.unassembled_ids.txt"
             output:
                 #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" if config["demultiplexing"]["bc_correction"]
                 #else "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.assigned.fna"  
                 "{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna"
             shell: 
                 "{config[qiime][path]}filter_fasta.py -f {input.fasta} -s {input.ids} -n -o {output}"
        rule remove_unassembled_rv:
             input:
                 fasta="{PROJECT}/samples/{sample}/splitLibsRC/seqs.no_unassigneds.fna",
                 ids="{PROJECT}/samples/{sample}/seqs.unassembled_ids.txt"
             output:
                 "{PROJECT}/samples/{sample}/splitLibsRC/seqs.assigned.fna"
             shell: 
                 "{config[qiime][path]}filter_fasta.py -f {input.fasta} -s {input.ids} -n -o {output}"
    else:
         rule rename_assembled_fw:
             input:
                 "{PROJECT}/samples/{sample}/splitLibs/seqs.no_unassigneds.fna"
             output:
                 "{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna"
             shell: 
                 "mv {input} {output}"
         rule rename_assembled_rv:
             input:
                 "{PROJECT}/samples/{sample}/splitLibsRC/seqs.no_unassigneds.fna"
             output:
                 "{PROJECT}/samples/{sample}/splitLibsRC/seqs.assigned.fna"
             shell: 
                 "mv {input} {output}"
    #IF ELSE HERE CHANGE validateDemuxScript, see what is going on for single end!!!
    #same here, different validations, but now always demux is in one step. 
  #  if config["demultiplexing"]["bc_mismatch"]:
    rule validateDemultiplex:
        input:
            split="{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna", 
            logSplit="{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt",
            logSplitDemux="{PROJECT}/samples/{sample}/barcodes/demux.log",
            allreads="{PROJECT}/samples/{sample}/barcodes/reads.fastq",
            unassigned="{PROJECT}/samples/{sample}/splitLibs/seqs.unassigned.fna",
            #heatMap="{PROJECT}/runs/{run}/{sample}_data/barcodes/demux.log"
            hmap="{PROJECT}/samples/{sample}/barcodes/sample_matrix.txt" if config["demultiplexing"]["create_tag_pairs_heatmap"] != "T"
            else "{PROJECT}/samples/{sample}/barcodes/heat_map.png"
               
        output:
            "{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt.validation"
        params:
            "{PROJECT}/samples/{sample}/splitLibs"
        script:
            "Scripts/validateSplitDemux.py"
    rule combine_accepted_reads:
        input:
            seqs="{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna", #marc as tmp
            tmpFlow="{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt.validation"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/combine_seqs_fw_rev.benchmark"
        shell:
            "ln -s $PWD/{input.seqs} {output}"
    
                
                
  #  else:
  #      rule validateDemultiplex:
  #          input:
  #              split="{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.assigned.fna", 
  #              splitRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.assigned.fna", 
  #              logSplit="{PROJECT}/runs/{run}/{sample}_data/splitLibs/split_library_log.txt",
  #              logSplitRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/split_library_log.txt",
  #              allreads="{PROJECT}/runs/{run}/{sample}_data/barcodes/reads.fastq",
  #              unassigned="{PROJECT}/runs/{run}/{sample}_data/seqs.unassigned.fna"
  #          output:
  #              "{PROJECT}/runs/{run}/{sample}_data/splitLibs/split_library_log.txt.validation"
  #          params:
  #              "{PROJECT}/runs/{run}/{sample}_data/splitLibs",
  #              "{PROJECT}/runs/{run}/{sample}_data/splitLibsRC"
  #          script:
  #              "Scripts/validateSplitNew.py"
  #      rule combine_accepted_reads:
  #          input:
  #              seqs="{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.assigned.fna", #marc as tmp
  #              seqsRC="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.assigned.fna", #marc as
  #              tmpFlow="{PROJECT}/runs/{run}/{sample}_data/splitLibs/split_library_log.txt.validation"
  #          output:
  #              "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
  #          benchmark:
  #              "{PROJECT}/runs/{run}/{sample}_data/combine_seqs_fw_rev.benchmark"
  #          shell:
  #              "cat {input.seqs} {input.seqsRC}  > {output}"
    

    if (config["demultiplexing"]["create_fastq_files"] == "T" or config["ANALYSIS_TYPE"] == "ASV") and config["LIBRARY_LAYOUT"] != "SE":
        rule write_dmx_files_fw:
            '''
            summary_wf.txt (output) is added after the new demux, as we assign samples in only one pass
            we only have one output  (summary_fw.txt) but in the past we also had summary_rv.txt, so this 
            new output will help with the WF 
            '''
            input:
                #dmx= "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
                #else "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.assigned.fna",
                dmx="{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.ori.txt",
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
                tmpFlow="{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt.validation"
            params:
                outdir="{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                sum="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt",
                wf="{PROJECT}/samples/{sample}/demultiplexed/summary_wf.txt"  
            benchmark:
                "{PROJECT}/samples/{sample}/demultiplexed/demultiplex_fq.benchmark"
            shell:
                "{config[java][command]}  -cp Scripts DemultiplexQiime --over-write --txt -a fw -b {config[demultiplexing][remove_bc]}  -d {input.dmx} -o {params.outdir} "
                "-r1 {input.r1} -r2 {input.r2} {config[demultiplexing][dmx_params]} > {output.wf}"
                 

        rule write_dmx_files_rv:
            """ 
            We supply the input 'overw' because it force to run first the write_dmx_files_fw which is the 
            one with the --over-write flag, so if the rule is re run, the generated files do not duplicate
            entries.
            """
            input:
                #dmx="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
                #else "{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.assigned.fna",
                dmx="{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.rc.txt", #if demux order_by_strand else validate arriba 
                overw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt",
                r2="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r1="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
            params:
                outdir="{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt"
            benchmark:
                "{PROJECT}/samples/{sample}/demultiplexed/demultiplex_fq_rc.benchmark"
            shell:
                "{config[java][command]}  -cp Scripts DemultiplexQiime  --txt -a rv -b {config[demultiplexing][remove_bc]}  -d {input.dmx} -o {params.outdir} "
                "-r1 {input.r1} -r2 {input.r2} {config[demultiplexing][dmx_params]}"
                
        rule write_dmx_files_fw_unpaired:
            '''
            summary_wf.txt (output) is added after the new demux, as we assign samples in only one pass
            we only have one output  (summary_fw.txt) but in the past we also had summary_rv.txt, so this 
            new output will help with the WF 
            '''
            input:
                #dmx= "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
                #else "{PROJECT}/runs/{run}/{sample}_data/splitLibs/seqs.assigned.fna",
                dmx="{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.assigned.ori.txt",
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz",
                tmpFlow="{PROJECT}/samples/{sample}/splitLibs/split_library_log.txt.validation"
            params:
                outdir="{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                sum="{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                wf="{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_wf.txt"  
            benchmark:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/demultiplex_fq.benchmark"
            shell:
                "{config[java][command]}  -cp Scripts DemultiplexQiime --over-write --txt -a fw -b {config[demultiplexing][remove_bc]}  -d {input.dmx} -o {params.outdir} "
                "-r1 {input.r1} -r2 {input.r2} {config[demultiplexing][dmx_params]} > {output.wf}"
                 

        rule write_dmx_files_rv_unapired:
            """ 
            We supply the input 'overw' because it force to run first the write_dmx_files_fw which is the 
            one with the --over-write flag, so if the rule is re run, the generated files do not duplicate
            entries.
            """
            input:
                #dmx="{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.no_unassigneds.fna" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["demultiplexing"]["add_unpair"] == "T"
                #else "{PROJECT}/runs/{run}/{sample}_data/splitLibsRC/seqs.assigned.fna",
                dmx="{PROJECT}/samples/{sample}/splitLibs_unpaired/seqs.assigned.rc.txt", #if demux order_by_strand else validate arriba 
                overw="{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                r2="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r1="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
            params:
                outdir="{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
            benchmark:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/demultiplex_fq_rc.benchmark"
            shell:
                "{config[java][command]}  -cp Scripts DemultiplexQiime  --txt -a rv -b {config[demultiplexing][remove_bc]}  -d {input.dmx} -o {params.outdir} "
                "-r1 {input.r1} -r2 {input.r2} {config[demultiplexing][dmx_params]}"
                
        if config["primers"]["remove"].lower() == "cfg" and config["ANALYSIS_TYPE"] == "ASV":
            rule remove_primers_dmx_files_cfg:
                input:
                    #fw=expand("{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run),
                    #rv=expand("{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run) #if int(config["demultiplexing"]["bc_mismatch"]) < 1
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                    rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
                    #else "{PROJECT}/runs/{run}/{sample}_data/demultiplexed/summary_wf.txt"
                params:
                    #"{PROJECT}/runs/{run}/{sample}_data/demultiplexed",
                    "{PROJECT}/samples/{sample}/demultiplexed" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired" ,
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    "PE",
                    "{PROJECT}/runs/"+ run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/remove_primers_fq.benchmark" 
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
                script:
                    "Scripts/removePrimersDemultiplex_cfg.py"
                #shell:
                    #"Scripts/removePrimersDemultiplex.sh {params} fastq.gz {config[demultiplexing][primers][fw_primer]} {config[demultiplexing][primers][rv_primer]} {config[demultiplexing][primers][min_overlap]} {config[demultiplexing][primers][extra_params]}"
        elif config["primers"]["remove"].lower() == "metadata"  and config["ANALYSIS_TYPE"] == "ASV":
            rule remove_primers_dmx_files_metadata:
                input:
                    config="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt" if config["demultiplexing"]["demultiplex"] == "T"
                    else config["metadata"],
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                    rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
                params:
                    "{PROJECT}/samples/{sample}/demultiplexed"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired",
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    "PE",
                    "{PROJECT}/runs/"+ run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/remove_primers_fq.benchmark"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
                script:
                    "Scripts/removePrimersDemultiplex.py"
        else:
            rule skip_remove_primers_dmx_files_PE:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
                    rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                    else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt" 
                shell:
                    "cat {input.fw} {input.rv} > {output}"
    elif (config["demultiplexing"]["create_fastq_files"] == "T" or config["ANALYSIS_TYPE"] == "ASV") and config["LIBRARY_LAYOUT"] == "SE":
        rule write_dmx_files_fw_SE:
            input:
                dmx="{PROJECT}/samples/{sample}/splitLibs/seqs.assigned.fna",
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
            params:
                outdir="{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
            benchmark:
                "{PROJECT}/samples/{sample}/demultiplexed/demultiplex_fq.benchmark"
            shell:
                "{config[java][command]}  -cp Scripts DemultiplexQiime --over-write --fasta -a fw -b {config[demultiplexing][remove_bc]}  -d {input.dmx} -o {params.outdir} "
                "-r {input.r1}  {config[demultiplexing][dmx_params]}"
        if config["primers"]["remove"].lower() == "cfg"  and config["ANALYSIS_TYPE"] == "ASV":
            rule remove_primers_dmx_files_SE_cfg:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
                params:
                    "{PROJECT}/samples/{sample}/demultiplexed",
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    "SE",
                    "{PROJECT}/runs/"+run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
                script:
                    "Scripts/removePrimersDemultiplex_cfg.py"
                #shell:
                #    "Scripts/removePrimersDemultiplexSE.sh {params} fastq.gz {config[demultiplexing][primers][fw_primer]} {config[demultiplexing][primers][min_overlap]} {config[demultiplexing][primers][extra_params]}"
        elif config["primers"]["remove"].lower() == "metadata"  and config["ANALYSIS_TYPE"] == "ASV":
            rule remove_primers_dmx_files_SE_metadata:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
                params:
                    "{PROJECT}/samples/{sample}/demultiplexed",
                    config["primers"]["extra_params"],
                    "fastq.gz",
                    "SE",
                    "{PROJECT}/runs/"+run+"/report_files/cutadapt.{sample}.fastq_summary.tsv",
                    run
                benchmark:
                    "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
                    #"{PROJECT}/runs/{run}/{sample}_data/demultiplexed/no_primer/cutadapt.log"
                script:
                    "Scripts/removePrimersDemultiplex.py"
        else:
            rule skip_remove_primers_dmx_files_SE:
                input:
                    fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt"
                output:
                    "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
                shell:
                    "cat {input.fw} > {output}"


    else:
        rule skip_dmx_file_creation_OTU:
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/" if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            shell:
                "touch {output}"

#if demultiplex{}
else:
    rule skip_demultiplexing_fq2fasta_OTU:
        input:
            #extended_reads="{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["LIBRARY_LAYOUT"] != "SE" 
            #else "{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq" if config["UNPAIRED_DATA_PIPELINE"] == "T" and config["LIBRARY_LAYOUT"] != "SE"
            #else "{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
            extended_reads="{PROJECT}/samples/{sample}/peared/seqs.assembled.fastq" if config["UNPAIRED_DATA_PIPELINE"] != "T" and config["LIBRARY_LAYOUT"] != "SE"
            else "{PROJECT}/samples/{sample}/peared/seqs.assembled.UNPAIRED.fastq" if config["UNPAIRED_DATA_PIPELINE"] == "T" and config["LIBRARY_LAYOUT"] != "SE"
            else "{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"            
            #tmp_pear_validation="{PROJECT}/samples/{sample}/peared/pear.log.validation",
            #tmp_pear_fq_validation="{PROJECT}/runs/{run}/{sample}_data/peared/qc/fq_fw_internal_validation.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
        shell:
            #to change vsearch --fastq_filter N_debres_r1.txt -fastaout test.txt
            #"fq2fa {input.extended_reads} {output}"
            #"sed -n '1~4s/^@/>/p;2~4p' {input.extended_reads} >  {output}"
            "zcat {input.extended_reads} | sed -n '1~4s/^@/>/p;2~4p'  | "
            "awk  '{{if($0 ~ \"^>\"){{seq=seq+1;print \">{wildcards.sample}_\"seq\" \"substr($1,2)}}else{{print $0}}}}' > {output}"
            if config["gzip_input"] == "T" and config["LIBRARY_LAYOUT"] == "SE" else
            "sed -n '1~4s/^@/>/p;2~4p' {input.extended_reads} | "
            "awk  '{{if($0 ~ \"^>\"){{seq=seq+1;print \">{wildcards.sample}_\"seq\" \"substr($1,2)}}else{{print $0}}}}' > {output}"
            
    if config["ANALYSIS_TYPE"] == "ASV" and config["LIBRARY_LAYOUT"] != "SE" and config["UNPAIRED_DATA_PIPELINE"] != "T":
        rule skip_dmx_file_creation_ASV_PE:
            input:
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz",
                r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/rv.fastq.gz"
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.txt_tmp" 
                if config["primers"]["remove"].lower() != "f" else
                "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq "
                if config["gzip_input"] == "F" else
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq.gz "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq.gz "

 
    elif config["ANALYSIS_TYPE"] == "ASV" and config["LIBRARY_LAYOUT"] == "SE":
        rule skip_dmx_file_creation_ASV_SE:
            input:
                r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq" if config["gzip_input"] == "F" else "{PROJECT}/samples/{sample}/rawdata/fw.fastq.gz"
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.txt_tmp"
                if config["primers"]["remove"].lower() != "f" else
                "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                if config["gzip_input"].lower() == "f" else
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq.gz"

    elif config["ANALYSIS_TYPE"] == "ASV" and config["LIBRARY_LAYOUT"] != "SE" and config["UNPAIRED_DATA_PIPELINE"] == "T":
         rule skip_dmx_file_creation_ASV_Unpaired:
            input:
                r1="{PROJECT}/samples/{sample}/peared/seqs.unassembled.forward.fastq",
                r2="{PROJECT}/samples/{sample}/peared/seqs.unassembled.reverse.fastq" 
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.txt_tmp"
                if config["primers"]["remove"].lower() != "f" else
                "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            shell:
                "touch {output} &&  ln -s $PWD/{input.r1} {params}{wildcards.sample}_1.fastq "
                " && ln -s $PWD/{input.r2} {params}{wildcards.sample}_2.fastq "


    else:
        rule skip_dmx_file_creation_ASV:
            params:
                "{PROJECT}/samples/{sample}/demultiplexed/"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.txt_tmp"
                if config["primers"]["remove"].lower() != "f"  and  config["UNPAIRED_DATA_PIPELINE"] != "T"
                else "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"
                if config["primers"]["remove"].lower() == "f"  and  config["UNPAIRED_DATA_PIPELINE"] != "T"
                else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.txt_tmp"
                if config["primers"]["remove"].lower() != "f"  and  config["UNPAIRED_DATA_PIPELINE"] == "T"
                else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            shell:
                "touch {output}" 

    if config["primers"]["remove"].lower() == "cfg":
        rule remove_primers_fq_files_cfg:
            input:
                fw="{PROJECT}/samples/{sample}/demultiplexed/summary.txt_tmp"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.txt_tmp"
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
                "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/remove_primers_fq.benchmark"
            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            script:
                "Scripts/removePrimersDemultiplex_cfg.py"
    elif config["primers"]["remove"].lower() == "metadata":
        rule remove_primers_fq_files_metadata:
            input:
                config="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt" if config["demultiplexing"]["demultiplex"] == "T"
                else config["metadata"],
                fw="{PROJECT}/samples/{sample}/demultiplexed/summary.txt_tmp"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else  "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.txt_tmp"
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
                "{PROJECT}/samples/{sample}/demultiplexed/remove_primers_fq.benchmark"

            output:
                "{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
                else   "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt"
            script:
                "Scripts/removePrimersDemultiplex.py"

if config["LIBRARY_LAYOUT"] != "SE":
    rule fastqCheckSum_PE:
        input:
            fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt",
            rv="{PROJECT}/samples/{sample}/demultiplexed/summary_rv.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_rv.txt"
        params:
            "{PROJECT}/samples/{sample}/demultiplexed/" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
        shell:
            "for file in {params}*.fastq.gz; do"
            "  fq=$(echo $file | awk -F\"/\" '{{print $NF}}');" 
            "  check=$(md5sum $file | cut -f1 -d \" \");"
            "  echo -e $check\"\\t\"$fq; "
            "done > {output}"   
else:
    rule fastqCheckSum_SE:
        input:
            fw="{PROJECT}/samples/{sample}/demultiplexed/summary_fw.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary_fw.txt"
        params:
            "{PROJECT}/samples/{sample}/demultiplexed/" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/"
        output:
            "{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
        shell:
            "for file in {params}*.fastq.gz; do"
            "  fq=$(echo $file | awk -F\"/\" '{{print $NF}}');" 
            "  check=$(md5sum $file | cut -f1 -d \" \");"
            "  echo -e $check\"\\t\"$fq; "
            "done > {output}"   

if config["ANALYSIS_TYPE"] == "ASV":
    #if config["dada2_filter"]["generateQAplots"] == "T":
    #    rule dada2_QA_Plots:
    #           input:
    #               expand("{PROJECT}/runs/{run}/{sample}_data/demultiplexed/summary.txt", PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
    #           output:
    #               "{PROJECT}/runs/{run}/asv/fw_QA_plots.pdf"
    #           benchmark:
    #               "{PROJECT}/runs/{run}/asv/qa_plots.benchmark"
    #           params:
    #               "{PROJECT}/runs/{run}/asv/"
    #           shell:
    #               "{config[Rscript][command]} Scripts/asvFilter.R $PWD {params} {input}"
    #else: 
    #    rule skip_QA_Plots:
    #           output:
    #               "{PROJECT}/runs/{run}/asv/no_qa_plots.txt"
    #           shell:
    #               "touch {output}"
    #rule dada2TruncationValues:
    #    input:
    #        "{PROJECT}/runs/{run}/asv/fw_QA_plots.pdf" if config["dada2_filter"]["generateQAplots"] == "T" else "{PROJECT}/runs/{run}/asv/no_qa_plots.txt"
 
    rule dada2Filter:
        input:
            expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]) 
            
        output:
            "{PROJECT}/runs/{run}/asv/filter_summary.out"
        benchmark:
            "{PROJECT}/runs/{run}/asv/filter.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/asvFilter.R $PWD " + str(config["dada2_filter"]["generateQAplots"]) + " " + str(config["dada2_filter"]["truncFW"]) + " " + str(config["dada2_filter"]["truncRV"]) + " "+str(config["dada2_filter"]["maxEE_FW"]) + " "+str(config["dada2_filter"]["maxEE_RV"]) + " " +str(config["dada2_filter"]["cpus"]) + " \"" +str(config["dada2_filter"]["extra_params"]) + "\" " +str(config["interactive"])+ " {output} " +config["primers"]["remove"] +" {input} " 

    rule validate_dada2Filter:
        input:
            "{PROJECT}/runs/{run}/asv/filter_summary.out"
        output:
            "{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        script:
            "Scripts/validateFilterASV.py"
    rule subset_for_errors_fw:
        input:
            expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
             "{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        output:
            temp("{PROJECT}/runs/{run}/asv/subset.fastq") 
        params:
            dir="{PROJECT}/runs/{run}/asv/"
            #ssize=config["dada2_asv"]["subset_size"]
        benchmark:
            "{PROJECT}/runs/{run}/asv/dada2.benchmark"
        shell:
            "cat test_out_demux/samples/NIOZ102/demultiplexed/filtered/*1.fastq.gz "
            "| gzip -d |  vsearch --fastx_subsample - --sample_size 200000  --fastaout test_out_demux/samples/NIOZ102/demultiplexed/filtered/subsample.fastq.gz"
             
    rule learn_errors:
        input:
            reads=expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
             flow="{PROJECT}/runs/{run}/asv/filter_summary.validation.txt"
        output:
            "{PROJECT}/runs/{run}/asv/errors.RData"
        params:
            dir="{PROJECT}/runs/{run}/asv/"
        benchmark:
            "{PROJECT}/runs/{run}/asv/dada2_err.benchmark"
        shell:
            "{config[Rscript][command]} Scripts/errDada2.R $PWD " +str(config["dada2_asv"]["cpus"]) + " "
            " "+ str(config["dada2_asv"]["generateErrPlots"])  + " {params.dir} "  + " "
            " "+ str(config["dada2_asv"]["nbases"]) +" " + str(config["binned_q_scores"])  + " "
            + "{input.reads}"
            
    rule run_dada2:
        input:
            reads=expand("{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
            else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt" , PROJECT=config["PROJECT"],sample=config["LIBRARY"]),
            wf="{PROJECT}/runs/{run}/asv/filter_summary.validation.txt",
            err="{PROJECT}/runs/{run}/asv/errors.RData"
        output:
            "{PROJECT}/runs/{run}/asv/stats_dada2.txt",
            "{PROJECT}/runs/{run}/asv/representative_seq_set.fasta",
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
            temp("{PROJECT}/runs/{run}/asv/dada2_asv_table.txt") 
        params:
            dir="{PROJECT}/runs/{run}/asv/",
            script="Scripts/asvDada2_bd.R"  if config["big_data_wf"].lower()  == "t"
            else "Scripts/asvDada2.R"
        benchmark:
            "{PROJECT}/runs/{run}/asv/dada2.benchmark"
        shell:
            "{config[Rscript][command]} {params.script} $PWD " +str(config["dada2_asv"]["pool"]) + " "
            " "+ str(config["dada2_asv"]["cpus"])  + " "+str(config["dada2_asv"]["generateErrPlots"]) + " "
            " "+ str(config["dada2_asv"]["extra_params"]) + " {params.dir} "  + " "+str(config["rm_reads"]["shorts"])  + " "
            " "+ str(config["rm_reads"]["longs"]) + " "+str(config["rm_reads"]["offset"])  + " "+str(config["dada2_asv"]["chimeras"])  + " "
            " "+str(config["dada2_taxonomy"]["db"]) + " "+str(config["dada2_taxonomy"]["add_sps"]["db_sps"])  + " "
            " "+str(config["dada2_taxonomy"]["add_sps"]["add"]) + " \""+str(config["dada2_taxonomy"]["extra_params"]) + "\" "
            " "+str(config["dada2_merge"]["minOverlap"]) +" "+str(config["dada2_merge"]["maxMismatch"])+" "
            " "+str(config["UNPAIRED_DATA_PIPELINE"]) +" " + " \""+str(config["dada2_taxonomy"]["add_sps"]["extra_params"]) + "\" "
            " "+str(config["interactive"]) + " " +str(config["rm_reads"]["non_interactive_behaviour"]) + " "  
            + "{input.reads}"

    rule asv_table:
        input:
            "{PROJECT}/runs/{run}/asv/dada2_asv_table.txt" 
        output:
            "{PROJECT}/runs/{run}/asv/asv_table.txt"
            # "{}" 
        shell:
            "cat {input} | awk '{{if(NR==1){{header=\"#OTU_ID\";for(i=1;i<=NF;i++){{header=header\"\\t\"$i}};print header}}else{{print $0}}}}'|   awk '{{ for (i=1; i<=NF; i++){{ a[NR,i] = $i }} }} NF>p {{ p = NF }} END {{ for(j=1; j<=p; j++) {{ str=a[1,j]; for(i=2; i<=NR; i++){{ str=str\"\\t\"a[i,j]; }} print str }} }}' > {output}"

    rule asv_tax_table:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
            "{PROJECT}/runs/{run}/asv/asv_table.txt" 
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt"
        benchmark:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/dada2.table.benchmark"
        shell:
            "cat {input[0]} | awk  -F \"\\t\" 'NR==FNR{{if(NR>1){{tax=$2;for(i=3;i<=NF;i++){{tax=tax\";\"$i}};h[$1]=tax;}}next;}} {{if(FNR==1){{print $0\"\\ttaxonomy\"}}else{{print $0\"\\t\"h[$1]}}}}' -  {input[1]} > {output}" 

    rule asv_to_biom:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt" 
        output:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom"
        benchmark:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/dada2.biom.benchmark"

        shell:
            "{config[biom][command]} convert -i {input[0]} -o {output} --table-type \"OTU table\" --to-hdf5 --process-obs-metadata taxonomy "

if config["align_vs_reference"]["align"] == "T":
    rule align_vs_reference:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/"
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
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
            output:
                temp("{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.fna")
            shell:
                "vsearch --fastx_revcomp {input} --fastaout {output}"
        rule concat_seqs:
            input:
                sq="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna",
                rc="{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.fna"
            output:
                temp("{PROJECT}/runs/{run}/{sample}_data/seqs_revcomplement.concat.fna")
            shell:
                "cat {input.sq} {input.rc} > {output}"

        rule cutadapt:
            input:
                #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                #else 
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna",
                config["metadata"]
            output:
                out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"
                #log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
            benchmark:
                "{PROJECT}/runs/{run}/{sample}_data/cutadapt.benchmark"
            params:
                config["primers"]["extra_params"],
                "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt",
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                "{PROJECT}/runs/{run}/report_files/cutadapt.{sample}.summary.tsv",
                "{sample}",
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
            script:
                "Scripts/remove_adapters_by_sample.py"

    else:
        if config["primers"]["remove"].lower() == "metadata":
            rule cutadapt:
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna",
                    "{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt" if config["demultiplexing"]["demultiplex"] == "T" 
                    else config["metadata"]
                output:
                    out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna",
                   # log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                benchmark:
                    "{PROJECT}/runs/{run}/{sample}_data/cutadapt.benchmark"
                params:
                    config["primers"]["extra_params"],
                    "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt",
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                    "{PROJECT}/runs/{run}/report_files/cutadapt.{sample}.summary.tsv",
                    "{PROJECT}/runs/{run}/{sample}_data/cutadapt_tmp/",
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                script:
                     "Scripts/remove_adapters_v2.py" # && ln -s ../../report_files/cutadapt.{wildcards.sample}.summary.tsv {params[4]}   "
        elif config["primers"]["remove"].lower() == "cfg":
            rule cutadapt:
                input:
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["align_vs_reference"]["align"] == "T"
                    else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
                    
                output:
                    out="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"
                    #log="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"

                benchmark:
                    "{PROJECT}/runs/{run}/{sample}_data/cutadapt.benchmark"
                params:
                    config["primers"]["extra_params"],
                    "{PROJECT}/runs/{run}/report_files/primers.{sample}.txt", 
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_removed.fna",
                    "{PROJECT}/runs/{run}/report_files/cutadapt.{sample}.summary.tsv",
                    "{PROJECT}/runs/{run}/{sample}_data/cutadapt_tmp/",
                    "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.log"
                script:
                    "Scripts/remove_adapters_v2.py"          


if config["align_vs_reference"]["align"] == "T":
    rule degap_alignment:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align" if config["primers"]["remove"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna" if config["primers"]["remove"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna"
        shell:
            "degapseq {input} {output}"

    rule single_line_fasta:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna" if config["primers"]["remove"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.oneline.fna" if config["primers"]["remove"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.oneline.fna"
        shell:
            "{config[java][command]} -cp Scripts FastaOneLine -f {input} -m 1 --write-discarded -o {output}"


#Creates file with sequence length ditribution
rule histogram:
    input:
        #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_nc.fna",
        #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna" #for tmp
        *selectInput()
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_lengths.txt"),
        temp("{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_hist.txt")
    shell:
        "cat {input} | grep -v '^>' | awk '{{print length}}' > {output[0]} "
        "&&  sort -g {output[0]} | uniq -c > {output[1]}"

rule histogram_chart:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_hist.txt",
        "{PROJECT}/runs/{run}/{sample}_data/seqs_accepted_lengths.txt"
    output:
        "{PROJECT}/runs/{run}/report_files/seqs_dist_hist.{sample}.png",
        "{PROJECT}/runs/{run}/{sample}_data/seqs_statistics.txt"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/",
        "{PROJECT}/runs/{run}/report_files/",
        "{sample}"
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
        "{PROJECT}/runs/{run}/{sample}_data/filter.benchmark"
    script:
        "Scripts/rmShortLong.py"


if config["chimera"]["search"] == "T":
#check for chimeric sequences
    if config["chimera"]["method"] == "usearch61":
        rule search_chimera:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
            #if (config["align_vs_reference"]["align"] == "T" and config["cutAdapters"] == "T") "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.no_adapter.degapped.fna"
            #elif (config["align_vs_reference"]["align"] == "T" and config["cutAdapters"] != "T") "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.align.degapped.fna"
            #elif (config["align_vs_reference"]["align"] != "T" and config["cutAdapters"] == "T") "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted_no_adapters.fna"
            #else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_accepted.fna"
            #*selectInput()
            output:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.txt"
            params:
                "{PROJECT}/runs/{run}/{sample}_data/chimera"
            benchmark:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimera.benchmark"
            shell:
                "{config[qiime][path]}identify_chimeric_seqs.py -m {config[chimera][method]} -i {input}  -o {params} --threads {config[chimera][threads]} {config[chimera][extra_params]}"
            #remove chimeras
    else:
        rule search_chimera_vsearch:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
            output:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.summary.txt"
            benchmark:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimera.benchmark"
            shell:
                "vsearch --{config[chimera][method]} {input} --threads {config[chimera][threads]} {config[chimera][extra_params]} --uchimeout {output}"
        rule filter_chimera_vsearch:
            input:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.summary.txt"
            output:
                "{PROJECT}/runs/{run}/{sample}_data/chimera/chimeras.txt"
            shell:
                "cat {input} | awk '$18==\"Y\"{{print $2\"\\t\"$1}}' > {output}"
    rule remove_chimera:
        """
        This script counts the sequence in input[0], the chimeras in input[2] if it is interactive it
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
        script:
            "Scripts/remove_chimera.py"


if config["demultiplexing"]["demultiplex"] == "T" and  config["ANALYSIS_TYPE"] != "ASV":
    rule count_samples_final:
        input:
            fasta="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",
            metadata="{PROJECT}/metadata/sampleList_mergedBarcodes_{sample}.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
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
        shell:
            "cat {input[0]} | awk 'NR==FNR{{if(NR>1){{h[$1]=$2;}}next}}{{if(FNR>1){{print $1\"\t\"h[$1]}}}}' - {input[1]} > {output}"
elif config["ANALYSIS_TYPE"] == "ASV" and config["demultiplexing"]["demultiplex"] != "T":

    rule count_samples_final:
        input:
            "{PROJECT}/runs/{run}/asv/filter_summary.out",
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        shell:
            "cat {input} | awk 'NR>1{{print $1\"\\t\"$2}}' > {output}"


else:
   rule count_samples_final:
        input:
            fasta="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T"
            else "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
        shell:
            "echo {wildcards.sample}\\t 100 > {output}"


rule distribution_chart:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt"
    output:
        "{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered.{sample}.dist.png"
    params:
        "{PROJECT}/runs/{run}/report_files/"
    shell:
        #"Scripts/sampleDist.py"
        "python  Scripts/sampleDist.manual.py {input} {output}"

if config["ANALYSIS_TYPE"]=="ASV":
    rule create_sample_log:
        input:
            ff=expand("{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.dist.txt",PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run) 
        output:
            "{PROJECT}/runs/{run}/samples.log"
        script:
            "Scripts/combineAllReads_asv.py"

else:
    rule combine_filtered_samples:
        input:
            allFiltered = expand("{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered_nc.fasta" if config["chimera"]["search"] == "T" 
            else  "{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta",PROJECT=config["PROJECT"],sample=config["LIBRARY"], run=run)
        output:
            "{PROJECT}/runs/{run}/seqs_fw_rev_combined.fasta",
            "{PROJECT}/runs/{run}/samples.log"
        benchmark:
            "{PROJECT}/runs/{run}/combine_seqs_fw_rev.benchmark"
        script:
            "Scripts/combineAllReads.py"

#Dereplicate
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
            "{PROJECT}/runs/{run}/derep/derep.benchmark"
        shell:
            "{config[derep][vsearch_cmd]} --derep_fulllength {input} --output {output[0]} --uc {output[1]} --strand {config[derep][strand]} "
            "--fasta_width 0 --minuniquesize {config[derep][min_abundance]} --sizeout" if config["pickOTU"]["m"] == "swarm"
            else "{config[derep][vsearch_cmd]} --derep_fulllength {input} --output {output[0]} --uc {output[1]} --strand {config[derep][strand]} "
            "--fasta_width 0 --minuniquesize {config[derep][min_abundance]}" 

if  config["pickOTU"]["m"] == "swarm" :
    rule cluster_OTUs:
        input:
            "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.fasta" 
        output:
            swarms="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt",
           # rep_seqs="{PROJECT}/runs/{run}/otu/representative_seq_set_swarm.fasta", (in case we want to generate reps -w {output.rep_seqs})
            uc="{PROJECT}/runs/{run}/otu/swarms.uc" 
        params:
            otuDir="{PROJECT}/runs/{run}/otu/"
        benchmark:
            "{PROJECT}/runs/{run}/otu.benchmark"
        #-i  {params.otuDir}swarm.struct
        shell:
            "swarm  -s {params.otuDir}swarm.stats -d {config[pickOTU][s]} -z "
            "-o {output.swarms}  -u {output.uc}  -t {config[pickOTU][cpus]} "
            "{config[pickOTU][extra_params]} < {input} "

else:
#cluster OTUs
    rule cluster_OTUs:
        input:
            "{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.fasta" if config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "swarm"
            and config["pickOTU"]["m"] != "usearch" else "{PROJECT}/runs/{run}/seqs_fw_rev_combined.fasta"
            #"{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
        output:
            "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt" if config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "swarm"
            and config["pickOTU"]["m"] != "usearch" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
        params:
            trieDir="{PROJECT}/runs/{run}/otu/"
        benchmark:
            "{PROJECT}/runs/{run}/otu.benchmark"
        shell:
            "{config[qiime][path]}pick_otus.py -m {config[pickOTU][m]} -i {input} "
            "-o {params.trieDir}  -s {config[pickOTU][s]} --threads {config[pickOTU][cpus]} {config[pickOTU][extra_params]} "

if config["derep"]["dereplicate"] == "T"  and config["pickOTU"]["m"] != "usearch" and config["pickOTU"]["m"] != "swarm":
    rule remap_clusters:
        input:
            otu_txt="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt",
            uc_derep="{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.uc"
        output:
            map="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt",
            log="{PROJECT}/runs/{run}/otu/remap.log"
        benchmark:
            "{PROJECT}/runs/{run}/derep/remap.benchmark"
        shell:
            "{config[java][command]} -cp Scripts/ClusterMapper/build/classes clustermapper.ClusterMapper uc2otu "
            "-uc {input.uc_derep} -otu {input.otu_txt} -o {output.map} > {output.log}"
elif config["pickOTU"]["m"] == "swarm":
    rule remap_clusters:
        input:
            uc_swarm="{PROJECT}/runs/{run}/otu/swarms.uc",
            uc_derep="{PROJECT}/runs/{run}/derep/seqs_fw_rev_combined_derep.uc"
        output:
            map="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt",
            log="{PROJECT}/runs/{run}/otu/remap.log"
        benchmark:
            "{PROJECT}/runs/{run}/derep/remap.benchmark"
        shell:
            "{config[java][command]} -cp Scripts/ClusterMapper/build/classes clustermapper.ClusterMapper uc2uc "
            "-uc {input.uc_derep} -uc2 {input.uc_swarm} -o {output.map} --full-uc --relabel -l OTU -lidx 1 > {output.log}"


#pick representative OTUs
rule pick_representatives:
    input:
        otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt" 
        if (config["derep"]["dereplicate"] == "T"  and config["pickOTU"]["m"] != "usearch" and config["pickOTU"]["m"] != "swarm")
        else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt" if config["pickOTU"]["m"] == "swarm" 
        else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt",
        filtered="{PROJECT}/runs/{run}/seqs_fw_rev_combined.fasta"
        #filtered="{PROJECT}/runs/{run}/{sample}_data/seqs_fw_rev_filtered.fasta"
    output:
        reps="{PROJECT}/runs/{run}/otu/representative_seq_set.fasta",
        log="{PROJECT}/runs/{run}/otu/representative_seq_set.log"
    benchmark:
        "{PROJECT}/runs/{run}/pick_reps.benchmark"
    shell:
        "{config[qiime][path]}pick_rep_set.py -m {config[pickRep][m]} -i {input.otus} "
        "-f {input.filtered} -o {output.reps} --log_fp {output.log} {config[pickRep][extra_params]}"


#assign taxonomy to representative OTUs
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.benchmark"
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
            otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_derep_otus.txt" if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch")
            or config["pickOTU"]["m"] == "swarm" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.complete.out"
        benchmark:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/complte_blast.assign_taxa.benchmark"
        shell:
            "cat {input.blastout} | cut -f1 | sort | uniq | grep -v -w -F -f - {input.otus} "
            "| awk '{{print $1\"\\tUnassigned\\t-\\t-\\t-\\t-\"}}' | cat {input.blastout} - > {output}"

    rule prepare_blast_for_stampa:
        """
         Take completed blast output file and make some reformat in order to
         fulfill stampa format.
        """
        input:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_blastn.complete.out"
        output:
            temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/hits.blast.out")
        benchmark:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.hits.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt"
        benchmark:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.map.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.hits.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
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
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
        shell:
            "Scripts/stampa_merge.py {params} {config[assignTaxonomy][vsearch][taxo_separator]}"
    rule normalize_taxo_out:
        """
         Normalize the output in terms oif format and names in order to be able
         to continue with the pipeline
        """
        input:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.out" if config["assignTaxonomy"]["map_lca"].lower() == "t"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/results.vsearch.no_stampa.out"  
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt"
        benchmark:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/stampa.benchmark"
        shell:
            "cat {input} |  awk -F\"\\t\" '{{print $1\"\\t\"$4\"\\t\"$3\"\\t\"$5}}' | sed 's/N;o;_;h;i;t/Unassigned/' > {output}"

else:
    rule assign_taxonomy:
        input:
            "{PROJECT}/runs/{run}/otu/representative_seq_set.fasta"
        output:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt",
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.log"
        params:
            outdir="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/"
        benchmark:
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/assign_taxa.benchmark"
        shell:
            "{config[qiime][path]}parallel_assign_taxonomy_{config[assignTaxonomy][qiime][method]}.py -i {input} --id_to_taxonomy_fp {config[assignTaxonomy][qiime][mapFile]} "
            "{config[assignTaxonomy][qiime][dbType]} {config[assignTaxonomy][qiime][dbFile]} --jobs_to_start {config[assignTaxonomy][qiime][jobs]} "
            "--output_dir {params.outdir}  {config[assignTaxonomy][qiime][extra_params]}"

rule make_otu_table:
    input:
        tax="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt",
        otus="{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_remapped_otus.txt" if (config["derep"]["dereplicate"] == "T" and config["pickOTU"]["m"] != "usearch") 
        or config["pickOTU"]["m"] == "swarm" else "{PROJECT}/runs/{run}/otu/seqs_fw_rev_combined_otus.txt"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.tmp.biom" 
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom.benchmark"
    shell:
        "{config[qiime][path]}make_otu_table.py -i {input.otus} -t {input.tax} -o {output} {config[makeOtu][extra_params]}"

rule prapare_obs_metadata:
    '''
    Prepare observation metadata 
    '''
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_tax_assignments.txt"
    output:
        temp("{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otu.metadata.tsv")
    shell:
        "cat {input} | awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\";print \"#OTUID\\tIdentity\\tACCs\"}}{{print $1,$3,$4}}' > {output}" 

rule add_obs_metadata:
    input:
        metadata="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otu.metadata.tsv",
        otu="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.tmp.biom"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom"
    shell:
        "{config[biom][command]} add-metadata -i {input.otu} -o {output} --observation-metadata-fp {input.metadata} --float-fields Identity"

#filter OTU table
rule summarize_taxa:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/otuTable_L6.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/asvTable_L6.txt"
    params:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/" 
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/summarize_taxa.benchmark"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/summarize_taxa.benchmark"
    shell:
        "{config[qiime][path]}summarize_taxa.py -i {input} -o {params} -a  {config[summTaxa][extra_params]}"

#Converts otu table from biom format to tsv
rule convert_table:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom",
        #"{PROJECT}/runs/{run}/otu/taxa_"+config["assignTaxonomy"]["method"]+"/"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt"
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt.benchmark"
    shell:
        "{config[biom][command]} convert -i {input[0]} -o {output} {config[biom][tableType]} "
        "{config[biom][headerKey]} {config[biom][outFormat]} {config[biom][extra_params]}"

#filter OTU table
rule filter_otu:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.biom"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.biom"
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_nosingletons.bio.benchmark"
    shell:
        "{config[qiime][path]}filter_otus_from_otu_table.py -i {input} -o {output} -n {config[filterOtu][n]} {config[filterOtu][extra_params]}"

#Convert to/from the BIOM table format.
rule convert_filter_otu:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.biom"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt"
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_nosingletons.txt.benchmark"
    shell:
        "{config[biom][command]} convert -i {input} -o {output} {config[biom][tableType]} "
        "{config[biom][headerKey]} {config[biom][outFormat]} {config[biom][extra_params]}"

#filter asv table
rule filter_asv:
    input:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom"
    output:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.biom"
    benchmark:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_nosingletons.bio.benchmark"
    shell:
        "{config[qiime][path]}filter_otus_from_otu_table.py -i {input} -o {output} -n {config[filterOtu][n]} {config[filterOtu][extra_params]}"

#Convert to/from the BIOM table format.
rule convert_filter_asv:
    input:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.biom"
    output:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt"
    benchmark:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_nosingletons.txt.benchmark"
    shell:
        "{config[biom][command]} convert -i {input} -o {output} {config[biom][tableType]} "
        "{config[biom][headerKey]} {config[biom][outFormat]} {config[biom][extra_params]}"

#Summarize singletons
rule summarize_taxa_no_singletons:
    input:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.biom"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.biom"
    output:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/otuTable_noSingletons_L6.txt"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/asvTable_noSingletons_L6.txt"
    params:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/"
    benchmark:
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/summarize_taxa.benchmark"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/summarize_taxa.benchmark"
    shell:
        "{config[qiime][path]}summarize_taxa.py -i {input} -o {params} -a {config[summTaxa][extra_params]}"

if  config["krona"]["report"].casefold() == "t" or config["krona"]["report"].casefold() == "true":
    rule krona_report:
        input:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.txt" if config["ANALYSIS_TYPE"] == "ASV" and config["krona"]["otu_table"].casefold() != "singletons"  
            else "{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.txt" 
            if config["ANALYSIS_TYPE"] == "ASV"
            else
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt" 
            if config["krona"]["otu_table"].casefold() != "singletons"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt"
        output:
            "{PROJECT}/runs/{run}/report_files/krona_report.dada2.html"
            if config["ANALYSIS_TYPE"] == "ASV" else
            "{PROJECT}/runs/{run}/report_files/krona_report."+config["assignTaxonomy"]["tool"].lower()+".html"
        params:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/" 
            if config["ANALYSIS_TYPE"] == "ASV"  else 
            "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/" 
        benchmark:
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/krona_report.benchmark" if config["ANALYSIS_TYPE"] == "ASV"
            else  "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/krona_report.benchmark"
        script:
            "Scripts/otu2krona.py"
else:
    rule skip_krona:
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
        otuNoSingleton="{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable_noSingletons.biom"
        if config["ANALYSIS_TYPE"] == "ASV" else
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.biom"
    output:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta"
         if config["ANALYSIS_TYPE"] == "ASV" else
         "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_noSingletons.fasta"
    benchmark:
        "{PROJECT}/runs/{run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.benchmark"
        if config["ANALYSIS_TYPE"] == "ASV" else
        "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/representative_seq_set_noSingletons.benchmark"
    shell:
        "{config[qiime][path]}filter_fasta.py -f {input.fastaRep} -o {output} -b {input.otuNoSingleton} {config[filterFasta][extra_params]}"
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
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/align_rep_seqs.benchmark" if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/align_rep_seqs.benchmark"
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
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/align_rep_seqs.benchmark"  if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/align_rep_seqs.benchmark"
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
            "{PROJECT}/runs/{run}/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.benchmark"  if config["ANALYSIS_TYPE"] == "ASV"
            else "{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.benchmark"
        shell:
            "{config[qiime][path]}make_phylogeny.py -i {input} -o {output} -t {config[makeTree][method]} {config[makeTree][extra_params]}"

if config["ANALYSIS_TYPE"] != "ASV": 
    rule report_all:
       input:
            a="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable.txt",
#            if config["ANALYSIS_TYPE"] == "OTU" else
#            "{PROJECT}/runs/{run}/asv/asv_table.biom",
            b="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary/otuTable_L6.txt",
#            if config["ANALYSIS_TYPE"] == "OTU" else
#            "{PROJECT}/runs/{run}/asv/summary/asv_table_L6.txt",
            c="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/otuTable_noSingletons.txt",
            d="{PROJECT}/runs/{run}/otu/taxonomy_"+config["assignTaxonomy"]["tool"].lower()+"/summary_noSingletons/otuTable_noSingletons_L6.txt",
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
            g="{PROJECT}/runs/{run}/samples.log"
       output:
            temp("{PROJECT}/runs/{run}/reporttmp_all.html")
       benchmark:
            "{PROJECT}/runs/{run}/report_all.benchmark"
       script:
            "Scripts/report_all_v2.py"
else:
    rule report_all_asv:
       input:
            a="{PROJECT}/runs/{run}/asv/taxonomy_dada2/asvTable.biom",
            b="{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary/asvTable_L6.txt",
            d="{PROJECT}/runs/{run}/asv/taxonomy_dada2/summary_noSingletons/asvTable_noSingletons_L6.txt",
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
            g="{PROJECT}/runs/{run}/samples.log"
       output:
            temp("{PROJECT}/runs/{run}/reporttmp_all.html")
       benchmark:
            "{PROJECT}/runs/{run}/report_all.benchmark"
       script:
            "Scripts/report_all_asv.py"


rule tune_report_all:
    input:
        "{PROJECT}/runs/{run}/reporttmp_all.html"
    output:
        "{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv_report_dada2.html"
    script:
        "Scripts/tuneReport.py"
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
        shell:
            "{config[wkhtmltopdf_command]}  {input} {output}"

rule report:
    input:
        report_all="{PROJECT}/runs/{run}/otu_report_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/asv_report_dada2.html",
        dist_chart="{PROJECT}/runs/{run}/report_files/seqs_fw_rev_filtered.{sample}.dist.png",
        dmxFiles="{PROJECT}/samples/{sample}/demultiplexed/summary.pcr.txt" if config["UNPAIRED_DATA_PIPELINE"] != "T"
        else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/summary.pcr.txt",
        md="{PROJECT}/samples/{sample}/demultiplexed/md5sum.txt"  if config["UNPAIRED_DATA_PIPELINE"] != "T"
        else "{PROJECT}/samples/{sample}/demultiplexed/unpaired/md5sum.txt"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/report.html")
    script:
        "Scripts/report_v2.py" if config["ANALYSIS_TYPE"] == "OTU" and config["LIBRARY_LAYOUT"] != "SE" 
        else "Scripts/report_v2_SE.py" if config["ANALYSIS_TYPE"] == "OTU" and config["LIBRARY_LAYOUT"] == "SE"  
        else "Scripts/report_asv.py"
rule tune_report:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/report.html"
    output:
        "{PROJECT}/runs/{run}/report_{sample}_"+config["assignTaxonomy"]["tool"].lower()+".html"
        if config["ANALYSIS_TYPE"] == "OTU" else
        "{PROJECT}/runs/{run}/report_{sample}_dada2.html"
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
    shell:
        "zip -r {output} {params} {input[0]} {input[1]}"
