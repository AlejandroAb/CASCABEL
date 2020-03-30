## Assign taxonomy

In order to perform taxonomy assignment, **Cascabel** can work with different methods and tools.

### OTU Workflow

When the [analysis workflow](../../wiki#233-analysis-type) is OTU all the options for taxonomy assignment are given within the configuration file according to the _Assign taxonomy_ section. In such part, most of the methods can work with a fasta file and a corresponding mapping file between the accessions and taxonomy paths _tab_ separated.

In this sense, for environmental community analysis (Qiime's Silva database releases)[https://www.arb-silva.de/download/archive/qiime/] provides files ready to use within Cascabel. i.e. [Silva_132_release.zip](https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip)

For instance, within the previous example, you can find the following files:

* SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.fna _(fasta file)_
* SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/taxonomy_7_levels.txt _(taxonomy mapping file)_    

However, [Silva fasta files](https://ftp.arb-silva.de/) (without any post processing) need a switable mapping file. For this cases we supply two mapping files:

* taxmap_embl_ssu_ref_132_taxid.map.gz This corresponds to the fasta file from Silva 132 Small Sub Units (16S/18S) available for download [here](https://ftp.arb-silva.de/release_132/Exports/SILVA_132_SSURef_tax_silva.fasta.gz)
* taxmap_embl_lsu_ref_132_taxid.map.gz This corresponds to the fasta file from Silva 132 Long Sub Units (23S/28S) available for download [here](https://ftp.arb-silva.de/release_132/Exports/SILVA_132_LSURef_tax_silva.fasta.gz)  

_Please notice that these files are gzip and thus, need to be un compressed prior to its use._
_
As the mapping file format is pretty simple you can recreate these mapping files from the fasta files with the following command:

<pre><code class="text">
zcat SILVA_132_SSURef_tax_silva.fasta.gz | grep "^>" | sed 's/>// ; s/ /\t/' | awk -F"\t" '{print $1"\t"$2}'  > taxmap_embl_ssu_ref_132_taxid.map
</code></pre>

### ASV Workflow

When the analysis type is ASV, the taxonomy assignment is performed within dada2 library using the RDP classifier. In this case, the database to supply can be downloaded from this (DADA2-formatted references page)[https://benjjneb.github.io/dada2/training.html]
