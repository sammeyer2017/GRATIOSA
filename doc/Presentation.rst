Presentation
============

GRATIOSA
--------
Genome Regulation Analysis Tool Incorporating Organization and Spatial Architecture

GRATIOSA facilitates the integration, comparison and combined statistical analyses of diffent data types such as Chip-Seq data, RNA-Seq data and genomic annotations. It includes (1) the integration of various data types and formats in a unified framework, allowing direct comparison; (2) an implementation of spatial statistical analyses of the relationship between transcriptional regulation and chromosome organization, greatly reducing the required code length; and (3) an improvement in reproducibility, in particular across different species in order to identify common regulatory mechanisms.

GRATIOSA is written in Python, and is targeted to computational biologists. The automated import of data and standardized statistical tests and procedures are combined with the flexibility of Python for developing custom subsequent analyses, depending on the particular regulatory mechanism under investigation. It is primarily designed to treat expression data (RNA-Seq or microarrays), and ChIP-Seq data, but can be used with any type of continuous signals along the genome (other binding signals, mutation rates, binding prediction profiles…) or lists of discrete features (annotated or predicted protein binding sites, Hi-C topological domain borders, …). 

How to install the package? 
---------------------------
To install this package, the recommended method is to execute the following command in your terminal 
``pip install GRATIOSA``
For simple testing, you can just download this package manually (using Code / Download ZIP) and work locally. 

General presentation 
--------------------
We describe here the recommended installation procedure. For a local use (for simple testing), please follow instructions proposed in the provided tutorials (Jupyter notebooks). 

Before using the package, the user needs to prepare and organize their data into predefined directories: see examples in the provided data/ directory and their use in tutorials. Each new organism corresponds to a directory with the associated name.  This directory must contain the reference sequence in fasta format and an annotation file in gff format in the "annotation" subdirectory. This can be created automatically from the NCBI database (see documentation of the Genome class). The pre-processed and formatted experimental data, provided by the user, should be placed in appropriate directories named according to the data type (see examples in tutorials). These experimental data files should usually be accompanied by an "info" file that the user needs to complete with information about the file organization, following the provided template.

Typical data files supported by GRATIOSA are summarized in the following table. Please see the detailed documentation of each function to see if and how alternate formats can be used:
|Data type | Rows | Columns | Standard format | Commonly used packages |
|:-------------------:|:-------:|:---------------:|:-----:|:-----------:|
|RNA-Seq coverage|positions|coverage|bedgraph|Bedtools|
|RNA-Seq expression|genes|locus, expression level (RPKM or FPKM)|tab|Cufflinks, HTSeq|
|RNA-Seq differential expression|genes|locus, fold-change, p-value|csv|limma, DESeq2|
|Microarray differential expression|genes|locus, fold-change, p-value|csv|limma|
|ChIP-Seq coverage|positions|bin start, bin end, coverage|bedgraph|deepTools|
|ChIP-Seq peaks|sites|position, value|bed|MACS|
|list of genomic position data|sites|position, value|csv||

How to use GRATIOSA?
------------------------
The user must define the location of the database (containing all datafiles for the investigated species) as an environment variable \$GRATIOSA_DB_PATH (e.g. export GRATIOSA_DB_PATH = /home/usr/documents/GRATIOSA/). Thus, the files must be in \$GRATIOSA_DB_PATH + data/organisms. The most convenient way to run GRATIOSA is probably to install a virtual environment. See tutorials for examples of use. 

Data analysis using the package is performed through Python commands, typically in a Jupyter notebook, with three major steps. Firstly, as the package is an object-oriented framework, the objects (Genome, Transcriptome, ChIP-Seq, etc.) need to be initialized, and then the data can be loaded as attributes. 

The second step involves data processing. During this step, attributes associated with genomic positions can be scaled to the gene level for further analysis. Additionally, the loaded continuous signals can be binned, smoothed, or averaged. A verification and graphical exploration of the signals can be performed at the end of the first or second step using a graphical function that plots the signals on the annotated genome. To prepare for statistical analysis, quantitative data can also be classified.

The last step is the statistical analysis with enrichment or proportion tests (for qualitative attributes) and Student/Wilcoxon-Mann-Whitney tests (for quantitative comparisons). These functions can also be used to handle custom data imported manually by the user. Results are saved as tables (in csv format) and can be visualized as annotated bar plots created with graphical functions included in the package. 


Documentation
-------------
https://gratiosa.readthedocs.io/en/latest/Presentation.html
