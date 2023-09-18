GRATIOSA
--------
Genome Regulation Analysis Tool Incorporating Organization and Spatial Architecture

GRATIOSA facilitates the integration, comparison and combined statistical analyses of diffent data types such as Chip-Seq data, RNA-Seq data and genomic annotations. It includes (1) the integration of various data types and formats in a unified framework, allowing direct comparison; (2) an implementation of spatial statistical analyses of the relationship between transcriptional regulation and chromosome organization, greatly reducing the required code length; and (3) an improvement in reproducibility, in particular across different species in order to identify common regulatory mechanisms.

GRATIOSA is written in Python, and is targeted to computational biologists. The automated import of data and standardized statistical tests and procedures are combined with the flexibility of Python for developing custom subsequent analyses, depending on the particular regulatory mechanism under investigation. It is primarily designed to treat expression data (RNA-Seq or microarrays), and ChIP-Seq data, but can be used with any type of continuous signals along the genome (other binding signals, mutation rates, binding prediction profiles…) or lists of discrete features (annotated or predicted protein binding sites, Hi-C topological domain borders, …). 


HOW TO USE THE PACKAGE ? 
------------------------
Before using the package, the user needs to prepare and organize their data and database. For each new organism, the user must create a folder with the organism's name. Inside this folder, a subfolder named "annotation" should be created and the user should add the sequence data in fasta format and an annotation file in gff format in this folder. The pre-processed and formatted experimental data, should also be placed in an appropriate folder, named according to the data type. These experimental data files should be accompanied by an info file that the user needs to complete with information about the file organization, following the info file template for that data type.

The data analysis using the package is performed, by executing Python commands, in three major steps. Firstly, as the package is an object-oriented framework, the objects (Genome, Transcriptome, ChIP-Seq, etc.) need to be initialized, and then the data can be loaded as attributes. In this first step, the sequence and annotation can be added to the Genome object, signals and sites can be added to the ChIP-Seq object, omics data can be added to the Transcriptome object, and functional annotations can be loaded onto the GO object.

The second step involves data processing. During this step, attributes associated with genomic positions can be scaled to the gene level for further analysis. Additionally, the loaded continuous signals can be binned, smoothed, or averaged. A verification and graphical exploration of the signals can be performed at the end of the first or second step using a graphical function that plots the signals on the annotated genome. To prepare for statistical analysis, quantitative data can also be classified.

The last step is the statistical analysis with enrichment or proportion tests (for qualitative attributes) and student's t-tests (for quantitative comparisons). These functions can also be used to treat new data imported by the user, as long as they are formatted as Python dictionaries . Tests results are saved as tables (in csv format) and can be visualized as annotated bar plots created with graphical functions included in the package. 


WARNING 
-------
Before using, please set the environment variable \$GRATIOSA_PATH in order to have \$GRATIOSA_PATH + data/organisms (e.g. \$GRATIOSA_PATH = /home/usr/documents/GRATIOSA/). The most convenient way to run GRATIOSA is probably to install a virtual environment.
