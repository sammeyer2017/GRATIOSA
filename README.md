# topo_database
Before using, please set the environment variable $TOPO_DATABASE_PATH in order to have $TOPO_DATABASE_PATH + data/organisms (e.g. $TOPO_DATABASE_PATH = /home/usr/documents/topo_database/)

requirements.txt > created using pip freeze, contains all the required packages and their version. To be noted that Python 2.7.15rc1 version is required. The most convenient way to run topo DB is probably to install a virtual environment.

The topo database path contains one folder for each organism that we want to study (e.g. data/dickeya, data/ecoli). Inside each organism folder, there are additional folders which are called and required to load information of the organism using the topo DB (e.g. annotation folder which contains the annotation file of the organism, expression folder, TSS/TTS folder). Whatever the type of information that needs to be loaded, the organization remains the same in each folder: there is one .info file that summarizes where is the file in the current folder that contains the useful information and how it is organized (e.g. expression folder that allows to load genome-wide expression data, expression.info file where is indicated what is the filename containing expression data, where is the tag column, first expression column, whether it is log or not, what is the column separator).

The topo database itself contains several program files (i.e. .py files) that have been organized in an object-oriented manner. Most files define classes (e.g. genome class) together with their attributes (e.g. a genome object has some genes and a genomic sequence) and their methods (e.g. load genes to a genome from an annotation file in an annotation folder). In addition, useful_functions.py contains miscellaneous functions sometimes called by the topo DB (e.g. functions used to process .bam and other files from RNAseq data, automatically annotate pyplots with brackets for confidence intervals). btssfinder.py is called by a genome method which predicts promoters starting from TSS data. globvar.py contains all the global variables. plot_genome contains useful functions to plot informations of a genome. .gitignore contains the files that need to be ignored in git repository (e.g. files used for modeling, to perform analysis based on genomic objects, .pyc files, virtual environment folder venv).

Example:
%load_ext autoreload # automatically reload modules with %autoreload      
import sys  
sys.path.append("modeling_spacer") # add another path for interpreter to search
sys.path.append("raph")
import Genome # contain main functions
gen = Genome.Genome(name="dickeya") # define a genome object
gen.load_fc_pval() # load expression variation (FC) data
gen.genes["Dda3937_04419"].fc_pval # access FC data for Dda3937_04419 gene
gen.load_TSS() # load TSS data
gen.TSSs["TSS_tex"][3312932] # access TSS object at position 3312932 from TSS_tex list of TSSs

In the future, it would be useful to define an environment variable indicating the results folder, instead of exporting results within the database