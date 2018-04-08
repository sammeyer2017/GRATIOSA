import os

try:  
    basedir = os.environ['TOPO_DATABASE_PATH']
except KeyError: 
    print "Please set the environment variable TOPO_DATABASE_PATH"
    sys.exit()

# --------- rpkm
# multiplicative constant btw rpkm values and coverage values: use same unit as rpkm for coverage
multconst=1.
# maximal value for ratio of rpkms to be considered as equal
maxrat=1.5
# minimal rpkm for non-zero
minrpkm=1.
# tresholds for log2(FC) where we consider a gene activated / repressed
fc_treshold_pos = 0.263 # FC > 1.2
fc_treshold_neg = -0.263 # fc < 1/1.2
pval_treshold = 0.05
# figure var
# for density circles
fig_width = 5;fig_height = 5
center_x = fig_width/2;center_y = fig_height/2
radius = fig_width / 1.4 # radius circles

# parameters for bTSS finder
PROM_LENGTH = 200
TSS_DOWNSTREAM = 50
#correspondent taxon for btssfinder
taxon={"e_coli_v3":"e","salmonella":"c"}


# --------- TSS
# error factor in small window coverage calculation:
# an increase by less than this factor at a TSS is not considered as significant
minfact=1.5
# minimum starts in TSS identification
minstart=10
# abortive region size
size_abortive=20
# buffer size for kon calculation
size_buffer=100
# maximal distance between gene and promoter TSS
maxTSSdist=300

# --------- zero expression
# window for considering zero expression: essential for defining TUs
zerowinl=50
zerowinm=1.


# -------------
# plotting

# constants
w=3.5
#h={("annot",.5), ("exp",3), ("mod",3)}
h=[.5,1.5,1.5]
dh=.05
# tick distances
tickd=[100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
# approximate tick nb
ticknb=4

# extensions
exts=[".pdf"] #,".eps"]
