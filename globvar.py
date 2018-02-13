basedir="/home/kwl06/Documents/stage_GB5/data/data_leger/"

# try:  
#    basedir = os.environ['TOPO_DATABASE_PATH']
# except KeyError: 
#    print "Please set the environment variable TOPO_DATABASE_PATH"
#    sys.exit()

# --------- rpkm
# multiplicative constant btw rpkm values and coverage values: use same unit as rpkm for coverage
multconst=1.
# maximal value for ratio of rpkms to be considered as equal
maxrat=1.5
# minimal rpkm for non-zero
minrpkm=1.

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
