import os
import sys

# basedir is the path to the topo_database folder which contains the database.
try:
    basedir = os.environ['TOPO_DATABASE_PATH']
except KeyError:
    print("Please set the environment variable TOPO_DATABASE_PATH")
    sys.exit()

# resdir is the default output folder.
resdir = basedir + 'resdir/'
