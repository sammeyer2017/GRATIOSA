import os
import sys

# basedir is the path to the folder which contains the database.
try:
    basedir = os.environ['GRATIOSA_DB_PATH']
except KeyError:
    print("Please set the environment variable GRATIOSA_DB_PATH")
    sys.exit()

# resdir is the default output folder.
resdir = basedir + 'resdir/'
