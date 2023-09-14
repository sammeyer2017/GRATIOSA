import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GRATIOSA",
    version="0.0.2",
    author="Forquet Raphaël, Pineau Maïwenn, Meyer Sam",
    author_email="sam.meyer@insa-lyon.fr",
    description=" ",
    long_description=long_description,
    url="https://github.com/sammeyer2017/topo_database.git",
    packages=['GRATIOSA'],
    install_requires=[
    	'numpy',
        'matplotlib',
        'pandas',
        'pysam',
        'scipy',
        'statsmodels'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
)

print(setuptools.find_packages())
