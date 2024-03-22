from setuptools import setup, find_packages
from setuptools.command.install import install
import sysconfig
import os
import problin_libs 
from os import walk, listdir
from os.path import join, normpath, isfile

def recursive_list_dir(path):
    listing=[]
    for x in walk(path):
        if isfile(x[0]):
            listing.append(x[0].split(path+'/')[1])
        for y in listdir(x[0]):
            z = normpath(join(x[0],y))
            if isfile(z):
                listing.append(z.split(path+'/')[1])
    return listing


param = {
        'name': problin_libs.PROGRAM_NAME,
        'version': problin_libs.PROGRAM_VERSION,
        'description': problin_libs.PROGRAM_DESCRIPTION,
        'author': problin_libs.PROGRAM_AUTHOR,
        'license': problin_libs.PROGRAM_LICENSE,
        'packages': find_packages(),
        'include_package_data': True,
        'scripts': ['run_mollusc.py'],
        'zip_safe': True,
        'install_requires': ['scipy>=1.3.1','cvxpy>=1.4','treeswift>=1.1.39', 'Mosek>=10.1.16', 'setuptools'],        
        'keywords': 'Spatial Lineage Tracing, Cell motility, Cell division, Maximum likelihood, Phylogeography',
        'long_description': "Maximum Likelihood Estimation Of Lineage and Location Using Single Cell Spatial Lineage tracing Data",
        'classifiers': [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
        'entry_points': {
            'console_scripts': [
                'run_mollusc= run_mollusc:main',
             ],
        },
}

setup(**param)
