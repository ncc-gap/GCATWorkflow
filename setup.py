#!/usr/bin/env python

from setuptools import setup, find_packages
from scripts.gcat_workflow import __version__

#import sys
#sys.path.append('./tests')

setup(
    name = 'gcat_workflow',
    version = __version__,
    description = 'Python tools for running gcat workflow for cancer genome and transcriptome sequencing analysis',
    keywords = 'cloud bioinformatics',
    author = 'Kenichi Chiba, Ai Okada and Yuichi Shiraishi',
    author_email = 'genomon.devel@gmail.com',
    url = 'https://github.com/ncc-ccat-gap/GCATWorkflow.git',
    license = 'GPLv3',
    
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Genome Analysis :: RNA-seq',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    
    package_dir = {'': 'scripts'},
    packages = find_packages("scripts"),
    package_data = {'gcat_workflow': ['*/data/*']},
    
    scripts = ['gcat_workflow', 'gcat_runner'],
    include_package_data = True,
    zip_safe = False,
    install_requires = [
        'snakemake',
        'pyyml',
        'drmaa',
    ],
    #test_suite = 'unit_tests.suite'
    test_suite = 'tests'
)
