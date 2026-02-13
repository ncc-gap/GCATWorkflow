[![UnitTest](https://github.com/ncc-ccat-gap/GCATWorkflow/actions/workflows/UnitTest.yml/badge.svg)](https://github.com/ncc-ccat-gap/GCATWorkflow/actions/workflows/UnitTest.yml)
![Python](https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9%20%7C%203.10-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

G-CAT Workflow
===============
G-CAT Workflow is a cancer genome and RNA sequencing data analysis pipeline that efficiently detects genomic variants and transcriptomic changes. Users can run G-CAT Workflow with ease of use in the supercomputer and also can run it on other high-performance computers.

## Manual
https://github.com/ncc-ccat-gap/GCATWorkflow

For developers, https://github.com/ncc-ccat-gap/GCATWorkflow/wiki

## Setup

0. Precondition

Make DRMAA and singularity available beforehand.

1. Install

```
git clone https://github.com/ncc-ccat-gap/GCATWorkflow.git
cd GCATWorkflow
python setup.py install
```

2. Pull container images

For example,
```
singularity pull docker://genomon/bwa_alignment:0.2.0
```

3. Edit config file

Edit `image` options, to pulled `.sif`.  
And edit paths of reference files.
```
vi ./tests/gcat.cfg
```

4. Edit sample.csv

Edit paths of sequence files.
```
vi ./tests/sample.csv
```

## How to use

1. Configure

case with germline mode
```
gcat_workflow germline ./tests/sample.csv ${output_dir} ./tests/gcat.cfg
```


2. `snakemake`
```
cd ${output_dir}
snakemake --cores 4 -k
```

case, dry-run
```
snakemake -n
```

case, re-run (force all)
```
snakemake --forceall
```


## License

Starting from v3.3.0, this software is free for academic and non-commercial research use only.
Commercial users must obtain a commercial license. Please contact the author.

See LICENSE for details.