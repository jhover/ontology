#	Ontology Utilities
This module provides an easy interface to Gene Ontology manipulations, in combination with species-specific gene annotations. 

## Install software and dependencies

These instructions assume familiarity with running bioinformatics pipelines...

* Install Conda. 
[https://docs.conda.io/projects/miniconda/en/latest/index.html](https://docs.conda.io/projects/miniconda/en/latest/index.html)

* Create an environment for the module pipeline.

```
conda create -n ontology python==3.9 
```

* Activate the environment

```
conda activate ontology
```

* Install additional Conda repositories

```
conda config --add channels conda-forge
conda config --add channels bioconda
```
* Install dependencies, useful tools, confirm update

```
# For Linux and MacOS
	conda install -y pandas numpy scipy   
```

* Clone the onotology project from the repository to the standard location. (This assumes you already have git installed. If not, install it first). 

```
mkdir ~/git
git clone https://github.com/jhover/ontology.git
```
All the code is currently organized so it is run directly from the git directory (rather than installed). 

* Create a working directory for your experiment, and copy in a metadata file and the fastq sequencing data files, and the default configuration file. 

```
```

## Initial Configuration
By default the module and scripts will take their defaults from a single configuration file, included in the distribution ~/git/ontology/etc/ontology.conf.
Often, these can be overridden from the command line. You can also copy the default configuration file, edit it, and use it from each command with the -c <configfile> switch. 

For each command, there are typical arguments. Additional arguments may be supported. Run the -h switch to see the full usage help. 

<IN PROGRESS>
