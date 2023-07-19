# geteilt
NARST workflow built in nextflow. Pulls paired short-read data from NCBI using SRRs or finds local reads, and runs through read metrics, coverage cutoff, assembly, and AR, plasmid, and point-mutation screening.

Please note that both the ARANDPLASMIDSCREEN and MUTATIONAL processes use custom databases only available to CDC users. If you are an external user, you should configure the workflow to run AMRFinderPlus instead of staramr. See [Parameters](#parameters) section for detailed instructions.

## Updates
7/19/2023:  
* Added print statement to encourage users to delete the work directory once finished
* Added bash script to create import sheets under scripts ⚠️Only for CDC users generating import sheets for NARMS database⚠️
* Option added for running workflow on local reads
* Option added for running NCBI's AMRFinderPlus instead of staramr
* ARANDPLASMIDSCREEN
	* Staramr updates databases to specific resfinder, pointfinder, and plasmidfinder commits
* COVERAGECUTOFF
	* Parameter added to account for variable read file naming formats

## Install
Navigate to your home directory and git clone the repository.
```bash
$ git clone https://github.com/ncezid-narst/geteilt.git
```
You will need to install Nextflow if you don't already have it: https://www.nextflow.io/docs/latest/getstarted.html

You will also need to install Singularity if you don't already have it: https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html

If you're working on CDC servers, run `module load nextflow/XX.XX.X` to load Nextflow, Singularity, and Java modules.

As of 6/12/2023, this workflow was developed on Nextflow version 22.10.6.

You will want an NCBI API Key: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

## Overview
READMETRICS: Calculates read metrics using CG-Pipeline  
COVERAGECUTOFF: Calculates minimum coverage cut-off based on read metrics  
ASSEMBLE: Assembles reads using Shovill  
ARANDPLASMIDSCREEN: Screens for resistance determinants, plasmid replicons, and point mutations using resfinder, plasmidfinder, and pointfinder databases respectively using Staramr  
AMRFINDER: Screens for resistance determinants, plasmid replicons, and point mutations using NCBI's AMRFinderPlus  
MUTATIONAL: Screens for additional point mutations for some genus/species combinations using Ariba  

## Parameters
Parameters for each process can be changed in `geteilt.config` under the first bracketed section `params`. Check out Check out [Resources](#resources) for links to each process's main github page to learn more about process-specific parameters.

Prior to running geteilt, make sure the INITIAL PARAMETERS are set accurately:
```java
	import_sheet = 'import.txt'
	ncbi_api_key = ''
	outdir = 'geteilt'
	local_reads = false
	reads = '' //reads/*_{1,2}.fastq.gz //reads/*_L001_R{1,2}_001.fastq.gz
	amrfinder = false
```
`import_sheet`: Tab-delimited text file wtih sample information. When pulling reads from NCBI, format should be as follows:
```bash
SRR     WGS     Genus   Species
SRR1910447      2013AY-0074     Campylobacter   coli
SRR8497481      2013AY-1307     Campylobacter   jejuni
SRR6219901      CFSAN070926     Campylobacter   upsaliensis
SRR24829758     PNUSAE138155    Escherichia     coli
SRR23509899     2022V-1177      Vibrio  cholerae
SRR7366937      PNUSAV000187    Vibrio  parahaemolyticus
SRR14250202     PNUSAY000051    Yersinia        enterocolitica
```
* SRR: SRA run number as listed on NCBI
* WGS: Sample ID (will be final name of sample)
* Genus: Sample genus
* Species: Sample species

When running the workflow using local reads, format should be as follows:
```bash
Read    WGS     Genus   Species
2013K-1828-M3235-22-039_S7      2013K-1828      Salmonella      enterica
2014C-3598-M347-22-044_S2       2014C-3598      Escherichia     coli
2014C-3599-M347-22-044_S3       2014C-3599      Escherichia     coli
2014C-3857-M347-23-003_S14      2014C-3857      Escherichia     coli
2014C-3946-M3235-17-036_S9      2014C-3946      Escherichia     coli
```
* Read: Name of read file. This name must match the filename of the read **minus** the string denoted under the `reads` parameter. For example, if the full name of the read file is `2013K-1828-M3235-22-039_S7_L001_R1_001.fastq.gz` and `reads = 'reads/*_L001_R{1,2}_001.fastq.gz'`, then the value under Read must be `2013K-1828-M3235-22-039_S7`.
* WGS: Sample ID (will be final name of sample)
* Genus: Sample genus
* Species: Sample species
 
`ncbi_api_key`: NCBI API Key to increase level of bandwidth/access when pulling reads. Necessary unless you'd like to edit the script directly in `geteilt.nf` under `Workflow` to not use an API key.  
`outdir`: Name of geteilt output directory. Default name is set to `geteilt`.  
`local_reads`: When `true`, the workflow will look for local reads using the `reads` parameter value. Default value is set to `false`.  
`reads`: Path to local directory with paired-end reads. May have to create custom string pattern depending on format of reads. Two possible formats are provided.  
`amrfinder`: When `true`, the workflow will run NCBI's AMRFinderPlus instead of staramr.  

## Processes
Directives for each process can be changed in `geteilt.config` under the second bracketed section `process`. This is where you can update the containers used in the processes that use one. Check out [Resources](#resources) to see a full list of all the containers and the tools' githubs.

## Profiles
Configuration settings for each profile can be changed in `geteilt.config` under the third bracketed section `profiles`. This is where you can update or create profiles that will dictate where and how each process is run. By default, there are two main profiles and three auxiliary profiles:

* `standard`: Will execute geteilt using the 'local' executor, running processes on the computer where Nextflow is launched. 
* `sge`: Will execute geteilt using a Sun Grid Engine cluster, running processes on the HPC (qsub).
* `short`: Auxiliary profile to change the sge default queue to short queue
* `gpu`: Auxiliary profile to chage the sge default queue to gpu queue
* `highmem`: Auxiliary profile to change the sge default queue to highmem queue

You can see how profiles are used in the next section **Usage**.

NOTE: The default profile settings were mostly pulled from recommendations made by CDC Scicomp in their Nextflow training called 'Reproducible, scalable, and shareable analysis workflows with Nextflow'. There is a good chance you will have to create/modify your own profile to run geteilt using your institution's computing environment. Check out [Resources](#resources) to learn more about creating profiles.

## Usage
Once you've made the necessary changes to the configuration file to run the workflow on your computing environment and have set up inital parameters, you can run stylo just as you would any nextflow workflow:
```bash
nextflow run /path/to/geteilt/schtappe/geteilt.nf -c /path/to/geteilt/config/geteilt.config
```
Nextflow is picky about single-hyphen flags vs. double-hyphen flags. Single-hyphens affect the nextflow command while double-hyphens affect the parameters in the configuration file. For example, to change the initial parameters without directly editing `geteilt.config`:
```bash
nextflow run /path/to/nanoporeWorkflow/schtappe/stylo.nf -c /path/to/nanoporeWorkflow/config/stylo.config \
  --import_sheet path/to/your/importsheet.txt \
  --ncbi_api_key ########## \
  --outdir youroutputdirectory \
```

By default, nextflow will run locally. If you want to specify a profile, use the `-profile` flag. For example, to qsub geteilt's processes:
```bash
nextflow run /.../geteilt/schtappe/geteilt.nf -c /.../geteilt/config/geteilt.config -profile sge
```

You can change the queue by adding the auxiliary profile name, separated by a comma:
```bash
nextflow run /.../geteilt/schtappe/geteilt.nf -c /.../geteilt/config/geteilt.config -profile sge,highmem
```
Run `nextflow help` or `nextflow run -help` for more information on nextflow flags.

## Output
Here's what geteilt output looks like per sample (directories only):
```bash
geteilt/
├── 2013AY-0074
│   ├── mutational
│   │   ├── mutational_23S
│   │   │   └── ariba.tmp.td0d420g
│   │   │       └── cluster
│   │   └── mutational_gyrA
│   │       └── ariba.tmp.10b2vw84
│   │           └── cluster
│   ├── readmetrics
│   ├── shovill
│   └── staramr
│       └── hits
```
Note: The `staramr` directory will be replaced by `amrfinder` when running NCBI's AMRFinderPlus over staramr.

## Resources
### Containers:
* LYVE-SET (CG-Pipeline):
  * https://hub.docker.com/r/staphb/lyveset
  * https://github.com/lskatz/lyve-SET
* Shovill:
  * https://hub.docker.com/r/staphb/shovill
  * https://github.com/tseemann/shovill
* Staramr:
  * https://hub.docker.com/r/staphb/staramr
  * https://github.com/phac-nml/staramr
* AMRFinderPlus:
  *  https://hub.docker.com/r/ncbi/amr
  *  https://github.com/ncbi/amr
* Ariba:
  * https://hub.docker.com/r/staphb/ariba
  * https://github.com/sanger-pathogens/ariba
 
