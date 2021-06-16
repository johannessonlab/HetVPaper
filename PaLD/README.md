# PaLD: A pipeline to estimate linkage disequilibrium in *Podospora anserina*

This pipeline is meant to be run using the vcf file with no missing data produced by `SNPpop.smk`. The output are:

- LD decay per chromosome of the entire Wageningen Collection (Figure S3 in the paper), plus some extras 
- LD heatmap-like plot inspired by [Hench et al. (2019)](https://www.nature.com/articles/s41559-019-0814-5) (Figure S10)

In the code, I ofter refer to the "V" and "A" mating groups. These are equivalent to the V and V1 alleles of *het-v*.

I ran the pipeline in a CentOS Linux environment in the slurm cluster [Uppmax](https://uppmax.uu.se/). It should work in a Mac OS too.

## Building the environment

First, to update conda:

    $ conda update -n base conda

The pipeline relies on the `PaLD` environment, which looks like this:

    $ conda create -n PaLD -c bioconda snakemake-minimal=5.4.4 vcftools=0.1.16 bcftools=1.9

To compress vcf files and make indexes:
    
    $ conda activate PaLD    
    $ conda install -c bioconda htslib=1.9

I also have it in a yaml file in case I need to build it again:

    $ cat envs/PaLD.yaml
```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
  - r
dependencies:
  - snakemake-minimal=5.4.4
  - vcftools=0.1.16
  - htslib=1.9
```

The pipeline also depends on a small yaml for the R environment. This one has to be present for the pipeline to run (using `--use-conda`).

    $ cat envs/PaLDPlot.yaml
```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
  - r
dependencies:
  - r-dplyr=0.8.3
  - r-cowplot=1.0.0 # It includes ggplot2
  - r-stringr=1.4.0
  - r-hexbin=1.27.3
  - r-gridextra=2.3 # It includes grid
  - r-vcfr=1.10.0
  - r-vegan=2.5_6 # Otherwise vcfR doesn't work
```

## Configuration file

This pipeline depends on a given configuration file including the samples, the path to the data, and the reference. The file looks like such:

```yaml
# Samples NOT in Wageningen (to be removed)
notWa: ["PaSp", "PaZp", "PaYp", "CBS433.50p", "CBS455.64m", "PaTgp"]

# Samples divided in mating groups
groupV: ["PaWa1p", "PaWa2m", "PaWa3m", "PaWa4p", "PaWa10p", "PaWa11m", "PaWa12p", "PaWa13m", "PaWa14p", "PaWa16p", "PaWa19m", "PaWa22m", "PaWa23p", "PaWa24m", "PaWa29p", "PaWa32p", "PaWa33m", "PaWa36p", "PaWa37m", "PaWa39m", "PaWa40m", "PaWa41p", "PaWa42m", "PaWa43p", "PaWa44m", "PaWa52p", "PaWa53m", "PaWa54m", "PaWa55p", "PaWa56m", "PaWa58m", "PaWa59m", "PaWa62p", "PaWa66m", "PaWa67p", "PaWa68m", "PaWa69p", "PaWa70m", "PaWa81p", "PaWa83m", "PaWa100p", "PaWa85p", "PaWa87p", "PaWa92p", "PaWa94p", "PaWa95p", "PaWa98m", "PaWa99p", "PaWa102p", "PaWa106p", "PaWa117m", "PaWa122m", "PaWa123p", "PaWa125m", "PaWa127m", "PaWa129p", "PaWa138m"]
groupA: ["PaWa7m", "PaWa8p", "PaWa9m", "PaWa15m", "PaWa17m", "PaWa18p", "PaWa21m", "PaWa25p", "PaWa26m", "PaWa27m", "PaWa28m", "PaWa38p", "PaWa45p", "PaWa46p", "PaWa47m", "PaWa49m", "PaWa57p", "PaWa60p", "PaWa61m", "PaWa63p", "PaWa64m", "PaWa71p", "PaWa72m", "PaWa76p", "PaWa77m", "PaWa78p", "PaWa79m", "PaWa86m", "PaWa88p", "PaWa89p", "PaWa91p", "PaWa96m", "PaWa97p", "PaWa101m", "PaWa103m", "PaWa104m", "PaWa105p", "PaWa107m", "PaWa108m", "PaWa109p", "PaWa115m", "PaWa116p", "PaWa118p", "PaWa124p", "PaWa126p", "PaWa128p", "PaWa137m", "PaWa142p", "PaWa143m"]

# The variants file (from the `SNPpop.smk` pipeline)
vcf: "../SNPpop/results/PodoPop-snps-NoTEs-gatkPASS-NoSibs-miss1.vcf.gz"

# The metadata with the het-v genotypes
PopData: "../SNPpop/data/SNPpopMetadata.csv"

# The R environment
envR: "envs/PaLDPlot.yaml"

# Number of chromosomes (assuming vcf has contigs names like "chromosome_1")
NCHR: 7

# Number of randomly sampled windows from each chromosome
NWIN: 30

# Length of LD windows for decay
WINLEN: 50000

# Thining of the SNPs for the intra and interchromosomal LD calculation
THIN: 1000

# Thining for het-v and het-r areas
THINhet: 1000

# MAF for LD heatmaps
MAF: 0.02
MAFhet: 0.02

# Define the het-v and het-r areas
hetvSTART: 725000
hetvEND: 1725000

hetrSTART: 2250000
hetrEND: 3250000

# Manipulate vcf files
vcfconcat4LD: "scripts/vcfconcat4LD.py"

# Plotting scripts
PaLDPlot: "scripts/Rscripts/PaLDPlot-nls-grs.R"
PaLDheatmap: "scripts/Rscripts/PaLDheatmap.R"
PaLD_SNPvsHetV: "scripts/Rscripts/PaLD_SNPvsHetV.R"

```
## Run pipeline in Uppmax

First, to get an idea of how the pipeline looks like we can make a rulegraph:
    
    $ conda activate PaLD
    $ snakemake --snakefile PaLD.smk --configfile PaLD_config.yaml --rulegraph | dot -Tpng > rulegraph.png

![rulegraph](rulegraph.png "rulegraph of PaLD.smk")

For testing:

    $ snakemake --snakefile PaLD.smk --configfile PaLD_config.yaml -pn

Run the pipeline:

    $ screen -R PaLD
    $ conda activate PaLD
    $ snakemake --snakefile PaLD.smk --configfile PaLD_config.yaml -p --cluster "sbatch -A snicXXXX-X-XXX -p core -n {params.threads} -t {params.time} --mail-user xxxxxxxxxxxx@xxxxx.xx --mail-type=ALL" -j 20 --keep-going --use-conda &> PaLD.log &

## Run pipeline locally

The same, but do:

    $ snakemake --snakefile PaLD.smk --configfile PaLD_config.yaml -p -j 20 --keep-going --use-conda &> PaLD.log &
