# Scripts and pipelines for the het-v paper
Scripts and Snakemake pipelines associated with the paper of Ament-Velásquez et al. "Allorecognition genes drive reproductive isolation in *Podospora anserina*", in preparation.

There are three main [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines. The first one, `SNPpop.smk`, must be run first. Some of the outputs are used by the other two pipelines. So, the order goes:

    1 - SNPpop.smk
    2 - DiversityStats.smk
    3 - PaLD.smk

Each pipeline in turn depends on a small conda environment(s), as explained in their own README files.

As output of all the pipelines there is a `results` folder that usually contains the figures used in the paper. 

I ran the pipelines in Uppsala University's supercomputer [UPPMAX](https://uppmax.uu.se/), which has a CentOS Linux operating system with a slurm scheduler. However, they should work fine also in other unix environments.

---

In addition to the bioinformatic pipelines, the folder `SLiMsimulations` contains the code and simulations using [SLiM 3](https://academic.oup.com/mbe/article/36/3/632/5229931) written by Ivain Martinossi.
  