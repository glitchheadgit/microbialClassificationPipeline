# PIPE_001

## Info

Pipeline consists of 3 chief scripts:

1. bbmap_index.nf - Makes indexes for your references and gives them appropriate identificator
2. main.nf - includes following submodules:
    * preprocessing
    * mapping
    * classify_reads
    * assembly
    * binning
3. additional.nf:
    * binrefinement
    * bins classification

All processes and their parameters can be found in subscripts/process.nf \
All information about submodules - in according nf scripts. Each submodule can be launched on its own as separate script!

## Requirements

### Conda

Install [Conda](https://conda.io/projects/conda/en/latest/user-guide/install)/[Miniconda](https://docs.anaconda.com/miniconda/)

### Nextflow

Create new environment with installed Conda/Miniconda:
```bash
conda install --name nextflow --yes -c bioconda nextflow==22.05
```

### MultiQC

Create new environment for MultiQC
```bash
conda install -p /path/to/pipeline/multiqc --yes -c bioconda multiqc==1.24.1
```

### GTDB-Tk

1. GTDB-Tk requires ~84G of external data that needs to be downloaded and unarchived:

```bash
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz # mirror: https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```
2. Link GTDB to the "gtdbtk_data" folder in pipeline:

```bash
ln -s /path/to/gtdbtk /path/to/pipeline/gtdbtk_data
```

### CheckM

1. CheckM DB requires ~1.4G of external data that needs to be downloaded and unarchived:

```bash
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xvzf checkm_data_2015_01_16.tar.gz
```
2. Link CheckM database to the "checkm_data" folder in pipeline:

```bash
ln -s /path/to/checkm_data /path/to/pipeline/checkm_data
```

### Kraken 2 or MetaPhlAn

1. Kraken requires ~16Gb of external data that needs to be downloaded and unarchived:

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20231009.tar.gz
tar xvzf k2_standard_16gb_20231009.tar.gz
```
2. Link Kraken2 db to the "kraken2_data" folder in pipeline:

```bash
ln -s /path/to/kraken2_data /path/to/pipeline/kraken2_data
```

---

1. Install MetaPhlAn in the pipeline directory:

```bash
mamba create -p /path/to/pipeline/metaphlan -c bioconda metaphlan
```
2. It requires ~15Gb of external data that needs to be downloaded:

```bash
mamba run -p /path/to/pipeline/metaphlan metaphlan --install 
```


### Installing metaWRAP

1. Install it in the pipeline directory:

```bash
mamba create -p /path/to/pipeline/mw-env -c ursky metawrap-mg=1.3.2
```

2. Make corrections to metaWRAP bin_refinement.sh script:

```bash
sed -i 's/(( $SIZE > 50000)) &&//g' /path/to/pipeline/mw-env/bin/metawrap-modules/bin_refinement.sh 
```

## How to start pipeline

### Make reference index

```bash
nextflow run bbmap_index.nf \
-with-report \
-with-conda \
--ref host.fna,human.fna
```

### Launch main.nf

```bash
nextflow run main.nf \
-with-conda \
-with-report \
--ref host.fna,human.fna \
--ref_number 1,2 \
--reads 'pipeline/fq_check/*{1,2}.fq' \
--maxbin2 --concoct --metabat2 \ # specify at least one binner
--output pipeline_results
```

Optional parameters of main.nf with their default options:

1. fastp options
    * --compress_level 7
    * --poly_g_min_len 5
    * --poly_x_min_len 10
    * --cut_mean_quality 27
    * --n_base_limit 0
    * --average_qual 25
    * --length_required 80
    * --working_thread_n 12
2. kraken2, bracken options
    * --classification_lvl "S"
3. bbmap options
    * --host_dir "host"
    * --human_dir "human"
    * --bbmap_memory false
    * --bbmap_cpus 4
4. metawrap binning options
    * --concoct false
    * --metabat2 false
    * --maxbin2 false
5. metawrap bin refinement options
    * metawrap_completion 50
    * metawrap_contamination 30
