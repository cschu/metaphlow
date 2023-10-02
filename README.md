# Metaphlow

Metaphlow is a nextflow workflow for running Metaphlan4. The workflow includes optional read preprocessing and host/human decontamination steps. It furthermore supports downstream runs of HUMAnN3 and samestr as well as optional Metaphlan3 runs. 

## Prerequisites & Requirements

The easiest way to handle metaphlow's dependencies is via Singularity containers. Alternatively, conda environments, software module systems or native installations can be used.

### Preprocessing

Preprocessing and QA is done with bbmap, fastqc, and multiqc.

### Decontamination/Host removal

Decontamination is done with kraken2 and additionally requires seqtk. 

#### Kraken2 database

Host removal requires a kraken2 host database.

### Metaphlan Profiling

The default supported Metaphlan version is 4. Additionally, Metaphlan3 can be run. The downstream analyses are based purely on Metaphlan4, though.

#### CHOCOPhlAn database for Metaphlan4

Get the `mpa_vOct22_CHOCOPhlAnSGB_202212` database from [here](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar), unpack the tarball, and point the `--mp4_db` parameter to the database's root directory. 

In `params.yml`:

```
mp4_db: "/path/to/mpa_vOct22_CHOCOPhlAnSGB_202212/"
```

On the command line:

```
--mp4_db "/path/to/mpa_vOct22_CHOCOPhlAnSGB_202212/"
```

#### CHOCOPhlAn database for Metaphlan3 (optional)

Get the `mpa_v30_CHOCOPhlAn_201901` database from [here](), unpack the tarball. Then point the `--mp3_db` parameter to the parent directory of the database's root directory and the `--mp3_db_version` parameter to the basename of the database's root directory.

In `params.yml`:

```
mp3_db: "/path/to/mp3_db/"
mp3_db_version: "mpa_v30_CHOCOPhlAn_201901"
```

On the command line:

```
--mp3_db "/path/to/mp3_db/"
--mp3_db_version "mpa_v30_CHOCOPhlAn_201901"
```

### Downstream Analyses (Optional)

Metaphlow supports HUMAnN3 and samestr downstream analyses.

#### HUMAnN3 databases

  Get the annotated CHOCOPhlAn db from [here](http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.tar.gz) and the annotated uniref db from [here](http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz), unpack the tarballs and set the respective parameters.

  In `params.yml`:

```
humann_nuc_db: "/path/to/full_chocophlan_db/"
humann_prot_db: "/path/to/uniref90_annotated_v201901b_full/"
```

On the command line:

```
--humann_nuc_db "/path/to/full_chocophlan_db/"
--humann_prot_db "/path/to/uniref90_annotated_v201901b_full/"
```

#### SameStr databases

TBD


## Running Metaphlow

A metaphlow run is controlled by environment-specific parameters (s. `run.config`) and studiy-specific parameters (s. `params.yml`). The parameters in the `params.yml` can be specified on the command line as well.

You can either clone metaphlow from GitHub and run it as follows

```
git clone https://github.com/cschu/metaphlow.git
nextflow run /path/to/metaphlow [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

Or, you can have nextflow pull it from github and run it from the `$HOME/.nextflow` directory.

```
nextflow run cschu/metaphlow [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

### Input files

Metaphlow supports fastq files. These can be uncompressed (but shouldn't) or compressed with gzip or bzip2. Sample data must be either be arranged in one directory per sample or in one global directory containing all files.

#### Per-sample input directories

All files in a sample directory will be associated with name of the sample folder. Paired-end mate files need to have matching prefixes. Mates 1 and 2 can be specified with suffixes `_[12]`, `_R[12]`, `.[12]`, `.R[12]`. Lane IDs or other read id modifiers have to precede the mate identifier. Files with names not containing either of those patterns will be assigned to be single-ended. Metaphlow assumes samples that consist of both single and paired end files to be paired end with all single end files being orphans (quality control survivors). 

#### Global study input directory

Samples have to be pure single ends or pure paired end. Mixed samples or multiple lanes are not supported.










