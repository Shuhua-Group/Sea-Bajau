# Sea Bajau genomics

This repository contains analysis scripts for *Haplotype-resolved genomes illuminate demographic isolation and complex variation in the Sea Bajau*.

Version 1.0 is the first release of the code archive. We are still cleaning up the scripts and improving the documentation. Future updates will simplify configuration, clarify the order of analysis steps, and add examples where possible.

The scripts are grouped by analysis below. Some retain settings from our computing environment, so input paths, software locations, and cluster settings need to be adjusted before use. Example paths use `/path/to` as a placeholder.

## Repository layout

| Directory                               | Manuscript analysis                                                                   |
| --------------------------------------- | ------------------------------------------------------------------------------------- |
| `workflows/01_assembly_qc`            | HiFi processing and assembly quality summaries, including NucFlag and Flagger outputs |
| `workflows/02_variant_discovery`      | Small-variant and structural-variant processing and catalog comparisons               |
| `workflows/03_population_structure`   | Dataset merging, PCA and ADMIXTURE                                                    |
| `workflows/04_demographic_isolation`  | ROH, IBD, LD decay, Relate and recent effective population size                       |
| `workflows/05_selection`              | FST, CMS and local selection analyses                                                 |
| `workflows/06_muc5ac`                 | MUC5AC extraction, alignment and visualization                                        |
| `workflows/07_digestive_cnv`          | PRSS, PGA and AMY graph genotyping, copy number summaries and PGS models              |
| `workflows/08_uniparental_demography` | Deep demographic, mtDNA and Y-chromosome analyses                                     |
| `docs`                                | Manuscript-to-code map                                                                |

## Getting started

The [analysis map](docs/MANUSCRIPT_CODE_MAP.md) links the manuscript analyses to the corresponding directories. 

Software requirements vary by analysis. Refer to the manuscript methods and the commands or module declarations in each script for the tools used.

## Data policy

Raw reads, assemblies, VCFs, sample metadata, and intermediate results are not included. See the [data policy](docs/DATA_POLICY.md) for details.

## Citation and contact

Please cite the associated paper when using these scripts. We will add the full citation when it becomes available.

For questions or bug reports, open a GitHub issue or contact Ziyi Yang at youngxyaaa@gmail.com.
