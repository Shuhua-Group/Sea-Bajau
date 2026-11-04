# Manuscript analysis to code map

| Manuscript result or method                            | Repository location                                   | Scope                                                                    |
| ------------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------------------ |
| Long-read sequencing and haplotype-resolved assemblies | `workflows/01_assembly_qc/hifi_processing`          | HiFi read processing, assembly and primary QC                            |
| Assembly quality and assembly-related bias             | `workflows/01_assembly_qc/flagger`, `revision_qc` | NucFlag, Flagger and assembly/SV QC summaries                            |
| Genome-wide SNV, Indel and SV landscape                | `workflows/02_variant_discovery`                    | Variant processing, normalization, comparison and difficult-region masks |
| PCA, regional affinities and ADMIXTURE                 | `workflows/03_population_structure`                 | APC/reference merging, PCA projections and ADMIXTURE CV summaries        |
| SFS, ROH, IBD, LD decay and recent Ne                  | `workflows/04_demographic_isolation`                | Isolation and bottleneck analyses, Relate, hapIBD and hapNe/IBDNe        |
| Highly differentiated variants and RHEB                | `workflows/05_selection_rheb`                       | FST, composite selection signals and regional plots                      |
| MUC5AC structural haplotypes                           | `workflows/06_muc5ac`                               | Region extraction, sequence alignment and haplotype plotting             |
| PRSS, PGA and AMY copy-number architecture             | `workflows/07_digestive_cnv`                        | Pangenome graphs, genotyping, population summaries and PGS models        |
| Deep ancestry and sex-biased history                   | `workflows/08_uniparental_demography`               | Deep demographic analysis plus mitochondrial and Y lineages              |
