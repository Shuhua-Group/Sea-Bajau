# Revised APC/CPC2/HPRC2 comparison

The primary comparison is **catalog-level**. It retains `GRCh38`, T2T-YAO, and every
other graph donor because a VCF record denotes a graph resource's encoded allele, not
a population-frequency estimate. It never merges sample genotype columns and never
uses AN/AF filtering for Venn or catalog-membership statistics.

`run_catalog_preparation.sh` performs the required common preprocessing: CHM13
normalization, multiallelic splitting, autosome restriction, and difficult-region
exclusion. It then makes an exact normalized sequence-allele membership table using
disk-backed sorting and streaming, rather than loading the catalogs into memory.

Main results:
1. `04.exact_membership.detailed.tsv` and `.counts.tsv`;
2. `plot_catalog_overlap.R`: a three-set Venn and membership composition by class;
3. `catalog_class_enrichment_bootstrap.R`: 1-Mb block-bootstrap CIs for each
   catalog's strictly specific alleles.

Exact matching should be the primary definition for sequence-resolved alleles. State
this explicitly. A tolerant breakpoint/size SV clustering analysis should be added as
a sensitivity analysis, not silently mixed with the exact small-variant definition.

Population/AF analyses are a separate Supplementary workflow. They must use curated
population sample manifests, exclude graph reference-only samples from their
denominators, derive AN thresholds from those manifests, and remain descriptive.
