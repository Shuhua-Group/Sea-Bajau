# Data and confidentiality policy

This repository is code-only.

Do not commit:

- the manuscript or unpublished reviewer material;
- participant identifiers, sample keys or consent-related records;
- raw PacBio HiFi or short-read data;
- assemblies, BAM/CRAM files or VCF/BCF files;
- controlled-access reference-panel data;
- unpublished result tables or figures;

Before every public release, inspect the staged files with `git diff --cached` and scan for identifiers, credentials and absolute paths.
