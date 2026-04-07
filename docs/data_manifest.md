# Data Manifest

## Overview
This document describes all datasets used in the TCGA LUAD genomics pipeline.

---

## Expression Data
- Source: UCSC Xena (TCGA LUAD)
- Type: RNA-seq gene expression (STAR counts)
- Samples: ~574
- Genes: ~20,000

---

## Metadata
- Source: TCGA clinical annotations
- Contains:
  - Sample ID
  - Sample type (Tumor / Normal)

---

## Mutation Data
- Source: LinkedOmics
- Type: Gene-level mutation matrix
- Samples: ~533
- Values: Mutation presence/absence

---

## Processed Data
Located in `data/processed/`:
- `expression_filtered.csv`
- `metadata_matched.csv`

---

## Notes
- Sample IDs were standardised to TCGA format (15 characters)
- Only overlapping samples across datasets were retained