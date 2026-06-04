<div align="center">
  <img src="assets/simplseq-readme-banner.png" alt="SIMPLseq-nf malaria amplicon pipeline" width="620">
</div>

# **SIMPLseq-nf:** *A Malaria Amplicon Pipeline with Integrated Genomic Transmission Analysis*

SIMPLseq-nf App is a local browser app for SIMPLseq malaria amplicon sequencing runs. It scans paired FASTQ files, runs the SIMPLseq Nextflow workflow, and shows SIMPLseq, DINEMITES, and Dcifer results in your browser.

Your data stays on your computer.

## Install

Install on Linux, WSL, or macOS:

```bash
curl -fsSL https://github.com/a-nadeem9/simplseq-malaria-amplicon-pipeline-nf/releases/download/v2.2.1/install-simplseq.sh | bash
simplseq run
```

## Use The App

1. Open the app with `simplseq run`.
2. In **Inputs**, choose the folder with your paired FASTQ files.
3. Click **Scan folder** to create `samples.csv`.
4. In **Run**, click **Run runtime check**.
5. Click **Start run**.
6. Use **Results** to view reports and download output files.
7. Use **DINEMITES** for new-infection analysis.
8. Use **Dcifer** for pairwise relatedness analysis.

## FASTQ Names

Supported paired-read names include:

```text
*_R1.fastq.gz / *_R2.fastq.gz
*_R1_001.fastq.gz / *_R2_001.fastq.gz
*_R1.fq.gz / *_R2.fq.gz
*_R1_001.fq.gz / *_R2_001.fq.gz
```
The app writes a sample sheet named `samples.csv`. You can edit it before starting the run if sample names or dates need correction.

## Main Outputs

Each run writes a new output folder. Common files include:

| File | Description |
| --- | --- |
| `reports/run_summary.html` | Main run report |
| `run_dada2/seqtab_iseq.tsv` | ASV count table |
| `run_dada2/ASV_mapped_table.tsv` | ASVs mapped to SIMPLseq targets |
| `run_dada2/asv_to_cigar.tsv` | ASV to CIGAR haplotype map |
| `run_dada2/seqtab_cigar.tsv` | Final CIGAR count table |

## DINEMITES Outputs

DINEMITES results are written to:

```text
<results>/dinemites/
```

Common files include:

| File | Description |
| --- | --- |
| `dinemites_allele_probabilities.tsv` | Per-allele probabilities |
| `dinemites_allele_key.tsv` | Short allele IDs mapped to exact alleles |
| `dinemites_new_infections.tsv` | New-infection summaries |
| `dinemites_molfoi.tsv` | molFOI summaries |
| `dinemites_plots/` | Subject plots |

## Dcifer Outputs

Dcifer results are written to:

```text
<results>/dcifer/
```

Common files include:

| File | Description |
| --- | --- |
| `dcifer_coi.tsv` | Complexity-of-infection estimates |
| `dcifer_pairwise_relatedness.tsv` | Pairwise relatedness estimates |
| `dcifer_relatedness_matrix.tsv` | Relatedness matrix |
| `dcifer_pvalue_matrix.tsv` | Raw p-value matrix |
| `dcifer_plots/` | Relatedness heatmaps |

<div style="color: grey;">

## References

- Schwabl P, Amaya-Romero J-E, Neafsey DE, et al. SIMPLseq: a high-sensitivity *Plasmodium falciparum* genotyping and PCR contamination tracking tool. *Malaria Journal*. 2026. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC12958562/" style="color: grey;">https://pmc.ncbi.nlm.nih.gov/articles/PMC12958562/</a>

- Broad Institute. malaria-amplicon-pipeline. GitHub. <a href="https://github.com/broadinstitute/malaria-amplicon-pipeline" style="color: grey;">https://github.com/broadinstitute/malaria-amplicon-pipeline</a>

- Nickols WA, Schwabl P, Niangaly A, Murphy SC, Crompton PD, Neafsey DE. Distinguishing new from persistent infections at the strain level using longitudinal genotyping data. 2025. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11839113/" style="color: grey;">https://pmc.ncbi.nlm.nih.gov/articles/PMC11839113/</a>

- Gerlovina I, Gerlovin B, Rodríguez-Barraquer I, Greenhouse B. Dcifer: an IBD-based method to calculate genetic distance between polyclonal infections. *Genetics*. 2022;222(2):iyac126. <a href="https://academic.oup.com/genetics/article/222/2/iyac126/6674513" style="color: grey;">https://academic.oup.com/genetics/article/222/2/iyac126/6674513</a>

</div>
