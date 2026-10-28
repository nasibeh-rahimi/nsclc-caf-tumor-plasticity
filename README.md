# nsclc-caf-tumor-plasticity
Integrative transcriptomic analysis of CAF-associated gene programs and tumor-state signatures in non-small cell lung cancer (NSCLC), using public GEO datasets and patient-level statistical analysis.
## Project overview

This computational research project investigates whether a cancer-associated fibroblast (CAF)-associated transcriptional program is associated with tumor-state signatures in non-small cell lung cancer (NSCLC).

The analysis combines public gene-expression datasets, differential expression analysis, gene-signature scoring, patient-level aggregation, correlation analysis, histology-adjusted regression, and sensitivity analyses.

**Study type:** Observational, hypothesis-generating transcriptomic analysis.

**Key finding:** In a cohort of 20 independent patients, the CAF-associated transcriptional score was not significantly associated with the primary EMT-associated score. A secondary stemness-associated analysis also showed no statistically significant association.

The results do not establish that CAFs cause tumor plasticity, nor do they rule out a biological relationship. The study is limited by its small cohort, heterogeneous histology, cross-sectional design, and use of bulk transcriptomic data.

## Research question

Is a CAF-associated transcriptional program associated with tumor plasticity-related states in NSCLC?

The project was framed as an association study. The bulk transcriptomic design cannot establish that CAFs cause tumor plasticity, temporal state transitions, or cell-to-cell interactions.

## Datasets

### GSE149955: CAF program discovery

- Platform: Affymetrix GPL570.
- Samples: 25 matched normal fibroblast (NF) and 25 CAF samples.
- Analysis: Paired differential expression analysis.
- Result: 360 genes met the false-discovery-rate threshold of FDR < 0.05, including 145 upregulated and 215 downregulated genes.

A signed CAF-associated transcriptional score was constructed from this fixed 360-gene program. This score represents a transcriptional program; it should not be interpreted as a direct measurement of CAF abundance.

### GSE33532: Tumor-state analysis

- Samples: 20 patients, each with four tumor subsamples (A/B/C/D) and one normal sample.
- Tumor subsamples: 80 in total, representing 20 patients rather than 80 independent patients.
- Patient-level histology: 10 adenocarcinoma, 6 mixed, and 4 squamous.
- Primary inference: tumor subsamples were aggregated at the patient level, giving **n = 20 independent observations**.

The CAF-associated gene program was transferred to GSE33532. All 360 genes were represented in the target expression matrix. This is gene-set coverage and transfer, not independent validation of the CAF–tumor association.

## Analysis strategy

1. Identify differentially expressed genes between matched CAF and normal fibroblast samples in GSE149955.
2. Construct a signed CAF-associated transcriptional score from the fixed gene program.
3. Transfer the program to GSE33532 and assess its coverage.
4. Calculate tumor-state scores and aggregate repeated tumor subsamples at the patient level.
5. Test the primary association between the CAF-associated score and the Byers 76-gene NSCLC EMT-associated score.
6. Evaluate the secondary Wong/Liu broad stemness-associated score and an exploratory Zakaria-derived putative cancer-stem-cell-associated transcriptional state.
7. Perform histology-adjusted regression and influence/leave-one-out sensitivity analyses.

The Byers score is interpreted as an EMT-associated transcriptional signature, not as a comprehensive measure of all forms of tumor plasticity. The Wong/Liu score is described as stemness-associated, not as a direct measurement of cancer stem-cell identity. The Zakaria-derived analysis is exploratory.

## Main results

### Primary analysis: CAF-associated program versus Byers EMT-associated score

- Pearson correlation: **r = -0.036, p = 0.880**
- Spearman correlation: **rho = -0.017, p = 0.945**
- Histology-adjusted CAF coefficient: **+0.115, p = 0.929**

**Interpretation:** There was no detectable association between the CAF-associated transcriptional program and the Byers EMT-associated state in this cohort.

### Secondary analysis: CAF-associated program versus Wong stemness-associated score

- Pearson correlation: **r = -0.315, p = 0.176**
- Histology-adjusted CAF coefficient: **-0.059, p = 0.108**

**Interpretation:** The estimated direction was negative, but the association was not statistically significant. This is not evidence that stronger CAF-associated activity reduces stemness.

### Exploratory Zakaria-derived analysis

Patient-level Pearson correlations:

- Zakaria-derived score versus Byers EMT-associated score: **r = 0.297, p = 0.204**
- Zakaria-derived score versus Wong stemness-associated score: **r = -0.262, p = 0.264**
- Zakaria-derived score versus CAF-associated score: **r = -0.183, p = 0.440**

These exploratory results were not used to rescue the primary hypothesis.

## Sensitivity analyses

The repeated-sample structure was handled at the patient level. The primary CAF–EMT relationship remained essentially null in leave-one-out analysis. Sensitivity results should be interpreted in the context of the small sample size and the limited precision of estimates.

## Limitations

- The primary analysis included only 20 independent patients.
- The cohort included mixed histologies.
- The data are cross-sectional and observational.
- Bulk transcriptomic measurements cannot resolve cellular composition or directly establish cell-state transitions.
- The analyses do not establish causality or prove cell-to-cell interactions.
- The findings are cohort-specific and hypothesis-generating; they do not biologically refute a possible role for CAFs in tumor plasticity.

## Scientific strengths

- CAF-associated program discovery using matched samples.
- Explicit transfer of a fixed gene program to a separate tumor dataset.
- Patient-level inference that accounts for repeated tumor subsamples.
- Histology-adjusted analysis.
- Influence and leave-one-out sensitivity analyses.
- Clear separation of primary, secondary, and exploratory analyses.
- No post-hoc removal of observations to obtain a preferred result.
- No signature-shopping after the primary null result.
- Explicit distinction between association and causality.

## Reproducibility and data availability

The project uses publicly available GEO datasets:

- [GSE149955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149955)
- [GSE33532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33532)

Raw GEO data are not included in this repository. Obtain the source datasets from GEO and follow the analysis scripts in the order documented in the repository. The project was developed using R; the exact R and package versions should be recorded before claiming full environment-level reproducibility.

## Project contents

The public repository is intended to contain selected analysis scripts, publication-ready figures, and selected result tables. Raw expression files, large intermediate objects, session checkpoints, and temporary files should remain outside the public repository unless there is a clear, documented reason to include them.

## Final conclusion

In this small, heterogeneous NSCLC bulk transcriptomic cohort, a CAF-associated transcriptional program was not detectably associated with the Byers EMT-associated state at the patient level. The analysis provides a framework for testing this hypothesis while highlighting the need for larger cohorts and higher-resolution data to investigate cellular mechanisms underlying CAF–tumor plasticity relationships.
