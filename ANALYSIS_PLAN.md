# Locked analysis plan

## Primary questions

1. Does WGD increase absolute PIK3CA dosage?
2. Does any WGD association persist after normalization to tumor ploidy?
3. Is chromosome 3q gain associated with local PIK3CA total copy number beyond genome-wide ploidy/WGD context?
4. Does the 3q association persist for mutant-specific dosage endpoints?
5. Is there evidence of systematic preferential enrichment of mutant-bearing PIK3CA copies?

## Prespecified WGD endpoints

- Total CN
- Total CN / ploidy
- Mutant CN
- Mutant CN / ploidy
- FAM

Two-sided Wilcoxon rank-sum tests; BH correction across all five endpoints.

## Prespecified adjusted 3q endpoints

- Total CN: `ploidy + WGD + 3q`
- Mutant CN: `ploidy + WGD + 3q`
- Mutant CN / ploidy: `WGD + 3q`
- FAM: `ploidy + WGD + 3q`

BH correction across these four adjusted 3q effects.

## Sensitivity analysis

Cook's distance threshold: `4/n`.

Primary models retain all evaluable tumors. A sensitivity model excludes observations above the threshold.

## External validation

MSK is used to validate:
- PIK3CA mutation prevalence
- mutation spectrum
- histologic context
- discrete CNA state composition

It is not presented as technical replication of the ABSOLUTE-derived mutant-dosage estimator.
