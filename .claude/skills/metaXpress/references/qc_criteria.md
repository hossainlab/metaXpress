# metaXpress Study QC Criteria Reference

## Overview

The 10-point QC scoring system is applied by `mx_qc_study()` to each study.
Each criterion scores 1 point (pass) or 0 points (fail). Studies scoring
below `qc_threshold` (default: 7) are flagged by `mx_filter_studies()`.

Based on: Heberle et al. (2025) Alzheimer's & Dementia. doi:10.1002/alz.70025

---

## The 10 QC Criteria

### Criterion 1 — Minimum Sample Size
**Threshold:** n >= 3 biological replicates per condition group  
**Why:** Fewer than 3 replicates makes variance estimation unreliable in DESeq2/edgeR.  
**How to check:**
```r
min_n <- min(table(metadata$condition))
score[1] <- as.integer(min_n >= 3)
```

---

### Criterion 2 — Sequencing Depth
**Threshold:** Median library size >= 10 million reads  
**Why:** Low depth leads to high zero-inflation and unreliable count estimates.  
**How to check:**
```r
lib_sizes <- colSums(counts)
score[2] <- as.integer(median(lib_sizes) >= 10e6)
```
**Note:** If only normalized data is available (e.g., FPKM/TPM), this criterion
cannot be assessed — score 0 and flag in QC report.

---

### Criterion 3 — Alignment Rate
**Threshold:** Mean alignment rate >= 70%  
**Why:** Low alignment suggests contamination, wrong reference genome, or poor library quality.  
**How to check:** Parse from STAR/HISAT2 log files if available via GEO supplementary files.
```r
# If alignment stats are in metadata:
score[3] <- as.integer(mean(metadata$alignment_rate, na.rm = TRUE) >= 0.70)
# If unavailable: score[3] <- NA (not penalized, but flagged)
```

---

### Criterion 4 — rRNA Contamination
**Threshold:** rRNA reads < 10% of total reads  
**Why:** High rRNA contamination reduces effective sequencing depth for mRNA.  
**How to check:** From RSeQC or FastQC reports in GEO supplementary files.
```r
score[4] <- as.integer(metadata$rrna_fraction < 0.10)
# If unavailable: score[4] <- NA
```

---

### Criterion 5 — Duplicate Rate
**Threshold:** PCR duplicate rate < 50%  
**Why:** High duplication inflates apparent counts and reduces library complexity.  
**How to check:** From Picard MarkDuplicates or FastQC reports.
```r
score[5] <- as.integer(metadata$duplicate_rate < 0.50)
```

---

### Criterion 6 — Gene Detection Rate
**Threshold:** >= 15,000 genes detected (count > 0) in >= 50% of samples  
**Why:** Low gene detection indicates poor library quality or wrong organism annotation.  
**How to check:**
```r
n_detected <- apply(counts > 0, 1, sum)
genes_in_half <- sum(n_detected >= (ncol(counts) * 0.5))
score[6] <- as.integer(genes_in_half >= 15000)
```

---

### Criterion 7 — Metadata Completeness
**Threshold:** Required columns present and non-missing in >= 90% of samples  
**Required columns:** `condition`, `sample_id`  
**Recommended columns:** `age`, `sex`, `tissue`, `batch`  
**How to check:**
```r
required <- c("condition", "sample_id")
has_required <- all(required %in% colnames(metadata))
completeness <- mean(!is.na(metadata[, required]))
score[7] <- as.integer(has_required && completeness >= 0.90)
```

---

### Criterion 8 — Clear Case/Control Definition
**Threshold:** `condition` column has exactly 2 levels with unambiguous labels  
**Why:** Ambiguous group labels (e.g., "sample1", "sample2") prevent correct DE analysis.  
**How to check:**
```r
n_levels <- length(unique(metadata$condition))
# Acceptable labels (case-insensitive): control/case, normal/disease,
# healthy/patient, WT/KO, untreated/treated, etc.
score[8] <- as.integer(n_levels == 2)
# Flag if labels are numeric or non-descriptive
```

---

### Criterion 9 — No Batch-Condition Confounding
**Threshold:** Batch variable (if present) is not perfectly confounded with condition  
**Why:** Perfect confounding makes batch correction impossible and DE results uninterpretable.  
**How to check:**
```r
if ("batch" %in% colnames(metadata)) {
  # Check if batch and condition are perfectly confounded
  tab <- table(metadata$condition, metadata$batch)
  # Confounded if every batch contains only one condition
  confounded <- all(apply(tab, 2, function(x) sum(x > 0) == 1))
  score[9] <- as.integer(!confounded)
} else {
  score[9] <- 1  # No batch variable = no confounding risk
}
```

---

### Criterion 10 — Raw Counts Available
**Threshold:** Count matrix contains non-negative integers (not FPKM/TPM/log-transformed)  
**Why:** DESeq2 and edgeR require raw integer counts. Pre-normalized data cannot be used directly.  
**How to check:**
```r
is_integer_like <- all(counts == floor(counts), na.rm = TRUE)
has_no_negatives <- all(counts >= 0, na.rm = TRUE)
max_val <- max(counts, na.rm = TRUE)
# FPKM/TPM values are typically < 100,000 and non-integer
# Raw counts can be in the millions
score[10] <- as.integer(is_integer_like && has_no_negatives && max_val > 1000)
```

---

## mx_qc_study() Implementation Template

```r
mx_qc_study <- function(study) {
  counts   <- study@counts
  metadata <- study@metadata
  score    <- integer(10)
  details  <- character(10)

  # Criterion 1: Sample size
  min_n    <- min(table(metadata$condition))
  score[1] <- as.integer(min_n >= 3)
  details[1] <- sprintf("Min samples per group: %d (need >= 3)", min_n)

  # Criterion 2: Sequencing depth
  lib_sizes <- colSums(counts)
  med_depth <- median(lib_sizes)
  score[2]  <- as.integer(med_depth >= 10e6)
  details[2] <- sprintf("Median depth: %.1fM reads (need >= 10M)",
                         med_depth / 1e6)

  # Criterion 3: Alignment rate (if available)
  if ("alignment_rate" %in% colnames(metadata)) {
    mean_align <- mean(metadata$alignment_rate, na.rm = TRUE)
    score[3]   <- as.integer(mean_align >= 0.70)
    details[3] <- sprintf("Mean alignment: %.1f%% (need >= 70%%)",
                           mean_align * 100)
  } else {
    score[3]   <- NA_integer_
    details[3] <- "Alignment rate: not available"
  }

  # Criterion 4: rRNA contamination (if available)
  if ("rrna_fraction" %in% colnames(metadata)) {
    mean_rrna <- mean(metadata$rrna_fraction, na.rm = TRUE)
    score[4]  <- as.integer(mean_rrna < 0.10)
    details[4] <- sprintf("Mean rRNA: %.1f%% (need < 10%%)",
                           mean_rrna * 100)
  } else {
    score[4]  <- NA_integer_
    details[4] <- "rRNA fraction: not available"
  }

  # Criterion 5: Duplicate rate (if available)
  if ("duplicate_rate" %in% colnames(metadata)) {
    mean_dup <- mean(metadata$duplicate_rate, na.rm = TRUE)
    score[5] <- as.integer(mean_dup < 0.50)
    details[5] <- sprintf("Mean duplicate rate: %.1f%% (need < 50%%)",
                           mean_dup * 100)
  } else {
    score[5]  <- NA_integer_
    details[5] <- "Duplicate rate: not available"
  }

  # Criterion 6: Gene detection
  n_detected    <- apply(counts > 0, 1, sum)
  genes_in_half <- sum(n_detected >= (ncol(counts) * 0.5))
  score[6]      <- as.integer(genes_in_half >= 15000)
  details[6]    <- sprintf("Genes detected in >=50%% samples: %d (need >= 15000)",
                            genes_in_half)

  # Criterion 7: Metadata completeness
  required      <- c("condition", "sample_id")
  has_required  <- all(required %in% colnames(metadata))
  score[7]      <- as.integer(has_required)
  details[7]    <- sprintf("Required metadata columns present: %s",
                            ifelse(has_required, "yes", "no"))

  # Criterion 8: Clear case/control
  n_levels <- length(unique(metadata$condition))
  score[8] <- as.integer(n_levels == 2)
  details[8] <- sprintf("Condition levels: %d (need exactly 2)",
                         n_levels)

  # Criterion 9: Batch confounding
  if ("batch" %in% colnames(metadata)) {
    tab        <- table(metadata$condition, metadata$batch)
    confounded <- all(apply(tab, 2, function(x) sum(x > 0) == 1))
    score[9]   <- as.integer(!confounded)
    details[9] <- sprintf("Batch-condition confounding: %s",
                           ifelse(confounded, "YES (problem)", "no"))
  } else {
    score[9]   <- 1L
    details[9] <- "No batch variable present"
  }

  # Criterion 10: Raw counts
  is_int  <- all(counts == floor(counts), na.rm = TRUE)
  is_pos  <- all(counts >= 0, na.rm = TRUE)
  is_raw  <- is_int && is_pos && max(counts, na.rm = TRUE) > 1000
  score[10]   <- as.integer(is_raw)
  details[10] <- sprintf("Raw integer counts: %s", ifelse(is_raw, "yes", "no"))

  # Compute total score (treat NA as 0 for scoring, but flag separately)
  total_score <- sum(score, na.rm = TRUE)

  # Store results
  study@qc_score <- total_score
  attr(study, "qc_details") <- data.frame(
    criterion = seq_len(10),
    score     = score,
    details   = details,
    stringsAsFactors = FALSE
  )

  if (total_score < 7)
    warning(sprintf("Study %s has low QC score (%d/10). Consider excluding.",
                    study@accession, total_score))

  return(study)
}
```

---

## QC Score Interpretation

| Score | Interpretation | Recommendation |
|-------|---------------|----------------|
| 9–10 | Excellent | Include |
| 7–8 | Good | Include |
| 5–6 | Marginal | Include with caution; run sensitivity analysis |
| 3–4 | Poor | Exclude (default threshold) |
| 0–2 | Unacceptable | Exclude |

**Default `qc_threshold` in `mx_filter_studies()`: 7**

Users can lower to 5 for exploratory analyses or raise to 9 for strict analyses.
Always report the QC threshold used in methods sections.
