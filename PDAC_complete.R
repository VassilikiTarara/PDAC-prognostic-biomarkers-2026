# PDAC historical analysis reproduction
# English figures and tables; editable PPTX and DOCX.
# WARNING: reproduces the original INCORRECT external score calculation.
# Historical results are for audit/comparison, not corrected validation.
# Saved gene-level Cox estimates and reported Table 4 C-index CIs are imported,
# not refitted/rebootstrapped. GO/KEGG workbook exports are archival only.
# Independent glmnet refitting is optional and cannot guarantee old coefficients.
#
# Run in a new R session: source("path/to/PDAC_legacy_reproduction_complete.R")
# The original workbook is never modified.
# Outputs are written to results/legacy_reproduction and replaced on rerun.

project_dir <- "C:/Users/User/OneDrive/Desktop/PhD_Med_Auth/PDAC reproduction"
run_optional_lasso_refit <- FALSE  # Diagnostic only; exact historical coefficients are reused.
col_low <- "#183A5A"
col_high <- "#C47A36"
col_other <- "#8B9AA7"

# ---------- Packages ----------
cran_packages <- c("readxl", "ggplot2", "ggrepel", "survival", "survminer",
                   "patchwork", "officer", "rvg", "flextable")
missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                               quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

# ---------- Shared helpers ----------
read_sheet <- function(file, sheet) {
  as.data.frame(readxl::read_excel(file, sheet = sheet, na = c("", "NA")))
}
clean_id <- function(x) trimws(sub("\\.\\.\\.[0-9]+$", "", x))
num <- function(x) as.numeric(as.character(x))
assert <- function(test, message) if (!isTRUE(test)) stop(message, call. = FALSE)
num_matrix <- function(df) {
  raw <- as.matrix(df)
  x <- trimws(as.character(raw))
  missing <- is.na(x) | x %in% c("", "NA", "NaN")
  values <- suppressWarnings(as.numeric(x))
  assert(!any(!missing & is.na(values)), "Unexpected nonnumeric expression values.")
  matrix(values, nrow = nrow(raw), ncol = ncol(raw), dimnames = dimnames(raw))
}
fetch <- function(url, path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    options(timeout = max(600, getOption("timeout")))
    download.file(url, path, mode = "wb", method = "libcurl")
  }
  assert(file.exists(path) && file.info(path)$size > 0, paste("Download failed:", path))
}
logrank <- function(data, time, event, group) {
  f <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, event, group))
  z <- survival::survdiff(f, data = data)
  pchisq(z$chisq, length(z$n) - 1, lower.tail = FALSE)
}
cox_row <- function(fit, term, label, cohort) {
  s <- summary(fit)
  data.frame(Cohort = cohort, Model = label, Patients = fit$n, Deaths = fit$nevent,
    HR = unname(s$conf.int[term, "exp(coef)"]),
    CI_low = unname(s$conf.int[term, "lower .95"]),
    CI_high = unname(s$conf.int[term, "upper .95"]),
    P_value = unname(s$coefficients[term, "Pr(>|z|)"]),
    C_index = unname(s$concordance[1]), check.names = FALSE)
}
cox_terms <- function(fit, model) {
  s <- summary(fit)
  data.frame(Model = model, Term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"], CI_low = s$conf.int[, "lower .95"],
    CI_high = s$conf.int[, "upper .95"],
    P_value = s$coefficients[, "Pr(>|z|)"], row.names = NULL)
}
model_display <- function(x) {
  x$`HR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", x$HR, x$CI_low, x$CI_high)
  x$`P value` <- format.pval(x$P_value, digits = 4, eps = 0.001)
  x$`C-index` <- sprintf("%.3f", x$C_index)
  x[, c("Cohort", "Model", "Patients", "Deaths", "HR (95% CI)", "P value", "C-index")]
}
plot_theme <- function(size = 12) {
  ggplot2::theme_classic(base_size = size, base_family = "Arial") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", color = col_low),
      plot.subtitle = ggplot2::element_text(color = "#56616B", size = 10),
      axis.text = ggplot2::element_text(color = "#303840"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, color = "#56616B"))
}
km_plot <- function(data, time, event, group, title, subtitle, max_time, step, caption = NULL) {
  f <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, event, group))
  fit <- survival::survfit(f, data = data)
  # survminer evaluates the model formula from the fit call.
  fit$call$formula <- f
  group_cox <- survival::coxph(f, data = data, ties = "efron")
  ci <- summary(group_cox)$conf.int[1, ]
  counts <- table(data[[group]])
  p <- logrank(data, time, event, group)
  annotation <- sprintf("Log-rank p = %.5f\nHigh vs Low: HR %.2f (95%% CI %.2f–%.2f)\nPatients = %d; deaths = %d",
    p, ci["exp(coef)"], ci["lower .95"], ci["upper .95"], nrow(data), sum(data[[event]]))
  unit <- if (time == "OS.time") "days" else "months"
  parts <- survminer::ggsurvplot(fit, data = data, palette = c(col_low, col_high),
    conf.int = TRUE, conf.int.alpha = 0.12, censor = TRUE, censor.shape = 3,
    censor.size = 2.5, linewidth = 0.9, risk.table = TRUE, risk.table.col = "strata",
    risk.table.y.text = TRUE, risk.table.y.text.col = TRUE,
    legend = "top", legend.title = "Score group",
    legend.labs = c(sprintf("Low (n = %d)", counts["Low"]), sprintf("High (n = %d)", counts["High"])),
    xlim = c(0, max_time), break.time.by = step, ylim = c(0, 1),
    xlab = paste0("Follow-up time (", unit, ")"), ylab = "Overall survival probability",
    ggtheme = plot_theme(),
    tables.theme = survminer::theme_cleantable(base_size = 11, base_family = "Arial"))
  main <- parts$plot + ggplot2::annotate("text", x = max_time / 30, y = 0.05,
    label = annotation, hjust = 0, vjust = 0, size = 3.4, family = "Arial", color = col_low)
  risk <- parts$table + ggplot2::labs(title = "Number at risk", y = NULL) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold", color = col_low),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 15)))
  patchwork::wrap_plots(main, risk, ncol = 1, heights = c(3.8, 1.15)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
      caption = paste("Shaded bands: 95% confidence intervals. Tick marks: censored observations.",
                      "HR: hazard ratio; CI: confidence interval.", caption),
      theme = ggplot2::theme(text = ggplot2::element_text(family = "Arial"),
        plot.title = ggplot2::element_text(size = 16, face = "bold", color = col_low),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#56616B"),
        plot.caption = ggplot2::element_text(size = 9, hjust = 0, color = "#56616B")))
}
risk_plot <- function(data, score, title, historical = FALSE) {
  d <- data[order(data[[score]]), , drop = FALSE]
  d$Rank <- seq_len(nrow(d)); d$Score <- d[[score]]
  cutoff <- median(d$Score)
  d$Group <- factor(ifelse(d$Score >= cutoff, "High", "Low"), levels = c("Low", "High"))
  n <- table(d$Group)
  ggplot2::ggplot(d, ggplot2::aes(Rank, Score, fill = Group)) +
    ggplot2::geom_col(width = 0.90) +
    ggplot2::geom_hline(yintercept = 0, color = "#59636D", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = cutoff, color = "#59636D", linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c(Low = col_low, High = col_high),
      labels = c(sprintf("Low (n = %d)", n["Low"]), sprintf("High (n = %d)", n["High"]))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(title = title, subtitle = sprintf("%d patients; median cutoff = %.3f", nrow(d), cutoff),
      x = "Patient rank (ordered by increasing score)", y = "Five-gene score", fill = "Score group",
      caption = paste("Dashed line: cohort median.", if (historical) "Historical incorrect-score calculation." else "")) +
    plot_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#E8EBEE"))
}
add_plot <- function(ppt, graph) {
  ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
  sz <- officer::slide_size(ppt)
  officer::ph_with(ppt, rvg::dml(code = print(graph), editable = TRUE),
    location = officer::ph_location(left = 0.25, top = 0.20,
      width = sz$width - 0.5, height = sz$height - 0.4))
}
word_table <- function(doc, x, heading, note = NULL) {
  doc <- officer::body_add_par(doc, heading, style = "heading 2")
  ft <- flextable::flextable(x)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::bg(ft, bg = col_low, part = "header")
  ft <- flextable::color(ft, color = "white", part = "header")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
  doc <- flextable::body_add_flextable(doc, ft)
  if (!is.null(note)) doc <- officer::body_add_par(doc, note)
  doc
}

# IPCW Brier score for a Cox model refitted within the evaluation cohort.
# Apparent/refitted prediction accuracy, not frozen-model external accuracy.
ipcw_brier <- function(fit, data, time, event, times) {
  assert(fit$n == nrow(data), "Brier calculation requires complete model data.")
  cf <- survival::survfit(stats::as.formula(sprintf("survival::Surv(%s, 1 - %s) ~ 1", time, event)), data = data)
  G <- function(t, before = FALSE) vapply(t, function(tt) {
    idx <- which(if (before) cf$time < tt else cf$time <= tt)
    if (length(idx)) as.numeric(cf$surv[max(idx)]) else 1
  }, numeric(1))
  sf <- survival::survfit(fit, newdata = data, se.fit = FALSE)
  vapply(times, function(t) {
    pred <- as.numeric(summary(sf, times = t, extend = TRUE)$surv)
    assert(length(pred) == nrow(data), "Unexpected survival prediction dimensions.")
    dead <- data[[time]] <= t & data[[event]] == 1
    alive <- data[[time]] > t
    values <- numeric(nrow(data))
    if (any(dead)) {
      g <- G(data[[time]][dead], before = TRUE)
      assert(all(g > 0), "Zero IPCW censoring survival at an event time.")
      values[dead] <- pred[dead]^2 / g
    }
    if (any(alive)) {
      g <- G(t); assert(g > 0, "Zero IPCW censoring survival at evaluation time.")
      values[alive] <- (1 - pred[alive])^2 / g
    }
    mean(values)
  }, numeric(1))
}
ibs <- function(b, times) sum(diff(times) * (head(b, -1) + tail(b, -1)) / 2) / diff(range(times))

# ---------- Main workflow ----------
main <- function() {
  assert(dir.exists(project_dir), paste("Project directory not found:", project_dir))
  input_dir <- file.path(project_dir, "input")
  results_dir <- file.path(project_dir, "results", "legacy_reproduction")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(results_dir, "PDAC_legacy_run_log.txt")
  sink(logfile, split = TRUE)
  on.exit({
    capture.output(sessionInfo(), file = file.path(results_dir, "R_session_info.txt"))
    sink()
  }, add = TRUE)
  cat("Historical PDAC reproduction started:", format(Sys.time()), "\n")
  cat("IMPORTANT: external risk score deliberately reproduces the historical error.\n")
  wb <- file.path(input_dir, "126295-Supplementary Material.xlsm")
  if (!file.exists(wb)) wb <- file.path(input_dir, "126295-Supplementary Material(1).xlsm")
  assert(file.exists(wb), "Supplementary workbook not found in input.")
  read <- function(sheet) read_sheet(wb, sheet)
  save_csv <- function(x, name) write.csv(x, file.path(results_dir, name), row.names = FALSE)
  needed_sheets <- c("tumor", "normal", "CPTAC_DE_all", "CPTAC_DE_significant",
    "LASSO_patient_scores", "LASSO_selected_genes", "Initial_signature_scores",
    "TCGA_multivariable_Cox", "GSE71729_validation_dataset")
  assert(all(needed_sheets %in% readxl::excel_sheets(wb)), "Required workbook sheets are missing.")

  # 1. Reproduce paired CPTAC differential abundance.
  cat("\n1. CPTAC paired differential abundance\n")
  tumor <- read("tumor"); normal <- read("normal")
  tids <- clean_id(names(tumor)[-1]); nids <- clean_id(names(normal)[-1])
  unique_nids <- names(table(nids))[table(nids) == 1]
  pairs <- intersect(tids, unique_nids)
  pairs <- pairs[nzchar(pairs)]
  assert(!anyDuplicated(tids), "Duplicate tumor sample IDs.")
  assert(length(pairs) == 80, "Expected 80 uniquely matched pairs; inspect workbook.")
  tg <- trimws(as.character(tumor[[1]])); ng <- trimws(as.character(normal[[1]]))
  tr <- which(!is.na(tg) & nzchar(tg) & !duplicated(tg))
  nr <- which(!is.na(ng) & nzchar(ng) & !duplicated(ng))
  genes <- intersect(tg[tr], ng[nr])
  tm <- num_matrix(tumor[tr[match(genes, tg[tr])], match(pairs, tids) + 1L, drop = FALSE])
  nm <- num_matrix(normal[nr[match(genes, ng[nr])], match(pairs, nids) + 1L, drop = FALSE])
  rownames(tm) <- rownames(nm) <- genes
  colnames(tm) <- colnames(nm) <- pairs
  expr <- cbind(nm, tm)
  # Original final filter: >=70% observed in EACH tissue group (56/80).
  keep <- rowSums(!is.na(nm)) >= ceiling(0.7 * ncol(nm)) &
          rowSums(!is.na(tm)) >= ceiling(0.7 * ncol(tm))
  patient <- factor(rep(pairs, 2))
  group <- factor(rep(c("Normal", "Tumor"), each = length(pairs)), levels = c("Normal", "Tumor"))
  design <- model.matrix(~ patient + group)
  fit <- limma::eBayes(limma::lmFit(expr[keep, , drop = FALSE], design))
  res <- limma::topTable(fit, coef = "groupTumor", number = Inf, adjust.method = "BH", sort.by = "P")
  res$Gene <- rownames(res)
  sig <- res[is.finite(res$adj.P.Val) & res$adj.P.Val < 0.05 & abs(res$logFC) >= 1, , drop = FALSE]
  saved <- read("CPTAC_DE_all")
  saved <- saved[!is.na(saved$Gene), , drop = FALSE]
  idx <- match(saved$Gene, res$Gene)
  assert(nrow(res) == 8056 && nrow(sig) == 29 && !anyNA(idx), "CPTAC gene counts do not match historical results.")
  delta_fc <- max(abs(num(saved$logFC) - res$logFC[idx]))
  delta_p <- max(abs(num(saved$P.Value) - res$P.Value[idx]))
  assert(delta_fc < 1e-8 && delta_p < 1e-8, "CPTAC model differs from saved results.")
  cat("Pairs:", length(pairs), "; tested:", nrow(res), "; significant:", nrow(sig), "\n")
  cat("Maximum logFC difference:", delta_fc, "; P-value difference:", delta_p, "\n")
  save_csv(res, "CPTAC_DE_all_reproduced.csv"); save_csv(sig, "CPTAC_DE_significant_reproduced.csv")

  vd <- res
  vd$category <- factor(ifelse(vd$adj.P.Val < 0.05 & vd$logFC <= -1, "Lower in tumor",
    ifelse(vd$adj.P.Val < 0.05 & vd$logFC >= 1, "Higher in tumor", "Other")),
    levels = c("Lower in tumor", "Other", "Higher in tumor"))
  vd$neglog <- -log10(pmax(vd$adj.P.Val, .Machine$double.xmin))
  labelled <- vd[vd$Gene %in% c("FKBP11", "SLC43A1", "SLC16A10", "PRSS3", "AMY2B") & vd$category != "Other", ]
  p_volcano <- ggplot2::ggplot(vd, ggplot2::aes(logFC, neglog, color = category)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#8B9AA7") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#8B9AA7") +
    ggrepel::geom_text_repel(data = labelled, ggplot2::aes(label = Gene), seed = 123, size = 3) +
    ggplot2::scale_color_manual(values = c("Lower in tumor" = col_low, Other = col_other, "Higher in tumor" = col_high), drop = FALSE) +
    ggplot2::labs(title = "Paired CPTAC proteomic differential expression",
      subtitle = "80 matched pairs; 8,056 proteins tested; 29 proteins meeting both thresholds",
      x = "Log2 fold change (tumor vs adjacent normal)", y = "−Log10(adjusted p-value)", color = NULL) + plot_theme()
  heat_genes <- head(sig$Gene[order(sig$adj.P.Val)], 20)
  z <- t(scale(t(expr[heat_genes, , drop = FALSE])))
  # Keep the matrix dimensions: scalar-first pmax/pmin can discard attributes.
  z[z < -2] <- -2; z[z > 2] <- 2
  hd <- expand.grid(Gene = rownames(z), Sample = seq_len(ncol(z)), stringsAsFactors = FALSE)
  hd$Z <- as.vector(z); hd$Gene <- factor(hd$Gene, levels = rev(heat_genes))
  p_heat <- ggplot2::ggplot(hd, ggplot2::aes(Sample, Gene, fill = Z)) +
    ggplot2::geom_tile() + ggplot2::geom_vline(xintercept = 80.5, color = "white", linewidth = 0.8) +
    ggplot2::scale_fill_gradient2(low = col_low, mid = "#F7F4EE", high = col_high,
      limits = c(-2, 2), na.value = "#E2E5E8", name = "Protein-wise\nZ score") +
    ggplot2::scale_x_continuous(breaks = c(40.5, 120.5), labels = c("Normal (n = 80)", "Tumor (n = 80)"), expand = c(0, 0)) +
    ggplot2::labs(title = "Differentially expressed proteins in CPTAC-PDAC",
      subtitle = "20 significant proteins with the smallest adjusted p-values",
      x = "Paired CPTAC samples", y = NULL, caption = "Gray tiles indicate missing measurements. Displayed Z scores are capped at ±2.") +
    plot_theme(11) + ggplot2::theme(legend.position = "right", axis.text.y = ggplot2::element_text(size = 9))

  # Archive existing enrichment results without rerunning database-dependent analyses.
  for (sheet in intersect(c("GO_results", "GO_results_downregulated", "KEGG_results_all",
                            "KEGG_results_downregulated"), readxl::excel_sheets(wb))) {
    save_csv(read(sheet), paste0("Archived_", sheet, ".csv"))
  }

  # 2. TCGA RNA and saved historical cohort.
  cat("\n2. TCGA expression and historical coefficients\n")
  tcga_file <- file.path(input_dir, "TCGA_PAAD_HiSeqV2.gz")
  fetch("https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/HiSeqV2.gz", tcga_file)
  rna <- read.delim(gzfile(tcga_file), check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  assert(!anyDuplicated(rna[[1]]), "Duplicate TCGA gene names.")
  samples <- names(rna)[-1]; samples <- samples[substr(samples, 14, 15) == "01"]
  patients <- substr(samples, 1, 12)
  assert(!anyDuplicated(patients), "Duplicate TCGA primary-tumor patients.")
  clin <- read("LASSO_patient_scores")
  clin <- clin[!is.na(clin$patient_id), , drop = FALSE]
  ids <- as.character(clin$patient_id)
  assert(nrow(clin) == 177 && !anyDuplicated(ids) && all(ids %in% patients), "Historical TCGA cohort mismatch.")
  candidate_genes <- intersect(sig$Gene, as.character(rna[[1]]))
  x <- t(num_matrix(rna[match(candidate_genes, rna[[1]]), match(samples[match(ids, patients)], names(rna)), drop = FALSE]))
  rownames(x) <- ids; colnames(x) <- candidate_genes
  assert(ncol(x) == 28 && all(is.finite(x)), "Expected 28 complete TCGA candidate genes.")
  sc <- read("LASSO_selected_genes")
  sc <- sc[!is.na(sc$Gene) & !is.na(sc$Coefficient), c("Gene", "Coefficient")]
  beta <- setNames(num(sc$Coefficient), as.character(sc$Gene))
  assert(all(names(beta) %in% colnames(x)), "A signature gene is missing from TCGA.")
  score <- as.numeric(x[, names(beta), drop = FALSE] %*% beta)
  assert(max(abs(score - num(clin$lasso_score))) < 1e-8, "TCGA score reconstruction mismatch.")
  tcga <- clin; tcga$lasso_score <- score
  tcga$OS <- num(tcga$OS); tcga$OS.time <- num(tcga$OS.time); tcga$age <- num(tcga$age)
  assert(all(is.finite(tcga$OS.time)) && all(tcga$OS.time > 0) &&
         all(tcga$OS %in% c(0, 1)), "Invalid TCGA survival data.")
  tcga$sex <- factor(tcga$sex, levels = c("Female", "Male"))
  tcga$stage_simple <- factor(tcga$stage_simple, levels = paste("Stage", c("I", "II", "III", "IV")))
  tcga$risk_group <- factor(ifelse(score >= median(score), "High", "Low"), levels = c("Low", "High"))
  assert(all(as.character(tcga$risk_group) == as.character(clin$group)), "TCGA group labels mismatch.")
  tcga_uni <- survival::coxph(survival::Surv(OS.time, OS) ~ lasso_score, data = tcga, ties = "efron")
  tcga_adj <- survival::coxph(survival::Surv(OS.time, OS) ~ lasso_score + age + sex + stage_simple, data = tcga, ties = "efron")
  save_csv(tcga, "TCGA_historical_patient_scores.csv"); save_csv(sc, "Historical_LASSO_coefficients.csv")

  # Optional refit: software changes/folds can change the original coefficients.
  if (run_optional_lasso_refit) {
    if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet", repos = "https://cloud.r-project.org")
    set.seed(123)
    cv <- glmnet::cv.glmnet(x, survival::Surv(tcga$OS.time, tcga$OS),
      family = "cox", alpha = 1, nfolds = 5, cox.ties = "breslow")
    saveRDS(cv, file.path(results_dir, "Diagnostic_current_glmnet_refit.rds"))
    cc <- as.matrix(stats::coef(cv, s = "lambda.min"))
    save_csv(data.frame(Gene = rownames(cc), Coefficient = cc[, 1]), "Diagnostic_current_glmnet_coefficients.csv")
  }

  # 3. Initial biological signature: recompute row-wise Z scores and negative mean.
  cat("\n3. Initial signature and TCGA models\n")
  initial_saved <- read("Initial_signature_scores")
  initial_saved <- initial_saved[!is.na(initial_saved$patient_id), , drop = FALSE]
  assert(setequal(initial_saved$patient_id, ids), "Initial-signature cohort differs from TCGA cohort.")
  biological_genes <- intersect(sig$Gene[sig$logFC < 0], colnames(x))
  initial_score <- -rowMeans(scale(x[, biological_genes, drop = FALSE]), na.rm = TRUE)
  initial <- tcga
  initial$signature_score <- initial_score
  initial_delta <- max(abs(initial_score - num(initial_saved$signature_score[match(ids, initial_saved$patient_id)])))
  assert(initial_delta < 1e-8, "Recomputed initial signature differs from saved scores.")
  initial_uni <- survival::coxph(survival::Surv(OS.time, OS) ~ signature_score, data = initial, ties = "efron")
  initial_adj <- survival::coxph(survival::Surv(OS.time, OS) ~ signature_score + age + sex + stage_simple, data = initial, ties = "efron")
  save_csv(initial, "TCGA_initial_signature_recomputed.csv")
  p_tcga_km <- km_plot(tcga, "OS.time", "OS", "risk_group", "Overall survival in the TCGA-PAAD cohort",
    sprintf("Five-gene score; cohort median cutoff = %.3f", median(score)), 3000, 500)
  p_tcga_risk <- risk_plot(tcga, "lasso_score", "Risk-score distribution in TCGA-PAAD")
  ta <- cox_terms(tcga_adj, "Adjusted score model")
  ta$Label <- c("Five-gene score", "Age (per year)", "Male vs Female", "Stage II vs I", "Stage III vs I", "Stage IV vs I")[match(ta$Term,
    c("lasso_score", "age", "sexMale", "stage_simpleStage II", "stage_simpleStage III", "stage_simpleStage IV"))]
  ta$Label <- factor(ta$Label, levels = rev(ta$Label))
  p_tcga_forest <- ggplot2::ggplot(ta, ggplot2::aes(HR, Label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#8B9AA7") +
    ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = Label), color = col_high, linewidth = 0.8) +
    ggplot2::geom_point(color = col_low, size = 3) + ggplot2::scale_x_log10() +
    ggplot2::labs(title = "Multivariable Cox regression in TCGA-PAAD", subtitle = sprintf("Age-, sex- and stage-adjusted score model; n = %d", tcga_adj$n),
      x = "Hazard ratio (95% CI; log scale)", y = NULL) + plot_theme()

  # Saved gene-level forest: no assumption that all 28 genes entered one model.
  gene_cox <- read("TCGA_multivariable_Cox")
  gene_cox <- gene_cox[!is.na(gene_cox$Gene), c("Gene", "HR", "CI_low", "CI_high", "P.Value", "adj.P.Val")]
  gene_cox <- gene_cox[order(gene_cox$HR, decreasing = TRUE), ]
  gene_cox$Gene <- factor(gene_cox$Gene, levels = rev(gene_cox$Gene))
  p_gene_forest <- ggplot2::ggplot(gene_cox, ggplot2::aes(HR, Gene)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#8B9AA7") +
    ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = Gene), color = col_high, linewidth = 0.6) +
    ggplot2::geom_point(color = col_low, size = 2) +
    ggplot2::labs(title = "Multivariable Cox regression – PDAC", subtitle = "Historical gene-level results from TCGA_multivariable_Cox",
      x = "Hazard ratio (95% CI)", y = NULL) + plot_theme(10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
  save_csv(gene_cox, "TCGA_saved_gene_level_Cox.csv")

  # 4. GEO: retain the ENTIRE matrix sample order for the historical bug.
  cat("\n4. GSE71729 historical external score\n")
  external <- read("GSE71729_validation_dataset")
  external <- external[!is.na(external$sample_id), , drop = FALSE]
  assert(nrow(external) == 125 && !anyDuplicated(external$sample_id), "Expected 125 unique external patients.")
  geo_file <- file.path(input_dir, "GSE71729_series_matrix.txt.gz")
  fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71729/matrix/GSE71729_series_matrix.txt.gz", geo_file)
  con <- gzfile(geo_file, "rt")
  metadata <- character(); header <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) { close(con); stop("GEO expression table not found.") }
    if (identical(line, "!series_matrix_table_begin")) break
    metadata <- c(metadata, line)
  }
  header <- strsplit(gsub('"', '', readLines(con, n = 1, warn = FALSE), fixed = TRUE), "\t", fixed = TRUE)[[1]]
  original_beta <- c(PRSS3 = 0.107430789, SLC16A10 = 0.014923690,
                     FKBP11 = -0.132749545, SLC43A1 = -0.109866376, AMY2B = -0.012718478)
  selected <- character()
  repeat {
    chunk <- readLines(con, n = 1000, warn = FALSE)
    if (!length(chunk)) break
    first <- gsub('"', '', sub("\t.*", "", chunk), fixed = TRUE)
    selected <- c(selected, chunk[first %in% names(original_beta)])
    if (any(chunk == "!series_matrix_table_end")) break
  }
  close(con)
  assert(length(selected) == 5, "Expected one GEO row per signature gene.")
  geo_rows <- do.call(rbind, strsplit(gsub('"', '', selected, fixed = TRUE), "\t", fixed = TRUE))
  gm <- matrix(num(geo_rows[, -1]), nrow = 5, dimnames = list(geo_rows[, 1], header[-1]))
  gm <- gm[names(original_beta), , drop = FALSE]
  assert(all(is.finite(gm)) && !anyDuplicated(colnames(gm)), "Invalid GEO expression matrix.")
  writeLines(metadata, file.path(results_dir, "GSE71729_original_metadata.txt"))
  # DELIBERATE HISTORICAL ERROR: vector recycling across the transposed matrix.
  # DO NOT subset to 125 patients or reorder columns before this calculation.
  legacy_all <- colSums(t(t(gm) * unname(original_beta)))
  corrected_all <- colSums(sweep(gm, 1, original_beta, "*"))
  ext_ids <- as.character(external$sample_id)
  assert(all(ext_ids %in% colnames(gm)), "External samples missing from GEO.")
  legacy <- data.frame(sample_id = ext_ids, OS_time = num(external$OS_time),
    OS_event = num(external$OS_event), risk_score = unname(legacy_all[ext_ids]))
  assert(all(is.finite(legacy$OS_time)) && all(legacy$OS_time >= 0) && all(legacy$OS_event %in% c(0, 1)), "Invalid external survival data.")
  score_delta <- max(abs(legacy$risk_score - num(external$risk_score)))
  assert(score_delta < 1e-8, "Legacy external score differs from Excel: check full GEO sample order.")
  legacy$risk_group <- factor(ifelse(legacy$risk_score >= median(legacy$risk_score), "High", "Low"), levels = c("Low", "High"))
  assert(all(as.character(legacy$risk_group) == as.character(external$risk_group)), "Legacy external group mismatch.")
  # Keep numerical subtype codes for model fitting. Biological labels below
  # follow the historical workbook, not inference from survival differences.
  legacy$tumor_subtype <- factor(as.character(external$tumor_subtype), levels = c("1", "2"))
  legacy$stroma_subtype <- factor(as.character(external$stroma_subtype), levels = c("1", "2", "3"))
  assert(!anyNA(legacy), "Missing historical external data.")
  ext_uni <- survival::coxph(survival::Surv(OS_time, OS_event) ~ risk_score, data = legacy, ties = "efron")
  ext_adj <- survival::coxph(survival::Surv(OS_time, OS_event) ~ risk_score + tumor_subtype + stroma_subtype, data = legacy, ties = "efron")
  ext_group <- survival::coxph(survival::Surv(OS_time, OS_event) ~ risk_group, data = legacy, ties = "efron")
  ext_p <- logrank(legacy, "OS_time", "OS_event", "risk_group")
  cat("Score difference:", score_delta, "; historical log-rank p:", ext_p, "\n")
  assert(abs(ext_p - 0.0436192) < 1e-6, "Historical log-rank result mismatch.")
  save_csv(legacy, "GSE71729_legacy_incorrect_patient_scores.csv")
  save_csv(data.frame(sample_id = ext_ids, Legacy_score = legacy$risk_score,
    Correct_gene_wise_score = unname(corrected_all[ext_ids])), "GSE71729_score_calculation_audit.csv")
  p_ext_km <- km_plot(legacy, "OS_time", "OS_event", "risk_group", "Overall survival in the GSE71729 cohort",
    sprintf("Historical incorrect-score calculation; median cutoff = %.3f", median(legacy$risk_score)), 60, 12,
    "Historical reproduction of the incorrect score calculation.")
  p_ext_risk <- risk_plot(legacy, "risk_score", "Risk-score distribution in GSE71729", TRUE)

  eu <- cox_terms(ext_uni, "Univariable"); eg <- cox_terms(ext_group, "Univariable")
  ea <- cox_terms(ext_adj, "Subtype-adjusted")
  ef <- rbind(eu, eg, ea)
  lookup <- c(risk_score = "Continuous risk score", risk_groupHigh = "High vs Low",
    tumor_subtype2 = "Tumour subtype: Basal vs Classical",
    stroma_subtype2 = "Stroma subtype: Normal vs Low",
    stroma_subtype3 = "Stroma subtype: Activated vs Low")
  ef$Label <- factor(unname(lookup[ef$Term]), levels = rev(unname(lookup)))
  ef$Significance <- factor(ifelse(ef$P_value < 0.05, "p < 0.05", "p ≥ 0.05"), levels = c("p ≥ 0.05", "p < 0.05"))
  dodge <- ggplot2::position_dodge(width = 0.35)
  p_ext_forest <- ggplot2::ggplot(ef, ggplot2::aes(HR, Label, color = Model, group = Model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#8B9AA7") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.13, position = dodge) +
    ggplot2::geom_point(ggplot2::aes(shape = Significance), size = 3, position = dodge) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(Univariable = col_low, "Subtype-adjusted" = col_high)) +
    ggplot2::scale_shape_manual(values = c("p ≥ 0.05" = 16, "p < 0.05" = 17)) +
    ggplot2::labs(title = "External validation: hazard ratios", subtitle = "GSE71729: historical incorrect-score calculation",
      x = "Hazard ratio (95% CI; log scale)", y = NULL, color = "Model", shape = NULL) +
    plot_theme(11) + ggplot2::theme(legend.position = "bottom")
  save_csv(ef, "GSE71729_legacy_Cox_terms.csv")

  # 5. Computed performance and prediction accuracy.
  cat("\n5. Tables and prediction accuracy\n")
  models <- rbind(cox_row(initial_uni, "signature_score", "Initial signature: univariable", "TCGA-PAAD"),
    cox_row(initial_adj, "signature_score", "Initial signature: adjusted", "TCGA-PAAD"),
    cox_row(tcga_uni, "lasso_score", "Five-gene score: univariable", "TCGA-PAAD"),
    cox_row(tcga_adj, "lasso_score", "Five-gene score: adjusted", "TCGA-PAAD"),
    cox_row(ext_uni, "risk_score", "Historical external score: univariable", "GSE71729"),
    cox_row(ext_adj, "risk_score", "Historical external score: subtype-adjusted", "GSE71729"),
    cox_row(ext_group, "risk_groupHigh", "Historical High vs Low", "GSE71729"))
  save_csv(models, "Recomputed_Cox_model_performance.csv")
  fits <- list(tcga_uni, tcga_adj, initial_uni, initial_adj, ext_uni, ext_adj)
  datas <- list(tcga, tcga[stats::complete.cases(tcga[, c("OS.time", "OS", "lasso_score", "age", "sex", "stage_simple")]), ],
    initial, initial[stats::complete.cases(initial[, c("OS.time", "OS", "signature_score", "age", "sex", "stage_simple")]), ], legacy, legacy)
  terms <- c("lasso_score", "lasso_score", "signature_score", "signature_score", "risk_score", "risk_score")
  labels <- c("Five-gene score: univariable", "Five-gene score: adjusted", "Initial signature: univariable", "Initial signature: adjusted",
    "Historical external score: univariable", "Historical external score: subtype-adjusted")
  accuracy <- do.call(rbind, lapply(seq_along(fits), function(i) {
    ext <- i >= 5; tname <- if (ext) "OS_time" else "OS.time"; ename <- if (ext) "OS_event" else "OS"
    # TCGA horizons use 365 days per year; external horizons use 12 months.
    times <- if (ext) c(12, 24, 36) else c(365, 730, 1095)
    grid <- if (ext) 0:36 else seq(0, 1095, by = 1095 / 36)
    b <- ipcw_brier(fits[[i]], datas[[i]], tname, ename, times)
    ig <- ibs(ipcw_brier(fits[[i]], datas[[i]], tname, ename, grid), grid)
    mean_surv <- "Not calculated"
    if (i %in% c(1, 3, 5)) {
      nd <- setNames(data.frame(mean(datas[[i]][[terms[i]]])), terms[i])
      sf <- survival::survfit(fits[[i]], newdata = nd)
      mean_surv <- paste(sprintf("%.1f", 100 * as.numeric(summary(sf, times = times, extend = TRUE)$surv)), collapse = "/")
    }
    data.frame(Cohort = if (ext) "GSE71729" else "TCGA-PAAD", Model = labels[i], Patients = fits[[i]]$n,
      `Score coefficient (Cox refit)` = sprintf("%.3f", stats::coef(fits[[i]])[terms[i]]),
      `Brier at 1/2/3 years` = paste(sprintf("%.3f", b), collapse = "/"),
      `Integrated Brier (0–3 years)` = sprintf("%.3f", ig),
      `Survival at mean score, 1/2/3 years (%)` = mean_surv, check.names = FALSE)
  }))
  save_csv(accuracy, "Recomputed_refitted_prediction_accuracy.csv")

  # Historical published tables: preserve supplied values, do NOT silently substitute computed values.
  table4 <- data.frame(`Prognostic model` = c(rep("Initial biologically derived signature", 2), rep("LASSO-derived 5-gene signature", 2)),
    `Cox model` = rep(c("Univariable", "Multivariable"), 2), `Hazard ratio` = c("1.13", "1.08", "6.41", "6.13"),
    `P value` = c("0.239", "0.480", "<0.001", "<0.001"), `Harrell's C-index` = c("0.522", "0.563", "0.587", "0.599"),
    `C-index 95% CI` = c("0.472–0.579", "0.498–0.626", "0.526–0.634", "0.534–0.653"), check.names = FALSE)
  table5 <- data.frame(Cohort = c(rep("TCGA-PAAD (n = 177; 93 events)", 4), rep("GSE71729 (n = 125; 84 events)", 2)),
    Model = c("LASSO 5-gene: univariable", "LASSO 5-gene: multivariable", "Initial signature: univariable", "Initial signature: multivariable",
      "LASSO 5-gene: univariable", "LASSO 5-gene: subtype-adjusted"),
    `Calibration slope (historical heading)` = c("1.861", "1.817", "0.124", "–", "0.154", "0.223"),
    `Brier at 1/2/3 years` = c("0.183/0.203/0.171", "–", "0.185/0.240/0.225", "–", "0.234/0.213/0.163", "–"),
    `Integrated Brier (0–3 years)` = c("0.168", "0.163", "0.186", "0.173", "0.189", "0.177"),
    `Baseline survival at 1/2/3 years (%) (historical heading)` = c("76.9/42.0/35.4", "–", "75.5/40.5/35.2", "–", "64.0/33.3/22.2", "–"), check.names = FALSE)
  save_csv(table4, "Historical_reported_Table_4.csv"); save_csv(table5, "Historical_reported_Table_5.csv")

  # 6. Editable PowerPoint, all figures in one deck.
  cat("\n6. Editable PowerPoint and Word\n")
  figures <- list(p_volcano, p_heat, p_gene_forest, p_tcga_km, p_tcga_risk, p_tcga_forest, p_ext_km, p_ext_risk, p_ext_forest)
  ppt <- officer::read_pptx()
  for (graph in figures) ppt <- add_plot(ppt, graph)
  ppt_file <- file.path(results_dir, "PDAC_legacy_all_figures_editable.pptx")
  print(ppt, target = ppt_file)

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "PDAC historical reproduction", style = "heading 1")
  doc <- officer::body_add_par(doc, paste("The GSE71729 score deliberately reproduces the historical vector-recycling error.",
    "Historical reported values and recomputed results are presented separately. This document is an audit record, not corrected validation."))
  doc <- word_table(doc, model_display(models), "Recomputed Cox model performance",
    "HRs for continuous scores are per one-unit increase. C-indices refer to the entire fitted model. TCGA adjusted models exclude patients with missing covariates.")
  doc <- word_table(doc, sc, "Historical five-gene coefficients",
    "Original TCGA coefficients reused for exact score reconstruction. Current glmnet refitting may differ.")
  doc <- officer::body_add_break(doc)
  et <- ea; et$`HR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", et$HR, et$CI_low, et$CI_high)
  et$`P value` <- format.pval(et$P_value, digits = 4, eps = 0.001)
  doc <- word_table(doc, et[, c("Term", "HR (95% CI)", "P value")], "Historical external subtype-adjusted Cox model",
    sprintf("Subtype terms use workbook codes 1/2 and 1/2/3. Historical median-split log-rank p = %.5f.", ext_p))
  doc <- officer::body_add_break(doc)
  doc <- word_table(doc, table4, "Table 4. Historical reported prognostic performance",
    "C-index confidence intervals are retained from the original reported table, not re-estimated. Multivariable TCGA models used 175 patients. CI: confidence interval; LASSO: least absolute shrinkage and selection operator.")
  doc <- officer::body_add_break(doc)
  doc <- word_table(doc, table5, "Table 5. Historical reported calibration, prediction accuracy and survival",
    "Original values and historical headings are retained. Score coefficients from cohort-refitted Cox models are not established external calibration slopes. Survival at the cohort mean score is distinct from baseline survival at score zero. The historical GSE71729 score is incorrect. GSE71729: Gene Expression Omnibus dataset; TCGA-PAAD: The Cancer Genome Atlas pancreatic adenocarcinoma.")
  doc <- officer::body_add_break(doc)
  doc <- word_table(doc, accuracy, "Recomputed apparent prediction accuracy of cohort-refitted models",
    "IPCW: inverse probability of censoring weighting. Monthly grid and trapezoidal integration over 0–3 years. TCGA uses 365 days/year; GSE71729 uses 12 months/year. These are cohort-refitted apparent prediction results, not validation of frozen TCGA survival probabilities. Values can differ slightly from historical Table 5.")
  word_file <- file.path(results_dir, "PDAC_legacy_tables_editable.docx")
  print(doc, target = word_file)

  # Preserve model objects, input hashes and software versions.
  saveRDS(list(tcga = tcga, initial = initial, legacy_external = legacy,
    tcga_uni = tcga_uni, tcga_adj = tcga_adj, initial_uni = initial_uni, initial_adj = initial_adj,
    external_uni = ext_uni, external_adj = ext_adj, external_group = ext_group,
    historical_beta = beta, original_external_beta = original_beta), file.path(results_dir, "PDAC_legacy_checkpoint.rds"))
  save_csv(data.frame(File = c(wb, tcga_file, geo_file), MD5 = unname(tools::md5sum(c(wb, tcga_file, geo_file)))), "Input_file_checksums.csv")
  assert(file.exists(ppt_file) && file.info(ppt_file)$size > 0, "PowerPoint was not saved.")
  assert(file.exists(word_file) && file.info(word_file)$size > 0, "Word document was not saved.")
  cat("\nCompleted. Slides:", length(ppt), "\nPowerPoint:", ppt_file, "\nWord:", word_file, "\n")
  cat("Session info:", file.path(results_dir, "R_session_info.txt"), "\n")
}

main()
