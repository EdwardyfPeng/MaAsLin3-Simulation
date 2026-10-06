args <- commandArgs(trailingOnly = TRUE)

p           <- as.numeric(args[1])
signal_type <- as.character(args[2])
u_depth     <- as.numeric(args[3])
target_fdr  <- as.numeric(args[4])
sim         <- as.numeric(args[5])

load("GALAXYMicrobLiver_study.RData")
library(tidyverse)
library(maaslin3)

generate_data_AA <- function(p, signal_type = c("prevalence", "abundance", "both"),
                             u = 0, seed = 1234) {
  set.seed(seed)
  
  n <- 500
  n_signal <- 30
  delta_zero <- 0.20 # Prevalence effect: control has 20 percentage points more structural zeros than case
  Delta <- 3  # Abundance effect
  signal_type <- match.arg(signal_type)

  # Prevalence of each taxon
  norm_count <- AA_real / rowSums(AA_real)
  col_means <- colMeans(norm_count > 0)
  indices <- which(col_means > 0.15)
  sorted_indices <- indices[order(col_means[indices], decreasing = TRUE)]
  dcount <- AA_real[,sorted_indices[seq_len(p)],drop = FALSE] # Take the top p taxa

  # Randomly sample 500 subjects
  sel_index <- sort(sample(1:nrow(dcount), n))
  dcount <- dcount[sel_index, ]

  Pi <- sweep(dcount,1,rowSums(dcount),"/") # Baseline composition, keep zeros

  Y <- sample(rep(c(0, 1), each = n / 2))
  control_idx <- which(Y == 0)
  case_idx <- which(Y == 1)

  n_add_zero <- round(delta_zero * length(control_idx))

  signal_pool <- 1:min(200, p)
  
  if (signal_type == "prevalence") {
    prevalence_idx <- sample(signal_pool, n_signal)
    abundance_idx <- integer(0)
    } else if (signal_type == "abundance") {
      prevalence_idx <- integer(0)
      abundance_idx <- sample(signal_pool, n_signal)
    } else {
      both_idx <- sample(signal_pool, n_signal)
      prevalence_idx <- both_idx
      abundance_idx <- both_idx
    }
  
  Pi_new <- Pi

  # Prevalence signals: add structural zeros to 20% of controls
  for (j in prevalence_idx) {
    positive_idx <- control_idx[Pi_new[control_idx, j] > 0]
    zero_idx <- sample(
      positive_idx,
      n_add_zero
    )
    Pi_new[zero_idx, j] <- 0
  }

  # Abundance signals: case abundance multiplied by 1 + U(0,3)
  if (length(abundance_idx) > 0) {
    fold_changes <- runif(
      length(abundance_idx),
      min = 0,
      max = Delta
    )
    for (j in seq_along(abundance_idx)) {
      col_j <- abundance_idx[j]
      Pi_new[case_idx, col_j] <-
        Pi_new[case_idx, col_j] * (1 + fold_changes[j])
    }
  }

  # True structural zeros before multinomial sampling
  structural_zero <- Pi_new == 0

  # Renormalize
  Pi_new <- Pi_new / rowSums(Pi_new)

  # Sequencing depth
  raw_depth <- rowSums(AA_real)

  target_median <- 50000 ## can change to 1000
  scale_factor <- target_median / median(raw_depth)
  drawn_depths <- round(raw_depth[sel_index] * scale_factor)
  # drawn_depths <- sample(drawn_depths, length(drawn_depths), replace = TRUE)
  adjusted_depths <- drawn_depths
  adjusted_depths[case_idx] <-  round(drawn_depths[case_idx] * (1 + u))

  # Multinomial sampling
  sim_count <- matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {
    sim_count[i, ] <- rmultinom(1,
      size = adjusted_depths[i],
      prob = Pi_new[i, ]
    )
  }

  colnames(sim_count) <- colnames(Pi)

  # Sampling zeros
  sampling_zero <- (sim_count == 0) & (Pi_new > 0)
  signal_indices <- sort(unique(c(prevalence_idx, abundance_idx)))

  return(list(
    Y = Y,
    X = sim_count,
    signal_indices = signal_indices,
    prevalence_indices = prevalence_idx,
    abundance_indices = abundance_idx,
    structural_zero = structural_zero,
    sampling_zero = sampling_zero,
    Pi_true = Pi_new
  ))
}

###### Main Simulation Function ######
run_single_iteration_AA <- function(p, signal_type, u_depth, target_fdr, sim) {
  cat("Generating AA data...\n")
  data_sim <- generate_data_AA(
    p = p,
    signal_type = signal_type,
    u = u_depth,
    seed = sim
  )

  X_sim <- data_sim$X
  Y_sim <- data_sim$Y
  true_signal <- data_sim$signal_indices

  rownames(X_sim) <- paste0("Sample_", seq_len(nrow(X_sim)))
  if (is.null(colnames(X_sim))) {
    colnames(X_sim) <- paste0("Taxon_", seq_len(ncol(X_sim)))
  }

  taxa_names <- colnames(X_sim)

  metadata <- data.frame(
    group = factor(
      Y_sim,
      levels = c(0, 1),
      labels = c("control", "case")
    ),
    log_reads = log(rowSums(X_sim)),
    row.names = rownames(X_sim)
  )

  feature_table_t <- as.data.frame(t(X_sim))

  names_to_idx <- function(nm) {
    idx <- match(nm, taxa_names)
    unique(idx[!is.na(idx)])
  }

  results_list <- list()
  make_row <- function(selected, method) {
    FP <- sum(!(selected %in% true_signal))
    TP <- sum(selected %in% true_signal)

    data.frame(
      p = p,
      signal_type = signal_type,
      u = u_depth,
      sim_id = sim,
      method = method,
      target_fdr = target_fdr,
      empirical_fdr = ifelse(length(selected) > 0, FP / (FP + TP), 0),
      power = TP / 30,
      n_selected = length(selected)
    )
  }

  make_fail <- function(method) {
    data.frame(
      p = p,
      signal_type = signal_type,
      u = u_depth,
      sim_id = sim,
      method = method,
      target_fdr = target_fdr,
      empirical_fdr = NA,
      power = NA,
      n_selected = NA
    )
  }

  run_maaslin_group <- function(evaluate_only, adjust_reads, prefix){
    out_dir <- tempfile(prefix)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    maaslin3::maaslin3(
      input_data = as.data.frame(X_sim),
      input_metadata = metadata,
      output = out_dir,
      formula = if (adjust_reads) ~ group + log_reads else ~ group,
      normalization = "TSS", transform = "LOG",
      min_abundance = 0, min_prevalence = 0, min_variance = 0,
      correction = "BH", standardize = FALSE,
      median_comparison_abundance = TRUE,
      median_comparison_prevalence = FALSE,
      evaluate_only = evaluate_only,
      warn_prevalence = FALSE,
      plot_summary_plot = FALSE, plot_associations = FALSE
    )
    res <- read.delim(file.path(out_dir, "all_results.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
    res[res$metadata == "group", , drop = FALSE]
  }

  # ---------- MaAsLin3 Joint ----------
  results_list[["MaAsLin3-Joint-read"]] <- tryCatch({
    cat("=== MaAsLin3 Joint ===\n")
    grp <- run_maaslin_group(NULL, TRUE, "maaslin3_joint_")
    grp <- grp[!is.na(grp$pval_joint), , drop = FALSE]
    grp <- grp[!duplicated(grp$feature), , drop = FALSE]  
    grp$q_group <- p.adjust(grp$pval_joint, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin3-Joint-read")
  }, error = function(e) {
    warning(sprintf("MaAsLin3-Joint failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin3-Joint-read")
  })

  # ---------- MaAsLin3 Abundance ----------
  results_list[["MaAsLin-abundance-read"]] <- tryCatch({
    cat("=== MaAsLin-abundance ===\n")
    grp <- run_maaslin_group("abundance", TRUE, "maaslin3_abun_")
    grp <- grp[grp$model == "abundance" & !is.na(grp$pval_individual), , drop = FALSE]
    grp$q_group <- p.adjust(grp$pval_individual, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin-abundance-read")
  }, error = function(e) {
    warning(sprintf("MaAsLin-abundance failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin-abundance-read")
  })

  # ---------- MaAsLin3 Prevalence ----------
  results_list[["MaAsLin3-Prevalence-read"]] <- tryCatch({
    cat("=== MaAsLin3 Prevalence ===\n")
    grp <- run_maaslin_group("prevalence", TRUE, "maaslin3_prev_")
    grp <- grp[grp$model == "prevalence" & !is.na(grp$pval_individual), , drop = FALSE]
    grp$q_group <- p.adjust(grp$pval_individual, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin3-Prevalence-read")
  }, error = function(e) {
    warning(sprintf("MaAsLin3-Prevalence failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin3-Prevalence-read")
  })

  # ---------- MaAsLin3 Joint ----------
  results_list[["MaAsLin3-Joint-noread"]] <- tryCatch({
    cat("=== MaAsLin3 Joint ===\n")
    grp <- run_maaslin_group(NULL, FALSE, "maaslin3_joint_")
    grp <- grp[!is.na(grp$pval_joint), , drop = FALSE]
    grp <- grp[!duplicated(grp$feature), , drop = FALSE]  
    grp$q_group <- p.adjust(grp$pval_joint, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin3-Joint-noread")
  }, error = function(e) {
    warning(sprintf("MaAsLin3-Joint failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin3-Joint-noread")
  })

  # ---------- MaAsLin3 Abundance ----------
  results_list[["MaAsLin-abundance-noread"]] <- tryCatch({
    cat("=== MaAsLin-abundance ===\n")
    grp <- run_maaslin_group("abundance", FALSE, "maaslin3_abun_")
    grp <- grp[grp$model == "abundance" & !is.na(grp$pval_individual), , drop = FALSE]
    grp$q_group <- p.adjust(grp$pval_individual, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin-abundance-noread")
  }, error = function(e) {
    warning(sprintf("MaAsLin-abundance failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin-abundance-noread")
  })

  # ---------- MaAsLin3 Prevalence ----------
  results_list[["MaAsLin3-Prevalence-noread"]] <- tryCatch({
    cat("=== MaAsLin3 Prevalence ===\n")
    grp <- run_maaslin_group("prevalence", FALSE, "maaslin3_prev_")
    grp <- grp[grp$model == "prevalence" & !is.na(grp$pval_individual), , drop = FALSE]
    grp$q_group <- p.adjust(grp$pval_individual, method = "BH")
    make_row(names_to_idx(grp$feature[grp$q_group <= target_fdr]), "MaAsLin3-Prevalence-noread")
  }, error = function(e) {
    warning(sprintf("MaAsLin3-Prevalence failed for sim %d: %s", sim, e$message))
    make_fail("MaAsLin3-Prevalence-noread")
  })
}

cat("Starting AA simulation...\n")
results <- run_single_iteration_AA(
  p, signal_type, u_depth, target_fdr, sim
)

output_file <- sprintf(
  "new_smgd_AA_results_p%d_%s_u%.2f_fdr%.2f_sim%d.RData",
  p, signal_type, u_depth, target_fdr, sim
)
save(results, file = output_file)
