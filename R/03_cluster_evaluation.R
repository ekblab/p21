## p21 / 03 | Cluster resolution evaluation
## Standalone: source("R/03_cluster_evaluation.R")
## Requires exported response_matrix_DEG_union_ashr.csv + module_membership.csv.
## Does not fit DESeq2, run GO, or regenerate biological classification.

if (!requireNamespace("here", quietly = TRUE)) stop("Missing package: here")
source(here::here("R", "03_shared.R"), local = TRUE)
p21_check(c("ekbSeq", "ggplot2", "cluster", "scales"))

## Evaluation of independent k=4..12 partitions, silhouettes, ARI and
## gene-overlap flow; k=8 reproducibility checked against exported membership.
report_dir <- p21_evaluation_dir
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
input <- file.path(p21_result_root, "response_matrix_DEG_union_ashr.csv")
membership_file <- file.path(p21_module_dir, "module_membership.csv")
if (!file.exists(input)) stop("Run script 03 once to generate: ", input)
if (!file.exists(membership_file)) stop("Baseline module_membership.csv is required for auditing k = 8")

## Exactly the same DEG-union response matrix that script 03 used for k-means.
response <- utils::read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
conditions <- c("HT1080_wt_7gy", "HT1080_p21_neg_7gy",
                "HT1080_wt_UVA", "HT1080_p21_neg_UVA",
                "HT1080_wt_UVB", "HT1080_p21_neg_UVB",
                "MCF10A_wt_7gy", "MCF10A_p21_neg_7gy",
                "MCF10A_wt_UVA", "MCF10A_p21_neg_UVA",
                "MCF10A_wt_UVB", "MCF10A_p21_neg_UVB")
if (!all(c("ENSEMBL", "SYMBOL", conditions) %in% names(response)))
  stop("The response matrix lacks required ENSEMBL/SYMBOL/12 condition columns")
if (anyNA(response$ENSEMBL) || anyDuplicated(response$ENSEMBL))
  stop("Invalid ENSEMBL IDs in the response matrix")
original <- as.matrix(response[, conditions, drop = FALSE]); storage.mode(original) <- "numeric"
rownames(original) <- response$ENSEMBL
if (any(!is.finite(original))) stop("All 12 response estimates must be finite")
magnitude <- sqrt(rowSums(original^2))
if (any(!is.finite(magnitude) | magnitude <= 0)) stop("Zero or invalid gene response profiles")
normalized <- sweep(original, 1, magnitude, "/")
stopifnot(ncol(normalized) == 12L,
          max(abs(rowSums(normalized^2) - 1)) < 1e-8)
message("Cluster evaluation: ", nrow(normalized), " genes × 12 normalized responses")

k_range <- 4:12
n_starts <- 40L             # identical initialization budget for all candidate k
n_repeat <- 8L             # subsampling stability replicates (80% genes each)
sample_fraction <- 0.80
subsample_starts <- 10L    # cheaper fits for perturbation replicates (same across k)
silhouette_n <- min(1500L, nrow(normalized))
silhouette_repeats <- 3L

## Label-invariant adjusted Rand index; no N×N pairwise matrix required.
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  pair <- function(x) x * (x - 1) / 2
  if (n < 2) return(NA_real_)
  row_pairs <- sum(pair(rowSums(tab)))
  col_pairs <- sum(pair(colSums(tab)))
  expected <- row_pairs * col_pairs / pair(n)
  denom <- (row_pairs + col_pairs) / 2 - expected
  if (denom <= 0) return(NA_real_)
  (sum(pair(tab)) - expected) / denom
}

## The same subsamples/distances for all k: fair, memory-bounded comparison.
set.seed(214)
sil_ids <- replicate(silhouette_repeats,
                     sample.int(nrow(normalized), silhouette_n), simplify = FALSE)
sil_dist <- lapply(sil_ids, function(ii) stats::dist(normalized[ii, , drop = FALSE]))

## Factor labels ordered by the mean *original* log2FC, as in the primary k=8 analysis.
relabel_clusters <- function(fit) {
  mean_effect <- vapply(seq_len(nrow(fit$centers)), function(j)
    mean(original[fit$cluster == j, , drop = FALSE]), numeric(1))
  rank <- order(mean_effect, decreasing = TRUE)
  sprintf("C%02d", match(fit$cluster, rank))
}

fits <- vector("list", length(k_range)); names(fits) <- paste0("k", k_range)
rows <- vector("list", length(k_range))
reviews <- list()
assignments <- list()
sil_details <- list()
ari_details <- list()

for (i in seq_along(k_range)) {
  k_i <- k_range[i]
  message(sprintf("Cluster evaluation %d/%d | k=%d", i, length(k_range), k_i))
  ## k=8 uses precisely the original random seed + nstart/iter.max settings.
  set.seed(if (k_i == 8L) 21L else 210L + k_i)
  fit <- stats::kmeans(normalized, centers = k_i, nstart = n_starts, iter.max = 150)
  fits[[i]] <- fit
  label <- relabel_clusters(fit)
  names(label) <- rownames(normalized)
  assignments[[as.character(k_i)]] <- label
  counts <- tabulate(fit$cluster, nbins = k_i)

  silhouette <- vapply(seq_along(sil_ids), function(j) {
    use <- sil_ids[[j]]
    groups <- fit$cluster[use]
    if (length(unique(groups)) < 2) return(NA_real_)
    mean(cluster::silhouette(groups, sil_dist[[j]])[, "sil_width"])
  }, numeric(1))
  if (any(!is.finite(silhouette))) stop("Nonfinite silhouette for k = ", k_i)
  sil_details[[i]] <- data.frame(k = k_i, repeat_id = seq_along(silhouette),
                                  mean_silhouette = silhouette)

  stability <- numeric(n_repeat)
  for (j in seq_len(n_repeat)) {
    set.seed(40000L + j)
    index <- sample.int(nrow(normalized), floor(sample_fraction * nrow(normalized)))
    set.seed(50000L + 100L * k_i + j)
    subset_fit <- stats::kmeans(normalized[index, , drop = FALSE], centers = k_i,
                                nstart = subsample_starts, iter.max = 150)
    stability[j] <- adjusted_rand(fit$cluster[index], subset_fit$cluster)
  }
  if (any(!is.finite(stability))) stop("Nonfinite ARI for k = ", k_i)
  ari_details[[i]] <- data.frame(k = k_i, repeat_id = seq_along(stability),
                                  adjusted_rand_index = stability)

  ## Descriptive profiles for MANUAL interpretation, not an arbitrary GO surrogate score.
  for (cl in sprintf("C%02d", seq_len(k_i))) {
    index <- which(label == cl)
    raw_mean <- colMeans(original[index, , drop = FALSE])
    centroid <- colMeans(normalized[index, , drop = FALSE])
    centroid <- centroid / sqrt(sum(centroid^2))
    cos <- as.vector(normalized[index, , drop = FALSE] %*% centroid)
    top <- order(-cos, -magnitude[index])[seq_len(min(3L, length(index)))]
    selected <- response$SYMBOL[index[top]]
    selected[is.na(selected) | !nzchar(selected)] <- response$ENSEMBL[index[top]][is.na(selected) | !nzchar(selected)]
    reviews[[length(reviews) + 1L]] <- data.frame(
      k = k_i, module = cl, n_genes = length(index),
      dominant_condition = names(which.max(abs(raw_mean))),
      dominant_mean_log2FC = raw_mean[which.max(abs(raw_mean))],
      profile_examples = paste(selected, collapse = "; "),
      biological_pattern_label = "", GO_BP_annotation = "",
      distinct_or_redundant = "", interpretation_notes = "",
      as.list(setNames(as.numeric(raw_mean), paste0("mean_", conditions))),
      check.names = FALSE)
  }
  rows[[i]] <- data.frame(k = k_i, n_genes = nrow(normalized),
    silhouette_mean = mean(silhouette), silhouette_sd = stats::sd(silhouette),
    stability_ARI_mean = mean(stability), stability_ARI_sd = stats::sd(stability),
    min_module_size = min(counts), median_module_size = stats::median(counts),
    max_module_size = max(counts), size_ratio_max_min = max(counts)/min(counts),
    kmeans_tot_withinss = fit$tot.withinss)
  utils::write.csv(data.frame(ENSEMBL = response$ENSEMBL, SYMBOL = response$SYMBOL,
                               module = unname(label)),
                   file.path(report_dir, sprintf("k%02d_gene_membership.csv", k_i)), row.names = FALSE)
}

## Explicit reproducibility check: running k=8 must recover the existing gene assignment.
ref <- utils::read.csv(membership_file, stringsAsFactors = FALSE)
if (!all(c("ENSEMBL", "module") %in% names(ref)) || anyDuplicated(ref$ENSEMBL) ||
    !setequal(ref$ENSEMBL, response$ENSEMBL)) stop("Incompatible k=8 baseline membership")
ref_cluster <- ref$module[match(response$ENSEMBL, ref$ENSEMBL)]
ref_fit <- fits[[which(k_range == 8L)]]
ref_ari <- adjusted_rand(ref_cluster, ref_fit$cluster)
if (!isTRUE(all.equal(ref_ari, 1, tolerance = 1e-10)))
  stop(sprintf("k=8 changed (ARI %.5f). Check input matrix, seed and previous analysis.", ref_ari))
message("Existing k=8 assignments recovered (ARI=1.0)")

metrics <- do.call(rbind, rows)
utils::write.csv(metrics, file.path(report_dir, "k_selection_metrics.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, sil_details), file.path(report_dir, "silhouette_subsamples.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, ari_details), file.path(report_dir, "stability_gene_subsampling.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, reviews), file.path(report_dir, "module_biological_review.csv"), row.names = FALSE)
utils::write.csv(data.frame(k = k_range,
    biological_interpretability_assessment = "",
    nonredundant_programs = NA_integer_,
    reviewer_decision = "",
    reviewer_notes = ""),
    file.path(report_dir, "k_biological_summary_template.csv"), row.names = FALSE, na = "")
utils::write.csv(data.frame(parameter = c("input", "gene_selection", "k_range", "normalization",
    "kmeans_nstart", "kmeans_itermax", "kmeans_seed_k8", "silhouette_sample_n",
    "silhouette_repeats", "stability_gene_fraction", "stability_repeats", "stability_nstart",
    "stability_seed_rule", "k8_vs_baseline_ARI", "transition_scan",
    "split_definition", "dominant_split_rule", "split_cosine_definition"),
    value = c(input, "DEG union, complete nonzero 12-condition profiles", "4:12",
              "gene-wise L2 norm (not variance scaling)", n_starts, 150, 21,
              silhouette_n, silhouette_repeats, sample_fraction, n_repeat,
              subsample_starts, "gene subset: 40000 + repeat; k-means: 50000 + 100*k + repeat", ref_ari,
              "all consecutive k from 4 to 12",
              "new group assigned to source with largest overlapping gene count",
              "among source groups producing >=2 dominant children, largest smaller branch gene count",
              "cosine of two new-cluster centroids computed from L2-normalized 12D gene profiles")),
    file.path(report_dir, "evaluation_metadata.csv"), row.names = FALSE)

npg <- p21_module_colors
pc <- stats::prcomp(normalized, center = TRUE, scale. = FALSE)
var_pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
gene_pc <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                      module = factor(ref_cluster, levels = names(npg)))
centroid <- t(vapply(names(npg), function(m)
  colMeans(normalized[ref_cluster == m, , drop = FALSE]), numeric(ncol(normalized))))
center_pc <- sweep(centroid, 2, pc$center, "-") %*% pc$rotation[, 1:2, drop = FALSE]
center_data <- data.frame(PC1 = center_pc[, 1], PC2 = center_pc[, 2],
                          module = factor(names(npg), levels = names(npg)))
stopifnot(all(is.finite(as.matrix(center_pc))))
point_alpha <- suppressWarnings(as.numeric(Sys.getenv("P21_PCA_ALPHA", "0.2")))
if (!is.finite(point_alpha) || point_alpha <= 0 || point_alpha > 1)
  stop("P21_PCA_ALPHA must be in (0, 1]")
p_pca <- ggplot2::ggplot(gene_pc, ggplot2::aes(PC1, PC2, color = module)) +
  ggplot2::geom_point(alpha = point_alpha, size = 1.1, stroke = 0) +
  ggplot2::geom_point(data = center_data, size = 4.2, shape = 21, fill = "white",
                      stroke = 1.25, show.legend = FALSE) +
  ggplot2::geom_label(data = center_data, ggplot2::aes(label = module),
                      fill = "white", color = "black", size = 3.4, fontface = "bold",
                      label.size = .20, label.padding = grid::unit(.12, "lines"),
                      show.legend = FALSE) +
  ggplot2::scale_color_manual(values = npg, drop = FALSE,
     guide = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 4.2))) +
  ggplot2::coord_equal() +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(title = "Radiation-response modules in PCA space",
    subtitle = "Each dot represents one gene",
    x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]),
    caption = "Clustering was performed in 12D using L2-normalized irradiation-response profiles; PCA is a two-dimensional visualization.")
p21_export(file.path(report_dir, "k8_gene_profile_PCA"), p_pca,
                          width = 6, height = 6)
utils::write.csv(data.frame(ENSEMBL = response$ENSEMBL, SYMBOL = response$SYMBOL,
                            module = ref_cluster, gene_pc[, 1:2]),
                 file.path(report_dir, "k8_gene_profile_PCA_coordinates.csv"), row.names = FALSE)

## Compact two-metric comparison; biological interpretability is reviewed separately.
graph_metrics <- rbind(data.frame(k = metrics$k, value = metrics$silhouette_mean,
                                  sd = metrics$silhouette_sd, measure = "Separation (silhouette)"),
                       data.frame(k = metrics$k, value = metrics$stability_ARI_mean,
                                  sd = metrics$stability_ARI_sd, measure = "Stability (ARI)"))
p_metrics <- ggplot2::ggplot(graph_metrics, ggplot2::aes(x = k, y = value)) +
  ggplot2::geom_vline(xintercept = 8, color = "grey75", linetype = "dashed") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = value - sd, ymax = value + sd),
                         width = .14, linewidth = .45, color = "grey40") +
  ggplot2::geom_line(linewidth = .65, color = "#3C5488") +
  ggplot2::geom_point(size = 2.6, color = "#3C5488") +
  ggplot2::facet_wrap(~measure, scales = "free_y", ncol = 2) +
  ggplot2::scale_x_continuous(breaks = k_range) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(title = "A | Choosing the number of response modules",
    subtitle = "Separation (silhouette) and stability (ARI) across k = 4–12",
    x = "Number of modules (k)", y = "Mean score (± SD)",
    caption = "ARI measures whether genes retain their cluster memberships across repeated 80% gene subsamples (1 = identical grouping). Dashed line: selected k = 8.")
p21_export(file.path(report_dir, "k_selection_metrics"), p_metrics,
                          width = 12.0, height = 5.3)

## Compare independently fitted solutions for every consecutive k=4..12.
## Algorithm is independent of GO results and module labels:
## (1) each new cluster inherits its dominant source (largest gene overlap);
## (2) when one source dominates >=2 new clusters, consider it a split;
## (3) choose the split with the largest *smaller* branch (absolute gene count).
## All candidate splits are retained; the report visualizes k=6->7 and k=7->8.
transitions <- list()
for (k_from in k_range[-length(k_range)]) {
  k_to <- k_from + 1L
  old <- assignments[[as.character(k_from)]]
  new <- assignments[[as.character(k_to)]]
  if (!identical(names(old), names(new))) stop("Gene order differs across k")
  overlap <- table(old, new)
  parent_of_child <- setNames(rownames(overlap)[max.col(t(overlap), ties.method = "first")],
                               colnames(overlap))
  source_counts <- table(parent_of_child)
  candidate_parents <- names(source_counts)[as.integer(source_counts) >= 2L]
  if (!length(candidate_parents)) {
    ## An incremental non-hierarchical fit need not produce a clean split.
    transitions[[as.character(k_from)]] <- NULL
    next
  }
  for (parent in candidate_parents) {
    children <- names(parent_of_child)[parent_of_child == parent]
    branch_counts <- sort(setNames(as.integer(overlap[parent, children]), children),
                          decreasing = TRUE)
    branch_a <- names(branch_counts)[1]; branch_b <- names(branch_counts)[2]
    n_parent <- sum(overlap[parent, ])
    n_a <- as.integer(branch_counts[1]); n_b <- as.integer(branch_counts[2])
    n_other <- n_parent - n_a - n_b
    ## True k-means feature space: mean of gene-wise L2-normalized profiles.
    z_a <- colMeans(normalized[new == branch_a, , drop = FALSE])
    z_b <- colMeans(normalized[new == branch_b, , drop = FALSE])
    centroid_cosine <- sum(z_a*z_b) / sqrt(sum(z_a*z_a) * sum(z_b*z_b))
    ## Secondary descriptive diagnostic in original ashr-log2FC space.
    a_raw <- colMeans(original[new == branch_a, , drop = FALSE])
    b_raw <- colMeans(original[new == branch_b, , drop = FALSE])
    raw_cosine <- sum(a_raw*b_raw) / sqrt(sum(a_raw*a_raw) * sum(b_raw*b_raw))
    transitions[[length(transitions) + 1L]] <- data.frame(
      k_from = k_from, k_to = k_to, parent = parent,
      n_parent = as.integer(n_parent), branch_a = branch_a, branch_b = branch_b,
      n_branch_a_from_parent = n_a, n_branch_b_from_parent = n_b,
      n_other = as.integer(n_other),
      branch_a_share_parent = n_a / n_parent,
      branch_b_share_parent = n_b / n_parent,
      two_branch_retention = (n_a + n_b) / n_parent,
      smaller_branch_n = min(n_a, n_b),
      centroid_cosine_L2 = centroid_cosine,
      mean_ashr_log2FC_cosine = raw_cosine,
      stringsAsFactors = FALSE)
  }
}
if (!length(transitions)) stop("No dominant splits identified for consecutive k solutions")
transition_all <- do.call(rbind, transitions)
transition_all <- transition_all[order(transition_all$k_from,
                                       -transition_all$smaller_branch_n,
                                       -transition_all$two_branch_retention), , drop = FALSE]
transition_all$dominant_split <- !duplicated(transition_all$k_from)
utils::write.csv(transition_all, file.path(report_dir, "k_transition_splits.csv"), row.names = FALSE)

targeted <- transition_all[transition_all$k_from %in% c(6L, 7L) &
                           transition_all$dominant_split, , drop = FALSE]
if (nrow(targeted) != 2L || !setequal(targeted$k_from, c(6L, 7L))) {
  stop("Cannot identify both k=6->7 and k=7->8 dominant splits. See k_transition_splits.csv")
}
targeted <- targeted[order(targeted$k_from), , drop = FALSE]
message("Dominant successive splits (L2 centroid cosine):")
for (i in seq_len(nrow(targeted))) {
  t <- targeted[i, ]
  message(sprintf("  k=%d->%d: %s -> %s (%d genes) / %s (%d genes); cosine=%.3f; retention=%.1f%%",
    t$k_from, t$k_to, t$parent, t$branch_a, t$n_branch_a_from_parent,
    t$branch_b, t$n_branch_b_from_parent, t$centroid_cosine_L2,
    100*t$two_branch_retention))
}
## Figure 11, full cluster tree: all partitions k=4,...,12.
## This is a gene-overlap flow graph (not hierarchical k-means).
## Each cluster node exists only within one independent k-means solution.
## Edges show the number of shared genes, including cross-connections;
## small flows are omitted from the plot only, not from the exported audit.
tree_nodes <- do.call(rbind, lapply(k_range, function(k_i) {
  lab <- assignments[[as.character(k_i)]]
  units <- sort(unique(lab))
  data.frame(k = k_i, module = units,
    n_genes = as.integer(table(factor(lab, levels = units))),
    mean_PC1 = vapply(units, function(id) mean(gene_pc$PC1[lab == id]), numeric(1)),
    stringsAsFactors = FALSE)
}))
## Use shared PCA coordinates to order clusters consistently within each k.
## Equally spaced rows prevent the labels from colliding at high k.
tree_nodes$y <- NA_real_
for (k_i in k_range) {
  ix <- which(tree_nodes$k == k_i)
  ix <- ix[order(tree_nodes$mean_PC1[ix])]
  tree_nodes$y[ix] <- seq(.04, .96, length.out = length(ix))
}
tree_nodes$node_key <- paste(tree_nodes$k, tree_nodes$module, sep = "_")
tree_nodes$label <- ifelse(tree_nodes$k == 8L,
                           sub("^C", "M", tree_nodes$module),
                           sub("^C0*", "", tree_nodes$module))
utils::write.csv(tree_nodes, file.path(report_dir, "k_clustree_nodes.csv"), row.names = FALSE)

tree_edges <- do.call(rbind, lapply(k_range[-length(k_range)], function(k_i) {
  old <- assignments[[as.character(k_i)]]
  new <- assignments[[as.character(k_i + 1L)]]
  counts <- as.data.frame(table(source = old, destination = new), stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  counts$k_from <- k_i
  counts$k_to <- k_i + 1L
  counts$source_n <- as.integer(table(old)[as.character(counts$source)])
  counts$target_n <- as.integer(table(new)[as.character(counts$destination)])
  counts$source_fraction <- counts$Freq / counts$source_n
  counts$target_fraction <- counts$Freq / counts$target_n
  counts$source_key <- paste(counts$k_from, counts$source, sep = "_")
  counts$target_key <- paste(counts$k_to, counts$destination, sep = "_")
  counts
}))
utils::write.csv(tree_edges, file.path(report_dir, "k_clustree_edges_all.csv"), row.names = FALSE)
tree_edges$x <- tree_edges$k_from
tree_edges$xend <- tree_edges$k_to
tree_edges$y <- tree_nodes$y[match(tree_edges$source_key, tree_nodes$node_key)]
tree_edges$yend <- tree_nodes$y[match(tree_edges$target_key, tree_nodes$node_key)]
if (anyNA(tree_edges$y) || anyNA(tree_edges$yend))
  stop("Failed to map full cluster-tree edges to nodes")
## Keep all counts in CSV; hide only <1% source-flow edges in the figure.
tree_edges <- tree_edges[tree_edges$source_fraction >= 0.01, , drop = FALSE]
tree_edges$highlight <- FALSE
for (i in seq_len(nrow(targeted))) {
  t <- targeted[i, ]
  ix <- tree_edges$k_from == t$k_from &
        tree_edges$source == t$parent &
        as.character(tree_edges$destination) %in% c(t$branch_a, t$branch_b)
  tree_edges$highlight[ix] <- TRUE
}
## Cosine similarity compares the two emerging CHILD centroids, not the
## parent-to-child edge. We place one label in the fork between the branches.
## All major splits are annotated so the k=6->7 / k=7->8 scores have context.
fork <- transition_all[transition_all$dominant_split, , drop = FALSE]
fork$x <- fork$k_from + 0.50
fork$child_a_key <- paste(fork$k_to, fork$branch_a, sep = "_")
fork$child_b_key <- paste(fork$k_to, fork$branch_b, sep = "_")
ya <- tree_nodes$y[match(fork$child_a_key, tree_nodes$node_key)]
yb <- tree_nodes$y[match(fork$child_b_key, tree_nodes$node_key)]
if (anyNA(ya) || anyNA(yb)) stop("Cannot locate split children in cluster tree")
fork$y <- (ya + yb) / 2
fork$label <- sprintf("%.2f", fork$centroid_cosine_L2)
fork$selected <- fork$k_from %in% c(6L, 7L)
utils::write.csv(fork, file.path(report_dir, "k_clustree_dominant_split_labels.csv"), row.names = FALSE)
# Draw unhighlighted flows behind prominent selected split paths.
bg <- tree_edges[!tree_edges$highlight, , drop = FALSE]
hi <- tree_edges[tree_edges$highlight, , drop = FALSE]
tree_nodes$node_color <- ifelse(tree_nodes$k == 8L,
   unname(npg[sub("^C", "M", tree_nodes$module)]), "#D8DDE3")
plot_tree <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = bg,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                 linewidth = Freq, alpha = source_fraction),
    color = "#9BA5B1", lineend = "round", show.legend = FALSE) +
  ggplot2::geom_segment(data = hi,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                 linewidth = Freq),
    color = "#3C5488", lineend = "round", alpha = .9, show.legend = FALSE) +
  ggplot2::geom_label(data = fork[!fork$selected, , drop = FALSE],
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE, label.size = .12, label.padding = grid::unit(.11, "lines"),
    fill = "white", color = "#677483", alpha = .97, size = 2.4,
    show.legend = FALSE) +
  ggplot2::geom_label(data = fork[fork$selected, , drop = FALSE],
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE, label.size = .42, label.padding = grid::unit(.16, "lines"),
    fill = "#F1F5FA", color = "#244879", fontface = "bold", size = 3,
    show.legend = FALSE) +
  ggplot2::geom_point(data = tree_nodes,
    ggplot2::aes(x = k, y = y, size = n_genes),
    shape = 21, fill = tree_nodes$node_color, color = "#56616F",
    stroke = .4, show.legend = FALSE) +
  ggplot2::geom_text(data = tree_nodes,
    ggplot2::aes(x = k, y = y, label = label),
    size = 2.25, fontface = "bold", color = "#14202C", show.legend = FALSE) +
  ggplot2::scale_linewidth_continuous(range = c(.10, 3.0), guide = "none") +
  ggplot2::scale_alpha_continuous(range = c(.10, .72), guide = "none") +
  ggplot2::scale_size_area(max_size = 11, guide = "none") +
  ggplot2::scale_x_continuous(breaks = k_range,
     labels = paste0("k = ", k_range), expand = ggplot2::expansion(add = .5)) +
  ggplot2::coord_cartesian(ylim = c(-.02, 1.03), clip = "off") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "#EEF0F2", linewidth = .35),
    axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(title = "B | How gene modules resolve as k increases",
     subtitle = paste0("All independent k-means solutions (k = 4–12). Node area: number of genes; ",
        "edge width: shared gene count; labels: cosine between split-off child profiles."),
     caption = paste0("Cosine labels describe the two major child profiles at each branching event, not individual edges. ",
        "Selected splits at k = 6→7 and 7→8 are highlighted. ",
        "Node IDs are specific to each k; independent k-means fits are not hierarchical. ",
        "Flows <1% of the source cluster are omitted visually but retained in CSV."))
p21_export(file.path(report_dir, "k_clustree_all"), plot_tree,
                          width = 17, height = 10.5)
message("Complete k=4–12 gene-overlap tree exported to PNG/SVG: k_clustree_all")

message("Evaluation complete; existing modules unchanged. Review: ", report_dir)
