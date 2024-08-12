library(infercnv)
library(tidyverse)
library(data.table)
library(parallelDist)
library(futile.logger)


annotations_legend = as.matrix(cbind(names(alphabet),alphabet))
colnames(annotations_legend) = c("name_ref_groups", "")

plot_cnv2 <- function (infercnv_obj, out_dir = ".", title = "inferCNV", obs_title = "Observations (Cells)",  ref_title = "References (Cells)", cluster_by_groups = TRUE,   cluster_references = TRUE, plot_chr_scale = FALSE, chr_lengths = NULL, k_obs_groups = 1, contig_cex = 1, x.center = mean(infercnv_obj@expr.data), x.range = "auto", hclust_method = "ward.D", custom_color_pal = NULL,  color_safe_pal = FALSE, output_filename = "infercnv", output_format = "png", png_res = 300, dynamic_resize = 0, ref_contig = NULL, write_expr_matrix = FALSE, write_phylo = FALSE, useRaster = TRUE) {
  
  if (!file.exists(out_dir)) {
    dir.create(out_dir)
  }
  plot_data = infercnv_obj@expr.data
  flog.info(paste("::plot_cnv:Start", sep = ""))
  flog.info(paste("::plot_cnv:Current data dimensions (r,c)=", 
                  paste(dim(plot_data), collapse = ","), " Total=", sum(plot_data, 
                                                                        na.rm = TRUE), " Min=", min(plot_data, na.rm = TRUE), 
                  " Max=", max(plot_data, na.rm = TRUE), ".", sep = ""))
  flog.info(paste("::plot_cnv:Depending on the size of the matrix", 
                  " this may take a moment.", sep = ""))
  if (write_expr_matrix) {
    expr_dat_file <- paste(out_dir, paste("expr.", output_filename, 
                                          ".dat", sep = ""), sep = "/")
    if ("matrix" %in% is(plot_data)) {
      write.table(as.matrix(plot_data), file = expr_dat_file, 
                  quote = FALSE, sep = "\t")
    }
  }
  if (!any(is.na(x.range))) {
    if ((length(x.range) == 1) & (x.range[1] == "auto")) {
      quantiles = quantile(plot_data[plot_data != x.center], 
                           c(0.01, 0.99))
      delta = max(abs(c(x.center - quantiles[1], quantiles[2] - 
                          x.center)))
      low_threshold = x.center - delta
      high_threshold = x.center + delta
      x.range = c(low_threshold, high_threshold)
      flog.info(sprintf("plot_cnv(): auto thresholding at: (%f , %f)", 
                        low_threshold, high_threshold))
    }
    else {
      low_threshold = x.range[1]
      high_threshold = x.range[2]
      if (low_threshold > x.center | high_threshold < x.center | 
          low_threshold >= high_threshold) {
        stop(paste("Error, problem with relative values of x.range: ", 
                   x.range, ", and x.center: ", x.center))
      }
    }
    plot_data[plot_data < low_threshold] <- low_threshold
    plot_data[plot_data > high_threshold] <- high_threshold
    infercnv_obj@expr.data <- plot_data
  }
  contigs = infercnv_obj@gene_order[['chr']]
  unique_contigs <- unique(contigs)
  n_contig <- length(unique_contigs)
  ct.colors <- infercnv:::get_group_color_palette()(n_contig)
  names(ct.colors) <- unique_contigs
  if (!is.null(custom_color_pal)) {
    custom_pal = custom_color_pal
  }
  else if (color_safe_pal == FALSE) {
    custom_pal <- color.palette(c("darkblue", "white", "darkred"), 
                                c(2, 2))
  }
  else {
    custom_pal <- color.palette(c("purple3", "white", "darkorange2"), 
                                c(2, 2))
  }
  ref_idx <- NULL
  if (infercnv:::has_reference_cells(infercnv_obj)) {
    ref_idx <- unlist(infercnv_obj@reference_grouped_cell_indices)
    ref_idx = ref_idx[order(ref_idx)]
  }
  contig_tbl <- table(contigs)[unique_contigs]
  col_sep <- cumsum(contig_tbl)
  col_sep <- col_sep[-1 * length(col_sep)]
  contig_labels = names(contig_tbl)
  contig_names = unlist(lapply(contig_labels, function(contig_name) {
    rep(contig_name, contig_tbl[contig_name])
  }))
  grouping_key_coln <- c()
  obs_annotations_names <- names(infercnv_obj@observation_grouped_cell_indices)
  obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data)))
  names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
  obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
  counter <- 1
  for (obs_index_group in obs_index_groupings) {
    obs_annotations_groups[obs_index_group] <- counter
    counter <- counter + 1
  }
  if (!is.null(ref_idx)) {
    obs_annotations_groups <- obs_annotations_groups[-ref_idx]
  }
  if (is.null(dynamic_resize) | dynamic_resize < 0) {
    flog.warn(paste("invalid dynamic_resize value: ", dynamic_resize, 
                    sep = ""))
    dynamic_resize = 0
  }
  dynamic_extension = 0
  nobs = length(unlist(infercnv_obj@observation_grouped_cell_indices))
  if (nobs > 200) {
    dynamic_extension = dynamic_resize * 3.6 * (nobs - 200)/200
  }
  grouping_key_coln[1] <- floor(123/(max(nchar(obs_annotations_names)) + 
                                       6))
  if (grouping_key_coln[1] < 1) {
    grouping_key_coln[1] <- 1
  }
  name_ref_groups = names(infercnv_obj@reference_grouped_cell_indices)
  if (is.null(name_ref_groups)) {
    grouping_key_coln[2] = 1
  }
  else {
    grouping_key_coln[2] <- floor(123/(max(nchar(name_ref_groups)) + 
                                         6))
    if (grouping_key_coln[2] < 1) {
      grouping_key_coln[2] <- 1
    }
  }
  grouping_key_rown <- c()
  grouping_key_rown[1] <- ceiling(length(obs_annotations_names)/grouping_key_coln[1])
  grouping_key_rown[2] <- ceiling(length(name_ref_groups)/grouping_key_coln[2])
  grouping_key_height <- c((grouping_key_rown[2] + 2) * 0.175, 
                           (grouping_key_rown[1] + 3) * 0.175)
  if (!is.na(output_format)) {
    if (output_format == "pdf") {
      pdf(paste(out_dir, paste(output_filename, ".pdf", sep = ""), sep = "/"), useDingbats = FALSE, width = 15, height = (12.22 + sum(grouping_key_height)) + dynamic_extension, paper = "special")
      print((12.22 + sum(grouping_key_height)) + dynamic_extension) # 13.445
    }
    else if (output_format == "png") {
      png_height = 8.22 + sum(grouping_key_height) + dynamic_extension
      if (!is.null(getOption("bitmapType")) && (getOption("bitmapType") == 
                                                "cairo") & (png_height > 32768/png_res)) {
        png_height = round((32767/png_res) - 5 * 10^(-3), 
                           2)
        flog.warn(paste0("Requested PNG output height too big at the current resolution, ", 
                         "using the max height instead. (cairo seems to have a size limit of 32767 (2^15-1) pixels ", 
                         "per dimension and 49151 (2^15+2^14-1)pixels for the sum of dimensions)"))
      }
      png(paste(out_dir, paste(output_filename, ".png", 
                               sep = ""), sep = "/"), width = 10, height = png_height, 
          units = "in", res = png_res)
    }
  }
  obs_data <- infercnv_obj@expr.data
  if (!is.null(ref_idx)) {
    obs_data <- plot_data[, -ref_idx, drop = FALSE]
    if (ncol(obs_data) == 1) {
      plot_data <- cbind(obs_data, obs_data)
      names(obs_data) <- c("", names(obs_data)[1])
    }
  }
  obs_data <- t(obs_data)
  ref_data_t <- NULL
  updated_ref_groups <- list()
  current_ref_count <- 1
  current_grp_idx <- 1
  plot_data <- infercnv_obj@expr.data
  ref_groups = infercnv_obj@reference_grouped_cell_indices
  for (ref_grp in ref_groups) {
    ref_data_t <- cbind(ref_data_t, plot_data[, ref_grp, 
                                              drop = FALSE])
    updated_ref_groups[[current_grp_idx]] = seq(current_ref_count, 
                                                current_ref_count + length(ref_grp) - 1)
    current_ref_count <- current_ref_count + length(ref_grp)
    current_grp_idx <- current_grp_idx + 1
  }
  ref_groups <- updated_ref_groups
  nb_breaks <- 16
  breaksList_t <- seq(x.range[1], x.range[2], length.out = nb_breaks)
  gene_position_breaks = NULL
  if (plot_chr_scale) {
    chr_name_list = unique(infercnv_obj@gene_order[["chr"]])
    if (is.null(chr_lengths)) {
      chr_lengths = c()
      for (chr_name in chr_name_list) {
        chr_lengths = c(chr_lengths, max(infercnv_obj@gene_order$stop[which(infercnv_obj@gene_order$chr == 
                                                                              chr_name)]) + 10000)
      }
      names(chr_lengths) = chr_name_list
    }
    gene_position_breaks = vector(mode = "integer", length = (length(unlist(infercnv_obj@gene_order$chr)) + 
                                                                1))
    sum_previous_contigs = 0
    gene_position_breaks[1] = 1
    current_idx = 2
    col_sep_idx = 1
    for (chr_name in chr_name_list) {
      index_pos = which(infercnv_obj@gene_order$chr == 
                          chr_name)
      latest_position = 1
      if (length(index_pos) > 1) {
        for (i in index_pos[2:length(index_pos)]) {
          gene_position_breaks[current_idx] = sum_previous_contigs + 
            ((latest_position + infercnv_obj@gene_order$start[i])/2)
          latest_position = max(infercnv_obj@gene_order$stop[i], 
                                latest_position)
          current_idx = current_idx + 1
        }
      }
      gene_position_breaks[current_idx] = sum_previous_contigs + 
        chr_lengths[chr_name]
      current_idx = current_idx + 1
      sum_previous_contigs = sum_previous_contigs + chr_lengths[chr_name]
      if (col_sep_idx != length(chr_name_list)) {
        col_sep[col_sep_idx] = sum_previous_contigs
        col_sep_idx = col_sep_idx + 1
      }
    }
    overlap_pos = which(gene_position_breaks[1:(length(gene_position_breaks) - 
                                                  1)] == gene_position_breaks[2:length(gene_position_breaks)])
    if (any(overlap_pos)) {
      gene_position_breaks[overlap_pos + 1] = gene_position_breaks[overlap_pos] + 
        1
    }
  }
  force_layout <- infercnv:::.plot_observations_layout(grouping_key_height = grouping_key_height, 
                                                       dynamic_extension = dynamic_extension)
  plot_cnv_observations2(infercnv_obj = infercnv_obj, obs_data = obs_data, 
                         file_base_name = out_dir, do_plot = !is.na(output_format), 
                         write_expr_matrix = write_expr_matrix, write_phylo = write_phylo, 
                         output_filename_prefix = output_filename, cluster_contig = ref_contig, 
                         contigs = contigs, contig_colors = ct.colors[contigs], 
                         contig_labels = contig_labels, contig_names = contig_names, 
                         col_pal = custom_pal, contig_seps = col_sep, num_obs_groups = k_obs_groups, 
                         obs_annotations_groups = obs_annotations_groups, obs_annotations_names = obs_annotations_names, 
                         grouping_key_coln = grouping_key_coln[1], cluster_by_groups = cluster_by_groups, 
                         cnv_title = title, cnv_obs_title = obs_title, contig_lab_size = contig_cex, 
                         breaksList = breaksList_t, gene_position_breaks = gene_position_breaks, 
                         x.center = x.center, hclust_method = hclust_method, layout_lmat = force_layout[["lmat"]], 
                         layout_lhei = force_layout[["lhei"]], layout_lwid = force_layout[["lwid"]], 
                         useRaster = useRaster)
  obs_data <- NULL
  if (!is.null(ref_idx)) {
    plot_cnv_references2(infercnv_obj = infercnv_obj, ref_data = ref_data_t, 
                         ref_groups = ref_groups, name_ref_groups = name_ref_groups, 
                         cluster_references = cluster_references, hclust_method = hclust_method, 
                         grouping_key_coln = grouping_key_coln[2], col_pal = custom_pal, 
                         contig_seps = col_sep, file_base_name = out_dir, 
                         do_plot = !is.na(output_format), write_expr_matrix = write_expr_matrix, 
                         output_filename_prefix = output_filename, cnv_ref_title = ref_title, 
                         breaksList = breaksList_t, gene_position_breaks = gene_position_breaks, 
                         x.center = x.center, layout_add = TRUE, useRaster = useRaster)
  }
  if (!is.na(output_format)) {
    dev.off()
  }
  return(list(cluster_by_groups = cluster_by_groups, k_obs_groups = k_obs_groups, 
              contig_cex = contig_cex, x.center = x.center, x.range = x.range, 
              hclust_method = hclust_method, color_safe_pal = color_safe_pal, 
              output_format = output_format, png_res = png_res, dynamic_resize = dynamic_resize))
}


plot_cnv_observations2 <- function (infercnv_obj, obs_data, col_pal, contig_colors, contig_labels, contig_names, contig_seps, num_obs_groups, file_base_name, do_plot = TRUE, write_expr_matrix, write_phylo, output_filename_prefix, cnv_title, cnv_obs_title, contig_lab_size = 1, contigs, cluster_contig = NULL, obs_annotations_groups, obs_annotations_names, grouping_key_coln, cluster_by_groups, breaksList, gene_position_breaks, x.center, hclust_method = "ward.D", testing = FALSE, layout_lmat = NULL, layout_lhei = NULL, layout_lwid = NULL, useRaster = useRaster) {
  flog.info("plot_cnv_observation:Start")
  flog.info(paste("Observation data size: Cells=", nrow(obs_data), 
                  "Genes=", ncol(obs_data), sep = " "))
  observation_file_base <- paste(file_base_name, sprintf("%s.observations.txt", 
                                                         output_filename_prefix), sep = .Platform$file.sep)
  dendrogram_file_path = paste(file_base_name, sprintf("%s.observations_dendrogram.txt", 
                                                       output_filename_prefix), sep = .Platform$file.sep)
  hcl_desc <- "General"
  hcl_group_indices <- seq_len(ncol(obs_data))
  if (!is.null(cluster_contig)) {
    hcl_contig_indices <- which(contig_names %in% cluster_contig)
    if (length(hcl_contig_indices) > 0) {
      hcl_group_indices <- hcl_contig_indices
      hcl_desc <- paste(cluster_contig, collapse = "_")
      flog.info(paste("plot_cnv_observation:Clustering only by contig ", 
                      cluster_contig))
      infercnv_obj@tumor_subclusters = NULL
    }
    else {
      flog.warn(paste("plot_cnv_observations: Not able to cluster by", 
                      cluster_contig, "Clustering by all genomic locations.", 
                      "To cluster by local genomic location next time", 
                      "select from:", unique(contig_names), collapse = ",", 
                      sep = " "))
    }
  }
  obs_dendrogram <- list()
  ordered_names <- NULL
  isfirst <- TRUE
  hcl_obs_annotations_groups <- vector()
  obs_seps <- c()
  sub_obs_seps <- c()
  if (!is.null(infercnv_obj@tumor_subclusters)) {
    if (cluster_by_groups) {
      split_groups = vector()
      for (i in seq_along(obs_annotations_names)) {
        if (!is.null(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])) {
          obs_dendrogram[[i]] = as.dendrogram(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])
          ordered_names <- c(ordered_names, infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$labels[infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$order])
          obs_seps <- c(obs_seps, length(ordered_names))
          hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                          rep(i, length(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$order)))
          if (write_phylo) {
            if (isfirst) {
              write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]), 
                         file = dendrogram_file_path)
              isfirst <- FALSE
            }
            else {
              write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]), 
                         file = dendrogram_file_path, append = TRUE)
            }
          }
        }
        else {
          if ((length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
               2) || (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
                      1)) {
            if (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
                2) {
              obs_dendrogram[[i]] <- .pairwise_dendrogram(colnames(infercnv_obj@expr.data[, 
                                                                                          unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                                          drop = FALSE]))
            }
            else {
              obs_dendrogram[[i]] <- .single_element_dendrogram(colnames(infercnv_obj@expr.data[, 
                                                                                                unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                                                drop = FALSE]))
            }
            ordered_names <- c(ordered_names, colnames(infercnv_obj@expr.data[, 
                                                                              unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                              drop = FALSE]))
            obs_seps <- c(obs_seps, length(ordered_names))
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                            rep(i, length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]))))
          }
          else {
            flog.error("Unexpected error, should not happen.")
            stop("Error")
          }
        }
        for (subcluster in names(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) {
          tmp = infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]][[subcluster]]
          tmp[] = subcluster
          split_groups = c(split_groups, tmp)
        }
      }
      split_groups <- split_groups[ordered_names]
      if (length(obs_dendrogram) > 1) {
        obs_dendrogram <- do.call(merge, obs_dendrogram)
      }
      else {
        obs_dendrogram <- obs_dendrogram[[1]]
      }
      for (subtumor in infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]) {
        sub_obs_seps <- c(sub_obs_seps, (sub_obs_seps[length(sub_obs_seps)] + 
                                           length(subtumor)))
      }
    }
    else {
      obs_hcl <- infercnv_obj@tumor_subclusters$hc[["all_observations"]]
      if (write_phylo) {
        write.tree(as.phylo(obs_hcl), file = dendrogram_file_path)
      }
      obs_dendrogram <- as.dendrogram(obs_hcl)
      ordered_names <- obs_hcl$labels[obs_hcl$order]
      if (num_obs_groups > 1) {
        split_groups <- cutree(obs_hcl, k = num_obs_groups)
      }
      else {
        split_groups = vector()
        for (subcluster in names(infercnv_obj@tumor_subclusters$subclusters[["all_observations"]])) {
          tmp = infercnv_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster]]
          tmp[] = subcluster
          split_groups = c(split_groups, tmp)
        }
      }
      split_groups <- split_groups[ordered_names]
      hcl_obs_annotations_groups <- obs_annotations_groups[ordered_names]
      flog.info("plot_cnv_observation:Writing observations by grouping.")
      for (cut_group in unique(split_groups)) {
        group_memb <- names(split_groups)[which(split_groups == 
                                                  cut_group)]
        memb_file <- file(paste(file_base_name, paste(hcl_desc, 
                                                      "HCL", cut_group, "members.txt", sep = "_"), 
                                sep = .Platform$file.sep))
        write.table(as.matrix(obs_data[group_memb, ]), 
                    memb_file)
        ordered_memb <- which(ordered_names %in% group_memb)
        if (is.null(obs_seps)) {
          obs_seps <- c(length(ordered_memb))
        }
        else {
          obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + 
                                     length(ordered_memb)))
        }
      }
      obs_seps <- c(obs_seps, length(ordered_names))
      sub_obs_seps = obs_seps
    }
  }
  else if (cluster_by_groups) {
    flog.info(paste("clustering observations via method: ", 
                    hclust_method, sep = ""))
    for (i in seq_len(max(obs_annotations_groups))) {
      cell_indices_in_group <- which(obs_annotations_groups == 
                                       i)
      num_cells_in_group <- length(cell_indices_in_group)
      flog.info(sprintf("Number of cells in group(%d) is %d", 
                        i, num_cells_in_group))
      if (num_cells_in_group < 2) {
        flog.info(sprintf("Skipping group: %d, since less than 2 entries", 
                          i))
        ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups == 
                                                                     i), , drop = FALSE]))
        obs_dendrogram[[length(obs_dendrogram) + 1]] = .single_element_dendrogram(unique_label = row.names(obs_data[which(obs_annotations_groups == 
                                                                                                                            i), , drop = FALSE]))
        if (write_phylo) {
          if (isfirst) {
            write(row.names(obs_data[which(obs_annotations_groups == 
                                             i), ]), file = dendrogram_file_path)
            isfirst <- FALSE
          }
          else {
            write(row.names(obs_data[which(obs_annotations_groups == 
                                             i), ]), file = dendrogram_file_path, append = TRUE)
          }
        }
      }
      else {
        data_to_cluster <- obs_data[cell_indices_in_group, 
                                    hcl_group_indices, drop = FALSE]
        flog.info(paste("group size being clustered: ", 
                        paste(dim(data_to_cluster), collapse = ","), 
                        sep = " "))
        group_obs_hcl <- hclust(parallelDist(data_to_cluster, 
                                             threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), 
                                method = hclust_method)
        ordered_names <- c(ordered_names, group_obs_hcl$labels[group_obs_hcl$order])
        group_obs_dend <- as.dendrogram(group_obs_hcl)
        obs_dendrogram[[length(obs_dendrogram) + 1]] <- group_obs_dend
        if (write_phylo) {
          if (isfirst) {
            write.tree(as.phylo(group_obs_hcl), file = dendrogram_file_path)
            isfirst <- FALSE
          }
          else {
            write.tree(as.phylo(group_obs_hcl), file = dendrogram_file_path, 
                       append = TRUE)
          }
        }
      }
      hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                      rep(i, num_cells_in_group))
      obs_seps <- c(obs_seps, length(ordered_names))
    }
    if (length(obs_dendrogram) > 1) {
      obs_dendrogram <- do.call(merge, obs_dendrogram)
    }
    else {
      obs_dendrogram <- obs_dendrogram[[1]]
    }
    split_groups <- rep(1, dim(obs_data)[1])
    names(split_groups) <- ordered_names
    sub_obs_seps = obs_seps
  }
  else {
    flog.info(paste("clustering observations via method: ", 
                    hclust_method, sep = ""))
    if (nrow(obs_data) > 1) {
      obs_hcl <- hclust(parallelDist(obs_data[, hcl_group_indices], 
                                     threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), method = hclust_method)
      if (write_phylo) {
        write.tree(as.phylo(obs_hcl), file = dendrogram_file_path)
      }
      obs_dendrogram <- as.dendrogram(obs_hcl)
      ordered_names <- obs_hcl$labels[obs_hcl$order]
      split_groups <- cutree(obs_hcl, k = num_obs_groups)
      split_groups <- split_groups[ordered_names]
      flog.info("plot_cnv_observation:Writing observations by grouping.")
      for (cut_group in unique(split_groups)) {
        group_memb <- names(split_groups)[which(split_groups == 
                                                  cut_group)]
        memb_file <- file(paste(file_base_name, paste(hcl_desc, 
                                                      "HCL", cut_group, "members.txt", sep = "_"), 
                                sep = .Platform$file.sep))
        write.table(as.matrix(obs_data[group_memb, ]), 
                    memb_file)
        ordered_memb <- which(ordered_names %in% group_memb)
        if (is.null(obs_seps)) {
          obs_seps <- c(length(ordered_memb))
        }
        else {
          obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + 
                                     length(ordered_memb)))
        }
      }
      obs_seps <- c(obs_seps, length(ordered_names))
    }
    else {
      obs_dendrogram <- .single_element_dendrogram(row.names(obs_data))
      ordered_names <- row.names(obs_data)
      obs_seps <- c(1)
      split_groups <- c(1)
      names(split_groups) <- ordered_names
    }
    sub_obs_seps = obs_seps
    hcl_obs_annotations_groups <- obs_annotations_groups[ordered_names]
  }
  if (length(obs_seps) > 1) {
    obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) - 
                                                         1):1]
  }
  idx = 1
  split_groups_as_idx = split_groups
  for (grp in unique(split_groups)) {
    to_replace = which(split_groups == grp)
    split_groups_as_idx[to_replace] = idx
    idx = idx + 1
  }
  split_groups_as_idx = as.integer(split_groups_as_idx)
  row_groupings <- infercnv:::get_group_color_palette()(length(table(split_groups_as_idx)))[split_groups_as_idx]
  #row_groupings <- cbind(row_groupings, infercnv:::get_group_color_palette()(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])
  alphabet =  c('1' = '#FC6FCF', '2' = '#BC589B','3' = '#EB545C', '4' = '#FFEE00', '5' = '#00475F')
  row_groupings <- cbind(row_groupings, alphabet[hcl_obs_annotations_groups])

  alphabet =  c('0' = '#8497B0', '1' = '#878787', '2' = '#F0CE58', '3' = '#EB545C','4' = '#D7EF9B', '5' = '#DCEAF7', '6' = '#FFEE00','7' = '#00475F', '8' = '#F297A7','9' = '#ABDDDE','10' = 'lightskyblue','11' = '#B487B7','12' = '#1D65A6','13' = '#F2A104','14' = '#FC6FCF','15' = '#0FFFFF', '16' = '#289E92','17' = '#C3FF00','18' = '#BC589B','19' = '#8ED973')
  annotations_legend = as.matrix(cbind(names(alphabet),alphabet))
  #annotations_legend <- cbind(obs_annotations_names, infercnv:::get_group_color_palette()(length(table(hcl_obs_annotations_groups))))
  colnames(annotations_legend) = c("name_ref_groups", "")
  flog.info("plot_cnv_observation:Writing observation groupings/color.")
  groups_file_name <- file.path(file_base_name, sprintf("%s.observation_groupings.txt", 
                                                        output_filename_prefix))
  file_groups <- cbind(split_groups, row_groupings[, 1], hcl_obs_annotations_groups, 
                       row_groupings[, 2])
  colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color", 
                             "Annotation Group", "Annotation Color")
  write.table(file_groups, groups_file_name)
  flog.info("plot_cnv_observation:Done writing observation groupings/color.")
  contigSepList <- infercnv:::create_sep_list(row_count = nrow(obs_data), 
                                              col_count = ncol(obs_data), row_seps = obs_seps, col_seps = contig_seps)
  obs_data <- obs_data[ordered_names, ]
  orig_row_names <- row.names(obs_data)
  row.names(obs_data) <- rep("", nrow(obs_data))
  flog.info("plot_cnv_observation:Writing observation heatmap thresholds.")
  heatmap_thresholds_file_name <- file.path(file_base_name, 
                                            sprintf("%s.heatmap_thresholds.txt", output_filename_prefix))
  write.table(breaksList, heatmap_thresholds_file_name, row.names = FALSE, 
              col.names = FALSE)
  flog.info("plot_cnv_observation:Done writing observation heatmap thresholds.")
  if (do_plot) {
    data_observations <- infercnv:::heatmap.cnv(obs_data, Rowv = obs_dendrogram, 
                                                Colv = FALSE, cluster.by.row = TRUE, cluster.by.col = FALSE, 
                                                main = cnv_title, ylab = cnv_obs_title, margin.for.labCol = 2, 
                                                xlab = "Genomic Region", key = TRUE, labCol = contig_labels, 
                                                cexCol = contig_lab_size, cexAt = c(1, contig_seps), 
                                                notecol = "black", density.info = "histogram", denscol = "blue", 
                                                trace = "none", dendrogram = "row", cexRow = 0.8, 
                                                breaks = breaksList, gene_position_breaks = gene_position_breaks, 
                                                scale = "none", x.center = x.center, color.FUN = col_pal, 
                                                if.plot = !testing, sepList = contigSepList, sep.color = c("black", 
                                                                                                           "black"), sep.lty = 1, sep.lwd = 1, RowIndividualColors = row_groupings, 
                                                annotations_legend = annotations_legend, grouping_key_coln = grouping_key_coln, 
                                                ColIndividualColors = contig_colors, key.title = "Distribution of Expression", 
                                                key.xlab = "Modified Expression", key.ylab = "Count", 
                                                force_lmat = layout_lmat, force_lwid = layout_lwid, 
                                                force_lhei = layout_lhei, useRaster = useRaster)
  }
  if ("matrix" %in% is(obs_data)) {
    if (write_expr_matrix) {
      flog.info(paste("plot_cnv_observations:Writing observation data to", 
                      observation_file_base, sep = " "))
      row.names(obs_data) <- orig_row_names
      write.table(as.matrix(t(obs_data[data_observations$rowInd, 
                                       data_observations$colInd])), file = observation_file_base)
    }
  }
}


plot_cnv_references2 <- function (infercnv_obj, ref_data, ref_groups, name_ref_groups,  cluster_references, hclust_method, grouping_key_coln, col_pal, contig_seps, file_base_name, do_plot = TRUE, write_expr_matrix,  output_filename_prefix, cnv_ref_title, breaksList, gene_position_breaks,  x.center = x.center, layout_lmat = NULL, layout_lwid = NULL,  layout_lhei = NULL, layout_add = FALSE, testing = FALSE, useRaster = useRaster) {
  flog.info("plot_cnv_references:Start")
  flog.info(paste("Reference data size: Cells=", ncol(ref_data), 
                  "Genes=", nrow(ref_data), sep = " "))
  number_references <- ncol(ref_data)
  reference_ylab <- NA
  reference_data_file <- paste(file_base_name, sprintf("%s.references.txt", 
                                                       output_filename_prefix), sep = .Platform$file.sep)
  ref_seps <- c()
  ordered_names <- c()
  if (!is.null(infercnv_obj@tumor_subclusters$hc[["all_references"]])) {
    ordered_names <- infercnv_obj@tumor_subclusters$hc[["all_references"]]$labels[infercnv_obj@tumor_subclusters$hc[["all_references"]]$order]
    split_groups <- rep(1, length(ordered_names))
    ref_data <- ref_data[, ordered_names, drop = FALSE]
  }
  else if (all(name_ref_groups %in% infercnv_obj@tumor_subclusters$subclusters)) {
    if (cluster_references) {
      split_groups <- c()
      for (i in seq_along(name_ref_groups)) {
        if (!is.null(infercnv_obj@tumor_subclusters$hc[[name_ref_groups[i]]])) {
          ordered_names <- c(ordered_names, infercnv_obj@tumor_subclusters$hc[[name_ref_groups[i]]]$labels[infercnv_obj@tumor_subclusters$hc[[name_ref_groups[i]]]$order])
          ref_seps <- c(ref_seps, length(ordered_names))
          split_groups <- c(split_groups, rep(i, length(infercnv_obj@tumor_subclusters$hc[[name_ref_groups[i]]]$order)))
        }
        else {
          if ((length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])) == 
               2) || (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])) == 
                      1)) {
            ordered_names <- c(ordered_names, colnames(infercnv_obj@expr.data[, 
                                                                              unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]]), 
                                                                              drop = FALSE]))
            ref_seps <- c(ref_seps, length(ordered_names))
            split_groups <- c(split_groups, rep(i, length(length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])))))
          }
          else {
            flog.error("Unexpected error, should not happen.")
            stop("Error")
          }
        }
      }
    }
    ref_data <- ref_data[, ordered_names, drop = FALSE]
  }
  else {
    if (length(ref_groups) > 1) {
      if (cluster_references) {
        order_idx <- lapply(ref_groups, function(ref_grp) {
          if (cluster_references && length(ref_grp) > 
              2) {
            ref_hcl <- hclust(parallelDist(t(ref_data[, 
                                                      ref_grp]), threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), 
                              method = hclust_method)
            ref_grp <- ref_grp[ref_hcl$order]
          }
          ref_grp
        })
        ref_seps <- head(cumsum(lengths(order_idx)), 
                         -1)
        split_groups <- unlist(mapply(rep, seq_along(order_idx), 
                                      lengths(order_idx)))
        order_idx <- unlist(order_idx)
      }
      else {
        ref_seps <- head(cumsum(lengths(ref_groups)), 
                         -1)
        split_groups <- unlist(mapply(rep, seq_along(ref_groups), 
                                      lengths(ref_groups)))
        order_idx <- unlist(ref_groups)
      }
    }
    else {
      if (cluster_references) {
        ref_hcl <- hclust(parallelDist(t(ref_data), threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), 
                          method = hclust_method)
        order_idx = ref_hcl$order
      }
      else {
        order_idx = unlist(ref_groups)
      }
      split_groups <- rep(1, length(ref_groups[[1]]))
    }
    ref_data <- ref_data[, order_idx, drop = FALSE]
  }
  if (number_references == 1) {
    ref_data <- cbind(ref_data, ref_data)
    names(ref_data) <- c("", names(ref_data)[1])
    colnames(ref_data) <- c("", colnames(ref_data)[1])
    split_groups <- rep(1, 2)
  }
  flog.info(paste("plot_cnv_references:Number reference groups=", 
                  length(ref_groups), sep = " "))
  ref_data <- t(ref_data)
  ref_orig_names <- row.names(ref_data)
  reference_ylab <- cnv_ref_title
  if (number_references == 1) {
    row.names(ref_data) = rep("", 2)
  }
  else {
    row.names(ref_data) <- rep("", number_references)
  }
  alphabet =  c('0' = '#8497B0', '1' = '#878787', '2' = '#F0CE58','4' = '#D7EF9B', '5' = '#DCEAF7', '8' = '#F297A7','9' = '#ABDDDE','10' = 'lightskyblue','11' = '#B487B7','12' = '#1D65A6','13' = '#F2A104', '15' = '#0FFFFF', '16' = '#289E92','17' = '#C3FF00','19' = '#8ED973')
  #row_groupings <- as.matrix(get_group_color_palette()(length(table(split_groups)))[split_groups])
  row_groupings <- as.matrix(alphabet[split_groups])
  
  plot_cnv_observations2 <- function (infercnv_obj, obs_data, col_pal, contig_colors, contig_labels, 
                                      contig_names, contig_seps, num_obs_groups, file_base_name, 
                                      do_plot = TRUE, write_expr_matrix, write_phylo, output_filename_prefix, 
                                      cnv_title, cnv_obs_title, contig_lab_size = 1, contigs, cluster_contig = NULL, 
                                      obs_annotations_groups, obs_annotations_names, grouping_key_coln, 
                                      cluster_by_groups, breaksList, gene_position_breaks, x.center, 
                                      hclust_method = "ward.D", testing = FALSE, layout_lmat = NULL, 
                                      layout_lhei = NULL, layout_lwid = NULL, useRaster = useRaster) 
  {
    flog.info("plot_cnv_observation:Start")
    flog.info(paste("Observation data size: Cells=", nrow(obs_data), 
                    "Genes=", ncol(obs_data), sep = " "))
    observation_file_base <- paste(file_base_name, sprintf("%s.observations.txt", 
                                                           output_filename_prefix), sep = .Platform$file.sep)
    dendrogram_file_path = paste(file_base_name, sprintf("%s.observations_dendrogram.txt", 
                                                         output_filename_prefix), sep = .Platform$file.sep)
    hcl_desc <- "General"
    hcl_group_indices <- seq_len(ncol(obs_data))
    if (!is.null(cluster_contig)) {
      hcl_contig_indices <- which(contig_names %in% cluster_contig)
      if (length(hcl_contig_indices) > 0) {
        hcl_group_indices <- hcl_contig_indices
        hcl_desc <- paste(cluster_contig, collapse = "_")
        flog.info(paste("plot_cnv_observation:Clustering only by contig ", 
                        cluster_contig))
        infercnv_obj@tumor_subclusters = NULL
      }
      else {
        flog.warn(paste("plot_cnv_observations: Not able to cluster by", 
                        cluster_contig, "Clustering by all genomic locations.", 
                        "To cluster by local genomic location next time", 
                        "select from:", unique(contig_names), collapse = ",", 
                        sep = " "))
      }
    }
    obs_dendrogram <- list()
    ordered_names <- NULL
    isfirst <- TRUE
    hcl_obs_annotations_groups <- vector()
    obs_seps <- c()
    sub_obs_seps <- c()
    if (!is.null(infercnv_obj@tumor_subclusters)) {
      if (cluster_by_groups) {
        split_groups = vector()
        for (i in seq_along(obs_annotations_names)) {
          if (!is.null(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])) {
            obs_dendrogram[[i]] = as.dendrogram(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]])
            ordered_names <- c(ordered_names, infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$labels[infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$order])
            obs_seps <- c(obs_seps, length(ordered_names))
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                            rep(i, length(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]$order)))
            if (write_phylo) {
              if (isfirst) {
                write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]), 
                           file = dendrogram_file_path)
                isfirst <- FALSE
              }
              else {
                write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[obs_annotations_names[i]]]), 
                           file = dendrogram_file_path, append = TRUE)
              }
            }
          }
          else {
            if ((length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
                 2) || (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
                        1)) {
              if (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 
                  2) {
                obs_dendrogram[[i]] <- .pairwise_dendrogram(colnames(infercnv_obj@expr.data[, 
                                                                                            unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                                            drop = FALSE]))
              }
              else {
                obs_dendrogram[[i]] <- .single_element_dendrogram(colnames(infercnv_obj@expr.data[, 
                                                                                                  unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                                                  drop = FALSE]))
              }
              ordered_names <- c(ordered_names, colnames(infercnv_obj@expr.data[, 
                                                                                unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), 
                                                                                drop = FALSE]))
              obs_seps <- c(obs_seps, length(ordered_names))
              hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                              rep(i, length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]))))
            }
            else {
              flog.error("Unexpected error, should not happen.")
              stop("Error")
            }
          }
          for (subcluster in names(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) {
            tmp = infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]][[subcluster]]
            tmp[] = subcluster
            split_groups = c(split_groups, tmp)
          }
        }
        split_groups <- split_groups[ordered_names]
        if (length(obs_dendrogram) > 1) {
          obs_dendrogram <- do.call(merge, obs_dendrogram)
        }
        else {
          obs_dendrogram <- obs_dendrogram[[1]]
        }
        for (subtumor in infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]) {
          sub_obs_seps <- c(sub_obs_seps, (sub_obs_seps[length(sub_obs_seps)] + 
                                             length(subtumor)))
        }
      }
      else {
        obs_hcl <- infercnv_obj@tumor_subclusters$hc[["all_observations"]]
        if (write_phylo) {
          write.tree(as.phylo(obs_hcl), file = dendrogram_file_path)
        }
        obs_dendrogram <- as.dendrogram(obs_hcl)
        ordered_names <- obs_hcl$labels[obs_hcl$order]
        if (num_obs_groups > 1) {
          split_groups <- cutree(obs_hcl, k = num_obs_groups)
        }
        else {
          split_groups = vector()
          for (subcluster in names(infercnv_obj@tumor_subclusters$subclusters[["all_observations"]])) {
            tmp = infercnv_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster]]
            tmp[] = subcluster
            split_groups = c(split_groups, tmp)
          }
        }
        split_groups <- split_groups[ordered_names]
        hcl_obs_annotations_groups <- obs_annotations_groups[ordered_names]
        flog.info("plot_cnv_observation:Writing observations by grouping.")
        for (cut_group in unique(split_groups)) {
          group_memb <- names(split_groups)[which(split_groups == 
                                                    cut_group)]
          memb_file <- file(paste(file_base_name, paste(hcl_desc, 
                                                        "HCL", cut_group, "members.txt", sep = "_"), 
                                  sep = .Platform$file.sep))
          write.table(as.matrix(obs_data[group_memb, ]), 
                      memb_file)
          ordered_memb <- which(ordered_names %in% group_memb)
          if (is.null(obs_seps)) {
            obs_seps <- c(length(ordered_memb))
          }
          else {
            obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + 
                                       length(ordered_memb)))
          }
        }
        obs_seps <- c(obs_seps, length(ordered_names))
        sub_obs_seps = obs_seps
      }
    }
    else if (cluster_by_groups) {
      flog.info(paste("clustering observations via method: ", 
                      hclust_method, sep = ""))
      for (i in seq_len(max(obs_annotations_groups))) {
        cell_indices_in_group <- which(obs_annotations_groups == 
                                         i)
        num_cells_in_group <- length(cell_indices_in_group)
        flog.info(sprintf("Number of cells in group(%d) is %d", 
                          i, num_cells_in_group))
        if (num_cells_in_group < 2) {
          flog.info(sprintf("Skipping group: %d, since less than 2 entries", 
                            i))
          ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups == 
                                                                       i), , drop = FALSE]))
          obs_dendrogram[[length(obs_dendrogram) + 1]] = .single_element_dendrogram(unique_label = row.names(obs_data[which(obs_annotations_groups == 
                                                                                                                              i), , drop = FALSE]))
          if (write_phylo) {
            if (isfirst) {
              write(row.names(obs_data[which(obs_annotations_groups == 
                                               i), ]), file = dendrogram_file_path)
              isfirst <- FALSE
            }
            else {
              write(row.names(obs_data[which(obs_annotations_groups == 
                                               i), ]), file = dendrogram_file_path, append = TRUE)
            }
          }
        }
        else {
          data_to_cluster <- obs_data[cell_indices_in_group, 
                                      hcl_group_indices, drop = FALSE]
          flog.info(paste("group size being clustered: ", 
                          paste(dim(data_to_cluster), collapse = ","), 
                          sep = " "))
          group_obs_hcl <- hclust(parallelDist(data_to_cluster, 
                                               threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), 
                                  method = hclust_method)
          ordered_names <- c(ordered_names, group_obs_hcl$labels[group_obs_hcl$order])
          group_obs_dend <- as.dendrogram(group_obs_hcl)
          obs_dendrogram[[length(obs_dendrogram) + 1]] <- group_obs_dend
          if (write_phylo) {
            if (isfirst) {
              write.tree(as.phylo(group_obs_hcl), file = dendrogram_file_path)
              isfirst <- FALSE
            }
            else {
              write.tree(as.phylo(group_obs_hcl), file = dendrogram_file_path, 
                         append = TRUE)
            }
          }
        }
        hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, 
                                        rep(i, num_cells_in_group))
        obs_seps <- c(obs_seps, length(ordered_names))
      }
      if (length(obs_dendrogram) > 1) {
        obs_dendrogram <- do.call(merge, obs_dendrogram)
      }
      else {
        obs_dendrogram <- obs_dendrogram[[1]]
      }
      split_groups <- rep(1, dim(obs_data)[1])
      names(split_groups) <- ordered_names
      sub_obs_seps = obs_seps
    }
    else {
      flog.info(paste("clustering observations via method: ", 
                      hclust_method, sep = ""))
      if (nrow(obs_data) > 1) {
        obs_hcl <- hclust(parallelDist(obs_data[, hcl_group_indices], 
                                       threads = infercnv:::infercnv.env$GLOBAL_NUM_THREADS), method = hclust_method)
        if (write_phylo) {
          write.tree(as.phylo(obs_hcl), file = dendrogram_file_path)
        }
        obs_dendrogram <- as.dendrogram(obs_hcl)
        ordered_names <- obs_hcl$labels[obs_hcl$order]
        split_groups <- cutree(obs_hcl, k = num_obs_groups)
        split_groups <- split_groups[ordered_names]
        flog.info("plot_cnv_observation:Writing observations by grouping.")
        for (cut_group in unique(split_groups)) {
          group_memb <- names(split_groups)[which(split_groups == 
                                                    cut_group)]
          memb_file <- file(paste(file_base_name, paste(hcl_desc, 
                                                        "HCL", cut_group, "members.txt", sep = "_"), 
                                  sep = .Platform$file.sep))
          write.table(as.matrix(obs_data[group_memb, ]), 
                      memb_file)
          ordered_memb <- which(ordered_names %in% group_memb)
          if (is.null(obs_seps)) {
            obs_seps <- c(length(ordered_memb))
          }
          else {
            obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + 
                                       length(ordered_memb)))
          }
        }
        obs_seps <- c(obs_seps, length(ordered_names))
      }
      else {
        obs_dendrogram <- .single_element_dendrogram(row.names(obs_data))
        ordered_names <- row.names(obs_data)
        obs_seps <- c(1)
        split_groups <- c(1)
        names(split_groups) <- ordered_names
      }
      sub_obs_seps = obs_seps
      hcl_obs_annotations_groups <- obs_annotations_groups[ordered_names]
    }
    if (length(obs_seps) > 1) {
      obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) - 
                                                           1):1]
    }
    idx = 1
    split_groups_as_idx = split_groups
    for (grp in unique(split_groups)) {
      to_replace = which(split_groups == grp)
      split_groups_as_idx[to_replace] = idx
      idx = idx + 1
    }
    split_groups_as_idx = as.integer(split_groups_as_idx)
    row_groupings <- infercnv:::get_group_color_palette()(length(table(split_groups_as_idx)))[split_groups_as_idx]
    #row_groupings <- cbind(row_groupings, infercnv:::get_group_color_palette()(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])
    alphabet =  c('0' = '#8497B0', '1' = '#878787', '2' = '#F0CE58', '3' = '#EB545C','4' = '#D7EF9B', '5' = '#DCEAF7', '6' = '#FFEE00','7' = '#00475F', '8' = '#F297A7','9' = '#ABDDDE','10' = 'lightskyblue','11' = '#B487B7','12' = '#1D65A6','13' = '#F2A104','14' = '#FC6FCF','15' = '#0FFFFF', '16' = '#289E92','17' = '#C3FF00','18' = '#BC589B','19' = '#8ED973')
    row_groupings <- cbind(row_groupings, alphabet[hcl_obs_annotations_groups])
    
    annotations_legend = as.matrix(cbind(names(alphabet),alphabet))
    #annotations_legend <- cbind(obs_annotations_names, infercnv:::get_group_color_palette()(length(table(hcl_obs_annotations_groups))))
    colnames(annotations_legend) = c("name_ref_groups", "")
    flog.info("plot_cnv_observation:Writing observation groupings/color.")
    groups_file_name <- file.path(file_base_name, sprintf("%s.observation_groupings.txt", 
                                                          output_filename_prefix))
    file_groups <- cbind(split_groups, row_groupings[, 1], hcl_obs_annotations_groups, 
                         row_groupings[, 2])
    colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color", 
                               "Annotation Group", "Annotation Color")
    write.table(file_groups, groups_file_name)
    flog.info("plot_cnv_observation:Done writing observation groupings/color.")
    contigSepList <- infercnv:::create_sep_list(row_count = nrow(obs_data), 
                                                col_count = ncol(obs_data), row_seps = obs_seps, col_seps = contig_seps)
    obs_data <- obs_data[ordered_names, ]
    orig_row_names <- row.names(obs_data)
    row.names(obs_data) <- rep("", nrow(obs_data))
    flog.info("plot_cnv_observation:Writing observation heatmap thresholds.")
    heatmap_thresholds_file_name <- file.path(file_base_name, 
                                              sprintf("%s.heatmap_thresholds.txt", output_filename_prefix))
    write.table(breaksList, heatmap_thresholds_file_name, row.names = FALSE, 
                col.names = FALSE)
    flog.info("plot_cnv_observation:Done writing observation heatmap thresholds.")
    if (do_plot) {
      data_observations <- infercnv:::heatmap.cnv(obs_data, Rowv = obs_dendrogram, 
                                                  Colv = FALSE, cluster.by.row = TRUE, cluster.by.col = FALSE, 
                                                  main = cnv_title, ylab = cnv_obs_title, margin.for.labCol = 2, 
                                                  xlab = "Genomic Region", key = TRUE, labCol = contig_labels, 
                                                  cexCol = contig_lab_size, cexAt = c(1, contig_seps), 
                                                  notecol = "black", density.info = "histogram", denscol = "blue", 
                                                  trace = "none", dendrogram = "row", cexRow = 0.8, 
                                                  breaks = breaksList, gene_position_breaks = gene_position_breaks, 
                                                  scale = "none", x.center = x.center, color.FUN = col_pal, 
                                                  if.plot = !testing, sepList = contigSepList, sep.color = c("black", 
                                                                                                             "black"), sep.lty = 1, sep.lwd = 1, RowIndividualColors = row_groupings, 
                                                  annotations_legend = annotations_legend, grouping_key_coln = grouping_key_coln, 
                                                  ColIndividualColors = contig_colors, key.title = "Distribution of Expression", 
                                                  key.xlab = "Modified Expression", key.ylab = "Count", 
                                                  force_lmat = layout_lmat, force_lwid = layout_lwid, 
                                                  force_lhei = layout_lhei, useRaster = useRaster)
    }
    if ("matrix" %in% is(obs_data)) {
      if (write_expr_matrix) {
        flog.info(paste("plot_cnv_observations:Writing observation data to", 
                        observation_file_base, sep = " "))
        row.names(obs_data) <- orig_row_names
        write.table(as.matrix(t(obs_data[data_observations$rowInd, 
                                         data_observations$colInd])), file = observation_file_base)
      }
    }
  }
  contigSepList <- infercnv:::create_sep_list(row_count = nrow(ref_data), 
                                              col_count = ncol(ref_data), row_seps = ref_seps, col_seps = contig_seps)
  flog.info("plot_cnv_references:Plotting heatmap.")
  if (do_plot) {
    data_references <- infercnv:::heatmap.cnv(ref_data, main = NULL, 
                                              ylab = reference_ylab, xlab = NULL, key = FALSE, 
                                              labCol = rep("", nrow(ref_data)), notecol = "black", 
                                              trace = "none", dendrogram = "none", Colv = FALSE, 
                                              Rowv = FALSE, cexRow = 0.4, breaks = breaksList, 
                                              gene_position_breaks = gene_position_breaks, scale = "none", 
                                              x.center = x.center, color.FUN = col_pal, sepList = contigSepList, 
                                              RowIndividualColors = row_groupings, annotations_legend = annotations_legend, 
                                              grouping_key_coln = grouping_key_coln, sep.color = c("black", 
                                                                                                   "black"), sep.lty = 1, sep.lwd = 1, if.plot = !testing, 
                                              force_lmat = layout_lmat, force_lwid = layout_lwid, 
                                              force_lhei = layout_lhei, force_add = layout_add, 
                                              useRaster = useRaster)
  }
  if ("matrix" %in% is(ref_data)) {
    if (write_expr_matrix) {
      row.names(ref_data) <- ref_orig_names
      flog.info(paste("plot_cnv_references:Writing reference data to", 
                      reference_data_file, sep = " "))
      write.table(as.matrix(t(ref_data[data_references$rowInd, 
                                       data_references$colInd])), file = reference_data_file)
    }
  }
}


#创建inferCNV对象，直接给相应的文件路径即可
infercnv_obj = CreateInfercnvObject(delim = '\t',chr_exclude = c("chrX", "chrY", "chrM"), raw_counts_matrix = 'result/CNV/inferCNV/expmat4cnv.tsv', annotations_file = 'result/CNV/inferCNV/cluster_id.csv', gene_order_file = 'result/CNV/inferCNV/transloc.txt', ref_group_names = c("0", "1", "2", "4", "5", "8", "9", "10", "11", "12", "13", "15", "16", "17", "19"))


infercnv_obj = infercnv::run(infercnv_obj, num_threads = 16, cutoff=0.1, out_dir='result/CNV/inferCNV/result', cluster_by_groups=TRUE, denoise=TRUE, write_expr_matrix = TRUE, write_phylo = TRUE, HMM=TRUE, debug = TRUE)


aa = readRDS('result/CNV/inferCNV/result/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')

plot_cnv2(infercnv_obj=aa, k_obs_groups=1, cluster_by_groups=TRUE, cluster_references=TRUE, plot_chr_scale=FALSE, chr_lengths=FALSE, out_dir='result/', title='RNA-derived CNV', output_filename='used_this', output_format= "pdf", x.center=1, x.range=c(-1,3), png_res=300, useRaster=TRUE)

aa = readRDS('result/CNV/inferCNV/result/run.final.infercnv_obj')
plot_cnv2(infercnv_obj=aa, k_obs_groups=1, cluster_by_groups=TRUE, cluster_references=TRUE, plot_chr_scale=FALSE, chr_lengths=NULL, out_dir='./', title='RNA-derived CNV', output_filename='used_this2', output_format= "pdf", x.center=1, x.range= "auto", useRaster=TRUE)

