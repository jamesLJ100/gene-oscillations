library(dplyr)
propagate_module_assignments <- function(sim, context = NULL) {

  feature_net <- sim$model$feature_network
  
  tf_modules <- bind_rows(
    feature_net %>%
      filter(!is.na(from_module)) %>%
      dplyr::select(gene = from, module = from_module),
    feature_net %>%
      filter(!is.na(to_module)) %>%
      dplyr::select(gene = to, module = to_module)
  ) %>%
    distinct()
  
  target_edges <- feature_net %>%
    filter(grepl("^Target", to)) %>%
    dplyr::select(from, to)
  
  
  gene_module <- setNames(tf_modules$module, tf_modules$gene)
  hops <- setNames(rep(0L, nrow(tf_modules)), tf_modules$gene)
  
  repeat {
    unresolved <- target_edges$to[!(target_edges$to %in% names(gene_module))]
    if (length(unresolved) == 0) break
    
    newly_resolved <- target_edges %>%
      filter(to %in% unresolved, from %in% names(gene_module)) %>%
      mutate(
        module = gene_module[from],
        hops   = hops[from] + 1L
      )
    
    if (nrow(newly_resolved) == 0) {
      unresolved <- unique(unresolved)
      prefix <- if (!is.null(context)) sprintf("[%s] ", context) else ""
      max_show <- 20L
      shown <- head(unresolved, max_show)
      target_list <- paste(shown, collapse = ", ")
      if (length(unresolved) > max_show) {
        target_list <- sprintf("%s, ... and %d more",
                               target_list, length(unresolved) - max_show)
      }
      warning(sprintf(
        "%s%d target(s) not traceable to a TF (cycle with no TF ancestor); labelling module/hops as NA: %s",
        prefix,
        length(unresolved),
        target_list
      ))
      gene_module <- c(gene_module, setNames(rep(NA_character_, length(unresolved)), unresolved))
      hops        <- c(hops,        setNames(rep(NA_integer_,   length(unresolved)), unresolved))
      break
    }

    gene_module <- c(gene_module, setNames(newly_resolved$module, newly_resolved$to))
    hops        <- c(hops,        setNames(newly_resolved$hops,   newly_resolved$to))
  }
  
  tibble(
    gene   = names(gene_module),
    module = unname(gene_module),
    hops   = unname(hops)
  )
}
# sim_file <- file.path(proj_root, "synthetic/data/dyngen", paste0("c50g207_1", "_sim.rds"))
# sim <- readRDS(sim_file)
# result <- propagate_module_assignments(sim)