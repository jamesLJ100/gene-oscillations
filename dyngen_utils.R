propagate_module_assignments <- function(tf_modules, target_edges) {
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
      warning(sprintf(
        "Could not resolve %d target(s): %s",
        length(unresolved),
        paste(unresolved, collapse = ", ")
      ))
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