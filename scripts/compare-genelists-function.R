


plot_gene_sets = function(
    genelists, 
    mat, 
    mat_cols = 1:24,
    output = "plot",
    plotType = "venn",
    heatmapTitle = NULL,
    heatmapRowNames = FALSE
) {
  
  # binary membership matrix
  df = tibble(
    group = names(genelists),
    Gene = genelists
  ) %>%
    unnest(Gene) %>%
    mutate(value = 1) %>%
    distinct(Gene, group, .keep_all = TRUE) %>%
    pivot_wider(
      names_from = group,
      values_from = value,
      values_fill = 0
    ) %>%
    column_to_rownames("Gene")
  
  # membership codes
  codes = df %>%
    rownames_to_column("Gene") %>%
    mutate(
      code = apply(select(., -Gene), 1, paste0, collapse = "")
    ) %>%
    select(Gene, code)
  
  # all non-empty overlap patterns
  code_list = codes %>%
    count(code) %>%
    pull(code)
  
  if(output == "plot"){
    # make heatmaps
    hms = lapply(code_list, function(cd) {
      
      genes = codes %>%
        filter(code == cd) %>%
        pull(Gene)
      
      title = names(genelists)[as.logical(as.integer(strsplit(cd, "")[[1]]))] %>%
        paste(collapse = "\n∩\n")
      
      mat[genes, mat_cols, drop = FALSE] %>%
        Heatmap(
          .,
          cluster_columns = F,
          col = generic_l2f_heatmap_colors(.),
          show_row_names = heatmapRowNames,
          row_title = title,
          row_title_rot = 0,
          show_column_names = T,
          column_title = heatmapTitle,
          # column_split = colsplits,
          # column_title_rot = 45,
          column_names_rot = 45,
          row_names_gp = gpar(cex = 0.5)
        )
    })
    
    # Euler/Venn plot
    print(plot_euler(df, plotType = plotType, plotTitle = NULL))
    
    # stack heatmaps
    print(Reduce(`%v%`, hms))
  }
  
  else if (output == "codes"){
    return(codes)
  }
  
  else {
    print("Invalid group type -- returning nothing")
  }

}






