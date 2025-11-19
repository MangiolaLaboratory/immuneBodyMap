get_brms_tar <- 
  function(
    cell_type, 
    gene_id, 
    base_path = '/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE/'
  ){
    
    cli::cli_alert_info("Locating cell type `{cell_type}` under {base_path} ...")
    
    # 1. List all versioned directories
    ct_list <- list.dirs(base_path, recursive = F) %>% basename() 
    ct_list <- ct_list[ct_list %>% str_starts('V[0-9]_')]
    ct_name <- str_remove(ct_list, 'V[0-9]_')
    
    # 2. Check cell_type existence
    if (!make.names(cell_type) %in% ct_name){
      cli::cli_alert_danger("Cell type '{cell_type}' not found.")
      stop("Cell type not found.")
    }
    
    # 3. Choose latest version
    ct_list = ct_list[ct_list %>% str_detect(make.names(cell_type))]
    version_no = ct_list %>% 
      str_extract("^V\\d+") %>% 
      str_remove("^V") %>% 
      as.integer() %>% 
      max(na.rm = TRUE)
    
    targets_path = glue::glue(
      '{base_path}V{version_no}_{make.names(cell_type)}/_targets'
    )
    # mapping_file <- glue::glue('{targets_path}/{make.names(cell_type)}_brms_mapping_table')
    
    cli::cli_alert_info("Using version V{version_no} at {targets_path}")
    
    # 4. Load mapping table
    # if (!file.exists(mapping_file)) {
    #   cli::cli_alert_danger("Mapping file not found at {mapping_file}")
    #   stop("Mapping file missing.")
    # }
    # 
    # mapping_table <- qs2::qs_read(mapping_file)
    mapping_table <- 
      read_delim(glue::glue('{targets_path}/meta/meta'), delim = '|') %>% 
      filter(name  %>% str_starts('estimates_chunk')) %>% 
      filter(type == 'branch' & is.na(error)) %>% 
      mutate(
        gene = warnings %>% 
          str_extract("(?<=Gene:___)ENSG\\d+(?=___)")
      )
    
    if (!gene_id %in% mapping_table$gene) {
      cli::cli_alert_danger("Gene ID '{gene_id}' not found in mapping.")
      stop("Gene ID not found.")
    }
    
    # 5. Load target object
    target_name <- mapping_table %>%
      filter(gene == gene_id) %>%
      pull(name) 
    
    if (length(target_name) != 1) {
      cli::cli_alert_danger("Expected 1 target for gene '{gene_id}', found {length(target_name)}.")
      stop("Ambiguous or missing target.")
    }
    
    brms_path <- glue::glue("{targets_path}/objects/{target_name}")
    
    if (!file.exists(brms_path)) {
      cli::cli_alert_danger("Target file not found at {brms_path}")
      stop("Target file missing.")
    }
    
    cli::cli_alert_success("Successfully loaded target for gene '{gene_id}'.")
    
    return(qs2::qs_read(brms_path))
  }