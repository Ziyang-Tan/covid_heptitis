library(dplyr)
library(readr)
library(ggpubr)

# parameters

bulk_data_dir <- '/Users/tan/Library/CloudStorage/OneDrive-KI.SE/TCR_processed_data/bulk'
bulk_proj_list <- c('P27104')

# bulk TCR--------------------------------------------------------------------------------------------------

# read data

sample_info_bulk <- lapply(bulk_proj_list, function(proj_id){
  return(
    read_delim(
      Sys.glob(file.path(bulk_data_dir, proj_id, '*sample_info.txt')),
      delim = '\t',
      show_col_types = FALSE
    ) %>%
      tibble::add_column(proj_id = proj_id)
  )
}) %>% do.call(what = rbind)

data_bulk <- lapply(sample_info_bulk$`NGI ID`, function(ngi_id){
  return(
    read_delim(
      Sys.glob(file.path(bulk_data_dir, '*', 'data', paste0(ngi_id, '*TRB*'))),
      delim = '\t',
      show_col_types = FALSE
    ) %>%
      tibble::add_column(`NGI ID` = ngi_id)
  )
}) %>% do.call(what = rbind) %>% left_join(sample_info_bulk, by = 'NGI ID')


# write to file
write_csv(data_bulk, file = file.path('.', 'data', 'bulk_data_ctrl.csv.gz'))

