library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

raw <- read_csv('data/bulk_data.csv.gz')
raw_ctrl <- read_csv('data/bulk_data_ctrl.csv.gz')
raw_ISAC <- read_csv('data/bulk_data_ISAC_ref.csv.gz') %>% 
  rename(`NGI ID` = ngi_id) %>% 
  left_join(read_csv('data/bulk_info_ISAC_ref.csv.gz'), by = 'NGI ID')
public_ref <- readxl::read_excel('data/public kids control.xlsx')

process_mixcr_raw <- function(raw){
  raw <- raw[,colSums(!is.na(raw))!=0]
  data <- raw %>% 
    mutate(V = sub('\\*.*$', '', allVHitsWithScore),
           D = sub('\\*.*$', '', allDHitsWithScore),
           J = sub('\\*.*$', '', allJHitsWithScore),
           C = sub('\\*.*$', '', allCHitsWithScore)) %>%
    rename(CDR3 = nSeqCDR3,
           CDR3aa = aaSeqCDR3) %>%
    select(cloneCount, cloneFraction, V, D, J, C, CDR3, CDR3aa, `NGI ID`, `User ID`, Mreads, proj_id) %>%
    mutate(normCloneCount = round(cloneCount/Mreads))
  return(data)
}

count_V_fraction <- function(data) {
  res <- lapply(unique(data$`User ID`), function(cur_sample){
    sub_data <- data %>% filter(`User ID` == cur_sample)
    tmp <- sub_data %>% 
      group_by(V) %>% 
      summarise(freq = sum(cloneFraction))
    colnames(tmp) <- c('V', cur_sample)
    return(tmp)
  }) %>% Reduce(f=function(x,y) full_join(x,y,by='V'))
}

res <- process_mixcr_raw(raw) %>% count_V_fraction()
res_ctrl <- process_mixcr_raw(raw_ctrl) %>% count_V_fraction()
res_ISAC <- process_mixcr_raw(raw_ISAC) %>% count_V_fraction()

rowOrder <- (res %>% rowwise() %>% mutate(mean = mean(c_across(contains('sample')), na.rm = T)) %>% arrange(desc(mean)))$V

res.long <- res %>% 
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group='Hepatitis')

res_ctrl.long <- res_ctrl %>%
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res_ctrl), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group = if_else(subject == 'PBMC Ctrl', 'PBMC control', 'ISAC'))

res_ISAC.long <- res_ISAC %>%
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res_ISAC), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group = 'ISAC') %>%
  filter(freq<0.5)

res_public.long <- public_ref %>%
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(public_ref), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group = 'public children control')

ggplot(rbind(res.long, res_ctrl.long, res_ISAC.long, res_public.long), aes(x=V, y=freq, color=group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(group=group), position=position_jitterdodge(), size=0.3, alpha=0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = 'figures/V_distribution_with_ISAC_ref.pdf', width = 15, height = 8)
