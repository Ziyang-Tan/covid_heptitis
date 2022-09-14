library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

raw <- read_csv('data/bulk_data.csv.gz')
raw_ctrl <- read_csv('data/bulk_data_ctrl.csv.gz')

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
  }) %>% Reduce(f=function(x,y) left_join(x,y,by='V'))
}

res <- process_mixcr_raw(raw) %>% count_V_fraction()
res_ctrl <- process_mixcr_raw(raw_ctrl) %>% count_V_fraction()

rowOrder <- (res %>% rowwise() %>% mutate(mean = mean(c_across(contains('sample')), na.rm = T)) %>% arrange(desc(mean)))$V

res.long <- res %>% 
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group='Hepatitis')

res_ctrl.long <- res_ctrl %>%
  mutate(V = factor(V, levels = rowOrder)) %>%
  tidyr::pivot_longer(cols = setdiff(colnames(res_ctrl), c('V')), names_to = 'subject', values_to = 'freq') %>%
  mutate(group = if_else(subject == 'PBMC Ctrl', 'PBMC control', 'ISAC34'))

ggplot(rbind(res.long, res_ctrl.long) %>% filter(group != 'ISAC34'), aes(x=V, y=freq, color=group)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = 'figures/V_distribution.pdf')
