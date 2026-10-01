tpm = read_tsv("data/a549-tak1-ko-tpm.txt") %>%
  mutate(
    condition = factor(condition, levels = unique(condition)),
    rescue = factor(rescue, levels = unique(rescue)),
    treatment = factor(treatment, levels = unique(treatment))
  )

dea = read_tsv("data/a549-tak1-ko-dea.txt") %>%
  mutate(sig = case_when(
    log2fold >= 1 & FDR <= 0.05 ~ "up",
    log2fold <=-1 & FDR <= 0.05 ~ "dn",
    TRUE ~ "ns"
  )) %>%
  mutate(
    condition = factor(condition, levels = unique(condition)),
    rescue = factor(rescue, levels = unique(rescue)),
    treatment = factor(treatment, levels = unique(treatment))
  )





mat = tpm %>%
  filter(treatment == "IL1B") %>%
  select(Gene, rep, condition, rescue, log2fold) %>%
  arrange(Gene, condition, rescue, rep) %>%
  unite("sample", c(condition, rescue, rep), sep = " ") %>%
  pivot_wider(names_from = "sample", values_from = "log2fold") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat_ns = tpm %>%
  filter(treatment == "NS") %>%
  select(Gene, rep, condition, rescue, log2tpm) %>%
  arrange(Gene, condition, rescue, rep) %>%
  unite("sample", c(condition, rescue, rep), sep = " ") %>%
  pivot_wider(names_from = "sample", values_from = "log2tpm") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat_il1b = tpm %>%
  filter(treatment == "IL1B") %>%
  select(Gene, rep, condition, rescue, log2tpm) %>%
  arrange(Gene, condition, rescue, rep) %>%
  unite("sample", c(condition, rescue, rep), sep = " ") %>%
  pivot_wider(names_from = "sample", values_from = "log2tpm") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat_a = mat_ns
mat_b = mat_il1b

colnames(mat_a) = paste("(NS)", colnames(mat_a))
colnames(mat_b) = paste("(IL1B)", colnames(mat_b))

mat_combined = cbind(mat_a, mat_b) %>%
  as.data.frame() %>%
  select(
    contains("WT naive"),
    contains("sham naive"),
    contains("TAK1KO naive"),
    contains("TAK1KO lipid"),
    contains("TAK1KO TAK1-smRNA"),
    contains("TAK1KO mTAK1-smRNA")
  ) %>%
  as.matrix()

rm(mat_a, mat_b)








dea_ns = read_tsv("data/effect-of-group-on-NS_dea.txt") %>%
  mutate(sig = case_when(
    log2fold >= 1 & FDR <= 0.05 ~ "up",
    log2fold <=-1 & FDR <= 0.05 ~ "dn",
    TRUE ~ "ns"
  )) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  group_by(Gene) %>%
  filter(!all(FDR == 1 & log2fold == 0))


dea_il1b = read_tsv("data/effect-of-group-on-IL1B_dea.txt") %>%
  mutate(sig = case_when(
    log2fold >= 1 & FDR <= 0.05 ~ "up",
    log2fold <=-1 & FDR <= 0.05 ~ "dn",
    TRUE ~ "ns"
  )) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  group_by(Gene) %>%
  filter(!all(FDR == 1 & log2fold == 0))