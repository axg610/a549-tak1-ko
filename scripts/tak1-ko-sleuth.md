### === prep cluster ===

```sh
ssh alex.gao1@arc.ucalgary.ca

salloc -c 6 --mem 128GB --time 05:00:00
module load R/4.4.1

cd /work/newton_lab/ag_analysis/a549-tak1-ko

mkdir -p sleuth

R
```

# Main models and tpms

### === prep sleuth inputs ===

```r
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/a549-tak1-ko")

library(dplyr)
library(readr)
library(sleuth)
library(tidyr)

t2g <- read_tsv("/work/newton_lab/ag_analysis/ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt") %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

grouplevels = c(
    "WT_naive",
    "sham_naive",
    "TAK1KO_naive",
    "TAK1KO_lipid",
    "TAK1KO_TAK1-smRNA",
    "TAK1KO_mTAK1-smRNA"
)

s2c = read_tsv("tak1-ko-meta.txt") %>%
    mutate(condition = factor(condition, levels = c("WT", "sham", "TAK1KO"))) %>%
    mutate(treatment = factor(treatment, levels = c("NS", "IL1B"))) %>%
    mutate(rescue = factor(rescue, levels = c("naive", "lipid", "TAK1-smRNA", "mTAK1-smRNA"))) %>%

    mutate(
        group = paste0(condition, "_", rescue),
        group = factor(
            group,
            levels = grouplevels
        )
    ) %>%
    select(sample, rep, time, group, condition, rescue, treatment, path = kallisto_path)

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}
```

### === create sleuth objects ===

```r
# no modelling, just for tpms
so = sleuth_prep(
    sample_to_covariates = s2c,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

saveRDS(so, "sleuth/so.rds")

# individual objects for modelling
for(g in grouplevels){

    print(paste0("======= Processing ", g, " ======="))

    internal.s2c = filter(s2c, group == g)

    internal.so = sleuth_prep(
        sample_to_covariates = internal.s2c,
        target_mapping = t2g,
        gene_mode = TRUE,
        aggregation_column = "Gene",
        filter = new_filter,
        num_cores = 4
    )

    internal.objName = paste0("so_", g)
    internal.objPath = paste0("sleuth/", internal.objName, ".rds")

    assign(internal.objName, internal.so)
    saveRDS(internal.so, internal.objPath)

    print(paste0("======= Finished ", g, " ======="))

}
```

### === pull tpm table ===
```r
tpm = kallisto_table(so, use_filtered = FALSE)

tpm_clean = tpm %>%
    select(Gene = target_id, sample, rep, group, condition, rescue, treatment, time, tpm) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(tpm == 0)) %>%
    ungroup() %>%
    mutate(log2tpm = log2(tpm + 0.1)) %>%
    group_by(Gene, rep, time, group) %>%
    mutate(log2fold = log2tpm - log2tpm[treatment == "NS"]) %>%
    ungroup() %>%
    mutate(fold = 2^log2fold) %>%
    arrange(Gene, rep, group, treatment) %>%
    select(Gene, sample, rep, group, condition, rescue, treatment, time, tpm, log2tpm, fold, log2fold)

write_tsv(tpm_clean, "sleuth/a549-tak1-ko-tpm.txt")
```

### === model treatment ===
```r
do_sleuth_modelling = function(sleuth.obj){
    sleuth.obj %>%
        sleuth_fit(., ~treatment) %>%
        sleuth_wt(., "treatmentIL1B")
}

`so_WT_naive`             <- do_sleuth_modelling(`so_WT_naive`)
`so_sham_naive`           <- do_sleuth_modelling(`so_sham_naive`)
`so_TAK1KO_naive`         <- do_sleuth_modelling(`so_TAK1KO_naive`)
`so_TAK1KO_lipid`         <- do_sleuth_modelling(`so_TAK1KO_lipid`)
`so_TAK1KO_TAK1-smRNA`    <- do_sleuth_modelling(`so_TAK1KO_TAK1-smRNA`)
`so_TAK1KO_mTAK1-smRNA`   <- do_sleuth_modelling(`so_TAK1KO_mTAK1-smRNA`)

a <- sleuth_results(`so_WT_naive`, "treatmentIL1B") %>%
    mutate(group = "WT_naive",
           condition = "WT",
           rescue = "naive",
           treatment = "IL1B",
           time = 2)

b <- sleuth_results(`so_sham_naive`, "treatmentIL1B") %>%
    mutate(group = "sham_naive",
           condition = "sham",
           rescue = "naive",
           treatment = "IL1B",
           time = 2)

c <- sleuth_results(`so_TAK1KO_naive`, "treatmentIL1B") %>%
    mutate(group = "TAK1KO_naive",
           condition = "TAK1KO",
           rescue = "naive",
           treatment = "IL1B",
           time = 2)

d <- sleuth_results(`so_TAK1KO_lipid`, "treatmentIL1B") %>%
    mutate(group = "TAK1KO_lipid",
           condition = "TAK1KO",
           rescue = "lipid",
           treatment = "IL1B",
           time = 2)

e <- sleuth_results(`so_TAK1KO_TAK1-smRNA`, "treatmentIL1B") %>%
    mutate(group = "TAK1KO_TAK1-smRNA",
           condition = "TAK1KO",
           rescue = "TAK1-smRNA",
           treatment = "IL1B",
           time = 2)

f <- sleuth_results(`so_TAK1KO_mTAK1-smRNA`, "treatmentIL1B") %>%
    mutate(group = "TAK1KO_mTAK1-smRNA",
           condition = "TAK1KO",
           rescue = "mTAK1-smRNA",
           treatment = "IL1B",
           time = 2)

dea_clean = rbind(a, b, c, d, e, f) %>%
    select(Gene = target_id, group, condition, rescue, treatment, time, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(
        time = factor(time, levels = c(2)),
        treatment = factor(treatment, levels = c("IL1B")),
        condition = factor(condition, levels = c("WT", "sham", "TAK1KO")),
        rescue = factor(rescue, levels = c("naive", "lipid", "TAK1-smRNA", "mTAK1-smRNA")),
        group = factor(group, levels = grouplevels)
    ) %>%
    arrange(Gene, time, group) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(is.na(log2fold))) %>%
    ungroup() %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR)
        )

write_tsv(dea_clean, "sleuth/a549-tak1-ko-dea.txt")
```





# Interaction testing to ask specific questions

### === prep sleuth inputs ===

```r
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/a549-tak1-ko")

library(dplyr)
library(readr)
library(sleuth)
library(tidyr)

t2g <- read_tsv(
    "/work/newton_lab/ag_analysis/ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt"
    ) %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

grouplevels = c(
    "WT_naive",
    "sham_naive",
    "TAK1KO_naive",
    "TAK1KO_lipid",
    "TAK1KO_TAK1-smRNA",
    "TAK1KO_mTAK1-smRNA"
)

s2c = read_tsv("tak1-ko-meta.txt") %>%
    mutate(condition = factor(condition, levels = c("WT", "sham", "TAK1KO"))) %>%
    mutate(treatment = factor(treatment, levels = c("NS", "IL1B"))) %>%
    mutate(rescue = factor(rescue, levels = c("naive", "lipid", "TAK1-smRNA", "mTAK1-smRNA"))) %>%

    mutate(
        group = paste0(condition, "_", rescue),
        group = factor(
            group,
            levels = grouplevels
        )
    ) %>%
    select(sample, rep, time, group, condition, rescue, treatment, path = kallisto_path)

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}
```

### Q1B: effect of sham on WT IL1B response?

```r

s2c_q1b = s2c %>%
    filter(group %in% c("WT_naive", "sham_naive"))

so_q1b = sleuth_prep(
    sample_to_covariates = s2c_q1b,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q1b = sleuth_fit(
    so_q1b,
    ~treatment*condition
)

so_q1b = sleuth_wt(
    so_q1b,
    "treatmentIL1B:conditionsham"
)

results_q1b = sleuth_results(so_q1b, test = "treatmentIL1B:conditionsham") %>%
    as_tibble() %>%
    mutate(
        contrast = "treatmentIL1B:conditionsham"
    ) %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, contrast, log2diff = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2diff = if_else(is.na(log2diff), 0, log2diff),
        FDR = if_else(is.na(FDR), 1, FDR)
        ) %>%
    arrange(Gene)

write_tsv(results_q1b, "sleuth/sham_IL1B_interaction_test.txt")
```



### Q2B: effect of knockout on WT IL1B response?

```r

s2c_q2b = s2c %>%
    filter(group %in% c("WT_naive", "TAK1KO_naive"))

so_q2b = sleuth_prep(
    sample_to_covariates = s2c_q2b,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q2b = sleuth_fit(
    so_q2b,
    ~treatment*condition
)

so_q2b = sleuth_wt(
    so_q2b,
    "treatmentIL1B:conditionTAK1KO"
)

results_q2b = sleuth_results(so_q2b, test = "treatmentIL1B:conditionTAK1KO") %>%
    as_tibble() %>%
    mutate(
        contrast = "treatmentIL1B:conditionTAK1KO"
    ) %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, contrast, log2diff = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2diff = if_else(is.na(log2diff), 0, log2diff),
        FDR = if_else(is.na(FDR), 1, FDR)
        ) %>%
    arrange(Gene)

write_tsv(results_q2b, "sleuth/KO_IL1B_interaction_test.txt")
```

### Q4B: effect of introducing native TAK1 smRNA to TAK1KO?

```r

s2c_q4b = s2c %>%
    filter(group %in% c("TAK1KO_naive", "TAK1KO_TAK1-smRNA")) %>%
    mutate(
        rescue = factor(rescue, levels = c("naive", "TAK1-smRNA"))
    )

so_q4b = sleuth_prep(
    sample_to_covariates = s2c_q4b,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q4b = sleuth_fit(
    so_q4b,
    ~treatment*rescue
)

so_q4b = sleuth_wt(
    so_q4b,
    "treatmentIL1B:rescueTAK1-smRNA"
)

results_q4b = sleuth_results(so_q4b, test = "treatmentIL1B:rescueTAK1-smRNA") %>%
    as_tibble() %>%
    mutate(
        contrast = "treatmentIL1B:rescueTAK1-smRNA"
    ) %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, contrast, log2diff = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2diff = if_else(is.na(log2diff), 0, log2diff),
        FDR = if_else(is.na(FDR), 1, FDR)
        ) %>%
    arrange(Gene)

write_tsv(results_q4b, "sleuth/nativeRescue_interaction_test.txt")
```

### Q4C: compare WT and rescued IL1B response?

```r
s2c_q4c = s2c %>%
    filter(group %in% c("WT_naive", "TAK1KO_TAK1-smRNA"))

so_q4c = sleuth_prep(
    sample_to_covariates = s2c_q4c,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q4c = sleuth_fit(
    so_q4c,
    ~treatment*group
)

so_q4c = sleuth_wt(
    so_q4c,
    "treatmentIL1B:groupTAK1KO_TAK1-smRNA"
)

results_q4c = sleuth_results(so_q4c, test = "treatmentIL1B:groupTAK1KO_TAK1-smRNA") %>%
    as_tibble() %>%
    mutate(
        contrast = "treatmentIL1B:groupTAK1KO_TAK1-smRNA"
    ) %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, contrast, log2diff = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2diff = if_else(is.na(log2diff), 0, log2diff),
        FDR = if_else(is.na(FDR), 1, FDR)
        ) %>%
    arrange(Gene)

write_tsv(results_q4c, "sleuth/wt_vs_rescued_interaction_test.txt")

```


# Tests to address questions regarding clonal process

### === prep sleuth inputs ===

```r
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/a549-tak1-ko")

library(dplyr)
library(readr)
library(sleuth)
library(tidyr)

t2g <- read_tsv(
    "/work/newton_lab/ag_analysis/ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt"
    ) %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

grouplevels = c(
    "WT_naive",
    "sham_naive",
    "TAK1KO_naive",
    "TAK1KO_lipid",
    "TAK1KO_TAK1-smRNA",
    "TAK1KO_mTAK1-smRNA"
)

s2c = read_tsv("tak1-ko-meta.txt") %>%
    mutate(condition = factor(condition, levels = c("WT", "sham", "TAK1KO"))) %>%
    mutate(treatment = factor(treatment, levels = c("NS", "IL1B"))) %>%
    mutate(rescue = factor(rescue, levels = c("naive", "lipid", "TAK1-smRNA", "mTAK1-smRNA"))) %>%

    mutate(
        group = paste0(condition, "_", rescue),
        group = factor(
            group,
            levels = grouplevels
        )
    ) %>%
    select(sample, rep, time, group, condition, rescue, treatment, path = kallisto_path)

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}
```

### clonalQ1: effect of group on NS

```r
s2c_q1 = s2c %>%
    filter(treatment == "NS")

so_q1 = sleuth_prep(
    sample_to_covariates = s2c_q1,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q1 = sleuth_fit(
    so_q1,
    ~group
) %>%
    sleuth_wt(., "groupsham_naive") %>%
    sleuth_wt(., "groupTAK1KO_naive") %>%
    sleuth_wt(., "groupTAK1KO_lipid") %>%
    sleuth_wt(., "groupTAK1KO_TAK1-smRNA") %>%
    sleuth_wt(., "groupTAK1KO_mTAK1-smRNA")

results = data.frame()

for(test in c(
    "groupsham_naive",
    "groupTAK1KO_naive",
    "groupTAK1KO_lipid",
    "groupTAK1KO_TAK1-smRNA",
    "groupTAK1KO_mTAK1-smRNA"
)){

    gr = gsub("group", "", test)

    df = sleuth_results(so_q1, test) %>%
    mutate(group = gr, treatment = "NS")

    results = rbind(results, df)

}

results_clean = results %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, treatment, group, log2fold = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR)
    ) %>%
    mutate(group = factor(
        group,
        levels = c(
            "sham_naive",
            "TAK1KO_naive",
            "TAK1KO_lipid",
            "TAK1KO_TAK1-smRNA",
            "TAK1KO_mTAK1-smRNA"
        )
    )) %>%
    arrange(Gene, group)

write_tsv(results_clean, "sleuth/effect-of-group-on-NS_dea.txt")
```

### clonalQ2: effect of group on IL1B

```r
s2c_q2 = s2c %>%
    filter(treatment == "IL1B")

so_q2 = sleuth_prep(
    sample_to_covariates = s2c_q2,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_q2 = sleuth_fit(
    so_q2,
    ~group
) %>%
    sleuth_wt(., "groupsham_naive") %>%
    sleuth_wt(., "groupTAK1KO_naive") %>%
    sleuth_wt(., "groupTAK1KO_lipid") %>%
    sleuth_wt(., "groupTAK1KO_TAK1-smRNA") %>%
    sleuth_wt(., "groupTAK1KO_mTAK1-smRNA")

results = data.frame()

for(test in c(
    "groupsham_naive",
    "groupTAK1KO_naive",
    "groupTAK1KO_lipid",
    "groupTAK1KO_TAK1-smRNA",
    "groupTAK1KO_mTAK1-smRNA"
)){

    gr = gsub("group", "", test)

    df = sleuth_results(so_q2, test) %>%
    mutate(group = gr, treatment = "IL1B")

    results = rbind(results, df)

}

results_clean = results %>%
    mutate(b = b/log(2)) %>%
    select(Gene = target_id, treatment, group, log2fold = b, FDR = qval) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR)
    ) %>%
    mutate(group = factor(
        group,
        levels = c(
            "sham_naive",
            "TAK1KO_naive",
            "TAK1KO_lipid",
            "TAK1KO_TAK1-smRNA",
            "TAK1KO_mTAK1-smRNA"
        )
    )) %>%
    arrange(Gene, group)

write_tsv(results_clean, "sleuth/effect-of-group-on-IL1B_dea.txt")
```



