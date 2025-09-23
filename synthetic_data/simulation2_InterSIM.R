library(data.table)
library(tidyverse)


library(InterSIM)

set.seed(2023110400)
# 低信号 ---------------------------------------------------------------------
# 设置参数
prop <- c(0.36, 0.12, 0.18, 0.10, 0.24)

effect <- 0.6
effect2 <- 0.6
effect3 <- 0.3


for (i in 1:25) {
  # 生成模拟数据
  sim.data <- InterSIM(
    n.sample = 500,
    cluster.sample.prop = prop,
    delta.methyl = effect,
    delta.expr = effect2,
    delta.protein = effect3,
    p.DMP = 0.08,
    p.DEG = NULL,
    p.DEP = NULL,
    sigma.methyl = NULL,
    sigma.expr = NULL,
    sigma.protein = NULL,
    cor.methyl.expr = NULL,
    cor.expr.protein = NULL,
    do.plot = F,
    sample.cluster = TRUE,
    feature.cluster = TRUE
  )
  
  sim_meth <- sim.data[["dat.methyl"]] %>% round(., 4) %>% as_tibble()
  sim_exp <- sim.data[["dat.expr"]] %>% round(., 4) %>% as_tibble()
  sim_pro <- sim.data[["dat.protein"]] %>% round(., 4) %>% as_tibble()
  
  sim_label <- sim.data[["clustering.assignment"]][["cluster.id"]] %>% 
    as_tibble() %>% 
    `colnames<-`("label")
  
  sim_label <- sim_label - 1
  
  dataset_path <- str_glue(
    "F:/r-env/大课题/模拟数据-intersim/情景1/k-3/信号比例-高/信号水平-低/sim{i}/"
  )
  
  dir.create(dataset_path, recursive = TRUE, showWarnings = FALSE)
  
  path_mrna <- str_glue("s1-k3-high-low-batch{i}-mrna.csv")
  path_meth <- str_glue("s1-k3-high-low-batch{i}-meth.csv")
  path_protein <- str_glue("s1-k3-high-low-batch{i}-prot.csv")
  path_label <- str_glue("s1-k3-high-low-batch{i}-label.csv")
  
  fwrite(sim_exp, str_c(dataset_path, path_mrna))
  fwrite(sim_meth, str_c(dataset_path, path_meth))
  fwrite(sim_pro, str_c(dataset_path, path_protein))
  fwrite(sim_label, str_c(dataset_path, path_label))
  
  # 可选：打印进度
  cat("已完成第", i, "次模拟\n")
}

# 高信号 ---------------------------------------------------------------------
# 设置参数
prop <- c(0.36, 0.12, 0.18, 0.10, 0.24)

effect <- 0.8
effect2 <- 0.8
effect3 <- 0.4


for (i in 1:25) {
  # 生成模拟数据
  sim.data <- InterSIM(
    n.sample = 500,
    cluster.sample.prop = prop,
    delta.methyl = effect,
    delta.expr = effect2,
    delta.protein = effect3,
    p.DMP = 0.12,
    p.DEG = NULL,
    p.DEP = NULL,
    sigma.methyl = NULL,
    sigma.expr = NULL,
    sigma.protein = NULL,
    cor.methyl.expr = NULL,
    cor.expr.protein = NULL,
    do.plot = F,
    sample.cluster = TRUE,
    feature.cluster = TRUE
  )
  
  sim_meth <- sim.data[["dat.methyl"]] %>% round(., 4) %>% as_tibble()
  sim_exp <- sim.data[["dat.expr"]] %>% round(., 4) %>% as_tibble()
  sim_pro <- sim.data[["dat.protein"]] %>% round(., 4) %>% as_tibble()
  
  sim_label <- sim.data[["clustering.assignment"]][["cluster.id"]] %>% 
    as_tibble() %>% 
    `colnames<-`("label")
  
  sim_label <- sim_label - 1
  
  dataset_path <- str_glue(
    "F:/r-env/大课题/模拟数据-intersim/情景1/k-3/信号比例-高/信号水平-高/sim{i}/"
  )
  
  dir.create(dataset_path, recursive = TRUE, showWarnings = FALSE)
  
  path_mrna <- str_glue("s1-k3-high-high-batch{i}-mrna.csv")
  path_meth <- str_glue("s1-k3-high-high-batch{i}-meth.csv")
  path_protein <- str_glue("s1-k3-high-high-batch{i}-prot.csv")
  path_label <- str_glue("s1-k3-high-high-batch{i}-label.csv")
  
  fwrite(sim_exp, str_c(dataset_path, path_mrna))
  fwrite(sim_meth, str_c(dataset_path, path_meth))
  fwrite(sim_pro, str_c(dataset_path, path_protein))
  fwrite(sim_label, str_c(dataset_path, path_label))
  
  # 可选：打印进度
  cat("已完成第", i, "次模拟\n")
}

