
setwd("F:/r-env/中期/方法/模拟数据")



library(tidyverse)
library(magrittr)
library(M3JF)
library(data.table)



nclust <- 3 # 类别数量


n_byClust_norm <- c(130, 200, 170) # 类别分布
n_byClust_beta <- c(175, 135, 190) # 类别数量
n_byClust_bin <- c(100, 175, 225)

noises <- 1
p_noise <- 0.005

feature_nums <- c(1000, 1000, 1000)
props <- c(0.04, 0.02, 0.01)

# 高
# means <- c(1.2, 1, 0.8)
# sds <- c(0.6, 0.5, 0.4)

# 中
means <- c(1, 0.8, 0.6)
sds <- c(0.5, 0.4, 0.3)
# 
# means <- c(0.8, 0.6, 0.4)
# sds <- c(0.4, 0.3, 0.2)

params_norms <- mapply(function(m, sd) {
  c(mean = m, sd = sd)
}, means, sds, SIMPLIFY = FALSE)


# params_beta <- list(
#   c(mean1= -1.4, mean2 = 1.4, sd1 = 0.7, sd2 = 0.7),
#   c(mean1= -1.2, mean2 = 1.2, sd1 =0.6, sd2 = 0.6),
#   c(mean1= -1, mean2 = 1, sd1 = 0.5, sd2 = 0.5)
# )

params_beta <- list(
  c(mean1 = -1, mean2 = 1, sd1 = 0.5, sd2 = 0.5),
  c(mean1 = -0.8, mean2 = 0.8, sd1 = 0.4, sd2 = 0.4),
  c(mean1 = -0.6, mean2 = 0.6, sd1 =0.3, sd2 = 0.3)
)
# 
# params_beta <- list(
#   c(mean1 = -0.8, mean2 = 0.8, sd1 = 0.4, sd2 = 0.4),
#   c(mean1 = -0.6, mean2 = 0.6, sd1 = 0.3, sd2 = 0.3),
#   c(mean1 = -0.4, mean2 = 0.4, sd1 =0.2, sd2 = 0.2)
# )

# params_bin <- list(c(p = 0.1),c(p = 0.08),c(p = 0.06))
params_bin <- list(c(p = 0.08),c(p = 0.06),c(p = 0.04))
# params_bin <- list(c(p = 0.06),c(p = 0.04),c(p = 0.02))


# save(realdata_snv_sub, file = "realdata_snv_sub-高信号比例.rda")

params_norms
params_beta
params_bin


# 循环开始
for (i in 3:25) {
  
  # 生成文件夹路径
  dataset_scenario1 <- str_glue(
    "F:/r-env/中期/方法/模拟试验-模拟数据/情景1/k-3/信号比例-高/信号水平-中/sim{i}/"
  )
  dataset_scenario2 <- str_glue(
    "F:/r-env/中期/方法/模拟试验-模拟数据/情景2/k-3/信号比例-高/信号水平-中/sim{i}/"
  )
  
  # 创建目录（如果不存在）
  # dir.create(dataset_scenario1, recursive = TRUE, showWarnings = FALSE)
  # dir.create(dataset_scenario2, recursive = TRUE, showWarnings = FALSE)
  
  # 文件名
  s1_mrna_batch <- str_glue("s1-k3-high-mid-batch{i}-mrna.csv")
  s1_meth_batch <- str_glue("s1-k3-high-mid-batch{i}-meth.csv")
  s1_mutate_batch <- str_glue("s1-k3-high-mid-batch{i}-mutate.csv")
  s2_mrna_batch <- str_glue("s2-k3-high-mid-batch{i}-mrna.csv")
  s2_meth_batch <- str_glue("s2-k3-high-mid-batch{i}-meth.csv")
  s2_mutate_batch <- str_glue("s2-k3-high-mid-batch{i}-mutate.csv")
  
  # 正态分布数据
  sim_norm <- simulateY(
    nclust = nclust,        # 类别数
    n_byClust = n_byClust_norm,  # 类别分布
    J = feature_nums[1],    # 特征数量
    flavor = "normal",      # 分布类型
    prop = props[2],        # 标签相关变量（阳性信号）的比例
    params = params_norms,  # 不同标签类的分布参数
    noise = noises          # 噪声的比例?? noise背景噪声的标准差 包括标签相关变量和无关变量
  )
  
  simdata_mrna <- sim_norm$data %>% as.data.frame() %>% 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  fwrite(simdata_mrna, str_c(dataset_scenario1, s1_mrna_batch), col.names = FALSE)
  
  # SVD变换后的情景2-mrna
  svd_sim_mrna <- svd(t(simdata_mrna))
  
  v_sim_mrna <- svd_sim_mrna$v %>% t
  
  fusion_mrna <- u_mrna %*% s_mrna %*% v_sim_mrna
  
  fusion_mrna %<>% t %>% as.data.frame() %>% 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  fwrite(fusion_mrna, str_c(dataset_scenario2, s2_mrna_batch), col.names = FALSE)
  
  # Beta分布数据
  sim_beta <- simulateY(
    nclust = nclust, 
    n_byClust = n_byClust_beta, 
    J = feature_nums[2], 
    flavor = "beta", 
    params = params_beta,
    prop = props[2], 
    noise = noises
  )
  
  simdata_meth <- sim_beta$data %>% as.data.frame() %>% 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  
  
  fwrite(simdata_meth, str_c(dataset_scenario1, s1_meth_batch), col.names = FALSE)
  
  svd_sim_meth <- svd(t(simdata_meth))
  
  v_sim_meth <- svd_sim_meth$v %>% t
  
  fusion_meth <- u_meth %*% s_meth %*% v_sim_meth
  
  fusion_meth %<>% t %>% as.data.frame() %>% 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  fusion_meth %<>% 
    mutate(across(where(is.numeric), ~ (.x - min(.x)) / (max(.x) - min(.x))))
  
  fwrite(fusion_meth,str_c(dataset_scenario2, s2_meth_batch), col.names = F)
  
  
  # 二项分布数据
  sim_bin <- simulateY(
    nclust = nclust, 
    n_byClust = n_byClust_bin, 
    J = feature_nums[3], 
    flavor = "binary", 
    params = params_bin, 
    prop = props[1], 
    noise = p_noise
  )
  simdata_mut <- sim_bin$data %>% as.data.frame()
  fwrite(simdata_mut, str_c(dataset_scenario1, s1_mutate_batch), col.names = FALSE)
  
  # fusion_mutate：打乱列顺序
  label_related_feat_mutate <- unlist(sim_bin[["positive"]])
  
  fusion_mutate <- simdata_mut[, label_related_feat_mutate] %>% 
    bind_cols(realdata_snv_sub)
  
  # 生成随机列顺序
  shuffled_columns <- sample(1000)
  
  fusion_mutate <- fusion_mutate[, shuffled_columns] 
  
  fwrite(fusion_mutate, str_c(dataset_scenario2, s2_mutate_batch), col.names = FALSE)
}


label <- tibble(label = c(rep(0,150), rep(1,150), rep(2,200)))

fwrite(label,"F:/r-env/中期/方法/模拟试验-模拟数据/情景1/k-3/label-k3.csv", col.names = F)

