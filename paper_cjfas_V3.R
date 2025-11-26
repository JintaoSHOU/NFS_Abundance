library(mgcv)
library(AICcmodavg)
library(dplyr)
library(stringr)
library(ggplot2)
library(Hmisc)
library(car)
library(readxl)
library(patchwork)
library(forecast)

CPUE1<-read_excel("AnnualData.xlsx",sheet = "NFS")
CPUE1<-ts(CPUE1$CPUE,start=1994,end=2024,frequency = 1)
CPUE1<-window(CPUE1,start=1995,end=2022)
plot(CPUE1)


CPUE1%>%tseries::adf.test()
CPUE1%>%diff()%>%tseries::adf.test()
CPUE1%>%diff()%>%diff()%>%tseries::adf.test()

CPUE1%>%tseries::kpss.test()
CPUE1%>%diff()%>%tseries::kpss.test()

CPUE1%>%diff()%>%diff()%>%ggtsdisplay()

########################################################ACF PACF
{
  png("ACF_PACF_CPUE1.png", width = 10, height = 8, units = "in", res = 300)
  # 设置图形布局为三行两列，调整边距和标签大小
  par(mfrow = c(3, 2), mar = c(5, 5, 4, 2), cex.lab = 1.8)
  
  # 第一行：d = 0
  acf(CPUE1, main = "", ylab = "ACF")
  pacf(CPUE1, main = "", ylab = "PACF")
  mtext("d = 0", side = 3, outer = FALSE, line = 1, cex = 2, adj = -0.25)
  
  # 第二行：d = 1
  diff1 <- diff(CPUE1)
  acf(diff1, main = "", ylab = "ACF")
  pacf(diff1, main = "", ylab = "PACF")
  mtext("d = 1", side = 3, outer = FALSE, line = 1, cex = 2, adj = -0.25)
  
  # 第三行：d = 2
  diff2 <- diff(diff1)
  acf(diff2, main = "", ylab = "ACF")
  pacf(diff2, main = "", ylab = "PACF")
  mtext("d = 2", side = 3, outer = FALSE, line = 1, cex = 2, adj = -0.25)
  
  # 恢复默认布局
  par(mfrow = c(1, 1))
  dev.off()
}

################################################################################### 定义待评估的 ARIMA 模型组合
{
  orders <- list(
    c(0, 0, 1), c(1, 0, 0), c(1, 0, 1),
    c(0, 1, 0), c(0, 1, 1), c(1, 1, 0), c(1, 1, 1),
    c(0, 1, 2), c(2, 1, 0), c(1, 1, 2), c(2, 1, 1),
    c(0, 2, 0), c(0, 2, 1), c(1, 2, 0), c(0, 2, 2),
    c(2, 2, 0), c(1, 2, 2), c(2, 2, 1)
  )
  
  # 初始化结果数据框
  results <- data.frame(
    order = character(),
    mean_AICc = numeric(),
    mean_AdjR2 = numeric(),
    mean_LB_stat = numeric(),
    mean_LB_pvalue = numeric(),
    tsCV_RMSE = numeric(),
    tsCV_MdAE = numeric(),
    tsCV_MASE = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (ord in orders) {
    aic_values <- c()
    aicc_values <- c()
    adjr2_values <- c()
    lb_statistics <- c()
    lb_pvalues <- c()
    
    for (i in 18:(length(CPUE1) - 1)) {
      subdata <- CPUE1[1:i]
      fit <- tryCatch(Arima(subdata, order = ord), error = function(e) NULL)
      if (is.null(fit)) next
      
      res <- residuals(fit)
      meany <- mean(subdata, na.rm = TRUE)
      r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
      n <- length(subdata)
      k <- length(fit$coef)
      if (n <= k + 1) next
      
      adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
      aic_val <- AIC(fit)
      aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
      
      adjr2_values <- c(adjr2_values, adjR2)
      aic_values <- c(aic_values, aic_val)
      aicc_values <- c(aicc_values, aicc_val)
      
      lb <- tryCatch(Box.test(res, type = "Ljung-Box"), error = function(e) NULL)
      if (!is.null(lb)) {
        lb_statistics <- c(lb_statistics, lb$statistic)
        lb_pvalues <- c(lb_pvalues, lb$p.value)
      }
    }
    
    # tsCV RMSE、MdAE 和 MASE
    fun_ts <- function(x, h) forecast::forecast(Arima(x, order = ord), h = h)
    e <- tryCatch(tsCV(CPUE1, fun_ts, h = 1, initial = 17), error = function(e) rep(NA, length(CPUE1)))
    
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mdae <- median(abs(e), na.rm = TRUE)
    
    # MASE：分子是预测误差的 MAE，分母是 naive 模型的 MAE
    mae_model <- mean(abs(e), na.rm = TRUE)
    naive_errors <- diff(CPUE1)  # y_t - y_{t-1}
    mae_naive <- mean(abs(naive_errors), na.rm = TRUE)
    mase <- mae_model / mae_naive
    
    # 存储结果
    results <- rbind(results, data.frame(
      order = paste(ord, collapse = "-"),
      mean_AICc = mean(aicc_values, na.rm = TRUE),
      mean_AdjR2 = mean(adjr2_values, na.rm = TRUE),
      mean_LB_stat = mean(lb_statistics, na.rm = TRUE),
      mean_LB_pvalue = mean(lb_pvalues, na.rm = TRUE),
      tsCV_RMSE = rmse,
      tsCV_MdAE = mdae,
      tsCV_MASE = mase
    ))
  }
  
  results <- results[order(results$tsCV_RMSE), ]
  print(results)
}

##############################################################################################ARIMA  auto 单变量
{
  funauto <- function(x, h, order){
    forecast::forecast(auto.arima(x ), h=h)
  }
  e <- forecast::tsCV(CPUE1, funauto, h=1,initial=17)
  e
  cvauto <- c(ME=mean(e, na.rm=TRUE), RMSE=sqrt(mean(e^2, na.rm=TRUE)), MAE=mean(abs(e), na.rm=TRUE),MdAE = median(abs(e), na.rm = TRUE))
  cvauto#RMSE=342.18
  
  ###
  {
    aic_values <- c()
    aicc_values <- c()
    adjr2_values <- c()
    lb_statistics <- c()  # 新增向量，用于存储每次循环中 Ljung-Box 检验统计量的值
    lb_pvalues <- c()
    for (i in 18:27) {
      # 每次取不同的数据子集进行拟合
      subdata <-CPUE1[1:i]
      fit <- auto.arima(subdata)
      res <- resid(fit)
      meany <- mean(subdata, na.rm = TRUE)
      r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
      n <- length(subdata)
      k <- length(fit$coef)
      adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
      adjr2_values <- c(adjr2_values, adjR2)
      aic_value <- AIC(fit)
      aic_values <- c(aic_values, aic_value)
      # 通过 AIC 计算 AICc
      aicc_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
      aicc_values <- c(aicc_values, aicc_value)
      
      # 进行 Ljung-Box 检验并获取统计量，存储到 lb_statistics 向量中
      lb_result <- Box.test(res, type = "Ljung-Box")
      lb_statistics <- c(lb_statistics, lb_result$statistic)
      lb_pvalues <- c(lb_pvalues, lb_result$p.value)
    }
    
    # 计算平均 AIC、AICc、Adjusted R² 和平均 Ljung-Box 检验统计量
    mean_aic <- mean(aic_values)
    mean_aicc <- mean(aicc_values)
    mean_adjr2 <- mean(adjr2_values)
    mean_lb_statistic <- mean(lb_statistics)  # 计算平均 Ljung-Box 检验统计量
    mean_lb_pvalue <- mean(lb_pvalues) 
    
    cat("平均 AICc 值：", mean_aicc, "\n")
    cat("平均 Adjusted R² 值：", mean_adjr2, "\n")
    cat("平均 Ljung-Box 检验统计量值：", mean_lb_statistic, "\n")
    cat("平均 Ljung-Box 检验P值：", mean_lb_pvalue, "\n")
  }
  
}

##########################################################GAM  t-1  单变量
{
  CPUE <- as.numeric(CPUE1)  # 确保是数值向量
  n <- length(CPUE)
  h <- 1
  initial <- 18
  
  # 存储测试误差
  errors <- rep(NA, n)
  
  # 用于保存每次模型的评估指标
  metrics <- data.frame(
    time = integer(),
    Adj_R2 = numeric(),
    RMSE = numeric(),
    AICc = numeric()
  )
  
  for (i in initial:(n - h)) {
    train_xt <- CPUE[1:i]
    train_xt_lag <- CPUE[1:(i - 1)]
    
    train_data <- data.frame(
      xt = train_xt[2:length(train_xt)],
      xt_lag1 = train_xt_lag
    )
    
    # 拟合 GAM 模型
    gam_model <- gam(xt ~ s(xt_lag1), data = train_data,method = "REML")
    
    # 预测第 i+1 个值（测试集上的预测）
    new_data <- data.frame(xt_lag1 = CPUE[i])
    pred <- predict(gam_model, newdata = new_data)
    
    # 测试误差
    test_error <- CPUE[i + 1] - pred
    errors[i + 1] <- test_error
    
    # 计算 RMSE（单步）
    rmse_i <- sqrt(test_error^2)
    
    # 计算 Adjusted R²
    r2 <- summary(gam_model)$r.sq
    edf <- sum(summary(gam_model)$edf)  # 有效自由度
    n_train <- nrow(train_data)
    adj_r2 <- 1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
    
    # 计算 AICc
    aicc_val <- AICc(gam_model)
    
    # 保存指标
    metrics <- rbind(metrics, data.frame(
      time = i,
      Adj_R2 = adj_r2,
      RMSE = rmse_i,
      AICc = aicc_val
    ))
  }
  
  # 汇总平均值
  metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")])
  
  # 计算总体 RMSE（测试集所有误差的RMSE）
  overall_rmse <- sqrt(mean(na.omit(errors)^2))
  overall_mdae <- median(abs(na.omit(errors)))
  
  
  # 输出结果
  cat("===== 汇总指标 =====\n")
  cat("平均 Adj R²       ：", round(metrics_mean["Adj_R2"], 4), "\n")
  cat("平均 AICc         ：", round(metrics_mean["AICc"], 4), "\n")
  cat("总体 RMSE（测试集）：", round(overall_rmse, 4), "\n")
  cat("总体 MdAE（测试集）：", round(overall_mdae, 4), "\n")
  
}
##########################################################GAM  t-2  单变量
{
  library(mgcv)
  library(AICcmodavg)
  
  CPUE <- as.numeric(CPUE1)
  n <- length(CPUE)
  h <- 1
  initial <- 18  # 为了保证有 lag1 和 lag2
  
  # 存储测试误差
  errors <- rep(NA, n)
  
  # 用于保存每次模型的评估指标
  metrics <- data.frame(
    time = integer(),
    Adj_R2 = numeric(),
    RMSE = numeric(),
    AICc = numeric()
  )
  
  for (i in initial:(n - h)) {
    train_xt <- CPUE[1:i]
    train_xt_lag1 <- CPUE[1:(i - 1)]
    train_xt_lag2 <- CPUE[1:(i - 2)]
    
    # 构建训练数据（从第3个值开始，保证 lag1、lag2 都有值）
    train_data <- data.frame(
      xt = train_xt[3:length(train_xt)],
      xt_lag1 = train_xt_lag1[2:length(train_xt_lag1)],
      xt_lag2 = train_xt_lag2
    )
    
    # 拟合 GAM 模型（包括两个滞后项）
    gam_model <- gam(xt ~ s(xt_lag1) + s(xt_lag2), data = train_data)
    
    # 构建测试输入（当前 t 使用 CPUE[i]，t-1 使用 CPUE[i-1]）
    new_data <- data.frame(
      xt_lag1 = CPUE[i],
      xt_lag2 = CPUE[i - 1]
    )
    
    # 预测
    pred <- predict(gam_model, newdata = new_data)
    
    # 测试误差
    test_error <- CPUE[i + 1] - pred
    errors[i + 1] <- test_error
    rmse_i <- sqrt(test_error^2)
    
    # Adjusted R²
    r2 <- summary(gam_model)$r.sq
    edf <- sum(summary(gam_model)$edf)
    n_train <- nrow(train_data)
    adj_r2 <- 1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
    
    # AICc
    aicc_val <- AICc(gam_model)
    
    # 保存结果
    metrics <- rbind(metrics, data.frame(
      time = i,
      Adj_R2 = adj_r2,
      RMSE = rmse_i,
      AICc = aicc_val
    ))
  }
  
  # 汇总平均值
  metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")])
  
  # 总体 RMSE
  overall_rmse <- sqrt(mean(na.omit(errors)^2))
  overall_mdae <- median(abs(na.omit(errors)))
  # 输出结果
  cat("===== GAM（含 Y[t-1], Y[t-2]）汇总指标 =====\n")
  cat("平均 Adj R²       ：", round(metrics_mean["Adj_R2"], 4), "\n")
  cat("平均 AICc         ：", round(metrics_mean["AICc"], 4), "\n")
  cat("平均 RMSE（逐步） ：", round(metrics_mean["RMSE"], 4), "\n")
  cat("总体 RMSE（测试集）：", round(overall_rmse, 4), "\n")
  cat("总体 MdAE（测试集）：", round(overall_mdae, 4), "\n")
  
}

##########
#########################################################################加变量 初步筛选一个,arima1.0.0
{
  data <- read.csv("NFS_env.csv")
  data <- data[, -c(1, 2)]  # 删除前两列
  data_scaled <- scale(data)  # 标准化
  
  # 如果你想保留为 data 对象本身，可以这样写：
  data <- as.data.frame(scale(data))
  
  # 初始化结果列表
  results_list <- list()
  
  for (i in 1:ncol(data)) {
    env1 <- data[, i]
    
    # 滑动预测 RMSE
    fun <- function(x, h, xreg, newxreg) {
      fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
      forecast(fit, h = h, xreg = newxreg)
    }
    e <- tsCV(CPUE1, fun, h = 1, xreg = env1, initial = 17)
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mdae <- median(abs(e), na.rm = TRUE)
    
    # AICc、Adjusted R²、R² 计算
    aic_values <- c()
    aicc_values <- c()
    adjr2_values <- c()
    r2_values <- c()
    
    for (j in 18:27) {
      subdata <- CPUE1[1:j]
      subxreg <- env1[1:j]
      
      fit <- tryCatch({
        Arima(subdata, order = c(1, 0, 0), xreg = subxreg)
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        res <- resid(fit)
        meany <- mean(subdata, na.rm = TRUE)
        r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
        n <- length(subdata)
        k <- length(fit$coef)
        
        adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
        aic_value <- AIC(fit)
        aicc_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
        
        aic_values <- c(aic_values, aic_value)
        aicc_values <- c(aicc_values, aicc_value)
        adjr2_values <- c(adjr2_values, adjR2)
        r2_values <- c(r2_values, r2)
      }
    }
    
    mean_aic <- mean(aic_values)
    mean_aicc <- mean(aicc_values)
    mean_adjr2 <- mean(adjr2_values)
    mean_r2 <- mean(r2_values)
    
    results_list[[i]] <- data.frame(
      Variable = colnames(data)[i],
      RMSE = rmse,
      MdAE = mdae,  # 新增 MdAE
      Mean_AICc = mean_aicc,
      Mean_Adjusted_R2 = mean_adjr2,
      Mean_R2 = mean_r2
    )
    
  }
  
  # 合并结果并排序
  final_results <- bind_rows(results_list)
  final_results <- final_results[order(final_results$RMSE), ]
  
  # 保存
  write.csv(final_results, "one_variable_arima1_0_0.csv", row.names = FALSE)
  print(final_results)
  
}

########################################################################识别不同生活史阶段  画图 ARIMA 1.0.0
{
  final_results<-read.csv("one_variable_arima1_0_0.csv")
  
  final_results <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 查看并保存结果
  print(final_results)
  
  # 先基于条件筛选数据
  baseline_rmse <- 69.36143 
  baseline_r2 <- 0.19953578  
  baseline_AICc <- 268.9233
  baseline_mdae <- 55.97373
  
  final_results_filtered <- final_results %>%
    filter(MdAE <= baseline_mdae)%>%
    filter(RMSE <= baseline_rmse)%>%
    # filter(Mean_AICc <= baseline_AICc)%>%
    filter(Mean_Adjusted_R2 >= baseline_r2)
  
  
  # 添加 life_stage 列
  library(dplyr)
  library(stringr)
  
  final_results_filtered <- final_results_filtered %>%
    mutate(
      stage_letter = str_extract(Variable, "_[EPAJS](?=_)|_[EPAJS]$") %>% str_remove("_"),
      life_stage = case_when(
        stage_letter == "E" ~ "Egg",
        stage_letter == "P" ~ "Paralarvae",
        stage_letter == "J" ~ "Juvenile",
        stage_letter == "S" ~ "Subadult",
        stage_letter == "A" ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  
  # 添加 rmse_percent 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      rmse_percent = (baseline_rmse - RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered$life_stage <- factor(final_results_filtered$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered <- final_results_filtered %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
  
  final_results_filtered <- final_results_filtered %>%
    filter(rmse_percent > 0.1)
  
  ave<-mean(final_results_filtered$RMSE)
  
  # 绘图
  p <- ggplot(final_results_filtered, aes(x = rmse_percent, y = reorder(Variable, rmse_percent), fill = Type)) +
    geom_col() +
    geom_text(aes(label = Variable), 
              hjust = 0, vjust = 0.5, size = 4, color = "black", nudge_x = 0.1) +
    facet_wrap(~life_stage, ncol = 3, nrow = 2, scales = "free_y") +
    labs(x = "RMSE Reduction (%)", y = NULL, title = NULL, fill = "Covariate Type") +
    theme_test() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    scale_fill_manual(values = c("Food availability" = "#8491B4B2", "Environmental condition" = "#FFC0CB")) +
    scale_x_continuous(limits = c(0, 23))  # x轴范围
  
  p
  # 保存
  ggsave("life_stage_rmse.png", p, width = 14, height = 8, dpi = 500)
  
}

########################################################################### 加变量 初步筛选一个,arima0.1.2
{
  
  data <- read.csv("NFS_env.csv")
  data <- data[, -c(1,2)]
  # 初始化结果列表
  results_list <- list()
  
  for (i in 1:ncol(data)) {
    env1 <- data[, i]
    
    # 滑动预测 RMSE
    fun <- function(x, h, xreg, newxreg) {
      fit <- Arima(x, order = c(0, 2, 1), xreg = xreg)
      forecast(fit, h = h, xreg = newxreg)
    }
    e <- tsCV(CPUE1, fun, h = 1, xreg = env1, initial = 21)
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mdae <- median(abs(e), na.rm = TRUE)
    
    # AICc、Adjusted R²、R² 计算
    aic_values <- c()
    aicc_values <- c()
    adjr2_values <- c()
    r2_values <- c()
    
    for (j in 22:27) {
      subdata <- CPUE1[1:j]
      subxreg <- env1[1:j]
      
      fit <- tryCatch({
        Arima(subdata, order = c(0, 2, 1), xreg = subxreg)
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        res <- resid(fit)
        meany <- mean(subdata, na.rm = TRUE)
        r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
        n <- length(subdata)
        k <- length(fit$coef)
        
        adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
        aic_value <- AIC(fit)
        aicc_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
        
        aic_values <- c(aic_values, aic_value)
        aicc_values <- c(aicc_values, aicc_value)
        adjr2_values <- c(adjr2_values, adjR2)
        r2_values <- c(r2_values, r2)
      }
    }
    
    mean_aic <- mean(aic_values)
    mean_aicc <- mean(aicc_values)
    mean_adjr2 <- mean(adjr2_values)
    mean_r2 <- mean(r2_values)
    
    results_list[[i]] <- data.frame(
      Variable = colnames(data)[i],
      RMSE = rmse,
      MdAE = mdae,  # 新增 MdAE
      Mean_AICc = mean_aicc,
      Mean_Adjusted_R2 = mean_adjr2,
      Mean_R2 = mean_r2
    )
    
  }
  
  # 合并结果并排序
  final_results <- bind_rows(results_list)
  final_results <- final_results[order(final_results$RMSE), ]
  
  # 保存
  write.csv(final_results, "one_variable_arima0_1_2.csv", row.names = FALSE)
  print(final_results)
  
}
###########################################################################识别不同生活史阶段  画图 ARIMA 0.1.2
{
  final_results<-read.csv("one_variable_arima0_1_2.csv")
  
  final_results <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 查看并保存结果
  print(final_results)
  
  # 先基于条件筛选数据
  #baseline_rmse <- 62.64065 
  #baseline_r2 <- 0.22552 
  #baseline_AICc <- 280.6484
  #baseline_mdae <- 26.78345
  
  #ARIMA(1,2,2)
  baseline_rmse <- 47.56 
  baseline_r2 <- 0.1053 
  baseline_AICc <- 277.4117
  baseline_mdae <- 40.4365
  
  #ARIMA(2,2,1)
  baseline_rmse <- 58.506 
  baseline_r2 <- 0.1879
  baseline_AICc <- 276.2653
  baseline_mdae <- 15.387
  
  #ARIMA(2,1,1)
  baseline_rmse <- 66.72 
  baseline_r2 <- 0.2164
  baseline_AICc <- 282.6037
  baseline_mdae <- 25.269
  
  #ARIMA(1,1,1)
  baseline_rmse <- 49.286 
  baseline_r2 <- 0.1456
  baseline_AICc <- 282.84
  baseline_mdae <- 40.62
  
  #ARIMA(2,1,0)
  baseline_rmse <- 62.64065 
  baseline_r2 <- 0.2255
  baseline_AICc <- 280.64
  baseline_mdae <- 26.78
  
  #ARIMA(0,2,1)
  baseline_rmse <- 54.86
  baseline_r2 <- 0.063
  baseline_AICc <- 275.53
  baseline_mdae <- 40.277
  
  
  
  final_results_filtered <- final_results %>%
    #filter(MdAE <= baseline_mdae)%>%
    filter(RMSE <= baseline_rmse)%>%
    #filter(Mean_AICc <= baseline_AICc)
    filter(Mean_Adjusted_R2 >= baseline_r2)
  
  
  # 添加 life_stage 列
  library(dplyr)
  library(stringr)
  
  final_results_filtered <- final_results_filtered %>%
    mutate(
      stage_letter = str_extract(Variable, "_[EPAJS](?=_)|_[EPAJS]$") %>% str_remove("_"),
      life_stage = case_when(
        stage_letter == "E" ~ "Egg",
        stage_letter == "P" ~ "Paralarvae",
        stage_letter == "J" ~ "Juvenile",
        stage_letter == "S" ~ "Subadult",
        stage_letter == "A" ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  
  # 添加 rmse_percent 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      rmse_percent = (baseline_rmse - RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered$life_stage <- factor(final_results_filtered$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered <- final_results_filtered %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
  
  # 绘图
  p <- ggplot(final_results_filtered, aes(x = rmse_percent, y = reorder(Variable, rmse_percent), fill = Type)) +
    geom_col() +
    geom_text(aes(label = Variable), 
              hjust = 0, vjust = 0.5, size = 4, color = "black", nudge_x = 0.1) +
    facet_wrap(~life_stage, ncol = 3, nrow = 2, scales = "free_y") +
    labs(x = "RMSE Reduction (%)", y = NULL, title = NULL, fill = "Covariate Type") +
    theme_test() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    scale_fill_manual(values = c("Food availability" = "#8491B4B2", "Environmental condition" = "#FFC0CB")) +
    scale_x_continuous(limits = c(0, 18))  # x轴范围
  
  p
  # 保存
  ggsave("life_stage_rmse.png", p, width = 14, height = 8, dpi = 500)
  
}

##############################################################################加变量 初步筛选一个,auto.arima
{
  
  data <- read.csv("NFS_env.csv")
  data <- data[, -c(1, 2)]
  
  # 添加所有变量的滞后项（Xt-1），并重命名
  lagged_data <- data[-1, ]
  colnames(lagged_data) <- paste0(colnames(data), "_lag1")
  
  # CPUE 对应对齐
  CPUE_lag <- CPUE1[-1]
  
  # 原始环境变量对齐
  data_aligned <- data[-28, ]
  
  # 合并原始变量 + 滞后变量
  full_data <- cbind(data_aligned, lagged_data)
  
  results_list <- list()
  
  for (i in 1:ncol(full_data)) {
    
    env1 <- data[, i]
    
    # 先做滑动预测，计算 RMSE
    fun <- function(x, h, xreg, newxreg){
      fit <- auto.arima(x, xreg = xreg)
      forecast::forecast(fit, h = h, xreg = newxreg)
    }
    e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = env1, initial = 21 )
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mdae <- median(abs(e), na.rm = TRUE)
    
    # 再计算 AICc、Adjusted R2、R2
    aic_values <- c()
    aicc_values <- c()
    adjr2_values <- c()
    r2_values <- c()
    
    for (j in 22:27) {
      subdata <- CPUE1[1:j]
      subxreg <- env1[1:j]
      
      fit <- auto.arima(subdata, xreg = subxreg)
      res <- resid(fit)
      
      meany <- mean(subdata, na.rm = TRUE)
      r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
      
      n <- length(subdata)
      k <- length(fit$coef)
      
      adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
      
      aic_value <- AIC(fit)
      aicc_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
      
      aic_values <- c(aic_values, aic_value)
      aicc_values <- c(aicc_values, aicc_value)
      adjr2_values <- c(adjr2_values, adjR2)
      r2_values <- c(r2_values, r2)
    }
    
    mean_aic <- mean(aic_values)
    mean_aicc <- mean(aicc_values)
    mean_adjr2 <- mean(adjr2_values)
    mean_r2 <- mean(r2_values)
    
    # 保存这一列的结果
    results_list[[i]] <- data.frame(
      Variable = colnames(full_data)[i],
      RMSE = rmse,
      MdAE = mdae,
      Mean_AICc = mean_aicc,
      Mean_Adjusted_R2 = mean_adjr2,
      Mean_R2 = mean_r2
    )
  }
  
  final_results <- bind_rows(results_list)
  
  # 查看结果
  print(final_results)
  
  write.csv(final_results,file = "one_variable_auto.arima_withlag.csv")
  
}
###########################################################################识别不同生活史阶段  画图 auto.arima
{
  
  final_results<-read.csv("one_variable_auto.arima_withlag.csv")
  
  
  final_results <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 先基于条件筛选数据
  baseline_rmse <- 76.86919
  baseline_r2 <- 0.3309528
  baseline_mdae <- 71.80055
  
  final_results_filtered <- final_results %>%
    filter(RMSE <= baseline_rmse) %>% 
    filter(MdAE <= baseline_mdae) %>% 
    filter(Mean_Adjusted_R2 >= baseline_r2)              
  
  # 添加 life_stage 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      stage_num = as.numeric(str_sub(Variable, -1, -1)),
      life_stage = case_when(
        stage_num == 1 ~ "Egg",
        stage_num == 2 ~ "Paralarvae",
        stage_num == 3 ~ "Juvenile",
        stage_num == 4 ~ "Subadult",
        stage_num == 5 ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  # 添加 rmse_percent 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      rmse_percent = (baseline_rmse - RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered$life_stage <- factor(final_results_filtered$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered <- final_results_filtered %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
  
  # 绘图
  p <- ggplot(final_results_filtered, aes(x = rmse_percent, y = reorder(Variable, rmse_percent), fill = Type)) +
    geom_col() +
    geom_text(aes(label = Variable), 
              hjust = 0, vjust = 0.5, size = 4, color = "black", nudge_x = 0.1) +
    facet_wrap(~life_stage, ncol = 3, nrow = 2, scales = "free_y") +
    labs(x = "RMSE Reduction (%)", y = NULL, title = NULL, fill = "Covariate Type") +
    theme_test() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    scale_fill_manual(values = c("Food availability" = "#8491B4B2", "Environmental condition" = "#FFC0CB")) +
    scale_x_continuous(limits = c(0, 27))  # x轴范围
  
  p
  # 保存
  ggsave("life_stage_rmse.png", p, width = 14, height = 8, dpi = 500)
  
}

#######################################################################################加变量 初步筛选一个,GAM
{
  data <- read.csv("NFS_env.csv")
  data <- data[, -c(1, 2)]  # 删除前两列
  data_scaled <- scale(data)  # 标准化
  
  # 如果你想保留为 data 对象本身，可以这样写：
  data <- as.data.frame(scale(data))
  
  CPUE <- as.numeric(CPUE1)  # 确保是数值向量
  n <- length(CPUE)
  h <- 1
  initial <- 18
  
  results_list <- list()
  
  for (i in 1:ncol(data)) {
    env <- data[, i]
    errors <- rep(NA, n)
    
    metrics <- data.frame(
      time = integer(),
      Adj_R2 = numeric(),
      RMSE = numeric(),
      AICc = numeric()
    )
    
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- env[2:j]
      
      train_data <- data.frame(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        env = train_env
      )
      
      # GAM 带协变量
      gam_model <- gam(xt ~ s(xt_lag1) + s(env), data = train_data,method = "REML")
      
      # 预测第 j+1 个值
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        env = env[j + 1]
      )
      pred <- predict(gam_model, newdata = new_data)
      
      test_error <- CPUE[j + 1] - pred
      errors[j + 1] <- test_error
      rmse_j <- sqrt(test_error^2)
      
      r2 <- summary(gam_model)$r.sq
      edf <- sum(summary(gam_model)$edf)
      n_train <- nrow(train_data)
      adj_r2 <- 1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
      
      aicc_val <- AICc(gam_model)
      
      metrics <- rbind(metrics, data.frame(
        time = j,
        Adj_R2 = adj_r2,
        RMSE = rmse_j,
        AICc = aicc_val
      ))
    }
    
    # 平均指标
    metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")])
    overall_rmse <- sqrt(mean(na.omit(errors)^2))
    overall_mdae <- median(abs(na.omit(errors)))
    
    # 保存结果
    results_list[[i]] <- data.frame(
      Variable = colnames(data)[i],
      Mean_Adj_R2 = metrics_mean["Adj_R2"],
      Mean_AICc = metrics_mean["AICc"],
      Mean_Step_RMSE = metrics_mean["RMSE"],
      Overall_RMSE = overall_rmse,
      Overall_MdAE = overall_mdae 
    )
  }
  
  # 汇总所有变量结果
  final_results <- bind_rows(results_list)
  final_results <- final_results[order(final_results$Overall_RMSE), ]
  # 查看并保存结果
  print(final_results)
  write.csv(final_results, file = "gam_env_one_variable3.csv", row.names = FALSE)
  xdfg
}
###############################################################################识别不同生活史阶段  画图 GAM 18
{
  final_results<-read.csv("gam_env_one_variable3.csv")
  
  final_results <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 查看并保存结果
  print(final_results)
  # 先基于条件筛选数据
  baseline_rmse <- 67.4674
  baseline_r2 <- 0.1981
  baseline_AICc <-258.6851 
  baseline_mdae <- 53.8388 
  
  final_results_filtered <- final_results %>%
    filter(Overall_RMSE <= baseline_rmse) %>%
    filter(Overall_MdAE <= baseline_mdae) %>%
    #filter(Mean_AICc <= baseline_AICc)%>%
    filter(Mean_Adj_R2 >= baseline_r2)             
  
  # 添加 life_stage 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      stage_letter = str_extract(Variable, "_[EPAJS](?=_)|_[EPAJS]$") %>% str_remove("_"),
      life_stage = case_when(
        stage_letter == "E" ~ "Egg",
        stage_letter == "P" ~ "Paralarvae",
        stage_letter == "J" ~ "Juvenile",
        stage_letter == "S" ~ "Subadult",
        stage_letter == "A" ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  # 添加 rmse_percent 列
  final_results_filtered <- final_results_filtered %>%
    mutate(
      rmse_percent = (baseline_rmse - Overall_RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered$life_stage <- factor(final_results_filtered$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered <- final_results_filtered %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
  
  # 绘图
  p <- ggplot(final_results_filtered, aes(x = rmse_percent, y = reorder(Variable, rmse_percent), fill = Type)) +
    geom_col() +
    geom_text(aes(label = Variable), 
              hjust = 0, vjust = 0.5, size = 4, color = "black", nudge_x = 0.1) +
    facet_wrap(~life_stage, ncol = 3, nrow = 2, scales = "free_y") +
    labs(x = "RMSE Reduction (%)", y = NULL, title = NULL, fill = "Covariate Type") +
    theme_test() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    scale_fill_manual(values = c("Food availability" = "#8491B4B2", "Environmental condition" = "#FFC0CB")) +
    scale_x_continuous(limits = c(0, 100))  # x轴范围
  
  p
  # 保存
  ggsave("life_stage_rmse.png", p, width = 14, height = 8, dpi = 500)
  
  
}


##############################################################################################ARIMA-GAM 比较
{
  library(ggplot2)
  data1<-read.csv("one_variable_arima1_0_0.csv")
  data2<-read.csv("gam_env_one_variable3.csv")
  
  # 创建统一列名的副本
  data1_renamed <- data1 %>%
    select(Mean_Adjusted_R2, RMSE) %>%
    rename(Adj_R2 = Mean_Adjusted_R2, RMSE = RMSE) %>%
    mutate(source = "data1")
  
  data2_renamed <- data2 %>%
    select(Mean_Adj_R2, Overall_RMSE) %>%
    rename(Adj_R2 = Mean_Adj_R2, RMSE = Overall_RMSE) %>%
    mutate(source = "data2")
  
  # 合并数据
  combined_data <- bind_rows(data1_renamed, data2_renamed)
  
  p<-ggplot(combined_data, aes(x = Adj_R2, y = RMSE)) +
    geom_hline(yintercept = 69.36143, color = "#1B4F72", linetype = "solid", size = 0.6) +
    geom_hline(yintercept = 67.4674 , color = "#E67E22", linetype = "solid", size = 0.6) +
    geom_vline(xintercept = 0.19953578 , color = "#1B4F72", linetype = "solid", size = 0.6) +
    geom_vline(xintercept = 0.1981 , color = "#E67E22", linetype = "solid", size = 0.6) +
    # 第一层：实心圆
    geom_point(
      aes(color = source),
      shape = 16,
      size = 3.5
    ) +
    # 第二层：空心圆叠加在实心圆上
    geom_point(
      aes(color = source),
      shape = 21,
      size = 3,          # 稍小一点
      fill = NA,         # 空心
      stroke = 0.05,
      color = alpha("white", 0.6)  # 设置透明度为 60%
    ) +
    scale_color_manual(
      values = c("data1" = "#1B4F72", "data2" = "#E67E22"),
      labels = c("data1" = "ARIMA", "data2" = "GAM")
    ) +
    guides(color = guide_legend(
      override.aes = list(
        shape = 21,
        size = 3.5,
        stroke = 0.3,
        fill = c("#1B4F72", "#E67E22"),      # 与数据点填充一致
        color = alpha("white", 0.6)          # 空心边框半透明白色
      )
    )) +
    labs(
      x = "Adjusted R²",
      y = "tsCV-RMSE",
      color = "Model"
    ) +
    theme_bw(base_size = 18)
  
  p
  ggsave("ARIMA_GAM.png", p , width = 9, height =7,dpi = 400)
  
}

######################找到共同选中的变量，并画图
{

  
{
  final_results<-read.csv("one_variable_arima1_0_0.csv")
  
  final_results_arima <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 先基于条件筛选数据
  baseline_rmse <- 69.36143 
  baseline_r2 <- 0.19953578  
  baseline_AICc <- 268.9233
  baseline_mdae <- 55.97373
  
  final_results_filtered_arima <- final_results_arima %>%
    filter(MdAE <= baseline_mdae)%>%
    filter(RMSE <= baseline_rmse)%>%
    # filter(Mean_AICc <= baseline_AICc)%>%
    filter(Mean_Adjusted_R2 >= baseline_r2)
  
  
  # 添加 life_stage 列
  library(dplyr)
  library(stringr)
  
  final_results_filtered_arima <- final_results_filtered_arima %>%
    mutate(
      stage_letter = str_extract(Variable, "_[EPAJS](?=_)|_[EPAJS]$") %>% str_remove("_"),
      life_stage = case_when(
        stage_letter == "E" ~ "Egg",
        stage_letter == "P" ~ "Paralarvae",
        stage_letter == "J" ~ "Juvenile",
        stage_letter == "S" ~ "Subadult",
        stage_letter == "A" ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  
  # 添加 rmse_percent 列
  final_results_filtered_arima <- final_results_filtered_arima %>%
    mutate(
      rmse_percent = (baseline_rmse - RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered_arima$life_stage <- factor(final_results_filtered_arima$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered_arima <- final_results_filtered_arima %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
}


#######
{
  final_results<-read.csv("gam_env_one_variable3.csv")
  
  final_results_gam <- final_results %>%
    mutate(Type= if_else(
      str_detect(Variable, "^Chl|^Nppv|^Phyc"),
      "Food availability",
      "Environmental condition"
    ))
  
  # 查看并保存结果
  print(final_results)
  # 先基于条件筛选数据
  baseline_rmse <- 67.4674
  baseline_r2 <- 0.1981
  baseline_AICc <-258.6851 
  baseline_mdae <- 53.8388 
  
  final_results_filtered_gam <- final_results_gam %>%
    filter(Overall_RMSE <= baseline_rmse) %>%
    filter(Overall_MdAE <= baseline_mdae) %>%
    #filter(Mean_AICc <= baseline_AICc)%>%
    filter(Mean_Adj_R2 >= baseline_r2)             
  
  # 添加 life_stage 列
  final_results_filtered_gam <- final_results_filtered_gam %>%
    mutate(
      stage_letter = str_extract(Variable, "_[EPAJS](?=_)|_[EPAJS]$") %>% str_remove("_"),
      life_stage = case_when(
        stage_letter == "E" ~ "Egg",
        stage_letter == "P" ~ "Paralarvae",
        stage_letter == "J" ~ "Juvenile",
        stage_letter == "S" ~ "Subadult",
        stage_letter == "A" ~ "Adult",
        TRUE ~ NA_character_
      )
    )
  
  # 添加 rmse_percent 列
  final_results_filtered_gam <- final_results_filtered_gam %>%
    mutate(
      rmse_percent = (baseline_rmse - Overall_RMSE) / baseline_rmse * 100
    )
  
  
  # 指定life_stage顺序
  final_results_filtered_gam$life_stage <- factor(final_results_filtered_gam$life_stage, 
                                              levels = c("Egg", "Paralarvae", "Juvenile", "Subadult", "Adult"))
  
  # 排序
  final_results_filtered_gam <- final_results_filtered_gam %>%
    group_by(life_stage) %>%
    arrange(desc(rmse_percent), .by_group = TRUE) %>%
    ungroup()
  
}
# 1️⃣ 找出共同选中的变量
common_vars <- intersect(final_results_filtered_arima$Variable,
                         final_results_filtered_gam$Variable)

# 2️⃣ 分别筛选出共同变量的结果
common_results_arima <- final_results_filtered_arima %>%
  filter(Variable %in% common_vars) %>%
  select(Variable, rmse_percent, life_stage, Type) %>%
  rename(rmse_percent_arima = rmse_percent)

common_results_gam <- final_results_filtered_gam %>%
  filter(Variable %in% common_vars) %>%
  select(Variable, rmse_percent, life_stage, Type) %>%
  rename(rmse_percent_gam = rmse_percent)

# 3️⃣ 合并两个模型结果
common_results_combined <- left_join(common_results_arima, common_results_gam, by = c("Variable", "life_stage", "Type"))

# 4️⃣ 计算两个模型的 rmse_percent 平均值
common_results_combined <- common_results_combined %>%
  mutate(rmse_percent_mean = rowMeans(select(., rmse_percent_arima, rmse_percent_gam), na.rm = TRUE))

# 5️⃣ 按 life_stage 排序（可选）
common_results_combined <- common_results_combined %>%
  arrange(life_stage, desc(rmse_percent_mean))

# 查看结果
print(common_results_combined)

final_results_filtered<-common_results_combined


{
  # 添加新水平
  levels(final_results_filtered$life_stage) <- c(levels(final_results_filtered$life_stage), "Embryo")
  # 替换 Egg 为 Embryo
  final_results_filtered$life_stage[final_results_filtered$life_stage == "Egg"] <- "Embryo"
  # 示例颜色（浅色系）
  stage_colors <- c(
    "Embryo" = "#CBAACB",        # pastel purple
    "Paralarvae" = "#AFCBFF", # pastel blue
    "Juvenile" = "#B5EAD7",   # pastel mint green
    "Subadult" = "#E2F0CB",   # pastel green
    "Adult" = "#FFDAC1"       # pastel peach
  )
  
  # 按你希望的顺序设置 life_stage 的因子水平
  final_results_filtered$life_stage <- factor(
    final_results_filtered$life_stage,
    levels = c("Embryo", "Paralarvae", "Juvenile", "Subadult", "Adult")  # 你想要的顺序
  )
  
  # 定义形状符号
  type_shapes <- c("Food availability" = 8, "Environmental condition" = 4)
  
  # 绘图
  p <- ggplot(final_results_filtered, aes(
    x = rmse_percent_mean,
    y = reorder(Variable, rmse_percent_mean),
    fill = life_stage
  )) +
    geom_col(width = 0.8) +
    
    # 添加符号：根据变量类型区分形状
    geom_point(
      aes(shape = Type),
      x = final_results_filtered$rmse_percent_mean,
      size = 4,
      fill = "gray30",
      color = "gray30",
      stroke = 1.1
    ) +
    
    # 设置颜色和形状
    scale_shape_manual(values = type_shapes) +
    scale_fill_manual(values = stage_colors) +
    
    scale_x_continuous(limits = c(0, 18)) +
    
    labs(
      x = "tsCV-RMSE Reduction (%)",
      y = NULL,
      fill = "Life Stage",
      shape = "Covariate Type"
    ) +
    guides(
      fill = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    theme_bw() +
    theme(
      #plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # 缩小左边距
      axis.text.y = element_text(size = 12, color = "black"),  # 显示变量名
      axis.ticks.length = unit(3, "pt"),  # 3pt 长度，朝外
      #axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)
    )
  p
  
  ggsave("selected_variable_ARIMA_gam.png",p,width =12,height =7, dpi=400)
}



}

######共线性
{

col <- final_results_filtered[[1]]  # 或者 final_results_filtered$Variable

select_enc <- data[, col]

a <- Hmisc::redun(~ ., data = select_enc, nk = 0)

out_vars <- a$Out

select_enc_reduced <- select_enc[, !colnames(select_enc) %in% out_vars]


library(corrplot)

X <- as.matrix(select_enc)
png("correlation_plot_arima_gam.png", width = 7, height = 6, units = "in", res = 500)
corrplot::corrplot(
  cor(X),
  type = "lower",
  method = "color",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  number.cex = 0.6,
  col = colorRampPalette(c("#D73027", "#FEE090", "#91BFDB", "#4575B4"))(200),
  addgrid.col = rgb(255/255, 255/255, 255/255, alpha = 0.2),  # 白色网格线，透明度为0.2
  tl.cex = 0.7
)
dev.off()






library(car)
library(dplyr)

# 1️⃣ 从结果中提取出变量名
selected_vars <- unique(final_results_filtered$Variable)
print(selected_vars)

# 2️⃣ 在原始数据（假设名为 data）中提取这些变量
# ⚠️ 请确认 data 是包含所有环境变量的完整数据框
select_enc_full <- data %>% select(all_of(selected_vars))

##没去除
full <- lm(CPUE1 ~ ., data=select_enc_full)
car::vif(full)
olsrr::ols_vif_tol(full)

##去除了
full <- lm(CPUE1 ~ ., data=select_enc_reduced)
car::vif(full)
olsrr::ols_vif_tol(full)






}

###########正方形图
{

fitA0<-Arima(cpue, order = c(1, 0, 0))

predicted <- fitted(fitA0)  # ARIMA 模型拟合值
observed <- cpue           # 实际值
df_plot <- data.frame(
  Observed = observed,
  Predicted = predicted
)
max_val <- max(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
min_val <- min(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
# 4️⃣ 绘图
library(grid)

pA0 <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by ARIMA model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=69.36, Adj R²=0.2\nNo env covariates",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(A0)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")

pA0



xreg <- select_enc_reduced[[1]]
fitA1<-Arima(cpue, order = c(1, 0, 0), xreg = xreg)

predicted <- fitted(fitA1)  # ARIMA 模型拟合值
observed <- cpue           # 实际值
df_plot <- data.frame(
  Observed = observed,
  Predicted = predicted
)
max_val <- max(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
min_val <- min(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
# 4️⃣ 绘图
pA1 <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by ARIMA model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=63.77, Adj R²=0.234\nWith PDO_E",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(A2)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


pA1

#
xreg <- select_enc_reduced[[2]]
fitA2<-Arima(cpue, order = c(1, 0, 0), xreg = xreg)

predicted <- fitted(fitA2)  # ARIMA 模型拟合值
observed <- cpue           # 实际值
df_plot <- data.frame(
  Observed = observed,
  Predicted = predicted
)
max_val <- max(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
min_val <- min(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
# 4️⃣ 绘图
pA2 <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by ARIMA model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=65.97, Adj R²=0.208\nWith O2_0_P_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(A3)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


pA2

#
xreg <- select_enc_reduced[[3]]
fitA3<-Arima(cpue, order = c(1, 0, 0), xreg = xreg)

predicted <- fitted(fitA3)  # ARIMA 模型拟合值
observed <- cpue           # 实际值
df_plot <- data.frame(
  Observed = observed,
  Predicted = predicted
)
max_val <- max(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
min_val <- min(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
# 4️⃣ 绘图
pA3 <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by ARIMA model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=69.35, Adj R²=0.28\nWith Temp_200_A_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(A1)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


pA3

#
xreg <- select_enc_reduced[[4]]
fitA4<-Arima(cpue, order = c(1, 0, 0), xreg = xreg)

predicted <- fitted(fitA4)  # ARIMA 模型拟合值
observed <- cpue           # 实际值
df_plot <- data.frame(
  Observed = observed,
  Predicted = predicted
)
max_val <- max(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
min_val <- min(c(df_plot$Observed, df_plot$Predicted), na.rm = TRUE)
# 4️⃣ 绘图
pA4 <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by ARIMA model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=68.2, Adj R²=0.208\nWith Chl_300_A_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(A4)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


pA4


#######################
library(mgcv)
library(ggplot2)

cpue <- CPUE1                      # 因变量 (ts)
env_data <- select_enc_reduced     # 自变量矩阵（但这次不使用）
 # 当前环境变量名称（仅用于标题或标识）
cat("当前拟合变量（未使用）:", var_name, "\n")

# 构造训练数据（只使用自回归项）
train_xt <- cpue[2:28]
train_xt_lag <- cpue[1:27]

train_data <- data.frame(
  xt = train_xt,
  xt_lag1 = train_xt_lag
)

# 去除 NA
train_data <- na.omit(train_data)

# ❌ 不包含协变量，只建模 s(xt_lag1)
gam_formula <- as.formula("xt ~ s(xt_lag1)")

# 拟合 GAM 模型
gam_model <- gam(gam_formula, data = train_data, method = "REML")

# 获取拟合值
train_data$Predicted <- predict(gam_model)
train_data$Observed  <- train_data$xt

# 拟合可视化（时间序列对比）
plot(train_data$Observed, type = "l", col = "black", lwd = 2, ylab = "CPUE", xlab = "Time")
lines(train_data$Predicted, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), col = c("black", "red"), lty = 1)

# 绘制 Observed vs Predicted 图
min_val <- min(train_data$Observed, train_data$Predicted)
max_val <- max(train_data$Observed, train_data$Predicted)

pG0 <- ggplot(train_data, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by GAM model" 
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=67.47, Adj R²=0.198\nNo env covariates",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(G0)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


###
library(mgcv)

cpue <- CPUE1                      # 因变量 (ts)
env_data <- select_enc_reduced     # 自变量矩阵
var_name <- colnames(env_data)[1]  # 取第一个变量
cat("当前拟合变量:", var_name, "\n")

# 提取该变量
xreg <- env_data[[var_name]]

# 构造训练数据
train_xt <- cpue[2:28]
train_xt_lag <- cpue[1:27]
train_env <- xreg[2:28]

# 合并为数据框（记得给train_env命名）
train_data <- data.frame(
  xt = train_xt,
  xt_lag1 = train_xt_lag,
  env = train_env   # 命名为 env
)
train_data <- na.omit(train_data)
gam_formula <- as.formula("xt ~ s(xt_lag1) + s(env)")
gam_model <- gam(gam_formula, data = train_data, method = "REML")
fitted_values <- fitted(gam_model)
print(fitted_values)
plot(train_data$xt, type = "l", col = "black", lwd = 2, ylab = "CPUE", xlab = "Time")
lines(fitted_values, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), col = c("black", "red"), lty = 1)
train_data$Predicted <- predict(gam_model)
train_data$Observed  <- train_data$xt

min_val <- min(train_data$Observed, train_data$Predicted)
max_val <- max(train_data$Observed, train_data$Predicted)

pG1<-ggplot(train_data, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by GAM model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=56.5, Adj R²=0.285\nWith PDO_E",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(G2)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")

##
library(mgcv)

cpue <- CPUE1                      # 因变量 (ts)
env_data <- select_enc_reduced     # 自变量矩阵
var_name <- colnames(env_data)[2]  # 取第一个变量
cat("当前拟合变量:", var_name, "\n")

# 提取该变量
xreg <- env_data[[var_name]]

# 构造训练数据
train_xt <- cpue[2:28]
train_xt_lag <- cpue[1:27]
train_env <- xreg[2:28]

# 合并为数据框（记得给train_env命名）
train_data <- data.frame(
  xt = train_xt,
  xt_lag1 = train_xt_lag,
  env = train_env   # 命名为 env
)
train_data <- na.omit(train_data)
gam_formula <- as.formula("xt ~ s(xt_lag1) + s(env)")
gam_model <- gam(gam_formula, data = train_data, method = "REML")
fitted_values <- fitted(gam_model)
print(fitted_values)
plot(train_data$xt, type = "l", col = "black", lwd = 2, ylab = "CPUE", xlab = "Time")
lines(fitted_values, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), col = c("black", "red"), lty = 1)
train_data$Predicted <- predict(gam_model)
train_data$Observed  <- train_data$xt

min_val <- min(train_data$Observed, train_data$Predicted)
max_val <- max(train_data$Observed, train_data$Predicted)

pG2<-ggplot(train_data, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by GAM model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=55, Adj R²=0.25\nWith O2_0_P_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(G3)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")

##
library(mgcv)

cpue <- CPUE1                      # 因变量 (ts)
env_data <- select_enc_reduced     # 自变量矩阵
var_name <- colnames(env_data)[3]  # 取第一个变量
cat("当前拟合变量:", var_name, "\n")

# 提取该变量
xreg <- env_data[[var_name]]

# 构造训练数据
train_xt <- cpue[2:28]
train_xt_lag <- cpue[1:27]
train_env <- xreg[2:28]

# 合并为数据框（记得给train_env命名）
train_data <- data.frame(
  xt = train_xt,
  xt_lag1 = train_xt_lag,
  env = train_env   # 命名为 env
)
train_data <- na.omit(train_data)
gam_formula <- as.formula("xt ~ s(xt_lag1) + s(env)")
gam_model <- gam(gam_formula, data = train_data, method = "REML")
fitted_values <- fitted(gam_model)
print(fitted_values)
plot(train_data$xt, type = "l", col = "black", lwd = 2, ylab = "CPUE", xlab = "Time")
lines(fitted_values, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), col = c("black", "red"), lty = 1)
train_data$Predicted <- predict(gam_model)
train_data$Observed  <- train_data$xt

min_val <- min(train_data$Observed, train_data$Predicted)
max_val <- max(train_data$Observed, train_data$Predicted)

pG3<-ggplot(train_data, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by GAM model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=56.1, Adj R²=0.399\nWith Temp_200_A_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(G1)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")


##
library(mgcv)

cpue <- CPUE1                      # 因变量 (ts)
env_data <- select_enc_reduced     # 自变量矩阵
var_name <- colnames(env_data)[4]  # 取第一个变量
cat("当前拟合变量:", var_name, "\n")

# 提取该变量
xreg <- env_data[[var_name]]

# 构造训练数据
train_xt <- cpue[2:28]
train_xt_lag <- cpue[1:27]
train_env <- xreg[2:28]

# 合并为数据框（记得给train_env命名）
train_data <- data.frame(
  xt = train_xt,
  xt_lag1 = train_xt_lag,
  env = train_env   # 命名为 env
)
train_data <- na.omit(train_data)
gam_formula <- as.formula("xt ~ s(xt_lag1) + s(env)")
gam_model <- gam(gam_formula, data = train_data, method = "REML")
fitted_values <- fitted(gam_model)
print(fitted_values)
plot(train_data$xt, type = "l", col = "black", lwd = 2, ylab = "CPUE", xlab = "Time")
lines(fitted_values, col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), col = c("black", "red"), lty = 1)
train_data$Predicted <- predict(gam_model)
train_data$Observed  <- train_data$xt

min_val <- min(train_data$Observed, train_data$Predicted)
max_val <- max(train_data$Observed, train_data$Predicted)

pG4<-ggplot(train_data, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, color = "black", fill = "white", size = 3, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +
  coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE by GAM model"
  ) +
  theme_test(base_size = 16) +
  annotation_custom(
    grid::textGrob(
      label = "RMSE=57.3, Adj R²=0.308\nWith Chl_300_A_lag",
      x = unit(0.06, "npc"),   # 相对坐标横向位置
      y = unit(0.9, "npc"),   # 相对坐标纵向位置
      hjust = 0,              # ⭐ 左对齐
      gp = gpar(fontsize = 16)
    )
  )+
  coord_cartesian(clip="off", ylim=c(100,500),xlim=c(100,500))+
  annotate("text", x =30, y = 520, label = "(G4)", 
           vjust = 0, hjust =0, size = 7, fontface = "bold")



################################################################################################
combined_plot <- (pA0 + pA3) / (pG0 + pG3) +
  plot_annotation(title = "") 

ggsave("Observe_predic.png",combined_plot,width = 12, height =12,dpi = 500)


combined_plot <- (pA1 + pA2+ pA4) / (pG1 + pG2+pG4) +
  plot_annotation(title = "") 

ggsave("Observe_predic_附件.png",combined_plot,width = 18, height =12,dpi = 500)

}

############变量响应曲线
{
  library(mgcv)
  library(dplyr)
  library(ggplot2)
  library(gratia)
  library(grid)
  
  CPUE <- as.numeric(CPUE1)
  n <- length(CPUE)
  h <- 1
  initial <- 18
  
  plot_list <- list()
  
  for (i in 1:10) {
    df <- select_enc_reduced[, i, drop = FALSE]
    var_name <- colnames(df)
    
    if (grepl("O2", var_name)) {
      plot_title <- paste0(var_name, " (mmol/m³)")
    } else if (grepl("Temp", var_name)) {
      plot_title <- paste0(var_name, " (℃)")
    } else if (grepl("Chl", var_name)) {
      plot_title <- paste0(var_name, " (mg/m³)")
    } else if (grepl("PDO", var_name)) {
      plot_title <- paste0(var_name, " (m)")
    } else {
      plot_title <- var_name
    }
    
    
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    
    sm <- smooth_estimates(gam_model) |> 
      dplyr::filter(smooth == smooths(gam_model)[2])
    
    resid_data <- data.frame(
      x = gam_model$model[[var_name]],
      .resid = resid(gam_model)
    )
    
    p <- ggplot(sm, aes(x = !!sym(var_name), y = est)) +
      geom_hline(yintercept = 0, color = "grey50", linetype = "solid", size = 1) +
      geom_ribbon(aes(ymin = est - 2 * se, ymax = est + 2 * se),
                  fill = "grey80", alpha = 0.5) +
      geom_line(size = 1.2, color = "orange") +
      geom_point(data = resid_data, aes(x = x, y = .resid),
                 color = "black", size = 3) +
      geom_rug(data = resid_data,
               aes(x = x),
               sides = "b",
               color = "black",
               length = unit(0.03, "npc"),
               inherit.aes = FALSE) +
      labs(
        x = plot_title,
        y = "Effect Size"
      ) +
      ylim(-200, 200) +
      theme_test(base_size = 14)
    
    plot_list[[i]] <- p
  }
  
  library(gridExtra)
  
  # 将所有图排列成 3 行 2 列
  p<-grid.arrange(
    plot_list[[1]], plot_list[[2]],
    plot_list[[3]], plot_list[[4]],
    nrow = 2, ncol = 2
  )
  
  ggsave("变量响应曲线_gam.png", p, width = 10, height =6,dpi = 400)
  #ggsave("变量响应曲线_GAM.png", p, width = 12, height =8,dpi = 400)
}

#############################################################################################敏感性
{
  #####auto.arima
  {
    CPUE1 <- as.numeric(CPUE1)  # 确保为数值型
    
    results_auto <- data.frame(
      initial = integer(),
      RMSE = numeric(),
      Adjusted_R2 = numeric(),
      AICc = numeric(),
      Ljung_Box = numeric(),
      Ljung_Box_p = numeric()
    )
    
    for (initial in 9:26) {
      # Step 1: 计算 tsCV 的 RMSE
      funauto <- function(x, h) {
        forecast::forecast(auto.arima(x), h = h)
      }
      e <- forecast::tsCV(CPUE1, funauto, h = 1, initial = initial)
      rmse_val <- sqrt(mean(e^2, na.rm = TRUE))
      mdae_val <- median(abs(e), na.rm = TRUE)
      
      # Step 2: 初始化各指标
      aic_values <- c()
      aicc_values <- c()
      adjr2_values <- c()
      lb_statistics <- c()
      lb_pvalues <- c()
      
      for (i in (initial + 1):(length(CPUE1) - 1)) {
        subdata <- CPUE1[1:i]
        
        fit <- tryCatch({
          auto.arima(subdata)
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          res <- residuals(fit)
          meany <- mean(subdata, na.rm = TRUE)
          r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
          n <- length(subdata)
          k <- length(fit$coef)
          
          if ((n - k - 1) > 0) {
            adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
            adjr2_values <- c(adjr2_values, adjR2)
            
            aic_val <- AIC(fit)
            aic_values <- c(aic_values, aic_val)
            aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
            aicc_values <- c(aicc_values, aicc_val)
            
            lb <- Box.test(res, type = "Ljung-Box")
            lb_statistics <- c(lb_statistics, lb$statistic)
            lb_pvalues <- c(lb_pvalues, lb$p.value)
          }
        }
      }
      
      # Step 3: 存入结果数据框
      results_auto <- rbind(results_auto, data.frame(
        initial = initial + 1,
        RMSE = round(rmse_val, 4),
        MdAE = round(mdae_val, 4),
        Adjusted_R2 = round(mean(adjr2_values, na.rm = TRUE), 4),
        AICc = round(mean(aicc_values, na.rm = TRUE), 4)
      ))
      
    }
    
    # 查看结果
    print(results_auto)
    
    # 可选：保存为 CSV 文件
    # write.csv(results_auto, "auto_arima_initial_results.csv", row.names = FALSE)
  }
  
  results_auto$order <- "Auto.ARIMA"
  results_auto <- results_auto %>% mutate(order = "Auto.ARIMA")
  
  #####GAM 1
  {
    library(mgcv)
    library(AICcmodavg)
    
    CPUE <- as.numeric(CPUE1)  # 确保是数值型向量
    n <- length(CPUE)
    h <- 1
    
    results_gam <- data.frame(
      initial = integer(),
      RMSE = numeric(),
      Adjusted_R2 = numeric(),
      AICc = numeric()
    )
    
    for (initial in 10:27) {
      errors <- rep(NA, n)
      metrics <- data.frame(
        time = integer(),
        Adj_R2 = numeric(),
        RMSE = numeric(),
        AICc = numeric()
      )
      
      for (i in initial:(n - h)) {
        train_xt <- CPUE[1:i]
        train_xt_lag <- CPUE[1:(i - 1)]
        
        train_data <- data.frame(
          xt = train_xt[2:length(train_xt)],
          xt_lag1 = train_xt_lag
        )
        
        gam_model <- tryCatch({
          gam(xt ~ s(xt_lag1), data = train_data, method = "REML")
        }, error = function(e) NULL)
        
        if (!is.null(gam_model)) {
          new_data <- data.frame(xt_lag1 = CPUE[i])
          pred <- predict(gam_model, newdata = new_data)
          test_error <- CPUE[i + 1] - pred
          errors[i + 1] <- test_error
          
          rmse_i <- sqrt(test_error^2)
          r2 <- summary(gam_model)$r.sq
          edf <- sum(summary(gam_model)$edf)
          n_train <- nrow(train_data)
          
          adj_r2 <- if ((n_train - edf - 1) > 0) {
            1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
          } else {
            NA
          }
          
          aicc_val <- tryCatch({
            AICc(gam_model)
          }, error = function(e) NA)
          
          metrics <- rbind(metrics, data.frame(
            time = i,
            Adj_R2 = adj_r2,
            RMSE = rmse_i,
            AICc = aicc_val
          ))
        }
      }
      
      # 汇总结果
      metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")], na.rm = TRUE)
      overall_rmse <- sqrt(mean(na.omit(errors)^2))
      overall_mdae <- median(abs(na.omit(errors)))
      
      results_gam <- rbind(results_gam, data.frame(
        initial = initial,
        RMSE = round(overall_rmse, 4),
        MdAE = round(overall_mdae, 4),
        Adjusted_R2 = round(metrics_mean["Adj_R2"], 4),
        AICc = round(metrics_mean["AICc"], 4)
      ))
      
    }
    
    # 查看结果
    print(results_gam)
    
    # 可选：保存为 CSV
    # write.csv(results_gam, "gam_model_initial_results.csv", row.names = FALSE)
  }
  
  results_gam$order <- "GAM(t-1)"
  results_gam <- results_gam %>% mutate(order = "GAM(t-1)")
  
  #####GAM 2
  {
    library(mgcv)
    library(AICcmodavg)
    
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    
    # 用于保存所有 initial 的评估结果
    results_gam_lag2 <- data.frame(
      initial = integer(),
      RMSE = numeric(),
      Adjusted_R2 = numeric(),
      AICc = numeric()
    )
    
    for (initial in 10:27) {
      errors <- rep(NA, n)
      metrics <- data.frame(
        time = integer(),
        Adj_R2 = numeric(),
        RMSE = numeric(),
        AICc = numeric()
      )
      
      for (i in initial:(n - h)) {
        train_xt <- CPUE[1:i]
        train_xt_lag1 <- CPUE[1:(i - 1)]
        train_xt_lag2 <- CPUE[1:(i - 2)]
        
        # 构建训练数据
        train_data <- data.frame(
          xt = train_xt[3:length(train_xt)],
          xt_lag1 = train_xt_lag1[2:length(train_xt_lag1)],
          xt_lag2 = train_xt_lag2
        )
        
        gam_model <- tryCatch({
          gam(xt ~ s(xt_lag1) + s(xt_lag2), data = train_data, method = "REML")
        }, error = function(e) NULL)
        
        if (!is.null(gam_model)) {
          # 测试集输入
          new_data <- data.frame(
            xt_lag1 = CPUE[i],
            xt_lag2 = CPUE[i - 1]
          )
          pred <- predict(gam_model, newdata = new_data)
          test_error <- CPUE[i + 1] - pred
          errors[i + 1] <- test_error
          
          rmse_i <- sqrt(test_error^2)
          r2 <- summary(gam_model)$r.sq
          edf <- sum(summary(gam_model)$edf)
          n_train <- nrow(train_data)
          
          adj_r2 <- if ((n_train - edf - 1) > 0) {
            1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
          } else {
            NA
          }
          
          aicc_val <- tryCatch({
            AICc(gam_model)
          }, error = function(e) NA)
          
          metrics <- rbind(metrics, data.frame(
            time = i,
            Adj_R2 = adj_r2,
            RMSE = rmse_i,
            AICc = aicc_val
          ))
        }
      }
      
      # 汇总结果
      metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")], na.rm = TRUE)
      overall_rmse <- sqrt(mean(na.omit(errors)^2))
      overall_mdae <- median(abs(na.omit(errors)))
      
      
      results_gam_lag2 <- rbind(results_gam_lag2, data.frame(
        initial = initial,
        # Stepwise_RMSE = round(metrics_mean["RMSE"], 4),
        RMSE = round(overall_rmse, 4),
        MdAE = round(overall_mdae, 4),
        Adjusted_R2 = round(metrics_mean["Adj_R2"], 4),
        AICc = round(metrics_mean["AICc"], 4)
      ))
    }
    
    # 查看结果
    print(results_gam_lag2)
    
  }
  
  results_gam_lag2$order <- "GAM(t-2)"
  results_gam_lag2 <- results_gam_lag2 %>% mutate(order = "GAM(t-2)")
  ##########ARIMA ORDER
  {
    library(forecast)
    
    orders <- list(
      c(0, 0, 1), c(1, 0, 0), c(1, 0, 1),
      c(0, 1, 0), c(0, 1, 1), c(1, 1, 0), c(1, 1, 1),
      c(0, 1, 2), c(2, 1, 0), c(1, 1, 2), c(2, 1, 1),
      c(0, 2, 0), c(0, 2, 1), c(1, 2, 0), c(0, 2, 2),
      c(2, 2, 0), c(1, 2, 2), c(2, 2, 1)
    )
    
    # 初始化结果数据框
    results <- data.frame(
      order = character(),
      initial = integer(),
      AICc = numeric(),
      Adjusted_R2 = numeric(),
      RMSE = numeric(),
      stringsAsFactors = FALSE
    )
    
    # 遍历每个 ARIMA 模型
    for (ord in orders) {
      order_str <- paste(ord, collapse = "-")
      
      # 遍历 initial 值
      for (initial in 9:26) {
        aicc_values <- c()
        adjr2_values <- c()
        
        for (i in initial:(length(CPUE1) - 1)) {
          subdata <- CPUE1[1:i]
          
          fit <- tryCatch(Arima(subdata, order = ord), error = function(e) NULL)
          if (is.null(fit)) next
          
          res <- residuals(fit)
          meany <- mean(subdata, na.rm = TRUE)
          r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
          n <- length(subdata)
          k <- length(fit$coef)
          
          if (n <= k + 1) next  # 避免 AICc 异常
          adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
          aic_val <- AIC(fit)
          aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
          
          adjr2_values <- c(adjr2_values, adjR2)
          aicc_values <- c(aicc_values, aicc_val)
        }
        
        # tsCV 计算 RMSE
        fun_ts <- function(x, h) forecast(Arima(x, order = ord), h = h)
        e <- tryCatch(tsCV(CPUE1, fun_ts, h = 1, initial = initial - 1), error = function(e) rep(NA, length(CPUE1)))
        rmse <- sqrt(mean(e^2, na.rm = TRUE))
        mdae <- median(abs(e), na.rm = TRUE)
        
        # 存储结果
        results <- rbind(results, data.frame(
          order = order_str,
          initial = initial+1,
          AICc = mean(aicc_values, na.rm = TRUE),
          Adjusted_R2 = mean(adjr2_values, na.rm = TRUE),
          RMSE = rmse,
          MdAE = mdae
        ))
      }
    }
    
    # 排序并查看结果
    results <- results[order(results$order, results$initial), ]
    print(results)
    
    result_order<-results
    
    result_order <- results %>%
      mutate(order = paste0("ARIMA(", gsub("-", ",", order), ")"))
    
  }
  
  
  final_results <- bind_rows(result_order, results_gam,results_gam_lag2,results_auto)
  
  # 3. 排序并查看合并后的结果
  final_results <- final_results %>% arrange(order, initial)
  print(final_results)
  
  
  
  #####################################################################################RW
  {
    
    library(ggplot2)
    library(dplyr)
    # 假设 results_all 中有所有模型（包括 Random Walk）结果
    # Step 1: 提取 ARIMA(0,1,0) baseline
    rw_rmse <- final_results %>%
      filter(order == "ARIMA(0,1,0)") %>%
      select(initial, MdAE) %>%
      rename(MdAE_rw = MdAE)
    
    # Step 2: 计算提升幅度
    results_with_gain <- final_results %>%
      left_join(rw_rmse, by = "initial") %>%
      mutate(
        MdAE_gain = (MdAE_rw - MdAE) / MdAE_rw,
        MdAE_gain_level = case_when(
          MdAE_gain > 0.3 ~ ">30%",
          MdAE_gain > 0.15 ~ "15-30%",
          MdAE_gain > 0 ~ "0-15%",
          MdAE_gain <= 0 ~ "<0"
        )
      )
    
    # Step 3: 设置绘图模型顺序
    model_order <- unique(results_with_gain$order)
    
    # Step 4: 绘图
    p <- ggplot(results_with_gain, aes(
      x = initial,
      y = factor(order, levels = rev(model_order)),
      fill = MdAE_gain_level
    )) +
      geom_tile(color = "gray40", size = 0.3) +
      scale_fill_manual(
        values = c(
          ">30%" = "#d73027",
          "15-30%" = "#fc8d59",
          "0-15%" = "#fee08b",
          "<0" = "gray60"
        ),
        breaks = c(">30%", "15-30%", "0-15%", "<0"),
        guide = guide_legend(na.translate = FALSE)
      ) +
      scale_x_continuous(breaks = seq(10, 27, 1), expand = c(0, 0)) +
      labs(
        x = "Training Length",
        y = "Model",
        fill = "MdAE over RW"
      ) +
      coord_fixed(ratio = 0.9) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 24),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks.direction = "out"
      )
    
    
    p
    ggsave("MDAE_gain_vs_rw.png", p, width = 12, height = 10, dpi = 400)
  }
  
}

#############################################################################测试栖息地的偏移程度
{
  library(dplyr)
  PHY<-get(load("nfs_env_PHY_1994_2022.RData"))
  # 原始数据都已 load 进来，例如 Temp, Sali, SSHMLD, Chl, Phyc, Nppv, O2
  # Temp <- combined_data 等
  summary(PHY)
  
  
  # 定义偏移量范围（去掉 0,0 的重复）
  spawn_offsets <- -3:3
  feed_offsets  <- -3:3
  
  for (s in spawn_offsets) {
    for (f in feed_offsets) {
      
      # 跳过 S0_F0（原始情况）
      if (s == 0 & f == 0) next
      
      ## ----------- 阶段纬度范围定义 -----------
      # spawning stages (1-3)
      lat1_min <- 23 + s; lat1_max <- 28 + s
      lat2_min <- 25 + s; lat2_max <- 30 + s
      lat3_min <- 28 + s; lat3_max <- 32 + s
      
      # feeding stages (4-5)
      lat4_min <- 35 + f; lat4_max <- 40 + f
      lat5_min <- 38 + f; lat5_max <- 46 + f
      
      ## ----------- Temp 示例（其他变量同理） -----------
      Temp<-PHY
      
      Temp1 <- Temp %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat1_min, Lat <= lat1_max,
               Month == 3) %>%
        group_by(Year) %>%
        summarise(Temp0_1 = mean(Temp_0, na.rm = TRUE),
                  Temp50_1 = mean(Temp_47, na.rm = TRUE),
                  .groups = "drop")
      
      Temp2 <- Temp %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat2_min, Lat <= lat2_max,
               Month %in% 1:4) %>%
        group_by(Year) %>%
        summarise(Temp0_2 = mean(Temp_0, na.rm = TRUE),
                  Temp50_2 = mean(Temp_47, na.rm = TRUE),
                  .groups = "drop")
      
      Temp3 <- Temp %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat3_min, Lat <= lat3_max,
               Month %in% 3:4) %>%
        group_by(Year) %>%
        summarise(Temp0_3 = mean(Temp_0, na.rm = TRUE),
                  Temp50_3 = mean(Temp_47, na.rm = TRUE),
                  .groups = "drop")
      
      Temp4 <- Temp %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat4_min, Lat <= lat4_max,
               Month %in% 5:6) %>%
        group_by(Year) %>%
        summarise(Temp0_4 = mean(Temp_0, na.rm = TRUE),
                  Temp50_4 = mean(Temp_47, na.rm = TRUE),
                  Temp100_4 = mean(Temp_92, na.rm =TRUE),
                  Temp200_4 = mean(Temp_222, na.rm =TRUE),
                  Temp300_4 = mean(Temp_318, na.rm =TRUE),
                  .groups = "drop")
      
      Temp5 <- Temp %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat5_min, Lat <= lat5_max,
               Month %in% 7:10) %>%
        group_by(Year) %>%
        summarise(Temp0_5 = mean(Temp_0, na.rm = TRUE),
                  Temp50_5 = mean(Temp_47, na.rm = TRUE),
                  Temp100_5 = mean(Temp_92, na.rm =TRUE),
                  Temp200_5 = mean(Temp_222, na.rm =TRUE),
                  Temp300_5 = mean(Temp_318, na.rm =TRUE),
                  .groups = "drop")
      
      Temp_all <- Temp1 %>%
        inner_join(Temp2, by="Year") %>%
        inner_join(Temp3, by="Year") %>%
        inner_join(Temp4, by="Year") %>%
        inner_join(Temp5, by="Year")
      ##################################################
      SSHMLD<-PHY
      
      SSHMLD1 <- SSHMLD %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat1_min, Lat <= lat1_max,
               Month == 3) %>%
        group_by(Year) %>%
        summarise(SSH_1 = mean(SSH, na.rm = TRUE),
                  MLD_1 = mean(MLD, na.rm = TRUE),
                  .groups = "drop")
      
      SSHMLD2 <- SSHMLD %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat2_min, Lat <= lat2_max,
               Month %in% 1:4) %>%
        group_by(Year) %>%
        summarise(SSH_2 = mean(SSH, na.rm = TRUE),
                  MLD_2 = mean(MLD, na.rm = TRUE),
                  .groups = "drop")
      
      SSHMLD3 <- SSHMLD %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat3_min, Lat <= lat3_max,
               Month %in% 3:4) %>%
        group_by(Year) %>%
        summarise(SSH_3 = mean(SSH, na.rm = TRUE),
                  MLD_3 = mean(MLD, na.rm = TRUE),
                  .groups = "drop")
      
      SSHMLD4 <- SSHMLD %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat4_min, Lat <= lat4_max,
               Month %in% 5:6) %>%
        group_by(Year) %>%
        summarise(SSH_4 = mean(SSH, na.rm = TRUE),
                  MLD_4 = mean(MLD, na.rm = TRUE),
                  .groups = "drop")
      
      SSHMLD5 <- SSHMLD %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat5_min, Lat <= lat5_max,
               Month %in% 7:10) %>%
        group_by(Year) %>%
        summarise(SSH_5 = mean(SSH, na.rm = TRUE),
                  MLD_5 = mean(MLD, na.rm = TRUE),
                  .groups = "drop")
      
      SSHMLD_all <- SSHMLD1 %>%
        inner_join(SSHMLD2, by="Year") %>%
        inner_join(SSHMLD3, by="Year") %>%
        inner_join(SSHMLD4, by="Year") %>%
        inner_join(SSHMLD5, by="Year")
      
      #########################################
      BGC<-get(load("nfs_env_BGC_1994_2022.RData"))
      
      Chl<-BGC
      
      Chl1 <- Chl %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat1_min, Lat <= lat1_max,
               Month == 3) %>%
        group_by(Year) %>%
        summarise(Chl0_1 = mean(Chl_1, na.rm = TRUE),
                  Chl50_1 = mean(Chl_47, na.rm = TRUE),
                  .groups = "drop")
      
      Chl2 <- Chl %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat2_min, Lat <= lat2_max,
               Month %in% 1:4) %>%
        group_by(Year) %>%
        summarise(Chl0_2 = mean(Chl_1, na.rm = TRUE),
                  Chl50_2 = mean(Chl_47, na.rm = TRUE),
                  .groups = "drop")
      
      Chl3 <- Chl %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat3_min, Lat <= lat3_max,
               Month %in% 3:4) %>%
        group_by(Year) %>%
        summarise(Chl0_3 = mean(Chl_1, na.rm = TRUE),
                  Chl50_3 = mean(Chl_47, na.rm = TRUE),
                  .groups = "drop")
      
      Chl4 <- Chl %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat4_min, Lat <= lat4_max,
               Month %in% 5:6) %>%
        group_by(Year) %>%
        summarise(Chl0_4 = mean(Chl_1, na.rm = TRUE),
                  Chl50_4 = mean(Chl_47, na.rm = TRUE),
                  Chl100_4 = mean(Chl_97, na.rm =TRUE),
                  Chl200_4 = mean(Chl_200, na.rm =TRUE),
                  Chl300_4 = mean(Chl_301, na.rm =TRUE),
                  .groups = "drop")
      
      Chl5 <- Chl %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat5_min, Lat <= lat5_max,
               Month %in% 7:10) %>%
        group_by(Year) %>%
        summarise(Chl0_5 = mean(Chl_1, na.rm = TRUE),
                  Chl50_5 = mean(Chl_47, na.rm = TRUE),
                  Chl100_5 = mean(Chl_97, na.rm =TRUE),
                  Chl200_5 = mean(Chl_200, na.rm =TRUE),
                  Chl300_5 = mean(Chl_301, na.rm =TRUE),
                  .groups = "drop")
      
      Chl_all <- Chl1 %>%
        inner_join(Chl2, by="Year") %>%
        inner_join(Chl3, by="Year") %>%
        inner_join(Chl4, by="Year") %>%
        inner_join(Chl5, by="Year")
      ###########################################
      Nppv<-BGC
      
      Nppv1 <- Nppv %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat1_min, Lat <= lat1_max,
               Month == 3) %>%
        group_by(Year) %>%
        summarise(Nppv0_1 = mean(Nppv_1, na.rm = TRUE),
                  Nppv50_1 = mean(Nppv_47, na.rm = TRUE),
                  .groups = "drop")
      
      Nppv2 <- Nppv %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat2_min, Lat <= lat2_max,
               Month %in% 1:4) %>%
        group_by(Year) %>%
        summarise(Nppv0_2 = mean(Nppv_1, na.rm = TRUE),
                  Nppv50_2 = mean(Nppv_47, na.rm = TRUE),
                  .groups = "drop")
      
      Nppv3 <- Nppv %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat3_min, Lat <= lat3_max,
               Month %in% 3:4) %>%
        group_by(Year) %>%
        summarise(Nppv0_3 = mean(Nppv_1, na.rm = TRUE),
                  Nppv50_3 = mean(Nppv_47, na.rm = TRUE),
                  .groups = "drop")
      
      Nppv4 <- Nppv %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat4_min, Lat <= lat4_max,
               Month %in% 5:6) %>%
        group_by(Year) %>%
        summarise(Nppv0_4 = mean(Nppv_1, na.rm = TRUE),
                  Nppv50_4 = mean(Nppv_47, na.rm = TRUE),
                  Nppv100_4 = mean(Nppv_97, na.rm =TRUE),
                  Nppv200_4 = mean(Nppv_200, na.rm =TRUE),
                  Nppv300_4 = mean(Nppv_301, na.rm =TRUE),
                  .groups = "drop")
      
      Nppv5 <- Nppv %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat5_min, Lat <= lat5_max,
               Month %in% 7:10) %>%
        group_by(Year) %>%
        summarise(Nppv0_5 = mean(Nppv_1, na.rm = TRUE),
                  Nppv50_5 = mean(Nppv_47, na.rm = TRUE),
                  Nppv100_5 = mean(Nppv_97, na.rm =TRUE),
                  Nppv200_5 = mean(Nppv_200, na.rm =TRUE),
                  Nppv300_5 = mean(Nppv_301, na.rm =TRUE),
                  .groups = "drop")
      
      Nppv_all <- Nppv1 %>%
        inner_join(Nppv2, by="Year") %>%
        inner_join(Nppv3, by="Year") %>%
        inner_join(Nppv4, by="Year") %>%
        inner_join(Nppv5, by="Year")
      ############################################
      O2<-BGC
      
      O21 <- O2 %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat1_min, Lat <= lat1_max,
               Month == 3) %>%
        group_by(Year) %>%
        summarise(O20_1 = mean(O2_1, na.rm = TRUE),
                  O250_1 = mean(O2_47, na.rm = TRUE),
                  .groups = "drop")
      
      O22 <- O2 %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat2_min, Lat <= lat2_max,
               Month %in% 1:4) %>%
        group_by(Year) %>%
        summarise(O20_2 = mean(O2_1, na.rm = TRUE),
                  O250_2 = mean(O2_47, na.rm = TRUE),
                  .groups = "drop")
      
      O23 <- O2 %>%
        filter(Lon >= 130, Lon <= 170,
               Lat >= lat3_min, Lat <= lat3_max,
               Month %in% 3:4) %>%
        group_by(Year) %>%
        summarise(O20_3 = mean(O2_1, na.rm = TRUE),
                  O250_3 = mean(O2_47, na.rm = TRUE),
                  .groups = "drop")
      
      O24 <- O2 %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat4_min, Lat <= lat4_max,
               Month %in% 5:6) %>%
        group_by(Year) %>%
        summarise(O20_4 = mean(O2_1, na.rm = TRUE),
                  O250_4 = mean(O2_47, na.rm = TRUE),
                  O2100_4 = mean(O2_97, na.rm =TRUE),
                  O2200_4 = mean(O2_200, na.rm =TRUE),
                  O2300_4 = mean(O2_301, na.rm =TRUE),
                  .groups = "drop")
      
      O25 <- O2 %>%
        filter(Lon >= 150, Lon <= 165,
               Lat >= lat5_min, Lat <= lat5_max,
               Month %in% 7:10) %>%
        group_by(Year) %>%
        summarise(O20_5 = mean(O2_1, na.rm = TRUE),
                  O250_5 = mean(O2_47, na.rm = TRUE),
                  O2100_5 = mean(O2_97, na.rm =TRUE),
                  O2200_5 = mean(O2_200, na.rm =TRUE),
                  O2300_5 = mean(O2_301, na.rm =TRUE),
                  .groups = "drop")
      
      O2_all <- O21 %>%
        inner_join(O22, by="Year") %>%
        inner_join(O23, by="Year") %>%
        inner_join(O24, by="Year") %>%
        inner_join(O25, by="Year")
      ####################################################
      
      DMI<-read_excel("DMI.xlsx")
      
      DMI <- DMI %>%
        pivot_longer(cols = -Year, names_to = "Month", values_to = "Value")
      
      ONI<-read_excel("ONI.xlsx")
      
      ONI$Month <- as.numeric(ONI$Month)
      DMI$Month <- as.numeric(DMI$Month)
      
      
      large_env<-inner_join(DMI,ONI,by=c("Year","Month"))
      large_env = large_env %>% rename("DMI"="Value")
      ###
      
      ONI1 = large_env %>% filter( Month==3)
      
      ONI1 = ONI1  %>% 
        select(Year, Month,ONI,PDO,AMO,DMI) %>%
        group_by(Year) %>%
        summarise(ONI1 = mean(ONI, na.rm = TRUE),
                  PDO1 = mean(PDO, na.rm = TRUE),
                  AMO1 = mean(AMO, na.rm = TRUE),
                  DMI1 = mean(DMI, na.rm = TRUE)
        )
      
      ONI2 = large_env %>% filter( Month>=1, Month<=4)
      
      ONI2 = ONI2  %>% 
        select(Year, Month,ONI,PDO,AMO,DMI) %>%
        group_by(Year) %>%
        summarise(ONI2 = mean(ONI, na.rm = TRUE),
                  PDO2 = mean(PDO, na.rm = TRUE),
                  AMO2 = mean(AMO, na.rm = TRUE),
                  DMI2 = mean(DMI, na.rm = TRUE)
        )
      
      
      ONI3 = large_env %>% filter( Month>=3, Month<=4)
      ONI3 = ONI3  %>% 
        select(Year, Month,ONI,PDO,AMO,DMI) %>%
        group_by(Year) %>%
        summarise(ONI3 = mean(ONI, na.rm = TRUE),
                  PDO3 = mean(PDO, na.rm = TRUE),
                  AMO3 = mean(AMO, na.rm = TRUE),
                  DMI3 = mean(DMI, na.rm = TRUE)
        )
      
      ONI4 = large_env %>% filter( Month>=5, Month<=6)
      ONI4 = ONI4  %>% 
        select(Year, Month,ONI,PDO,AMO,DMI) %>%
        group_by(Year) %>%
        summarise(ONI4 = mean(ONI, na.rm = TRUE),
                  PDO4 = mean(PDO, na.rm = TRUE),
                  AMO4 = mean(AMO, na.rm = TRUE),
                  DMI4 = mean(DMI, na.rm = TRUE)
        )
      
      
      ONI5 = large_env %>% filter( Month>=7, Month<=10)
      ONI5 = ONI5  %>% 
        select(Year, Month,ONI,PDO,AMO,DMI) %>%
        group_by(Year) %>%
        summarise(ONI5 = mean(ONI, na.rm = TRUE),
                  PDO5 = mean(PDO, na.rm = TRUE),
                  AMO5 = mean(AMO, na.rm = TRUE),
                  DMI5 = mean(DMI, na.rm = TRUE)
        )
      
      
      ONI_all <- ONI1 %>%
        inner_join(ONI2, by="Year") %>%
        inner_join(ONI3, by="Year") %>%
        inner_join(ONI4, by="Year") %>%
        inner_join(ONI5, by="Year")
      
      
      #####
      library(dplyr)
      library(purrr)
      
      #############################################
      
      ## ----------- 合并所有变量 -----------
      NFS_env <- Temp_all %>%
        #inner_join(Sali_all,   by="Year") %>%
        inner_join(SSHMLD_all, by="Year") %>%
        inner_join(Chl_all,    by="Year") %>%
        # inner_join(Phyc_all,   by="Year") %>%
        inner_join(Nppv_all,   by="Year") %>%
        inner_join(ONI_all,   by="Year") %>%
        inner_join(O2_all,     by="Year")
      
      
      
      old_names <- colnames(NFS_env)
      
      # 阶段映射
      stage_map <- c("1" = "E", "2" = "P", "3" = "J", "4" = "S", "5" = "A")
      
      # 气候变量名集合
      climate_vars <- c("ONI", "PDO", "AMO", "DMI", "SSH", "MLD")
      
      # 特殊变量名（O2 要单独处理）
      special_vars <- c("O2")
      
      # 新版转换函数
      convert_name <- function(name) {
        # 处理大尺度气候变量（如 ONI_3 或 ONI3）
        for (clim in climate_vars) {
          if (grepl(paste0("^", clim, "_?[1-5]$"), name)) {
            stage_num <- gsub(paste0("^", clim, "_?"), "", name)
            stage_letter <- stage_map[stage_num]
            return(paste0(clim, "_", stage_letter))
          }
        }
        
        # 处理 O2 类型的变量：O2后接水层数字，再接下划线和阶段数字
        matches_o2 <- regmatches(name, regexec("^O2([0-9]+)_([1-5])$", name))[[1]]
        if (length(matches_o2) == 3) {
          depth <- matches_o2[2]
          stage_letter <- stage_map[matches_o2[3]]
          return(paste0("O2_", depth, "_", stage_letter))
        }
        
        # 处理一般变量：如 Chl100_2、Temp150_3
        matches <- regmatches(name, regexec("^([A-Za-z]+)([0-9]+)_([1-5])$", name))[[1]]
        if (length(matches) == 4) {
          var <- matches[2]
          depth <- matches[3]
          stage_num <- matches[4]
          stage_letter <- stage_map[stage_num]
          return(paste0(var, "_", depth, "_", stage_letter))
        }
        
        # 默认返回原始列名
        return(name)
      }
      
      old_names <- colnames(NFS_env)
      new_names <- sapply(old_names, convert_name)
      colnames(NFS_env) <- new_names
      
      
      
      env_full<-NFS_env
      
      # 假设 env_full 是包含 1994–2022 年的完整环境数据（29 行）
      # 创建当前年份数据（1995–2022）
      data1 <- env_full[2:29,-c(1) ]
      
      # 创建滞后一年数据（1994–2021），并重命名列名加 lag1
      data2 <- env_full[1:28, -c(1)]
      colnames(data2) <- paste0(colnames(data2), "_lag")
      
      # 合并当前年份与滞后数据（按行对齐）
      env_combined <- cbind(data1, data2)
      
      colnames(env_combined)
      # 可选：添加年份列
      years <- 1995:2022
      env_combined$Year <- years
      env_combined <- env_combined[, c(ncol(env_combined), 1:(ncol(env_combined) - 1))]  # 把Year放到第一列
      
      ## ----------- 输出文件 -----------
      suffix_s <- ifelse(s > 0, paste0("+", s), as.character(s))
      suffix_f <- ifelse(f > 0, paste0("+", f), as.character(f))
      filename <- paste0("NFS_env_S", suffix_s, "_F", suffix_f, ".csv")
      
      write.csv(env_combined, filename, row.names = FALSE)
      message("已输出: ", filename)
    }
  }
  
  
  
  
  
  library(forecast)
  library(dplyr)
  library(tools)
  
  
  
  CPUE1<-read_excel("AnnualData.xlsx",sheet = "NFS")
  CPUE1<-ts(CPUE1$CPUE,start=1995,end=2022,frequency = 1)
  CPUE1<-window(CPUE1,start=1995,end=2022)
  plot(CPUE1)
  
  
  CPUE1%>%tseries::adf.test()
  CPUE1%>%diff()%>%tseries::adf.test()
  CPUE1%>%diff()%>%diff()%>%tseries::adf.test()
  
  CPUE1%>%tseries::kpss.test()
  CPUE1%>%diff()%>%tseries::kpss.test()
  
  CPUE1%>%diff()%>%diff()%>%ggtsdisplay()
  
  
  
  
  
  # ===================== 用户设置 =====================
  env_folder <- "."   # 存放 NFS_env_S*_F*.csv 的目录
  cpue_varname <- "CPUE1"   # 这里改成你环境中定义的 CPUE 向量名
  output_dir <- file.path(env_folder, "arima_results_")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # ====================================================
  
  # 取出 CPUE 向量
  if (!exists(cpue_varname, envir = .GlobalEnv)) {
    stop("找不到 CPUE 向量: ", cpue_varname, "，请先在 R 里定义")
  }
  CPUE <- get(cpue_varname, envir = .GlobalEnv)
  
  # 列出所有偏移文件
  files <- list.files(env_folder, pattern = "^NFS_env_S.*_F.*\\.csv$", full.names = TRUE)
  
  process_one_file <- function(file, CPUE, outdir) {
    message("处理: ", basename(file))
    
    df <- read.csv(file)
    # 删除前两列（假设是 Year, Month）
    data <- as.data.frame(scale(df[, -c(1, 2)]))
    
    results_list <- list()
    for (i in 1:ncol(data)) {
      env1 <- data[, i]
      
      # tsCV RMSE / MdAE
      forecast_fun <- function(x, h, xreg, newxreg) {
        fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
        forecast(fit, h = h, xreg = newxreg)
      }
      e <- tsCV(CPUE, forecast_fun, h = 1, xreg = env1, initial = 17)
      rmse <- sqrt(mean(e^2, na.rm = TRUE))
      mdae <- median(abs(e), na.rm = TRUE)
      
      # AICc / R² / AdjR²
      aic_values <- c(); aicc_values <- c(); adjr2_values <- c(); r2_values <- c()
      for (j in 18:(length(CPUE)-1)) {
        subdata <- CPUE[1:j]
        subxreg <- env1[1:j]
        
        fit <- tryCatch(Arima(subdata, order = c(1, 0, 0), xreg = subxreg), error = function(e) NULL)
        if (!is.null(fit)) {
          res <- resid(fit)
          meany <- mean(subdata, na.rm = TRUE)
          r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((subdata - meany)^2, na.rm = TRUE)
          n <- length(subdata); k <- length(fit$coef)
          adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
          aic_value <- AIC(fit)
          aicc_value <- aic_value + (2 * k * (k + 1)) / (n - k - 1)
          
          aic_values <- c(aic_values, aic_value)
          aicc_values <- c(aicc_values, aicc_value)
          adjr2_values <- c(adjr2_values, adjR2)
          r2_values <- c(r2_values, r2)
        }
      }
      
      results_list[[i]] <- data.frame(
        Variable = colnames(data)[i],
        RMSE = rmse,
        MdAE = mdae,
        Mean_AICc = mean(aicc_values),
        Mean_Adjusted_R2 = mean(adjr2_values),
        Mean_R2 = mean(r2_values)
      )
    }
    
    final_results <- bind_rows(results_list)
    final_results <- final_results[order(final_results$RMSE), ]
    
    # 输出文件名：保留原偏移名，加后缀
    base <- file_path_sans_ext(basename(file))
    outname <- file.path(outdir, paste0(base, "_onevar_arima1_0_0.csv"))
    write.csv(final_results, outname, row.names = FALSE, fileEncoding = "UTF-8")
    message("  已保存: ", outname)
  }
  
  # 批量处理
  for (f in files) {
    process_one_file(f, CPUE, output_dir)
  }
  
  
  
  
  ######################################
  library(dplyr)
  library(stringr)
  library(purrr)
  
  baseline_rmse <- 69.36143 
  baseline_r2 <- 0.19953578  
  baseline_AICc <- 268.9233
  baseline_mdae <- 55.97373
  
  # ---------- 找到所有结果文件 ----------
  files <- list.files(pattern = "_onevar_arima1_0_0\\.csv$", full.names = TRUE)
  if(length(files) == 0) stop("未找到匹配文件")
  
  # ---------- 定义单文件处理函数 ----------
  process_one_result_file <- function(path) {
    df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) {
      warning("无法读取：", path); return(NULL)
    })
    if (is.null(df) || !"Variable" %in% colnames(df)) return(NULL)
    
    df2 <- df %>%
      filter(MdAE <= baseline_mdae,
             #RMSE <= baseline_rmse,
             #Mean_AICc <= baseline_AICc,
             Mean_Adjusted_R2 >= baseline_r2) %>%
      mutate(stage_letter = str_extract(Variable, "(?<=_)[EPAJS](?=_|$)"),
             life_stage = case_when(
               stage_letter == "E" ~ "Egg",
               stage_letter == "P" ~ "Paralarvae",
               stage_letter == "J" ~ "Juvenile",
               stage_letter == "S" ~ "Subadult",
               stage_letter == "A" ~ "Adult",
               TRUE ~ NA_character_
             ),
             File = basename(path)  # 文件名列作为索引
      ) %>%
      select(Variable, RMSE,MdAE, File)
    
    return(df2)
  }
  
  # ---------- 批量处理 ----------
  final_long <- map_dfr(files, process_one_result_file)
  
  # ---------- 保存 ----------
  write.csv(final_long, "variable_rmse_allsituation.csv", row.names = FALSE, fileEncoding = "UTF-8")
  
  message("完成：已输出 variable_rmse_allsituation.csv （列=Variable, RMSE, File）")
  
  #######################################
  # 对每个文件计算 rmse_percent 平均值
  file_avg <- final_long %>%
    group_by(File) %>%
    summarise(mean_rmse = mean(RMSE, na.rm = TRUE)) %>%
    arrange(desc(mean_rmse)) %>%
    mutate(rmse_percent_diff = (66.27513 - mean_rmse) / 66.27513)
  
  # 查看结果
  print(file_avg)
  # 保存
  write.csv(file_avg, "file_avg_rmse_percent.csv", row.names = FALSE, fileEncoding = "UTF-8")
  
  
  
  
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  
  # ---------- 提取偏移量 ----------
  file_avg <- file_avg %>%
    mutate(
      # 提取 Spawning 偏移 S
      S = str_extract(File, "(?<=_S)\\+?-?\\d+") %>% as.numeric(),
      # 提取 Feeding 偏移 F
      F = str_extract(File, "(?<=_F)\\+?-?\\d+") %>% as.numeric()
    )
  
  # ---------- 画图 ----------
  p <- ggplot(file_avg, aes(x = S, y = F, fill = rmse_percent_diff)) +
    geom_tile(color = "grey30") +
    scale_fill_gradient2(
      low = "grey50",    # 深灰
      mid = "white",     # 中心白色
      high = "grey80",   # 浅灰
      midpoint = 0,
      limits = c(min(file_avg$rmse_percent_diff, na.rm = TRUE),
                 max(file_avg$rmse_percent_diff, na.rm = TRUE)),
      name = "RMSE diff"
    ) +
    scale_x_continuous(
      breaks = -5:5,
      labels = paste0(abs(-5:5), ifelse(-5:5 < 0, "S", "N")),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = -5:5,
      labels = paste0(abs(-5:5), ifelse(-5:5 < 0, "S", "N")),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5)) +
    labs(x = "Spawning Ground offset", y = "Feeding Ground offset") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 24, lineheight = 1.5),
      legend.title = element_text(size = 26),
      legend.spacing.y = unit(0.6, "cm"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.direction = "out"
    )
  
  # ---------- 保存 ----------
  ggsave("rmse_offset.png", p, width = 8, height = 6, dpi = 400)
  
}

###########arima 多变量
{
  select_enc_reduced <- select_enc[, !colnames(select_enc) %in% out_vars]
  filtered_data <- select_enc_reduced
  library(forecast)
  
  # 确保响应变量为 CPUE1
  y <- as.numeric(CPUE1)
  
  # 滑动窗口设置
  for (win in 10) {
    results_df <- data.frame(
      cols_used = character(),
      Num_Vars = integer(),
      Variable_Indices = character(),
      AICc = numeric(),
      adjR2 = numeric(),
      R2 = numeric(),
      RMSE = numeric(),
      stringsAsFactors = FALSE
    )
    
    num_cols <- ncol(filtered_data)
    
    # 所有变量组合（1~8个变量）
    combinations <- unlist(
      lapply(1:4, function(x) combn(1:num_cols, x, simplify = FALSE)),
      recursive = FALSE
    )
    
    # 滑动预测函数
    forecast_fun <- function(x, h, xreg, newxreg) {
      fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
      forecast::forecast(fit, h = h, xreg = newxreg)
    }
    
    for (combo in combinations) {
      df <- filtered_data[, combo, drop = FALSE]
      df <- ts(df, start = 1995, end = 2022)
      
      # 滑动 RMSE
      e <- tryCatch({
        tsCV(y, forecast_fun, h = 1, xreg = df, initial = 17)
      }, error = function(err) rep(NA, length(y)))
      
      RMSE_val <- sqrt(mean(e^2, na.rm = TRUE))
      MdAE_val <- median(abs(e), na.rm = TRUE)  
      
      
      aic_values <- c()
      aicc_values <- c()
      adjr2_values <- c()
      r2_values <- c()
      
      for (i in 18:27) {
        sub_y <- y[1:i]
        sub_xreg <- df[1:i, , drop = FALSE]
        
        fit <- tryCatch({
          Arima(sub_y, order = c(1, 0, 0), xreg = sub_xreg)
        }, error = function(e) NULL)
        
        if (!is.null(fit)) {
          res <- resid(fit)
          mean_y <- mean(sub_y, na.rm = TRUE)
          r2 <- 1 - sum(res^2, na.rm = TRUE) / sum((sub_y - mean_y)^2, na.rm = TRUE)
          
          n <- length(sub_y)
          k <- length(fit$coef)
          
          if ((n - k - 1) > 0) {
            adjR2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
            aic_val <- AIC(fit)
            aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
            
            aic_values <- c(aic_values, aic_val)
            aicc_values <- c(aicc_values, aicc_val)
            adjr2_values <- c(adjr2_values, adjR2)
          }
          
          r2_values <- c(r2_values, r2)
        }
      }
      
      # 记录组合结果
      results_df <- rbind(results_df, data.frame(
        cols_used = paste(colnames(filtered_data)[combo], collapse = ","),
        Num_Vars = length(combo),
        Variable_Indices = paste0("(", paste(combo, collapse = ","), ")"),
        AICc = if (length(aicc_values) > 0) mean(aicc_values, na.rm = TRUE) else NA,
        adjR2 = if (length(adjr2_values) > 0) mean(adjr2_values, na.rm = TRUE) else NA,
        R2 = mean(r2_values, na.rm = TRUE),
        RMSE = RMSE_val,
        MdAE = MdAE_val,
        stringsAsFactors = FALSE
      ))
    }
    
    # 保存结果
    write.csv(results_df, paste0("arima_nfs", win, ".csv"), row.names = FALSE)
  }
  
  
}
#################gam 多变量
{
  select_enc_reduced <- select_enc[, !colnames(select_enc) %in% out_vars]
  filtered_data<-select_enc_reduced 
  
  # 基础数据
  CPUE <- as.numeric(CPUE1)
  n <- length(CPUE)
  h <- 1
  initial <- 18
  
  # 所有变量名
  filtered_vars <- colnames(filtered_data)
  
  # 所有非空组合（从1到全部）
  all_combos <- unlist(
    lapply(1:length(filtered_vars), function(k) {
      combn(filtered_vars, k, simplify = FALSE)
    }), recursive = FALSE
  )
  
  # 保存结果
  combo_results <- list()
  
  for (combo in all_combos) {
    env_subset <- filtered_data[, combo, drop = FALSE]
    errors <- rep(NA, n)
    metrics <- data.frame(time = integer(), Adj_R2 = numeric(), RMSE = numeric(), AICc = numeric())
    
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- env_subset[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      # GAM 模型公式：xt ~ s(xt_lag1) + s(var1) + s(var2) + ...
      gam_formula <- as.formula(
        paste("xt ~ s(xt_lag1)", paste0(" + s(", names(train_env), ")", collapse = ""))
      )
      
      gam_model <- tryCatch({
        gam(gam_formula, data = train_data, method = "REML")
      }, error = function(e) {
        message("模型拟合失败：", e$message)
        return(NULL)
      })
      
      if (is.null(gam_model)) next  # 跳过该组合
      
      
      # 预测下一个值
      new_data <- cbind(
        xt_lag1 = CPUE[j],
        env_subset[j + 1, , drop = FALSE]
      )
      
      pred <- predict(gam_model, newdata = new_data)
      test_error <- CPUE[j + 1] - pred
      errors[j + 1] <- test_error
      rmse_j <- sqrt(test_error^2)
      
      r2 <- summary(gam_model)$r.sq
      edf <- sum(summary(gam_model)$edf)
      n_train <- nrow(train_data)
      adj_r2 <- 1 - (1 - r2) * (n_train - 1) / (n_train - edf - 1)
      aicc_val <- AICc(gam_model)
      
      metrics <- rbind(metrics, data.frame(time = j, Adj_R2 = adj_r2, RMSE = rmse_j, AICc = aicc_val))
    }
    
    # 计算平均值
    metrics_mean <- colMeans(metrics[, c("Adj_R2", "RMSE", "AICc")])
    overall_rmse <- sqrt(mean(na.omit(errors)^2))
    overall_mdae <- median(abs(na.omit(errors)))
    
    # 保存该组合的结果
    combo_results[[paste(combo, collapse = "+")]] <- data.frame(
      Variables = paste(combo, collapse = "+"),
      Num_Vars = length(combo),
      Mean_Adj_R2 = metrics_mean["Adj_R2"],
      Mean_AICc = metrics_mean["AICc"],
      Mean_Step_RMSE = metrics_mean["RMSE"],
      Overall_RMSE = overall_rmse,
      Overall_MdAE = overall_mdae
    )
  }
  
  # 汇总所有组合结果
  final_combo_results <- bind_rows(combo_results)
  final_combo_results <- final_combo_results[order(final_combo_results$Overall_RMSE), ]
  
  # 创建一个查找向量：变量名 -> 索引
  var_index_lookup <- setNames(seq_along(filtered_vars), filtered_vars)
  
  # 为每个组合提取变量索引
  final_combo_results$Variable_Indices <- sapply(strsplit(final_combo_results$Variables, "\\+"), function(vars) {
    indices <- var_index_lookup[vars]
    paste(indices, collapse = ",")
  })
  
  
  write.csv(final_combo_results, "GAM_results_all_variable_combinations.csv", row.names = FALSE)
  
}

######################选中图
{
######GAM RMSE
{
  # 读取数据
  data <- read.csv("GAM_results_all_variable_combinations.csv")
  colnames(data)
  # 添加变量个数和新指标列
  data <- data %>%
    mutate(
      n_var = str_count(Variable_Indices, ",") + 1,
      RMSE_pct_dec = (67.4671 - Overall_RMSE) / 67.4671 ,
      MdAE_pct_dec = (53.8388 - Overall_MdAE) / 53.8388 ,
      AdjR2_pct_inc = (Mean_Adj_R2 - 0.1981) / 0.1981,
      AICc_diff = Mean_AICc - 258.6844
    )
  
  best_models <- data %>%
    filter(n_var >= 1 & n_var <= 4, AdjR2_pct_inc > 0) %>%
    group_by(n_var) %>%
    slice_max(order_by = RMSE_pct_dec, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 查看结果
  best_models %>% select(n_var, Variable_Indices, RMSE_pct_dec, AdjR2_pct_inc, AICc_diff)
  
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  
  
  # 所有变量名来自 select_enc_reduced
  all_variables <- colnames(select_enc_reduced)
  
  # 拆分已选模型数据，标注颜色
  plot_data <- best_models %>%
    mutate(
      var_index_list = str_split(Variable_Indices, ","),   # list of character vectors
      var_name_list = purrr::map(var_index_list, ~ all_variables[as.integer(.x)])
    ) %>%
    unnest(var_name_list) %>%
    rename(cols_used_name = var_name_list) %>%
    mutate(
      fill_color = case_when(
        RMSE_pct_dec > 0.3 ~ ">30%",
        RMSE_pct_dec > 0.1 ~ "10-30%",
        RMSE_pct_dec > 0 ~ "0-10%",
        RMSE_pct_dec < 0 ~ "<0",
        TRUE ~ "RMSE <= 0"
      )
    )
  
  # 创建完整13x13格子框架
  full_grid <- expand_grid(
    n_var = 1:4,
    cols_used_name = all_variables
  )
  
  # 合并格子数据和绘图数据
  full_plot_data <- full_grid %>%
    left_join(plot_data, by = c("n_var", "cols_used_name"))
  
  # 绘图
  p1<-ggplot(full_plot_data, aes(x = cols_used_name, y = factor(n_var, levels = 6:1), fill = fill_color)) +
    geom_tile(color = "gray50", size = 0.3) +  # 所有格子都有灰色边框
    scale_fill_manual(values = c(
      ">30%"   = "#d73027",   # 深红（显著提升）
      "10-30%" = "#fc8d59",   # 橙红
      "0-10%" = "#fee08b",   # 明黄
      "<0"     = "gray70",  # 柔和草绿（提升较小但不刺眼）
      "RMSE <= 0" = "white"
    ),
    breaks = c(
      ">30%",
      "10-30%",
      "0-10%",
      "<0"
    ),  
    na.value = "white",
    guide = guide_legend(na.translate = FALSE)
    ) + 
    scale_x_discrete(limits = all_variables) +
    labs(
      x = NULL,
      y =NULL,
      title = "GAM",
      fill = "Del RMSE"
    ) +
    coord_fixed(ratio = 1) +  # 关键：压缩横向比例（数值越小横向越紧凑）
    theme_minimal() +
    theme(
      plot.title = element_text(size = 22,hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
  p1
  ggsave("RMSE_GAM_t-1.png",p1,width = 6, height =8,dpi = 300)
}

######GAM ADJ R2
{
  # 读取数据
  data <- read.csv("GAM_results_all_variable_combinations.csv")
  colnames(data)
  # 添加变量个数和新指标列
  data <- data %>%
    mutate(
      n_var = str_count(Variable_Indices, ",") + 1,
      RMSE_pct_dec = (67.4671 - Overall_RMSE) / 67.4671 ,
      MdAE_pct_dec = (53.8388 - Overall_MdAE) / 53.8388 ,
      AdjR2_pct_inc = (Mean_Adj_R2 - 0.1981) / 0.1981,
      AICc_diff = Mean_AICc - 258.6844
    )
  
  best_models <- data %>%
    filter(n_var >= 1 & n_var <= 4, RMSE_pct_dec > 0) %>%
    group_by(n_var) %>%
    slice_max(order_by =AdjR2_pct_inc , n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 查看结果
  best_models %>% select(n_var, Variable_Indices, RMSE_pct_dec, AdjR2_pct_inc, AICc_diff)
  
  # 所有变量名来自 select_enc_reduced
  all_variables <- colnames(select_enc_reduced)
  
  # 拆分已选模型数据，标注颜色
  plot_data <- best_models %>%
    mutate(
      var_index_list = str_split(Variable_Indices, ","),   # list of character vectors
      var_name_list = purrr::map(var_index_list, ~ all_variables[as.integer(.x)])
    ) %>%
    unnest(var_name_list) %>%
    rename(cols_used_name = var_name_list) %>%
    mutate(
      fill_color = case_when(
        AdjR2_pct_inc > 0.3 ~ ">30%",
        AdjR2_pct_inc > 0.1 ~ "10-30%",
        AdjR2_pct_inc > 0 ~ "0-10%",
        AdjR2_pct_inc < 0 ~ "<0",
        TRUE ~ "RMSE <= 0"
      )
    )
  
  
  # 创建完整13x13格子框架
  full_grid <- expand_grid(
    n_var = 1:4,
    cols_used_name = all_variables
  )
  
  # 合并格子数据和绘图数据
  full_plot_data <- full_grid %>%
    left_join(plot_data, by = c("n_var", "cols_used_name"))
  
  # 绘图
  p3<-ggplot(full_plot_data, aes(x = cols_used_name, y = factor(n_var, levels = 6:1), fill = fill_color)) +
    geom_tile(color = "gray50", size = 0.3) +  # 所有格子都有灰色边框
    scale_fill_manual(values = c(
      ">30%"   = "#d73027",   # 深红（显著提升）
      "10-30%" = "#fc8d59",   # 橙红
      "0-10%" = "#fee08b",   # 明黄
      "<0"     = "gray70",  # 柔和草绿（提升较小但不刺眼）
      "RMSE <= 0" = "white"
    ),
    breaks = c(
      ">30%",
      "10-30%",
      "0-10%",
      "<0"
    ),  
    na.value = "white",
    guide = guide_legend(na.translate = FALSE)
    ) + 
    scale_x_discrete(limits = all_variables) +
    labs(
      x = NULL,
      y = NULL,
      fill = expression(DelAdj~R^2)
    ) +
    coord_fixed(ratio = 1) +  # 关键：压缩横向比例（数值越小横向越紧凑）
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
  p3
  ggsave("ADJ_r2_GAM_t-1_选中.png",p3,width = 6, height =8,dpi = 300)
}

######ARIMA RMSE
{
  # 读取数据
  data <- read.csv("arima_nfs10.csv")
  colnames(data)
  # 添加变量个数和新指标列
  data <- data %>%
    mutate(
      n_var = str_count(Variable_Indices, ",") + 1,
      RMSE_pct_dec = (69.36143 - RMSE) / 69.36143 ,
      AdjR2_pct_inc = (adjR2 - 0.19953578) / 0.19953578,
      AICc_diff = AICc - 268.9233
    )
  
  best_models <- data %>%
    filter(n_var >= 1 & n_var <= 4) %>%
    group_by(n_var) %>%
    slice_max(order_by =RMSE_pct_dec, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 查看结果
  best_models %>% select(n_var, Variable_Indices, RMSE_pct_dec, AdjR2_pct_inc, AICc_diff)
  
  # 所有变量名来自 select_enc_reduced
  all_variables <- colnames(select_enc_reduced)
  
  # 拆分已选模型数据，标注颜色
  plot_data <- best_models %>%
    mutate(
      # 去掉括号，拆分成字符向量
      var_index_list = str_remove_all(Variable_Indices, "[\\(\\)]") %>%
        str_split(","),
      # 清理空格并转成整数索引
      var_index_list = purrr::map(var_index_list, ~ as.integer(str_trim(.x))),
      # 排除无效索引
      var_index_list = purrr::map(var_index_list, ~ .x[!is.na(.x) & .x > 0 & .x <= length(all_variables)]),
      # 索引转变量名
      var_name_list = purrr::map(var_index_list, ~ all_variables[.x])
    ) %>%
    unnest(var_name_list) %>%
    rename(cols_used_name = var_name_list) %>%
    mutate(
      fill_color = case_when(
        RMSE_pct_dec > 0.3 ~ ">30%",
        RMSE_pct_dec > 0.1 ~ "10-30%",
        RMSE_pct_dec > 0 ~ "0-10%",
        RMSE_pct_dec < 0 ~ "<0",
        TRUE ~ "RMSE <= 0"
      )
    )
  # 创建完整13x13格子框架
  full_grid <- expand_grid(
    n_var = 1:4,
    cols_used_name = all_variables
  )
  
  # 合并格子数据和绘图数据
  full_plot_data <- full_grid %>%
    left_join(plot_data, by = c("n_var", "cols_used_name"))
  
  # 绘图
  p2<-ggplot(full_plot_data, aes(x = cols_used_name, y = factor(n_var, levels = 10:1), fill = fill_color)) +
    geom_tile(color = "gray50", size = 0.3) +  # 所有格子都有灰色边框
    scale_fill_manual(values = c(
      ">30%"   = "#d73027",   # 深红（显著提升）
      "10-30%" = "#fc8d59",   # 橙红
      "0-10%" = "#fee08b",   # 明黄
      "<0"     = "gray70",  # 柔和草绿（提升较小但不刺眼）
      "RMSE <= 0" = "white"
    ),
    breaks = c(
      ">30%",
      "10-30%",
      "0-10%",
      "<0"
    ),  
    na.value = "white",
    guide = guide_legend(na.translate = FALSE)
    ) + 
    scale_x_discrete(limits = all_variables) +
    labs(
      x = NULL,
      y = NULL,
      title = "ARIMA",
      fill = "Del RMSE"
    ) +
    coord_fixed(ratio = 1) +  # 关键：压缩横向比例（数值越小横向越紧凑）
    theme_minimal() +
    theme(
      plot.title = element_text(size =22,hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
  p2
  ggsave("RMSE_ARIMA1_0_0_选中.png",p2,width = 6, height =8,dpi = 300)
}

######ARIMA ADJ R2
{
  # 读取数据
  data <- read.csv("arima_nfs10.csv")
  colnames(data)
  # 添加变量个数和新指标列
  data <- data %>%
    mutate(
      n_var = str_count(Variable_Indices, ",") + 1,
      RMSE_pct_dec = (69.36143 - RMSE) / 69.36143 ,
      AdjR2_pct_inc = (adjR2 - 0.19953578) / 0.19953578,
      AICc_diff = AICc - 268.9233
    )
  
  best_models <- data %>%
    filter(n_var >= 1 & n_var <= 4) %>%
    group_by(n_var) %>%
    slice_max(order_by = AdjR2_pct_inc, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ###补齐，方便画图
  best_models <- best_models %>%
    add_row(
      cols_used= "1,2,3,4",
      Num_Vars = 4,
      adjR2 = 0,
      AICc = 0,
      RMSE = 0,
      R2 = 0,
      Variable_Indices = "(1,2,3,4)",
      n_var = 4,
      RMSE_pct_dec = 0.1,
      AdjR2_pct_inc = -0.1,
      AICc_diff = 0
    )%>%
    add_row(
      cols_used= "1,2,3,4",
      Num_Vars = 4,
      adjR2 = 0,
      AICc = 0,
      RMSE = 0,
      R2 = 0,
      Variable_Indices = "(1,2,3,4)",
      n_var = 4,
      RMSE_pct_dec = 0.1,
      AdjR2_pct_inc = -0.1,
      AICc_diff = 0
    )
  
  # 查看结果
  best_models %>% select(n_var, Variable_Indices, RMSE_pct_dec, AdjR2_pct_inc, AICc_diff)
  
  # 所有变量名来自 select_enc_reduced
  all_variables <- colnames(select_enc_reduced)
  
  # 拆分已选模型数据，标注颜色
  plot_data <- best_models %>%
    mutate(
      # 去掉括号，拆分成字符向量
      var_index_list = str_remove_all(Variable_Indices, "[\\(\\)]") %>%
        str_split(","),
      # 清理空格并转成整数索引
      var_index_list = purrr::map(var_index_list, ~ as.integer(str_trim(.x))),
      # 排除无效索引
      var_index_list = purrr::map(var_index_list, ~ .x[!is.na(.x) & .x > 0 & .x <= length(all_variables)]),
      # 索引转变量名
      var_name_list = purrr::map(var_index_list, ~ all_variables[.x])
    ) %>%
    unnest(var_name_list) %>%
    rename(cols_used_name = var_name_list) %>%
    mutate(
      fill_color = case_when(
        AdjR2_pct_inc > 0.3 ~ ">30%",
        AdjR2_pct_inc > 0.1 ~ "10-30%",
        AdjR2_pct_inc > 0 ~ "0-10%",
        AdjR2_pct_inc < 0 ~ "<0",
        TRUE ~ "RMSE <= 0"
      )
    )
  # 创建完整13x13格子框架
  full_grid <- expand_grid(
    n_var = 1:4,
    cols_used_name = all_variables
  )
  
  # 合并格子数据和绘图数据
  full_plot_data <- full_grid %>%
    left_join(plot_data, by = c("n_var", "cols_used_name"))
  
  # 绘图
  p4<-ggplot(full_plot_data, aes(x = cols_used_name, y = factor(n_var, levels = 10:1), fill = fill_color)) +
    geom_tile(color = "gray50", size = 0.3) +  # 所有格子都有灰色边框
    scale_fill_manual(values = c(
      ">30%"   = "#d73027",   # 深红（显著提升）
      "10-30%" = "#fc8d59",   # 橙红
      "0-10%" = "#fee08b",   # 明黄
      "<0"     = "gray70",  # 柔和草绿（提升较小但不刺眼）
      "RMSE <= 0" = "white"
    ),
    breaks = c(
      ">30%",
      "10-30%",
      "0-10%",
      "<0"
    ),  
    na.value = "white",
    guide = guide_legend(na.translate = FALSE)
    ) + 
    scale_x_discrete(limits = all_variables) +
    labs(
      x = NULL,
      y = NULL,
      fill = expression(DelAdj~R^2)
    ) +
    coord_fixed(ratio = 1) +  # 关键：压缩横向比例（数值越小横向越紧凑）
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20)
    )
  p4
  ggsave("Adj_R2_ARIMA1_0_0_选中.png",p4,width = 6, height =8,dpi = 300)
}


library(cowplot)
library(ggplot2)

# 横轴标题
x_axis_title <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "Selected Variables"), size = 9) +
  theme_void() +
  theme(plot.margin = margin(-5, 0, 0, 0))  # 上、右、下、左 调整

# 纵轴标题
y_axis_title <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "Number of Variables"), angle = 90, size = 9) +
  theme_void() +
  theme(plot.margin = margin(0, -5, 0, 0))  # 上、右、下、左 调整

# 拼接主图
top_bottom <- plot_grid(p, x_axis_title, nrow = 2, rel_heights = c(1, 0.05))
final_plot <- plot_grid(y_axis_title, top_bottom, ncol = 2, rel_widths = c(0.05, 1))

final_plot


ggsave("all_metrics_选中.png", final_plot , width = 12, height =10,dpi = 300)

}

##############GAM 和ARIMA散点展示图
{
  data <- read.csv("GAM_results_all_variable_combinations.csv")  # 你的文件名
  
  # 提取变量个数
  data$num_vars <- sapply(strsplit(as.character(data$Variable_Indices), ","), length)
  
  # 画图
  library(ggplot2)
  
  p1 <- ggplot(data, aes(x = num_vars, y = Overall_RMSE)) +
    geom_hline(yintercept = 67.4671, color = "grey50", linetype = "dashed", size = 1.6) +
    geom_point(
      shape = 21,
      size = 2.5,
      fill = "black",
      color = "black",
      stroke = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line",
      size = 1.5,
      color = "orange"
    ) +
    
    #geom_smooth(
    #method = "lm",
    #aes(group = Habitats, color = Habitats),
    #se = TRUE,
    #linewidth = 1.2,
    #color = "orange"
    #) +
    labs(
      x = NULL,
      y = "RMSE",
      title = "GAM"
    ) +
    theme_classic() +
    scale_x_continuous(breaks = unique(data$num_vars)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(size = 16,hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none"  # 因为只有一组，不需要legend
    )
  p1
  
  p3 <- ggplot(data, aes(x = num_vars, y = Mean_Adj_R2)) +
    geom_hline(yintercept = 0.1981 , color = "grey50", linetype = "dashed", size = 1.6) +
    geom_point(
      shape = 21,
      size = 2.5,
      fill = "black",
      color = "black",
      stroke = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line",
      size = 1.5,
      color = "orange"
    ) +
    
    #geom_smooth(
    #method = "lm",
    #aes(group = Habitats, color = Habitats),
    #se = TRUE,
    #linewidth = 1.2,
    #color = "orange"
    #) +
    labs(
      x = NULL,
      y = "AdjR square",
      title = NULL
    ) +
    theme_classic() +
    scale_x_continuous(breaks = unique(data$num_vars)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none"  # 因为只有一组，不需要legend
    )
  
  ####################################################################################################
  data <- read.csv("arima_nfs10.csv")  # 你的文件名
  
  # 提取变量个数
  data$num_vars <- sapply(strsplit(as.character(data$Variable_Indices), ","), length)
  
  # 画图
  library(ggplot2)
  
  p2 <- ggplot(data, aes(x = Num_Vars, y = RMSE)) +
    geom_hline(yintercept = 69.36143, color = "grey50", linetype = "dashed", size = 1.6) +
    geom_point(
      shape = 21,
      size =2.5,
      fill = "black",
      color = "black",
      stroke = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line",
      size = 1.5,
      color = "orange"
    ) +
    
    #geom_smooth(
    #method = "lm",
    #aes(group = Habitats, color = Habitats),
    #se = TRUE,
    #linewidth = 1.2,
    #color = "orange"
    #) +
    labs(
      x = NULL,
      y = NULL,
      title ="ARIMA"
    ) +
    theme_classic() +
    scale_x_continuous(breaks = unique(data$num_vars)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(size = 16,hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none"  # 因为只有一组，不需要legend
    )
  
  p4 <- ggplot(data, aes(x = Num_Vars, y = adjR2)) +
    geom_hline(yintercept = 0.19953578 , color = "grey50", linetype = "dashed", size = 1.6) +
    geom_point(
      shape = 21,
      size = 2.5,
      fill = "black",
      color = "black",
      stroke = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line",
      size = 1.5,
      color = "orange"
    ) +
    
    #geom_smooth(
    #method = "lm",
    #aes(group = Habitats, color = Habitats),
    #se = TRUE,
    #linewidth = 1.2,
    #color = "orange"
    #) +
    labs(
      x = NULL,
      y = NULL,
      title = NULL
    ) +
    theme_classic() +
    scale_x_continuous(breaks = unique(data$num_vars)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none"  # 因为只有一组，不需要legend
    )
  
  
  
  library(patchwork)
  
  p <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")& 
    theme(legend.text = element_text(size = 26),  # 调整图例字体大小
          legend.title = element_text(size = 28)) # 可选，调整图例标题字体
  
  p <- p + plot_annotation(
    theme = theme(
      plot.caption = element_text(size = 20, hjust = 0.52, margin = margin(t =10))
    ),
    caption = "Number of variables"  # 作为横坐标统一标题
  )
  
  p
  
  ggsave("RMSE_R2_distribute.png",p,width =10,height =6, dpi=400)
  
}

###########
##########################################################################################ARIMA 交叉验证图
#在跑之前要先加载VIF,筛选出选出的关键变量
{
  data <- read.csv("arima_nfs10.csv")
  # 过滤 adjR2 大于 0.3816971
  filtered_data <- data %>%
    filter(adjR2 > 0.2)%>%
    filter(RMSE < 69.36)%>%
    filter(MdAE < 55.97)
    #filter(AICc<268.92)
  
  
  # 按 RMSE 升序排列
  ranked_models <- filtered_data %>%
    arrange(RMSE)
  
  # 选出前10个最优模型
  top10_models <- ranked_models %>%
    slice(1:10)
  
  # 查看结果
  print(top10_models)
  
  nfs_test<-window(CPUE1,start=2013,end=2022)
  
  # 拆分字符串为数字向量
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[1])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_1 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[2])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_2 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[3])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_3 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[4])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_4 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[5])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_5 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[6])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_6 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[7])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_7 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[8])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_8 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[9])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_9 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[10])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  fun <- function(x, h, xreg, newxreg) {
    fit <- Arima(x, order = c(1, 0, 0), xreg = xreg)
    forecast(fit, h = h, xreg = newxreg)
  }
  e <- forecast::tsCV(CPUE1, fun, h = 1, xreg = df, initial = 17)
  e_nfs_10 <- (nfs_test%>%as.numeric())-((e)%>%na.omit())
  
  
  nfs_test = nfs_test%>%na.omit()%>%as.numeric()
  Year<-2013:2022
  rmse_tets_nfs<-cbind(Year,e_nfs_1,e_nfs_2,e_nfs_3,e_nfs_4,e_nfs_5,
                       e_nfs_6,e_nfs_7,e_nfs_8,e_nfs_9,e_nfs_10,nfs_test)
  
  write.csv(rmse_tets_nfs,file="rmse_tets_nfs_ARIMA1_0_0.csv")
  
  rmse_tets_nfs<-read.csv("rmse_tets_nfs_ARIMA1_0_0.csv")
  
  
  p<-ggplot(rmse_tets_nfs) +
    # 散点图 (黑色点)
    
    # 各线条的绘制并映射到颜色
    geom_line(aes(x = Year, y = e_nfs_1, color = "1"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_2, color = "2"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_3, color = "3"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_4, color = "4"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_5, color = "5"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_6, color = "6"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_7, color = "7"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_8, color = "8"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_9, color = "9"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_10, color = "10"), size = 0.8, na.rm = TRUE) +
    geom_point(aes(x = Year, y = nfs_test, color = "Observed Value"), size = 3) +
    theme_bw() +
    labs(x = "Test Year" , y = "CPUE (metric tons/y/v)",title = "ARIMA", color = "Model") + # 添加图例标题
    # 自定义颜色
    scale_color_manual(
      values = c(
        "Observed Value" = "black",
        "1" = "#e0e0f8",  # 淡紫灰
        "2" = "#bdbdf6",  # 柔淡紫
        "3" = "#9e9ee5",  # 稍深冷紫
        
        # 绿系（模型 4~6）
        "4" = "#d5f5e3",  # 浅薄荷绿
        "5" = "#a2e3cb",  # 中绿蓝
        "6" = "#69c9b2",  # 稍深绿蓝
        
        # 蓝系（模型 7~10）
        "7" = "#d2edf2",  # 浅灰蓝
        "8" = "#93cfdc",  # 中蓝绿
        "9" = "#66b2cc",  # 稍深青蓝
        "10" = "#3b97bb"  # 深湖蓝
      ),
      breaks = c( "1", "2", "3", "4", "5", "6", "7", "8", "9","10") # 控制图例项的顺序
    ) +
    scale_x_continuous(breaks = c(2013, 2015, 2017, 2019, 2021)) +
    #scale_y_continuous(breaks = c(100, 150, 200, 250, 300, 350)) +
    theme(
      legend.position = "right", # 图例在右侧
      legend.text = element_text(size = 13),
      legend.title  = element_text(size = 16),
      legend.key.height = unit(0.7, "cm"),
      legend.key.width = unit(0.8, "cm"),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      panel.border = element_rect(linewidth = 0.8, color = "black"),
      plot.title = element_text(size=20)
    )
  p
  
  ggsave("10cross_test_nfs_ARIMA1_0_0.png", p, width = 9, height = 6, dpi = 400)
  
}

###################################################################################################最后预测 ARIMA
{
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[1])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit1<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[2])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit2<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[3])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit3<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[4])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit4<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[5])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit5<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[6])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit6<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[7])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit7<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[8])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit8<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[9])
  cols <- as.numeric(unlist(strsplit(cols_str, ",")))
  df <- select_enc_reduced[, cols, drop = FALSE]
  df<-ts(df,start=1995,end=2022)
  fit9<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
  
  
  ############################
  fitted_values <- fitted(fit1)
  residuals_values <- residuals(fit1)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci1 <- fitted_values - z * residual_sd
  upper_ci1 <- fitted_values + z * residual_sd
  
  #2
  fitted_values <- fitted(fit2)
  residuals_values <- residuals(fit2)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci2 <- fitted_values - z * residual_sd
  upper_ci2 <- fitted_values + z * residual_sd
  
  #3
  fitted_values <- fitted(fit3)
  residuals_values <- residuals(fit3)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci3 <- fitted_values - z * residual_sd
  upper_ci3 <- fitted_values + z * residual_sd
  
  #41
  fitted_values <- fitted(fit4)
  residuals_values <- residuals(fit4)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci4 <- fitted_values - z * residual_sd
  upper_ci4 <- fitted_values + z * residual_sd
  
  #42
  fitted_values <- fitted(fit5)
  residuals_values <- residuals(fit5)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci5 <- fitted_values - z * residual_sd
  upper_ci5<- fitted_values + z * residual_sd
  
  #5
  fitted_values <- fitted(fit6)
  residuals_values <- residuals(fit6)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci6<- fitted_values - z * residual_sd
  upper_ci6 <- fitted_values + z * residual_sd
  
  #6
  fitted_values <- fitted(fit7)
  residuals_values <- residuals(fit7)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci7 <- fitted_values - z * residual_sd
  upper_ci7 <- fitted_values + z * residual_sd
  
  #71
  fitted_values <- fitted(fit8)
  residuals_values <- residuals(fit8)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci8 <- fitted_values - z * residual_sd
  upper_ci8 <- fitted_values + z * residual_sd
  
  #72
  fitted_values <- fitted(fit9)
  residuals_values <- residuals(fit9)
  alpha <- 0.05  # 置信水平为 95%
  z <- qnorm(1 - alpha / 2)  # z 值
  residual_sd <- sd(residuals_values)
  # 计算上下界
  lower_ci9 <- fitted_values - z * residual_sd
  upper_ci9 <- fitted_values + z * residual_sd
  
  
  alllowerci<-cbind(lower_ci1,lower_ci2,lower_ci3,lower_ci4,
                    lower_ci5,lower_ci6,lower_ci7,lower_ci8,lower_ci9)
  
  allupperci<-cbind(upper_ci1,upper_ci2,upper_ci3,upper_ci4,upper_ci5,
                    upper_ci6,upper_ci7,upper_ci8,upper_ci9)
  
  
  
  # 找出每一行的最小值
  # 计算每一行的均值
  mean_lower_ci <- apply(alllowerci, 1, mean, na.rm = TRUE)
  mean_upper_ci <- apply(allupperci, 1, mean, na.rm = TRUE)
  
  
  ci_data <- data.frame(
    time = time(CPUE1),
    lower = mean_lower_ci,
    upper = mean_upper_ci
  )
  
  
  library(dplyr)
  
  
  #####此处只展示了1996-2022年的结果，以匹配GAM 模型
  ci_data <- ci_data %>%
    mutate(
      lower = pmax(lower, 0),
      upper = pmax(upper, 0)
    ) %>%
    filter(time >= 1996)
  
  
  p <- ggplot() +
    geom_ribbon(data = ci_data, aes(x = time, ymin = lower, ymax = upper), 
                fill = "gray80", alpha = 0.4) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit1)), t >= 1996),
              aes(x = t, y = y, color = "1"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit2)), t >= 1996),
              aes(x = t, y = y, color = "2"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit3)), t >= 1996),
              aes(x = t, y = y, color = "3"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit4)), t >= 1996),
              aes(x = t, y = y, color = "4"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit5)), t >= 1996),
              aes(x = t, y = y, color = "5"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit6)), t >= 1996),
              aes(x = t, y = y, color = "6"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit7)), t >= 1996),
              aes(x = t, y = y, color = "7"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit8)), t >= 1996),
              aes(x = t, y = y, color = "8"), size = 0.8) +
    geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit9)), t >= 1996),
              aes(x = t, y = y, color = "9"), size = 0.8) +
    geom_point(data = subset(data.frame(t = time(CPUE1), y = CPUE1), t >= 1996),
               aes(x = t, y = y, color = "Observed Value"), size = 4) +
    theme_bw() +
    labs(x = "Time", y = "CPUE (metric tons/y/v)", color = "Model") +
    scale_color_manual(
      values = c(
        "Observed Value" = "black",
        "1" = "#c1b4f0",  # 柔和紫
        "2" = "#a89de0",  # 中紫
        "3" = "#907fcf",  # 稍深一点但不暗
        
        # 绿系（模型 4~6）
        "4" = "#aedeb2",  # 柔中绿
        "5" = "#87d6c0",  # 浅绿蓝
        "6" = "#63cdb1",  # 稍深但柔和绿蓝
        
        # 蓝系（模型 7~10）
        "7" = "#a2d6f9",  # 柔蓝
        "8" = "#7cc2f3",  # 稍亮蓝
        "9" = "#58b0e6"  # 中蓝
        #"10" = "#469cd4"
      ),
      breaks = c( "1", "2", "3", "4", "5", "6", "7", "8", "9") # 控制图例项的顺序
    ) +
    scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      legend.position = "right", # 图例在右侧
      legend.text = element_text(size = 17),
      legend.title  = element_text(size = 20),
      legend.key.height = unit(1.1, "cm"),
      legend.key.width = unit(1.2, "cm"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 23),
      axis.title.y = element_text(size = 23),
      panel.border = element_rect(linewidth = 1, color = "black")
    )
  
  p
  
  ggsave("nfs_forecast_1995_2022_ARIMA_1_0_0.png",p,width = 12, height =8,dpi = 400)
  
}


##########################################################################################GAM 画交叉验证图
#在跑之前要先加载VIF,筛选出选出的关键变量
{
  data <- read.csv("GAM_results_all_variable_combinations.csv")
  # 过滤 adjR2 大于 0.3816971
  filtered_data <- data %>%
    filter(Mean_Adj_R2 > 0.1981)%>%
    filter(Overall_RMSE < 67.47)%>%
    filter(Overall_MdAE < 53.84)
  
  # 按 RMSE 升序排列
  ranked_models <- filtered_data %>%
    arrange(Overall_RMSE)
  
  # 选出前10个最优模型
  top10_models <- ranked_models %>%
    slice(1:10)
  
  # 查看结果
  print(top10_models)
  
  
  
  nfs_test<-window(CPUE1,start=2013,end=2022)
  
  # 拆分字符串为数字向量
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[1], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_1 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[2], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_2 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[3], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_3 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[4], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_4 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[5], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_5 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[6], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_6 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[7], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_7 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[8], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_8 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[9], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_9 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[10], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    CPUE <- as.numeric(CPUE1)
    n <- length(CPUE)
    h <- 1
    initial <- 18
    errors <- rep(NA, n)
    for (j in initial:(n - h)) {
      train_xt <- CPUE[1:j]
      train_xt_lag <- CPUE[1:(j - 1)]
      train_env <- df[2:j, , drop = FALSE]
      
      train_data <- cbind(
        xt = train_xt[2:length(train_xt)],
        xt_lag1 = train_xt_lag,
        train_env
      )
      
      env_terms <- paste0("s(", colnames(train_env), ")", collapse = " + ")
      formula_str <- paste0("xt ~ s(xt_lag1) + ", env_terms)
      gam_formula <- as.formula(formula_str)
      
      gam_model <- gam(gam_formula, data = train_data, method = "REML")
      
      new_data <- data.frame(
        xt_lag1 = CPUE[j],
        df[j + 1, , drop = FALSE]
      )
      pred <- predict(gam_model, newdata = new_data)
      errors[j + 1] <- CPUE[j + 1] - pred
    }
    print(errors)
    e_nfs_10 <- (nfs_test%>%as.numeric())-((errors)%>%na.omit())
  }
  
  nfs_test = nfs_test%>%na.omit()%>%as.numeric()
  Year<-2013:2022
  rmse_tets_nfs<-cbind(Year,e_nfs_1,e_nfs_2,e_nfs_3,e_nfs_4,e_nfs_5,
                       e_nfs_6,e_nfs_7,e_nfs_8,e_nfs_9,e_nfs_10,nfs_test)
  
  write.csv(rmse_tets_nfs,file="rmse_tets_nfs.csv")
  
  rmse_tets_nfs<-read.csv("rmse_tets_nfs.csv")
  
  
  p<-ggplot(rmse_tets_nfs) +
    # 散点图 (黑色点)
    
    # 各线条的绘制并映射到颜色
    geom_line(aes(x = Year, y = e_nfs_1, color = "1"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_2, color = "2"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_3, color = "3"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_4, color = "4"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_5, color = "5"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_6, color = "6"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_7, color = "7"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_8, color = "8"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_9, color = "9"), size = 0.8, na.rm = TRUE) +
    geom_line(aes(x = Year, y = e_nfs_9, color = "10"), size = 0.8, na.rm = TRUE) +
    
    geom_point(aes(x = Year, y = nfs_test, color = "Observed Value"), size = 3) +
    theme_bw() +
    labs(x = "Test Year" , y = "CPUE (metric tons/y/v)",title = "GAM", color = "Model") + # 添加图例标题
    # 自定义颜色
    scale_color_manual(
      values = c(
        "Observed Value" = "black",
        "1" = "#e0e0f8",  # 淡紫灰
        "2" = "#bdbdf6",  # 柔淡紫
        "3" = "#9e9ee5",  # 稍深冷紫
        
        # 绿系（模型 4~6）
        "4" = "#d5f5e3",  # 浅薄荷绿
        "5" = "#a2e3cb",  # 中绿蓝
        "6" = "#69c9b2",  # 稍深绿蓝
        
        # 蓝系（模型 7~10）
        "7" = "#d2edf2",  # 浅灰蓝
        "8" = "#93cfdc",  # 中蓝绿
        "9" = "#66b2cc",  # 稍深青蓝
        "10" = "#3b97bb"  # 深湖蓝
        
      ),
      breaks = c( "1", "2", "3", "4", "5", "6", "7", "8", "9","10") # 控制图例项的顺序
    ) +
    scale_x_continuous(breaks = c(2013, 2015, 2017, 2019, 2021)) +
    #scale_y_continuous(breaks = c(100, 150, 200, 250, 300, 350)) +
    theme(
      legend.position = "right", # 图例在右侧
      legend.text = element_text(size = 13),
      legend.title  = element_text(size = 16),
      legend.key.height = unit(0.7, "cm"),
      legend.key.width = unit(0.8, "cm"),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      panel.border = element_rect(linewidth = 0.8, color = "black"),
      plot.title = element_text(size=20)
    )
  
  p
  ggsave("10cross_test_nfs_GAM.png", p, width = 9, height = 6, dpi = 600)
  
}

###################################################################################################最后预测 GAM
{
  ####1
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[1], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model1 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[2], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model2 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[3], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model3 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[4], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model4 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[5], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model5 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[6], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model6 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[7], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model7 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[8], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model8 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  
  {
    cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[9], ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    CPUE <- as.numeric(CPUE1)
    data_fit <- data.frame(
      CPUE = CPUE,
      xt_lag1 = c(NA, CPUE[-length(CPUE)]),
      df
    )
    data_fit <- na.omit(data_fit)
    env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
    gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
    gam_model9 <- gam(gam_formula, data = data_fit, method = "REML")
  }
  #########
  {
    #1
    fitted_values <- fitted(gam_model1)
    residuals_values <- residuals(gam_model1)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci1 <- fitted_values - z * residual_sd
    upper_ci1 <- fitted_values + z * residual_sd
    
    #2
    fitted_values <- fitted(gam_model2)
    residuals_values <- residuals(gam_model2)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci2 <- fitted_values - z * residual_sd
    upper_ci2 <- fitted_values + z * residual_sd
    
    #3
    fitted_values <- fitted(gam_model3)
    residuals_values <- residuals(gam_model3)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci3 <- fitted_values - z * residual_sd
    upper_ci3 <- fitted_values + z * residual_sd
    
    #4
    fitted_values <- fitted(gam_model4)
    residuals_values <- residuals(gam_model4)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci4 <- fitted_values - z * residual_sd
    upper_ci4<- fitted_values + z * residual_sd
    
    #5
    fitted_values <- fitted(gam_model5)
    residuals_values <- residuals(gam_model5)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci5<- fitted_values - z * residual_sd
    upper_ci5 <- fitted_values + z * residual_sd
    
    #6
    fitted_values <- fitted(gam_model6)
    residuals_values <- residuals(gam_model6)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci6 <- fitted_values - z * residual_sd
    upper_ci6 <- fitted_values + z * residual_sd
    
    #71
    fitted_values <- fitted(gam_model7)
    residuals_values <- residuals(gam_model7)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci7 <- fitted_values - z * residual_sd
    upper_ci7 <- fitted_values + z * residual_sd
    
    #8
    fitted_values <- fitted(gam_model8)
    residuals_values <- residuals(gam_model8)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci8 <- fitted_values - z * residual_sd
    upper_ci8 <- fitted_values + z * residual_sd
    
    #9
    fitted_values <- fitted(gam_model9)
    residuals_values <- residuals(gam_model9)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci9 <- fitted_values - z * residual_sd
    upper_ci9 <- fitted_values + z * residual_sd
    
    
    alllowerci<-cbind(lower_ci1,lower_ci2,lower_ci3,lower_ci4,
                      lower_ci5,lower_ci6,lower_ci7,lower_ci8,lower_ci9)
    
    allupperci<-cbind(upper_ci1,upper_ci2,upper_ci3,upper_ci4,upper_ci5,
                      upper_ci6,upper_ci7,upper_ci8,upper_ci9)
    
    
    
    # 找出每一行的最小值
    # 计算每一行的均值
    mean_lower_ci <- apply(alllowerci, 1, mean, na.rm = TRUE)
    mean_upper_ci <- apply(allupperci, 1, mean, na.rm = TRUE)
  }
  
  CPUElag1<-window(CPUE1,start=1996,end=2022)
  
  ci_data <- data.frame(
    time = time(CPUElag1),
    lower = mean_lower_ci,
    upper = mean_upper_ci
  )
  
  
  library(dplyr)
  
  ci_data <- ci_data %>%
    mutate(
      lower = pmax(lower, 0),
      upper = pmax(upper, 0)
    )
  
  p<-ggplot() +
    geom_ribbon(data = ci_data, aes(x = time, ymin = lower, ymax = upper), 
                fill = "gray80", alpha = 0.4) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model1), color = "1"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model2), color = "2"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model3), color = "3"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model4), color = "4"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model5), color = "5"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model6), color = "6"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model7), color = "7"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model8), color = "8"), size = 0.8) +
    geom_line(aes(x = time(CPUElag1), y = fitted(gam_model9), color = "9"), size = 0.8) +
    
    geom_point(aes(x = time(CPUElag1), y = CPUElag1, color = "Observed Value"), size =4) +
    labs(x = "Time", y = "CPUE (metric tons/y/v)", color = "Model") +
    scale_color_manual(
      values = c(
        "Observed Value" = "black",
        "1" = "#c1b4f0",  # 柔和紫
        "2" = "#a89de0",  # 中紫
        "3" = "#907fcf",  # 稍深一点但不暗
        
        # 绿系（模型 4~6）
        "4" = "#aedeb2",  # 柔中绿
        "5" = "#87d6c0",  # 浅绿蓝
        "6" = "#63cdb1",  # 稍深但柔和绿蓝
        
        # 蓝系（模型 7~10）
        "7" = "#a2d6f9",  # 柔蓝
        "8" = "#7cc2f3",  # 稍亮蓝
        "9" = "#58b0e6"
        
      ),
      breaks = c( "1", "2", "3", "4", "5", "6", "7", "8", "9" ) # 控制图例项的顺序
    ) +
    theme_bw() +
    scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      legend.position = "right", # 图例在右侧
      legend.text = element_text(size = 17),
      legend.title  = element_text(size = 20),
      legend.key.height = unit(1.1, "cm"),
      legend.key.width = unit(1.2, "cm"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 23),
      axis.title.y = element_text(size = 23),
      panel.border = element_rect(linewidth = 1, color = "black")
    )
  
  p
  ggsave("nfs_forecast_1995_2022_GAM.png",p,width = 12, height =8,dpi = 400)
  
  
}


#########################################################################################model forecast and test

{
  {
    cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[1])
    cols <- as.numeric(unlist(strsplit(cols_str, ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    fit3<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
    
    cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[2])
    cols <- as.numeric(unlist(strsplit(cols_str, ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    fit4<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
    
    cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[3])
    cols <- as.numeric(unlist(strsplit(cols_str, ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    fit5<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
    
    
    cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[4])
    cols <- as.numeric(unlist(strsplit(cols_str, ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    fit6<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
    
    cols_str <- gsub("[()]", "", ranked_models$Variable_Indices[5])
    cols <- as.numeric(unlist(strsplit(cols_str, ",")))
    df <- select_enc_reduced[, cols, drop = FALSE]
    df<-ts(df,start=1995,end=2022)
    fit7<-Arima(CPUE1,xreg = df,order = c(1, 0, 0))
    

    ###########################
    
    #3
    fitted_values <- fitted(fit3)
    residuals_values <- residuals(fit3)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci3 <- fitted_values - z * residual_sd
    upper_ci3 <- fitted_values + z * residual_sd
    
    #41
    fitted_values <- fitted(fit4)
    residuals_values <- residuals(fit4)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci4 <- fitted_values - z * residual_sd
    upper_ci4 <- fitted_values + z * residual_sd
    
    #42
    fitted_values <- fitted(fit5)
    residuals_values <- residuals(fit5)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci5 <- fitted_values - z * residual_sd
    upper_ci5<- fitted_values + z * residual_sd
    
    #5
    fitted_values <- fitted(fit6)
    residuals_values <- residuals(fit6)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci6<- fitted_values - z * residual_sd
    upper_ci6 <- fitted_values + z * residual_sd
    
    #6
    fitted_values <- fitted(fit7)
    residuals_values <- residuals(fit7)
    alpha <- 0.05  # 置信水平为 95%
    z <- qnorm(1 - alpha / 2)  # z 值
    residual_sd <- sd(residuals_values)
    # 计算上下界
    lower_ci7 <- fitted_values - z * residual_sd
    upper_ci7 <- fitted_values + z * residual_sd
    
    
    
    alllowerci<-cbind(lower_ci3,lower_ci4,
                      lower_ci5,lower_ci6,lower_ci7)
    
    allupperci<-cbind(upper_ci3,upper_ci4,upper_ci5,
                      upper_ci6,upper_ci7)
    
    
    
    # 找出每一行的最小值
    # 计算每一行的均值
    mean_lower_ci <- apply(alllowerci, 1, mean, na.rm = TRUE)
    mean_upper_ci <- apply(allupperci, 1, mean, na.rm = TRUE)
    
    
    ci_data <- data.frame(
      time = time(CPUE1),
      lower = mean_lower_ci,
      upper = mean_upper_ci
    )
    
    
    library(dplyr)
    
    
    #####此处只展示了1996-2022年的结果，以匹配GAM 模型
    ci_data <- ci_data %>%
      mutate(
        lower = pmax(lower, 0),
        upper = pmax(upper, 0)
      ) %>%
      filter(time >= 1996)
    
    
    p <- ggplot() +
      geom_ribbon(data = ci_data, aes(x = time, ymin = lower, ymax = upper), 
                  fill = "gray80", alpha = 0.4) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit3)), t >= 1996),
                aes(x = t, y = y, color = "3"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit4)), t >= 1996),
                aes(x = t, y = y, color = "4"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit5)), t >= 1996),
                aes(x = t, y = y, color = "5"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit6)), t >= 1996),
                aes(x = t, y = y, color = "6"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit7)), t >= 1996),
                aes(x = t, y = y, color = "7"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model1), color = "1"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model2), color = "2"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model3), color = "3"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model4), color = "4"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model5), color = "5"), size = 0.8) +
      
      geom_point(data = subset(data.frame(t = time(CPUE1), y = CPUE1), t >= 1996),
                 aes(x = t, y = y, color = "Observed Value"), size = 4) +
      theme_bw() +
      labs(x = "Time", y = "CPUE (metric tons/y/v)", color = "Model") +
      scale_color_manual(
        values = c(
          "Observed Value" = "black",
          
          "3" = "#907fcf",  # 稍深一点但不暗
          
          # 绿系（模型 4~6）
          "4" = "#aedeb2",  # 柔中绿
          "5" = "#87d6c0",  # 浅绿蓝
          "6" = "#63cdb1",  # 稍深但柔和绿蓝
          
          # 蓝系（模型 7~10）
          "7" = "#a2d6f9" 
          
          #"10" = "#469cd4"
        ),
        breaks = c(  "3", "4", "5", "6", "7") # 控制图例项的顺序
      ) +
      scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
        legend.position = "right", # 图例在右侧
        legend.text = element_text(size = 17),
        legend.title  = element_text(size = 20),
        legend.key.height = unit(1.1, "cm"),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),
        panel.border = element_rect(linewidth = 1, color = "black")
      )
    
    
    
  }
  
  
  {
    ####1
    {
      cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[1], ",")))
      df <- select_enc_reduced[, cols, drop = FALSE]
      df<-ts(df,start=1995,end=2022)
      CPUE <- as.numeric(CPUE1)
      data_fit <- data.frame(
        CPUE = CPUE,
        xt_lag1 = c(NA, CPUE[-length(CPUE)]),
        df
      )
      data_fit <- na.omit(data_fit)
      env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
      gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
      gam_model1 <- gam(gam_formula, data = data_fit, method = "REML")
    }
    
    
    
    {
      cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[2], ",")))
      df <- select_enc_reduced[, cols, drop = FALSE]
      df<-ts(df,start=1995,end=2022)
      CPUE <- as.numeric(CPUE1)
      data_fit <- data.frame(
        CPUE = CPUE,
        xt_lag1 = c(NA, CPUE[-length(CPUE)]),
        df
      )
      data_fit <- na.omit(data_fit)
      env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
      gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
      gam_model2 <- gam(gam_formula, data = data_fit, method = "REML")
    }
    
    {
      cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[3], ",")))
      df <- select_enc_reduced[, cols, drop = FALSE]
      df<-ts(df,start=1995,end=2022)
      CPUE <- as.numeric(CPUE1)
      data_fit <- data.frame(
        CPUE = CPUE,
        xt_lag1 = c(NA, CPUE[-length(CPUE)]),
        df
      )
      data_fit <- na.omit(data_fit)
      env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
      gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
      gam_model3 <- gam(gam_formula, data = data_fit, method = "REML")
    }
    
    {
      cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[4], ",")))
      df <- select_enc_reduced[, cols, drop = FALSE]
      df<-ts(df,start=1995,end=2022)
      CPUE <- as.numeric(CPUE1)
      data_fit <- data.frame(
        CPUE = CPUE,
        xt_lag1 = c(NA, CPUE[-length(CPUE)]),
        df
      )
      data_fit <- na.omit(data_fit)
      env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
      gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
      gam_model4 <- gam(gam_formula, data = data_fit, method = "REML")
    }
    
    {
      cols <- as.numeric(unlist(strsplit(ranked_models$Variable_Indices[5], ",")))
      df <- select_enc_reduced[, cols, drop = FALSE]
      df<-ts(df,start=1995,end=2022)
      CPUE <- as.numeric(CPUE1)
      data_fit <- data.frame(
        CPUE = CPUE,
        xt_lag1 = c(NA, CPUE[-length(CPUE)]),
        df
      )
      data_fit <- na.omit(data_fit)
      env_terms <- paste0("s(", colnames(df), ")", collapse = " + ")
      gam_formula <- as.formula(paste0("CPUE ~ s(xt_lag1) + ", env_terms))
      gam_model5 <- gam(gam_formula, data = data_fit, method = "REML")
    }
    
    
    {
      #1
      fitted_values <- fitted(gam_model1)
      residuals_values <- residuals(gam_model1)
      alpha <- 0.05  # 置信水平为 95%
      z <- qnorm(1 - alpha / 2)  # z 值
      residual_sd <- sd(residuals_values)
      # 计算上下界
      lower_ci1 <- fitted_values - z * residual_sd
      upper_ci1 <- fitted_values + z * residual_sd
      
      #2
      fitted_values <- fitted(gam_model2)
      residuals_values <- residuals(gam_model2)
      alpha <- 0.05  # 置信水平为 95%
      z <- qnorm(1 - alpha / 2)  # z 值
      residual_sd <- sd(residuals_values)
      # 计算上下界
      lower_ci2 <- fitted_values - z * residual_sd
      upper_ci2 <- fitted_values + z * residual_sd
      
      #3
      fitted_values <- fitted(gam_model3)
      residuals_values <- residuals(gam_model3)
      alpha <- 0.05  # 置信水平为 95%
      z <- qnorm(1 - alpha / 2)  # z 值
      residual_sd <- sd(residuals_values)
      # 计算上下界
      lower_ci3 <- fitted_values - z * residual_sd
      upper_ci3 <- fitted_values + z * residual_sd
      
      #4
      fitted_values <- fitted(gam_model4)
      residuals_values <- residuals(gam_model4)
      alpha <- 0.05  # 置信水平为 95%
      z <- qnorm(1 - alpha / 2)  # z 值
      residual_sd <- sd(residuals_values)
      # 计算上下界
      lower_ci4 <- fitted_values - z * residual_sd
      upper_ci4<- fitted_values + z * residual_sd
      
      #5
      fitted_values <- fitted(gam_model5)
      residuals_values <- residuals(gam_model5)
      alpha <- 0.05  # 置信水平为 95%
      z <- qnorm(1 - alpha / 2)  # z 值
      residual_sd <- sd(residuals_values)
      # 计算上下界
      lower_ci5<- fitted_values - z * residual_sd
      upper_ci5 <- fitted_values + z * residual_sd
      
      
      
      
      alllowerci<-cbind(lower_ci1,lower_ci2,lower_ci3,lower_ci4,
                        lower_ci5)
      
      allupperci<-cbind(upper_ci1,upper_ci2,upper_ci3,upper_ci4,upper_ci5
      )
      
      
      
      # 找出每一行的最小值
      # 计算每一行的均值
      mean_lower_ci <- apply(alllowerci, 1, mean, na.rm = TRUE)
      mean_upper_ci <- apply(allupperci, 1, mean, na.rm = TRUE)
    }
    
    CPUElag1<-window(CPUE1,start=1996,end=2022)
    
    ci_data_gam <- data.frame(
      time = time(CPUElag1),
      lower = mean_lower_ci,
      upper = mean_upper_ci
    )
    
    
    library(dplyr)
    
    ci_data_gam <- ci_data_gam %>%
      mutate(
        lower = pmax(lower, 0),
        upper = pmax(upper, 0)
      )
    
    ci_data_combined <- ci_data %>%
      inner_join(ci_data_gam, by = "time", suffix = c("_arima", "_gam")) %>%
      transmute(
        time = time,
        lower = (lower_arima + lower_gam) / 2,
        upper = (upper_arima + upper_gam) / 2
      )
    
    
    
    p_forecast<-ggplot() +
      geom_ribbon(data = ci_data_combined, aes(x = time, ymin = lower, ymax = upper), 
                  fill = "gray70", alpha = 0.4) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit3)), t >= 1996),
                aes(x = t, y = y, color = "ARIMA1"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit4)), t >= 1996),
                aes(x = t, y = y, color = "ARIMA2"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit5)), t >= 1996),
                aes(x = t, y = y, color = "ARIMA3"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit6)), t >= 1996),
                aes(x = t, y = y, color = "ARIMA4"), size = 0.8) +
      geom_line(data = subset(data.frame(t = time(CPUE1), y = fitted(fit7)), t >= 1996),
                aes(x = t, y = y, color = "ARIMA5"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model1), color = "GAM1"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model2), color = "GAM2"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model3), color = "GAM3"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model4), color = "GAM4"), size = 0.8) +
      geom_line(aes(x = time(CPUElag1), y = fitted(gam_model5), color = "GAM5"), size = 0.8) +
      
      
      geom_point(aes(x = time(CPUElag1), y = CPUElag1, color = "Observed Value"), size =4) +
      labs(x = "Time", y = "CPUE (metric tons/y/v)", color = "Model") +
      scale_color_manual(
        values = c(
          "Observed Value" = "black",
          "ARIMA1" = "#c1b4f0",  # 紫 - 浅
          "ARIMA2" = "#a89de0",
          "ARIMA3" = "#907fcf",
          "ARIMA4" = "#7c68c2",
          "ARIMA5" = "#6a55b5",  # 紫 - 深
          "GAM1"   = "#a2d6f9",  # 蓝 - 浅
          "GAM2"   = "#7cc2f3",
          "GAM3"   = "#58b0e6",
          "GAM4"   = "#469cd4",
          "GAM5"   = "#2e83c2"   # 蓝 - 深
        ),
        breaks = c( "ARIMA1", "ARIMA2", "ARIMA3", "ARIMA4", "ARIMA5","GAM1", "GAM2", "GAM3", "GAM4", "GAM5") # 控制图例项的顺序
      ) +
      theme_bw() +
      scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
        legend.position = "right", # 图例在右侧
        legend.text = element_text(size =14),
        plot.title = element_text(size = 22),
        legend.title  = element_text(size = 18),
        legend.key.height = unit(1.1, "cm"),
        legend.key.width = unit(1.1, "cm"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),
        panel.border = element_rect(linewidth = 1, color = "black")
      )
    
    p_forecast
    
    ggsave("final_forecast.png",p_forecast,width =12,height =6, dpi=400)
    
    
    
  }
}












