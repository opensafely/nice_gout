# install.packages("doBy")
# install.packages("astsa")
# install.packages("ggplot2")
# install.packages("sweep")
# install.packages("forecast")
# install.packages("timetk")
# install.packages("recipes")
# install.packages("fs")
# install.packages("here")
# install.packages("svglite")
# install.packages("gridExtra")
# install.packages("stringr")
# install.packages("lubridate")
# install.packages("strucchange")
# install.packages("prophet")

library(dplyr)
library(ggplot2)
library(fs)
library(zoo)
library(svglite)
library(here)
library(doBy)
library(astsa)
library(sweep)
library(forecast)
library(timetk)
library(gridExtra)
library(stringr)
library(lubridate)

# For running locally
#setwd("C:/Users/k1754142/OneDrive/PhD Project/OpenSAFELY NICE/nice_gout/")

# Create directories if needed
dir_create(here::here("output/figures"), showWarnings = FALSE, recurse = TRUE)
dir_create(here::here("output/tables"), showWarnings = FALSE, recurse = TRUE)

# For storing log
sink("logs/sarima_log.txt")

# Import incidence data from previous step
df <-read.csv("output/tables/arima_standardised.csv")

# Rename variables in the data and ensure dates in correct format
df<- df %>% select(disease, disease_full, mo_year_diagn, incidence, numerator)
#df<- df %>% select(disease, disease, mo_year_diagn, incidence, numerator) 
df$month <- substr(df$mo_year_diagn, 1, 3)
df$mo_year_diagn <- gsub("-", " ", df$mo_year_diagn)
df$mon_year <- df$mo_year_diagn
df$mo_year_diagn <- as.Date(paste0("01 ", df$mo_year_diagn), format = "%d %b %Y")
month_lab <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Extract list of diseases from data
disease_list <- unique(df$disease)

# Define start and end dates from the data
start_date <- min(df$mo_year_diagn, na.rm = TRUE)
end_date   <- max(df$mo_year_diagn, na.rm = TRUE)

# Generate vectors from the above
start <- c(year(start_date), month(start_date))
end <- c(year(end_date), month(end_date))

# Set intervention date of interest (passed from yaml)
args <- commandArgs(trailingOnly = TRUE)
intervention_raw <- args[1]
intervention_date <- as.Date(intervention_raw, format = "%Y-%m-%d")
intervention <- c(year(intervention_date), month(intervention_date))
print(intervention_date)

## Define number of months before intervention
n_preintervention <- (year(intervention_date) - year(start_date)) * 12 + (month(intervention_date) - month(start_date))
print(n_preintervention)

## Define max number of months and years (rounded up) in series
max_index <- (year(end_date) - year(start_date)) * 12 + (month(end_date) - month(start_date)) + 1
print(max_index)
max_years <- ceiling(max_index / 12)

# Define the variables to loop over (can also do for counts)
variables <- c("incidence")

# For writing output tables
first_write <- TRUE

# Loop through diseases
for (j in 1:length(disease_list)) {
  
  dis <- disease_list[j]
  df_dis <- df[df$disease == dis, ]
  df_dis <- df_dis %>%  mutate(index=1:n()) #create an index variable
  disease_full <- unique(df_dis$disease_full)
  #disease_full <- unique(df_dis$disease)

  # Label axes
  y_label <- "Monthly incidence rate per 100,000 population"
  x_label <- "Year"
  
  # Keep data from before intervention date
  df_obs <- df_dis[which(df_dis$index<=n_preintervention),]
  
  # Loop through incidence (+/- counts if needed)
  for (i in 1:length(variables)) {
    var <- variables[i]
    
    dis_title <- str_to_title(str_replace_all(dis, "_", " "))
    
    # Check pre-intervention series is complete
    series_values <- as.numeric(df_obs[[var]])
    
    if (length(series_values) == 0 || any(!is.finite(series_values))) {
      
      message(
        "Skipping ", dis, " / ", var,
        ": no complete pre-intervention series."
      )
      
      png(
        paste0("output/figures/obs_pred_", dis, "_", var, ".png"),
        width = 1200,
        height = 800,
        res = 150
      )
      plot.new()
      text(
        0.5, 0.5,
        "Insufficient data for SARIMA model",
        cex = 1.3
      )
      dev.off()
      
      next
    }
    
    # Convert to time series object
    df_obs_rate <- ts(series_values, frequency = 12, start = start)
    assign(paste0("ts_", var), df_obs_rate)
    
    # Convert to time series object
    df_obs_rate <- ts(series_values, frequency = 12, start=start)
    assign(paste0("ts_", var), df_obs_rate)
    
    # Plot time series for raw data; 1st order difference; 1st order seasonal difference
    svg(filename = paste0("output/figures/raw_pre_covid_", var, "_", dis, ".svg"), width = 8, height = 6)
    plot(df_obs_rate, ylim=c(), type='l', col="blue", xlab="Year", ylab=y_label)
    dev.off()
    svg(filename = paste0("output/figures/differenced_pre_covid_", var, "_", dis, ".svg"), width = 8, height = 6)
    plot(diff(df_obs_rate),type = "l");abline(h=0,col = "red")
    dev.off()
    svg(filename = paste0("output/figures/seasonal_pre_covid_", var, "_", dis, ".svg"), width = 8, height = 6)
    plot(diff(diff(df_obs_rate),12),type = "l");abline(h=0,col = "red")
    dev.off()
    
    # Use auto.arima to fit SARIMA model (identifying terms that optimise BIC/AIC); Nb. for models with poor fit on visual inspection/diagnostics, have explored different models
    suggested.rate<- auto.arima(df_obs_rate, max.p = 5, max.q = 5,  max.P = 2,  max.Q = 2, stepwise=FALSE, trace=TRUE)
    # if (dis == "") {
    #   suggested.rate <- forecast::Arima(df_obs_rate, order = c(0,0,0), seasonal = c(0,1,1))
    # }
    
    print(dis)  
    print(suggested.rate)
    
    # Non-seasonal differences suggested by KPSS; then explore diagnostics to see if this is appropriate; also explore alternative models that include drift constants
    d <- ndiffs(df_obs_rate) 
    print(paste("Non-seasonal differences:", d))
    
    # AIC
    aic_value <- AIC(suggested.rate)
    print(paste("AIC:", round(aic_value, 2)))
    
    # BIC
    bic_value <- BIC(suggested.rate)
    print(paste("BIC:", round(bic_value, 2)))
    
    # Sigma2 (residual variance)
    sigma2_value <- suggested.rate$sigma2
    print(paste("Sigma^2:", round(sigma2_value, 2)))
    
    # Log-likelihood
    loglik_value <- logLik(suggested.rate)
    print(paste("Log Likelihood:", round(as.numeric(loglik_value), 2)))
    
    m1.rate <-suggested.rate
    m1.rate
    
    # Check residual diagnostics
    res <- residuals(m1.rate)
    
    par(mfrow = c(2,2))
    plot(res, main = "Residuals", ylab = "Residuals")
    acf(res, main = "Residuals ACF", ylab = "ACF")
    pacf(res, main = "Residuals PACF", ylab = "PACF")
    hist(res, main = "Residuals Histogram", xlab = "Residuals", ylab = "Density", col = "lightgray")
    par(mfrow = c(1,1))
    dev.off()
    
    theme_centered <- theme(plot.title = element_text(hjust = 0.5))
    
    p1 <- autoplot(ts(res)) +
      ggtitle("Residuals") +
      xlab("Time (months)") + ylab("Residuals") +
      theme_centered
    
    p2 <- ggAcf(res, lag.max = 36) +
      ggtitle("Residuals ACF") + ylab("ACF") +
      theme_centered
    
    p3 <- ggPacf(res, lag.max = 36) +
      ggtitle("Residuals PACF") + ylab("PACF") +
      theme_centered
    
    p4 <- ggplot(data.frame(res = res), aes(x = res)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgray", color = "black") +
      stat_function(fun = dnorm, args = list(mean = mean(res), sd = sd(res)), color = "red", size = 0.5) +
      ggtitle("Residuals Histogram") +
      xlab("Residuals") + ylab("Density") +
      theme_centered
    
    n <- length(res)
    k <- length(coef(m1.rate))
    lags <- unique(pmin(c(12, 24), n - 1))
    
    k_eff <- function(m) max(0, min(k, m - 1))
    
    mk_lb_str <- function(m) {
      lb <- Box.test(res, lag = m, type = "Ljung-Box", fitdf = k_eff(m))
      sprintf("Q(%d): p = %s", m, formatC(lb$p.value, format = "f", digits = 2))
    }
    
    # RMSE
    y <- as.numeric(df_obs_rate)
    y_hat <- as.numeric(fitted(m1.rate))
    
    # RMSE (absolute)
    rmse <- sqrt(mean((y - y_hat)^2, na.rm = TRUE))
    
    # Normalized RMSE as % of the mean
    nrmse <- rmse / mean(y, na.rm = TRUE) * 100
    
    # Bai–Perron and CUSUM tests
    bp_str <- NULL
    cusum_str <- NULL
    
    if (requireNamespace("strucchange", quietly = TRUE)) {
      library(strucchange)
      
      # Bai–Perron test on the series with a minimum segment size
      bp_full <- breakpoints(df_obs_rate ~ 1, h = max(12, frequency(df_obs_rate)))
      
      # Selects optimal number of structural breaks by BIC
      bic_vals <- BIC(bp_full)
      k_grid <- 0:(length(bic_vals) - 1)
      k_bp <- k_grid[which.min(bic_vals)]
      
      # Extract break dates (if any)
      if (k_bp > 0) {
        bp_k <- breakpoints(bp_full, breaks = k_bp)
        tt <- time(df_obs_rate)
        bd_times <- breakdates(bp_k)
        bd_row   <- sapply(bd_times, function(bt) which.min(abs(tt - bt)))
        bd_dates <- df_obs$mo_year_diagn[bd_row]
        bp_str <- sprintf("Bai–Perron: k=%d; breaks=%s",
                          k_bp, paste(as.character(bd_dates), collapse = ", "))
      } else {
        bp_str <- "Bai–Perron: k=0 (no breaks)"
      }
      
      # CUSUM test on residuals
      cusum_fit <- sctest(res ~ 1, type = "OLS-CUSUM")
      cusum_str <- sprintf("CUSUM: p = %s", formatC(cusum_fit$p.value, format = "f", digits = 2))
    }
    
    cap <- paste0(
      "RMSE = ", formatC(rmse, format = "f", digits = 2),
      " (", formatC(nrmse, format = "f", digits = 1), "% of mean) | ",
      "Ljung–Box test results: ",
      paste(vapply(lags, mk_lb_str, character(1)), collapse = " | ")
    )
    
    # Append Bai–Perron and CUSUM strings when computed
    bp_cusum <- paste(na.omit(c(bp_str, cusum_str)), collapse = " | ")
    if (nzchar(bp_cusum)) {
      cap <- paste(cap, bp_cusum, sep = "\n")
    }
    
    caption_grob <- grid::textGrob(
      cap, x = 0.5, hjust = 0.5,
      gp = grid::gpar(fontsize = 12, lineheight = 1.1)
    )
    
    g <- gridExtra::arrangeGrob(
      p1, p2, p3, p4, ncol = 2,
      bottom = caption_grob
    )
    
    #ggsave(sprintf("output/figures/auto_residuals_%s_%s.svg", var, as.character(dis)[1]), plot = g, width = 8, height = 6, device = "svg")
    ggsave(sprintf("output/figures/auto_residuals_%s_%s.png", var, as.character(dis)[1]), plot = g, width = 8, height = 6, device = "png")
    
    # Forecast from intervention date and convert to time series object
    fc.rate  <- forecast(m1.rate, h = (max_index - n_preintervention), level = 95, bootstrap=TRUE, npaths=10000)
    
    # Forecasted rates 
    fc.ratemean <- ts(as.numeric(fc.rate$mean), start=intervention, frequency=12)
    fc.ratelower <- ts(as.numeric(fc.rate$lower), start=intervention, frequency=12) # Lower 95% prediction interval
    fc.rateupper <- ts(as.numeric(fc.rate$upper), start=intervention, frequency=12) # Upper 95% prediction interval
    
    # Flatten the matrix into a vector
    fc_rate<- data.frame(
      YearMonth = as.character(as.yearmon(time(fc.rate$mean))),
      mean = as.numeric(as.matrix(fc.rate$mean)),
      lower = as.numeric(as.matrix(fc.rate$lower)),
      upper = as.numeric(as.matrix(fc.rate$upper))
    )
    
    # Merge forecasted and observed data
    df_new <- df_dis %>% left_join(fc_rate, by = c("mon_year" = "YearMonth"))
    df_new$mean <- ifelse(is.na(df_new$mean), df_new[[var]], df_new$mean) #If NA (i.e. pre-forecast), replace as observed incidence
    df_new$lower <- ifelse(is.na(df_new$lower), df_new[[var]], df_new$lower) #If NA (i.e. pre-forecast), replace as observed incidence
    df_new$upper <- ifelse(is.na(df_new$upper), df_new[[var]], df_new$upper) #If NA (i.e. pre-forecast), replace as observed incidence
    
    df_new <- df_new %>%
      arrange(index) %>%
      mutate(moving_average = rollmean(get(var), k = 3, fill = NA, align = "center"))
    
    df_new <- df_new %>%
      arrange(index) %>%
      mutate(mean_ma = rollmean(mean, k = 3, fill = NA, align = "center"))
    
    # Plot observed and expected graphs
    c1<- 
      ggplot(data = df_new,aes(x = mo_year_diagn))+
      geom_point(aes(y = .data[[var]]), color="#5E716A", alpha = 0.20, size=1.5)+
      geom_line(aes(y = moving_average), color = "#5E716A", linetype = "solid", size=0.50)+
      geom_point(data = df_new %>% filter(mo_year_diagn >= intervention_date), aes(y = mean), color="orange", alpha = 0.20, size=1.5)+
      geom_line(data = df_new %>% filter(mo_year_diagn >= intervention_date), aes(y = mean_ma), color = "orange", linetype = "solid", size=0.50)+
      geom_ribbon(data = df_new %>% filter(mo_year_diagn >= intervention_date), aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey")+
      geom_vline(xintercept = as.numeric(intervention_date), linetype = "dashed", color = "grey")+
      scale_x_date(limits = c(as.Date(sprintf("%d-01-01", start[1])), as.Date(sprintf("%d-01-01", end[1] + 1))),
      breaks = seq(as.Date(sprintf("%d-01-01", start[1])), as.Date(sprintf("%d-01-01", end[1] + 1)), by = "1 year"), date_labels = "%Y")+
      theme_minimal(base_size = 12)+
      xlab(x_label)+
      ylab(y_label)+
      theme(
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "grey"),
        axis.ticks = element_line(color = "grey"),
        axis.text = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11, margin = margin(t = 10)),
        axis.title.y = element_text(size = 11, margin = margin(r = 10, l = 10)), 
        plot.title = element_text(size = 11, hjust = 0.5, face = "plain"))+
      annotate("text", x = intervention_date, y = Inf, label = "COVID-19", vjust = -0.5, hjust = 0.5, size = 3.0, color = "navy", inherit.aes = FALSE)+        
        coord_cartesian(clip = "off")+
        theme(plot.margin = margin(t = 20))
      #ggtitle(paste("Observed vs. expected", tolower(disease_full), "diagnoses"))
    
    saveRDS(c1, file = paste0("output/figures/obs_pred_", var, "_", dis, ".rds"))
    #ggsave(filename = paste0("output/figures/obs_pred_", var, "_", dis, ".svg"), plot = c1, width = 1200, height = 600, units = "px", device = svglite::svglite)
    ggsave(filename = paste0("output/figures/obs_pred_", var, "_", dis, ".png"), plot = c1, width = 1200, height = 600, units = "px", dpi = 144, bg = "white", device = "png", type="cairo")
    
    print(c1)
    
    # Store y-axis values for Prophet sensitivity analyses
    gb <- ggplot_build(c1)
    pp <- gb$layout$panel_params[[1]]
    
    y_limits <- if (!is.null(pp$y.range)) pp$y.range else pp$y$range$range
    y_breaks <- if (!is.null(pp$y.major)) pp$y.major else pp$y$breaks
    
    # Calculate absolute and relative differences between observed and expected values at different time intervals
    a<- c(n_preintervention, (n_preintervention + 12), (n_preintervention + 24), (n_preintervention + 36), n_preintervention)
    b<- c((n_preintervention + 12), (n_preintervention + 24), (n_preintervention + 36), max_index, max_index)
    
    # Number of intervals
    n <- length(a)
    results_list <- vector("list", n)
    
    # Table suffixes
    suffix <- as.character(year(intervention_date) + seq_len(n) - 1L)
    suffix[n] <- "total"
    
    for (i in seq_len(n)) {
      
      observed_val <- df_new %>%
        filter(index > a[i] & index <= b[i]) %>%
        summarise(observed = sum(get(var))) %>%
        select(observed)
      
      predicted_val <- df_new %>%
        filter(index > a[i] & index <= b[i]) %>%
        mutate(se = (upper - lower) / (2 * 1.96)) %>%
        summarise(sum_pred_rate = sum(mean), total_var = sum(se^2)) %>%
        mutate(
          sum_pred_rate_l = (sum_pred_rate - (1.96 * sqrt(total_var))),
          sum_pred_rate_u = (sum_pred_rate + (1.96 * sqrt(total_var)))
        ) %>%
        select(sum_pred_rate, sum_pred_rate_l, sum_pred_rate_u)
      
      observed_predicted_rate <- merge(observed_val, predicted_val) %>%
        mutate(change_rate = observed - sum_pred_rate,
               change_ratelow = observed - sum_pred_rate_l,
               change_ratehigh = observed - sum_pred_rate_u, 
               change_rate_per = (observed - sum_pred_rate) * 100 / sum_pred_rate,
               change_rate_per_low = (observed - sum_pred_rate_l) * 100 / sum_pred_rate_l, 
               change_rate_per_high = (observed - sum_pred_rate_u) * 100 / sum_pred_rate_u
        ) %>%
        mutate(observed = round(observed, 2),
               sum_pred_rate = round(sum_pred_rate, 2),
               sum_pred_rate_l = round(sum_pred_rate_l, 2),
               sum_pred_rate_u = round(sum_pred_rate_u, 2),
               change_rate = round(change_rate, 2),
               change_ratelow = round(change_ratelow, 2),
               change_ratehigh = round(change_ratehigh, 2),
               change_rate_per = round(change_rate_per, 2),
               change_rate_per_low = round(change_rate_per_low, 2),
               change_rate_per_high = round(change_rate_per_high, 2)
        )
      
      current_summary <- observed_predicted_rate %>%
        select(observed, sum_pred_rate, sum_pred_rate_l, sum_pred_rate_u, 
               change_rate, change_ratelow, change_ratehigh,
               change_rate_per, change_rate_per_low, change_rate_per_high)
      
      colnames(current_summary)[1:10] <- paste0(colnames(current_summary)[1:10], ".", suffix[i])
      
      if (i == 1) {
        rates.summary <- current_summary
      } else {
        rates.summary <- merge(rates.summary, current_summary, by = "row.names", all = TRUE) %>%
          select(-Row.names)
      }
    }
    
    rates.summary <- rates.summary %>%
      mutate(disease = dis_title) %>% 
      mutate(measure = var) %>% 
      select(measure, everything()) %>% 
      select(disease, everything())
    
    # Output to csv
    new_row <- rates.summary
    file_name <- "output/tables/change_incidence_byyear.csv"
    
    if (file.exists(file_name)) {
      existing_data <- read.csv(file_name)
      updated_data <- rbind(existing_data, new_row)
      write.csv(updated_data, file_name, row.names = FALSE)
    } else {
      write.csv(new_row, file_name, row.names = FALSE)
    }
    
    # Sensitivity analysis using Prophet forecasting method
    if (requireNamespace("prophet", quietly = TRUE)) {
      
      library(prophet)
      
      df_prop_train <- df_obs %>%
        dplyr::transmute(ds = as.Date(mo_year_diagn), y = .data[[var]])
      
      m.prophet <- prophet(
        df_prop_train,
        yearly.seasonality = TRUE,
        weekly.seasonality = FALSE,
        daily.seasonality = FALSE,
        seasonality.mode = "additive",
        interval.width = 0.95
      )
      
      future <- make_future_dataframe(m.prophet, periods = (max_index - n_preintervention), freq = "month")
      
      cut <- as.POSIXct(intervention_date, tz = "UTC")
      
      fc_prophet <- predict(m.prophet, future) %>%
        dplyr::mutate(ds = as.POSIXct(ds, tz = "UTC")) %>%
        dplyr::filter(ds >= cut) %>%
        dplyr::transmute(
          YearMonth = format(ds, "%b %Y"),
          mean_prophet = yhat,
          lower_prophet = yhat_lower,
          upper_prophet = yhat_upper
        )
      
      df_new2 <- df_dis %>%
        dplyr::left_join(fc_prophet, by = c("mon_year" = "YearMonth")) %>%
        dplyr::mutate(
          mean = ifelse(is.na(mean_prophet), .data[[var]], mean_prophet),
          lower = ifelse(is.na(lower_prophet), .data[[var]], lower_prophet),
          upper = ifelse(is.na(upper_prophet), .data[[var]], upper_prophet)
        ) %>%
        dplyr::arrange(index) %>%
        dplyr::mutate(
          moving_average = zoo::rollmean(.data[[var]], k = 3, fill = NA, align = "center"),
          mean_ma = zoo::rollmean(mean, k = 3, fill = NA, align = "center")
        )
      
      c_prophet <-
        ggplot(df_new2, aes(x = mo_year_diagn)) +
        geom_point(aes(y = .data[[var]]), color = "#5E716A", alpha = 0.25, size = 1.5)+
        geom_line(aes(y = moving_average), color = "#5E716A", linewidth = 0.7)+
        geom_ribbon(data = df_new2 %>% dplyr::filter(mo_year_diagn >= intervention_date), aes(ymin = lower, ymax = upper), alpha = 0.18, fill = "#2c7fb8")+
        geom_point(data = df_new2 %>% dplyr::filter(mo_year_diagn >= intervention_date), aes(y = mean), color = "#2c7fb8", alpha = 0.25, size = 1.2)+
        geom_line(data = df_new2 %>% dplyr::filter(mo_year_diagn >= intervention_date), aes(y = mean_ma), color = "#2c7fb8", linewidth = 0.7)+
        geom_vline(xintercept = as.numeric(intervention_date), linetype = "dashed", color = "grey")+
        scale_x_date(breaks = seq(as.Date(paste0(start[1], "-01-01")), as.Date(paste0(end[1] + 1, "-01-01")), by = "2 years"), date_labels = "%Y")+
        coord_cartesian(ylim = y_limits) +
        scale_y_continuous(breaks = y_breaks) +
        scale_color_manual(values = c("Observed" = "#5E716A", "Expected" = "#2c7fb8")) +
        theme_minimal() +
        xlab("") + ylab("") +
        theme(
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "grey"),
          axis.ticks = element_line(color = "grey"),
          axis.text  = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 10, margin = margin(t = 10)),
          axis.title.y = element_text(size = 10, margin = margin(r = 10)),
          plot.title   = element_text(size = 14, hjust = 0.5, face = "plain")
        )
        #ggtitle(paste0(disease_full))
      
      saveRDS(c_prophet, file = paste0("output/figures/prophet_", var, "_", dis, ".rds"))
      ggsave(filename = paste0("output/figures/prophet_", var, "_", dis, ".svg"),
             plot = c_prophet, width = 8, height = 6, device = "svg")
      
      print(c_prophet)
      
      # Calculate absolute and relative differences between observed and expected values for Prophet at different time intervals
      a<- c(n_preintervention, (n_preintervention + 12), (n_preintervention + 24), (n_preintervention + 36), n_preintervention)
      b<- c((n_preintervention + 12), (n_preintervention + 24), (n_preintervention + 36), max_index, max_index)
      
      # Number of intervals
      n <- length(a)
      results_list <- vector("list", n)
    
      # Table suffixes
      suffix <- as.character(year(intervention_date) + seq_len(n) - 1L)
      suffix[n] <- "total"
      
      for (i in seq_len(n)) {
        
        observed_val <- df_new2 %>%
          filter(index > a[i] & index <= b[i]) %>%
          summarise(observed = sum(get(var))) %>%
          select(observed)
        
        predicted_val <- df_new2 %>%
          filter(index > a[i] & index <= b[i]) %>%
          mutate(se = (upper - lower) / (2 * 1.96)) %>%
          summarise(sum_pred_rate = sum(mean), total_var = sum(se^2)) %>%
          mutate(
            sum_pred_rate_l = (sum_pred_rate - (1.96 * sqrt(total_var))),
            sum_pred_rate_u = (sum_pred_rate + (1.96 * sqrt(total_var)))
          ) %>%
          select(sum_pred_rate, sum_pred_rate_l, sum_pred_rate_u)
        
        observed_predicted_rate <- merge(observed_val, predicted_val) %>%
          mutate(change_rate = observed - sum_pred_rate,
                 change_ratelow = observed - sum_pred_rate_l,
                 change_ratehigh = observed - sum_pred_rate_u, 
                 change_rate_per = (observed - sum_pred_rate) * 100 / sum_pred_rate,
                 change_rate_per_low = (observed - sum_pred_rate_l) * 100 / sum_pred_rate_l, 
                 change_rate_per_high = (observed - sum_pred_rate_u) * 100 / sum_pred_rate_u
          ) %>%
          mutate(observed = round(observed, 2),
            sum_pred_rate = round(sum_pred_rate, 2),
            sum_pred_rate_l = round(sum_pred_rate_l, 2),
            sum_pred_rate_u = round(sum_pred_rate_u, 2),
            change_rate = round(change_rate, 2),
            change_ratelow = round(change_ratelow, 2),
            change_ratehigh = round(change_ratehigh, 2),
            change_rate_per = round(change_rate_per, 2),
            change_rate_per_low = round(change_rate_per_low, 2),
            change_rate_per_high = round(change_rate_per_high, 2)
          )
        
        current_summary <- observed_predicted_rate %>%
          select(observed, sum_pred_rate, sum_pred_rate_l, sum_pred_rate_u, 
                 change_rate, change_ratelow, change_ratehigh,
                 change_rate_per, change_rate_per_low, change_rate_per_high)
        
        colnames(current_summary)[1:10] <- paste0(colnames(current_summary)[1:10], ".", suffix[i])
        
        if (i == 1) {
          rates.summary <- current_summary
        } else {
          rates.summary <- merge(rates.summary, current_summary, by = "row.names", all = TRUE) %>%
            select(-Row.names)
        }
      }
      
      rates.summary <- rates.summary %>%
        mutate(disease = dis_title) %>% 
        mutate(measure = var) %>% 
        select(measure, everything()) %>% 
        select(disease, everything())
      
      new_row   <- rates.summary
      file_name <- "output/tables/change_incidence_byyear_prophet.csv"

      if (file.exists(file_name)) {
        existing_data <- read.csv(file_name)
        updated_data <- rbind(existing_data, new_row)
        write.csv(updated_data, file_name, row.names = FALSE)
      } else {
        write.csv(new_row, file_name, row.names = FALSE)
      }
      
    } else {
      message("Prophet package not installed, skipping forecast.")
    }
  }
}    

graphics.off()
sink()

# # Manual checks for diseases with poor fitting on visual inspection####################
# 
#   # Incidence data - use age and sex-standardised rates for incidence rates and unadjusted for counts
#   df <-read.csv("output/tables/incidence_rates_rounded.csv")
# 
#   #Rename variables in the datafile
#   df<- df %>% select(disease, year, mo_year_diagn, incidence, numerator)
#   df$month <- substr(df$mo_year_diagn, 1, 3)
#   df$mo_year_diagn <- gsub("-", " ", df$mo_year_diagn)
#   df$mon_year <- df$mo_year_diagn
#   df$mo_year_diagn <- as.Date(paste0("01 ", df$mo_year_diagn), format = "%d %b %Y")
#   month_lab <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# 
#   disease_list <- unique(df$disease[df$disease == "Vasc"])
# 
#   # Define the variables to loop over
#   variables <- c("incidence")
#   y_labels <- c("Incidence per 100,000 population")
# 
#   # Loop through diseases
#   for (j in 1:length(disease_list)) {
# 
#     dis <- disease_list[j]
#     df_dis <- df[df$disease == dis, ]
#     df_dis <- df_dis %>%  mutate(index=1:n()) #create an index variable 1,2,3...
# 
#     # Manually set titles based on the disease 
#     max_index <- max(df_dis$index)
# 
#     #Keep only data from before March 2020 and save to separate df
#     df_obs <- df_dis[which(df_dis$index<48),]
# 
#     # Loop through incidence
#     for (i in 1:length(variables)) {
#       var <- variables[i]
#       y_label <- y_labels[i]
# 
#       #Convert to time series object
#       df_obs_rate <- ts(df_obs[[var]], frequency=12, start=c(2016,4))
#       assign(paste0("ts_", var), df_obs_rate)
# 
#       #Use auto.arima to fit SARIMA model - will identify p/q parameters that optimise AIC - for models with poor fit on visual inspection, explore different models
#       if (dis == "Vasc") {
#         #suggested.rate<- auto.arima(df_obs_rate, max.p = 5, max.q = 5,  max.P = 2,  max.Q = 2, stepwise=FALSE, trace=TRUE)
#         suggested.rate<- auto.arima(df_obs_rate, D=1, max.p = 5, max.q = 5,  max.P = 2,  max.Q = 2, stepwise=FALSE, trace=TRUE)
#         #suggested.rate<- arima(df_obs_rate, order=c(1,1,3), seasonal=list(order=c(1,1,0),period=12))
#         #suggested.rate<- auto.arima(df_obs_rate, max.p = 5, max.q = 5,  max.P = 2,  max.Q = 2, stepwise=FALSE, trace=TRUE)
#         } else {
#         suggested.rate<- auto.arima(df_obs_rate, max.p = 5, max.q = 5,  max.P = 2,  max.Q = 2, stepwise=FALSE, trace=TRUE)
#       }
#       print(suggested.rate)
#       aic_value <- AIC(suggested.rate)
#       print(paste("AIC:", aic_value))
#       bic_value <- BIC(suggested.rate)
#       print(paste("BIC:", bic_value))
#       m1.rate <-suggested.rate
#       m1.rate
# 
#       #Monthly prediction intervals from the model
#       confint(m1.rate)
# 
#       #check residuals
#       svg(filename = paste0("output/figures/auto_residuals_", var, "_", dis, ".svg"), width = 8, height = 6)
#       checkresiduals(suggested.rate)
#       dev.off()
#       Box.test(suggested.rate$residuals, lag = max_index - 48 + 1, type = "Ljung-Box")
# 
#       #Forecast from March 2020
#       fc.rate  <- forecast(m1.rate, h = max_index - 48)
# 
#       #Forecasted rates
#       fc.ratemean <- ts(as.numeric(fc.rate$mean), start=c(2020,3), frequency=12)
#       fc.ratelower <- ts(as.numeric(fc.rate$lower[,2]), start=c(2020,3), frequency=12) #lower 95% CI
#       fc.rateupper <- ts(as.numeric(fc.rate$upper[,2]), start=c(2020,3), frequency=12) #upper 95% CI
# 
#       fc_rate<- data.frame(
#         YearMonth = as.character(as.yearmon(time(fc.rate$mean))), # Year and month
#         mean = as.numeric(as.matrix(fc.rate$mean)),
#         lower = as.numeric(as.matrix(fc.rate$lower[, 2])),
#         upper = as.numeric(as.matrix(fc.rate$upper[, 2]))
#         # Flatten the matrix into a vector
#       )
# 
#       df_new <- df_dis %>% left_join(fc_rate, by = c("mon_year" = "YearMonth"))
#       df_new$mean <- ifelse(is.na(df_new$mean), df_new[[var]], df_new$mean) #If NA (i.e. pre-forecast), replace as = observed incidence
#       df_new$lower <- ifelse(is.na(df_new$lower), df_new[[var]], df_new$lower) #If NA (i.e. pre-forecast), replace as = observed incidence
#       df_new$upper <- ifelse(is.na(df_new$upper), df_new[[var]], df_new$upper) #If NA (i.e. pre-forecast), replace as = observed incidence
# 
#       df_new <- df_new %>%
#         arrange(index) %>%
#         mutate(moving_average = rollmean(get(var), k = 3, fill = NA, align = "center"))
# 
#       df_new <- df_new %>%
#         arrange(index) %>%
#         mutate(mean_ma = rollmean(mean, k = 3, fill = NA, align = "center"))
# 
#       #observed and predicted graphs to check model predictions against observed values
#       c1<-
#         ggplot(data = df_new,aes(x = mo_year_diagn))+
#         geom_point(aes(y = .data[[var]]), color="#5E716A", alpha = 0.25, size=1.5)+
#         geom_line(aes(y = moving_average), color = "#5E716A", linetype = "solid", size=0.70)+
#         geom_point(data = df_new %>% filter(mo_year_diagn > as.Date("2020-02-01")), aes(y = mean), color="orange", alpha = 0.25, size=1.5)+
#         geom_line(data = df_new %>% filter(mo_year_diagn > as.Date("2020-02-01")), aes(y = mean_ma), color = "orange", linetype = "solid", size=0.65)+
#         #geom_ribbon(data = df_new %>% filter(mo_year_diagn > as.Date("2020-02-01")), aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "grey")+
#         geom_segment(x = as.Date("2020-03-01"),
#                      xend = as.Date("2020-03-01"),
#                      y = min(df_new[[var]], na.rm = TRUE) * 0.85,
#                      yend = max(df_new[[var]], na.rm = TRUE) * 1.15,
#                      linetype = "dashed",
#                      color = "grey")+
#         scale_x_date(breaks = seq(as.Date("2016-01-01"), as.Date("2026-01-01"), by = "2 years"),
#                      date_labels = "%Y")+
#         scale_y_continuous(limits = c(min(df_new[[var]], na.rm = TRUE) * 0.85,
#                                       max(df_new[[var]], na.rm = TRUE) * 1.15),
#                            breaks = pretty(df_new[[var]], n = 4),
#                            expand = expansion(mult = c(0.05, 0.05)))+
#         theme_minimal()+
#         #xlab("Year of diagnosis")+
#         #ylab(y_label)+
#         xlab("")+
#         ylab("")+
#         theme(
#           legend.title = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(color = "grey"),
#           axis.ticks = element_line(color = "grey"),
#           axis.text = element_text(size = 10, color = "black"),
#           axis.title.x = element_text(size = 10, margin = margin(t = 10)),
#           axis.title.y = element_text(size = 10, margin = margin(r = 10)),
#           plot.title = element_text(size = 14, hjust = 0.5, face = "plain")
#         ) +
#         ggtitle(disease_full)
# 
#       saveRDS(c1, file = paste0("output/figures/test_", var, "_", dis, ".rds"))
#       ggsave(filename = paste0("output/figures/test_", var, "_", dis, ".svg"), plot = c1, width = 8, height = 6, device = "svg")
# 
#       print(c1)
#     }
#   }
# 
# dev.off()
# graphics.off()
# sink()