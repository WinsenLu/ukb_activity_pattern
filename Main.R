########################################################
# 0. Environment preparation & Function package loading
########################################################
# If not installed, please run it first：
# install.packages(c("dplyr","lubridate","tidyverse","survival","readr"))

library(dplyr)
library(lubridate)
library(tidyverse)
library(survival)
library(readr)

########################################################
# Part 1. Data quality control 
########################################################
## 1.1 Read the original data
ms_data  <- read.csv("MS_original.csv",
                     stringsAsFactors = FALSE, strip.white = TRUE)
withdraw <- read.csv("Withdraw.csv",
                     stringsAsFactors = FALSE)

n_cal_no   <- sum(ms_data$Data_calibration == "No", na.rm = TRUE)
n_wear_no  <- sum(ms_data$Data_wear_time == "No", na.rm = TRUE)
n_acc_hi   <- sum(ms_data$Overall_acceleration_average > 100, na.rm = TRUE)
n_withdraw <- sum(ms_data$Participant_ID %in% withdraw$Participant_ID, na.rm = TRUE)
n_mv_blank <- sum(ms_data$Moderate_Vigorous_Day_average == "", na.rm = TRUE)
n_wear_short<- sum(ms_data$Wear_duration_overall < 4.6, na.rm = TRUE)

cat("Data_calibration = No :", n_cal_no, "\n")
cat("Data_wear_time = No   :", n_wear_no, "\n")
cat("Overall_acceleration_average > 100 :", n_acc_hi, "\n")
cat("Withdraw participants :", n_withdraw, "\n")
cat("MV_Day_average blank  :", n_mv_blank, "\n")
cat("Wear_duration < 4.6   :", n_wear_short, "\n")

## 1.2 One-time elimination
ms_clean <- ms_data %>%
  filter(Data_calibration != "No",
         Data_wear_time != "No",
         Overall_acceleration_average <= 100,
         !Participant_ID %in% withdraw$Participant_ID,
         Moderate_Vigorous_Day_average != "",
         Wear_duration_overall >= 4.6)

cat("After QC removal remaining :", nrow(ms_clean), "\n")

## 1.3 Handling of missing BMI 
original_count <- nrow(ms_clean)
ms_clean <- ms_clean %>% filter(!is.na(BMI))
cat("Remove NA in BMI :", original_count - nrow(ms_clean),
    "| Final QC sample :", nrow(ms_clean), "\n")

########################################################
# Part 2. Covariate preprocessing
########################################################
## 2.1 Read in the data after quality control
df0 <- ms_clean    

## 2.2 The four classifications of races
df0 <- df0 %>%
  mutate(ethnic_low = tolower(Ethnic_background),
         Race = case_when(
           str_detect(ethnic_low, "and|mixed|multiple") ~ "Other",
           ethnic_low == "african" ~ "Black",
           str_detect(ethnic_low, "british|irish|any other white") ~ "White",
           str_detect(ethnic_low, "caribbean|any other black") ~ "Black",
           str_detect(ethnic_low, "indian|pakistani|bangladeshi|chinese|asian") ~ "Asian",
           TRUE ~ "Other")) %>%
  select(-ethnic_low)

## 2.3 Define the diet scoring function
score_tertile <- function(x){
  x_num <- parse_number(as.character(x))
  x_num[is.na(x_num) & !is.na(x)] <- 1
  terc <- quantile(x_num, probs = c(1/3, 2/3), na.rm = TRUE)
  case_when(is.na(x_num) ~ NA_real_,
            x_num <= terc[1] ~ 1,
            x_num <= terc[2] ~ 2,
            TRUE ~ 3)
}
meat_score <- function(x){
  case_when(x == "Never" ~ 0,
            x == "Less than once a week" ~ 1,
            x == "Once a week" ~ 2,
            x == "2-4 times a week" ~ 3,
            x == "5-6 times a week" ~ 4,
            x == "Once or more daily" ~ 5,
            TRUE ~ NA_real_)
}
salt_score <- function(x){
  case_when(x == "Never/rarely" ~ 1,
            x == "Sometimes" ~ 2,
            x == "Usually" ~ 3,
            x == "Always" ~ 4,
            TRUE ~ NA_real_)
}

## 2.4 Calculate the dietary score
df_scored <- df0 %>%
  mutate(
    cook_v  = score_tertile(Cooked_vegetable_intake),
    raw_v   = score_tertile(Salad_raw_vegetable_intake),
    fresh_f = score_tertile(Fresh_fruit_intake),
    dried_f = score_tertile(Dried_fruit_intake),
    fruit_veg_tot = rowSums(select(., cook_v, raw_v, fresh_f, dried_f), na.rm = FALSE),
    across(c(Beef_intake, Lamb_mutton_intake, Pork_intake, Processed_meat_intake),
           ~ meat_score(.), .names = "{.col}_s"),
    meat_tot = rowSums(select(., ends_with("_s")), na.rm = FALSE),
    salt_s = salt_score(Salt_added_to_food)
  )

med_fv  <- median(df_scored$fruit_veg_tot, na.rm = TRUE)
med_meat<- median(df_scored$meat_tot, na.rm = TRUE)
med_salt<- median(df_scored$salt_s, na.rm = TRUE)

df_final <- df_scored %>%
  mutate(
    fv_group   = ifelse(fruit_veg_tot > med_fv, "High_fv", "Low_fv"),
    meat_group = ifelse(meat_tot > med_meat, "High_meat", "Low_meat"),
    salt_group = ifelse(salt_s > med_salt, "High_salt", "Low_salt"),
    Diet_quality = case_when(
      fv_group == "Low_fv" & meat_group == "High_meat" ~ "Unhealthy",
      fv_group == "High_fv" & meat_group == "Low_meat" & salt_group == "Low_salt" ~ "Healthy",
      TRUE ~ "Intermediate"
    )
  )

## 2.5 Define activity pattern
# 2.5.1. Wearing time in minutes (Columns 6-12)
wear_cols <- 6:12
df_final[wear_cols] <- df_final[wear_cols] * 60

# 2.5.2. Split Moderate_Vigorous_Day_average (column 5) into a 7-day matrix
mv_matrix <- df_final %>%
  pull(`Moderate_Vigorous_Day_average`) %>%
  strsplit(",") %>%
  lapply(function(x) as.numeric(trimws(x))) %>%
  do.call(rbind, .)   
colnames(mv_matrix) <- paste0("MV_", c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))

# 3. Calculate the number of MVPA minutes per day
mv_minutes <- mv_matrix * as.matrix(df_final[wear_cols])

# 4. Weekly total and maximum two-day total
total_mv  <- rowSums(mv_minutes, na.rm = TRUE)
top2_mv   <- apply(mv_minutes, 1, function(x){
  sum(tail(sort(x, na.last = NA), 2), na.rm = TRUE)
})

df_final <- df_final %>%
  mutate(
    Total_MV_Minutes      = total_mv,
    Top2_Days_MV_Minutes  = top2_mv
  )

# 6. Three classifications of movement patterns: 112.1(25%), 150(guidelines), 226.6(50%), and 388.8(75%)
df_final <- df_final %>%
  mutate(
    Activity_Type = case_when(
      Total_MV_Minutes < 388.8                    ~ "Inactive",
      Top2_Days_MV_Minutes / Total_MV_Minutes >= 0.5 ~ "Active WW",
      TRUE                                       ~ "Active regular"
    )
  )

########################################################
# Part 3. Survival analysis/Cox model
########################################################
## 3.1 Date variable reading (Forced character → Date)
dat <- df_final %>% 
  mutate(
    Date_G35_first_reported = as.Date(Date_G35_first_reported, format = "%Y/%m/%d"),
    Date_of_death           = as.Date(Date_of_death,           format = "%Y/%m/%d"),
    baseline_date           = as.Date(substr(End_time_of_wear, 1, 10)),
    end_fu                  = as.Date("2024-08-31")
  )

## 3.2 Construct survival variables
dat <- dat %>%
  mutate(
    event_date = pmin(Date_G35_first_reported, Date_of_death, end_fu, na.rm = TRUE),
    time       = as.numeric(event_date - baseline_date),
    status     = ifelse(!is.na(Date_G35_first_reported) &
                          Date_G35_first_reported <= event_date, 1, 0),
    has_G35    = !is.na(Date_G35_first_reported)
  ) %>%
  filter(time > 0)

## 3.3 Factorized covariates
dat <- dat %>%
  mutate(
    Activity_Type  = factor(Activity_Type,  levels = c("Inactive", "Active regular", "Active WW")),
    Smoking_status = factor(Smoking_status, levels = c("Never", "Previous", "Current")),
    Race           = factor(Race,           levels = c("White", "Black", "Asian", "Other")),
    Sex            = factor(Sex,            levels = c("Female", "Male")),
    Diet_quality   = factor(Diet_quality,   levels = c("Intermediate", "Healthy", "Unhealthy"))
  )

## 3.4 Summary of event numbers
g35_by_activity <- dat %>%
  group_by(Activity_Type) %>%
  summarise(n_total = n(),
            n_g35   = sum(has_G35),
            pct_g35 = round(100 * n_g35 / n_total, 2))
print(g35_by_activity)

## 3.5 Cox proportional hazards model
cox_model <- coxph(Surv(time, status) ~ Activity_Type + Age + Sex +
                     Townsend_deprivation_index + Smoking_status +
                     Vitamin_D + BMI + Race + Diet_quality,
                   data = dat)

summary(cox_model)

# 4. Proportional hazard hypothesis testing
ph_check <- cox.zph(cox_model)
print(ph_check)           
plot(ph_check)           

## 4.1. Proportional hazard test
ph_check <- cox.zph(cox_model)

## 4.2. Drawing: One for each non-reference level
cn <- colnames(ph_check$y)
lev_plot <- cn[grepl("Activity_Type", cn)]   
par(mfrow = c(1, length(lev_plot)))          
for(col in lev_plot){
  t  <- ph_check$time
  rs <- ph_check$y[, col]
  plot(t, rs, xlab = "Event time",
       ylab = "Schoenfeld residual",
       main  = col, col = "steelblue", pch = 19)
  lines(lowess(t, rs), col = "red", lwd = 2)
  abline(h = 0, lty = 2)
}
par(mfrow = c(1,1))


