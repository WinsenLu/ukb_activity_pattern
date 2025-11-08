# 0. Nagvis -------------------------------------------------------------
library(readr)
library(dplyr)
library(survival)

# 1. up load  ----------------------------------
dat <- read_csv(
  "MS_Subgroup.csv",
  col_types = cols(
    Date_G35_first_reported = col_character(),
    Date_of_death           = col_character(),
    baseline_date           = col_character()
  )
)

# 2. Manual Date conversion ------------------------------------
dat <- dat %>%
  mutate(
    Date_G35_first_reported = as.Date(Date_G35_first_reported, format = "%Y/%m/%d"),
    Date_of_death           = as.Date(Date_of_death,           format = "%Y/%m/%d"),
    baseline_date           = as.Date(baseline_date,           format = "%Y/%m/%d"),
    end_fu                  = as.Date("2024-08-31")
  )

# 3. Survival variable: Take G35 as the event ------------------------------------------------
dat <- dat %>%
  mutate(
    event_date = pmin(Date_G35_first_reported, Date_of_death, end_fu, na.rm = TRUE),
    time       = as.numeric(event_date - baseline_date),
    status     = ifelse(!is.na(Date_G35_first_reported) &
                          Date_G35_first_reported <= event_date, 1, 0),
    has_G35    = !is.na(Date_G35_first_reported)   
  ) %>%
  filter(time > 0)   

# 4. factorization ---------------------------------------------------------------
dat <- dat %>%
  mutate(
    Activity_Type  = factor(Activity_Type,  levels = c("Inactive", "Active regular", "Active WW")),
    Smoking_status = factor(Smoking_status, levels = c("Never", "Previous", "Current")),
    Race           = factor(Race,           levels = c("White", "Black", "Asian", "Other")),
    Sex            = factor(Sex,            levels = c("Female", "Male")),
    Diet_quality   = factor(Diet_quality,   levels = c("Intermediate", "Healthy", "Unhealthy"))
  )

# 5. Summary of the number of patients for each Activity_Type --------------------------------------
g35_by_activity <- dat %>%
  group_by(Activity_Type) %>%
  summarise(
    n_total   = n(),
    n_g35     = sum(has_G35),
    pct_g35   = round(100 * n_g35 / n_total, 2)
  ) %>%
  ungroup()

print(g35_by_activity)

# 6. Cox ---------------------------------------
cox_model <- coxph(Surv(time, status) ~ Activity_Type + Age + Sex +
                     Townsend_deprivation_index + Smoking_status +
                     Vitamin_D + BMI + Race + Diet_quality,
                   data = dat)
summary(cox_model)

# 7. Subgroup analysis + multiplicative interaction test
## 7.1 Make sure Age_cat is a factor and the horizontal order is clear
dat <- dat %>%
  mutate(
    Age_cat = factor(Age_cat, levels = c("≤65", ">65"))
  )

## 7.2 Full population model: Add interaction terms (multiplicative scale)
cox_int <- coxph(
  Surv(time, status) ~ Activity_Type * Age_cat + Sex + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat
)
summary(cox_int)
## 7.3 Wald test for one-time output of interaction items
print("Global interaction verification（Activity_Type × Age_cat）:")
anova(cox_int, test = "Chisq") %>%               
  slice_tail(n = 1) %>%                         
  select(-loglik) %>% print()

## 7.4 Subgroup models: Model One was proposed respectively among people aged ≤65 and >65
cox_sub_65 <- coxph(
  Surv(time, status) ~ Activity_Type + Sex + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Age_cat == "≤65")
)
summary(cox_sub_65)

cox_sub_65plus <- coxph(
  Surv(time, status) ~ Activity_Type + Sex + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Age_cat == ">65")
)
summary(cox_sub_65plus)

## 8.1 Full population model: Add the Activity_Type * Sex interaction item
cox_int_sex <- coxph(
  Surv(time, status) ~ Activity_Type * Sex + Age + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat
)

summary(cox_int_sex)

## 8.2 Global interaction Wald test（Activity_Type × Sex）
print("Global interaction Wald test（Activity_Type × Sex）:")
anova(cox_int_sex, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 8.3 Subgroup model: Model one is respectively drafted in Female and Male
cox_sub_female <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Sex == "Female")
)
summary(cox_sub_female)

cox_sub_male <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Race +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Sex == "Male")
)
summary(cox_sub_male)

# 9. Activity_Type × Race_cat multiplicative interaction
## 9.1 Make sure that Race_cat is a factor and the horizontal order is clear
dat <- dat %>%
  mutate(
    Race_cat = factor(Race_cat, levels = c("White", "Non-white"))   
  )

## 9.2 Full population model: Add the Activity_Type * Race_cat interaction item
cox_int_race <- coxph(
  Surv(time, status) ~ Activity_Type * Race_cat + Age + Sex +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat
)

summary(cox_int_race)

## 9.3 Global interaction Wald test（Activity_Type × Race_cat）
print("Global interaction verification（Activity_Type × Race_cat）:")
anova(cox_int_race, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 9.4 Subgroup model: Model One is drafted respectively in White and Non-white
cox_sub_white <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Race_cat == "White")
)
summary(cox_sub_white)

cox_sub_nonwhite <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex +
    Townsend_deprivation_index + Smoking_status +
    Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Race_cat == "Non-white")
)
summary(cox_sub_nonwhite)

# 10. Activity_Type × Townsend_deprivation_index_cat multiplication interaction
## 10.1 Make sure that Townsend_deprivation_index_cat is a factor and the order is clear
dat <- dat %>%
  mutate(
    Townsend_deprivation_index_cat = factor(Townsend_deprivation_index_cat,
                                            levels = c("Low", "Medium", "High"))
  )

## 10.2 Full population model: Add interaction items
cox_int_town <- coxph(
  Surv(time, status) ~ Activity_Type * Townsend_deprivation_index_cat + Age + Sex + Race +
    Smoking_status + Vitamin_D + BMI + Diet_quality,
  data = dat
)

summary(cox_int_town)

## 10.3 Global interaction Wald test
print("Global interaction verification（Activity_Type × Townsend_deprivation_index_cat）:")
anova(cox_int_town, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 10.4 Subgroup model: Fitted respectively in the Low/Medium/High subgroups
cox_sub_low <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Smoking_status + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Townsend_deprivation_index_cat == "Low")
)
summary(cox_sub_low)

cox_sub_medium <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Smoking_status + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Townsend_deprivation_index_cat == "Medium")
)
summary(cox_sub_medium)

cox_sub_high <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Smoking_status + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Townsend_deprivation_index_cat == "High")
)
summary(cox_sub_high)

# 11. Activity_Type × Smoking_status multiplicative interaction
## 11.1 Make sure the order of the Smoking_status factor
# dat$Smoking_status <- factor(dat$Smoking_status, levels = c("Never","Previous","Current"))

## 11.2 Full population model: Add interaction items
cox_int_smoke <- coxph(
  Surv(time, status) ~ Activity_Type * Smoking_status + Age + Sex + Race +
    Townsend_deprivation_index + Vitamin_D + BMI + Diet_quality,
  data = dat
)

summary(cox_int_smoke)

## 11.3 Global interaction Wald test (Activity_Type × Smoking_status)
print("Global interaction verification（Activity_Type × Smoking_status）:")
anova(cox_int_smoke, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 11.4 Subgroup models: Fit respectively in Never/Previous/Current
cox_sub_never <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Smoking_status == "Never")
)
summary(cox_sub_never)

cox_sub_previous <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Smoking_status == "Previous")
)
summary(cox_sub_previous)

cox_sub_current <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Vitamin_D + BMI + Diet_quality,
  data = dat %>% filter(Smoking_status == "Current")
)
summary(cox_sub_current)

# 12. Activity_Type × Diet_quality multiplicative interaction
## 12.2 Full population model: Add interaction items
cox_int_diet <- coxph(
  Surv(time, status) ~ Activity_Type * Diet_quality + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + BMI,
  data = dat
)

summary(cox_int_diet)

## 12.3 Global interaction Wald test (Activity_Type × Diet_quality)
print("Global interaction verification（Activity_Type × Diet_quality）:")
anova(cox_int_diet, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 12.4 Subgroup models: Fit respectively in Intermediate/Healthy/Unhealthy
cox_sub_intermediate <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + BMI,
  data = dat %>% filter(Diet_quality == "Intermediate")
)
summary(cox_sub_intermediate)

cox_sub_healthy <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + BMI,
  data = dat %>% filter(Diet_quality == "Healthy")
)
summary(cox_sub_healthy)

cox_sub_unhealthy <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + BMI,
  data = dat %>% filter(Diet_quality == "Unhealthy")
)
summary(cox_sub_unhealthy)

# 13. Activity_Type × Vitamin_D_cat multiplicative interaction
## 13.1 Make sure Vitamin_D_cat is a factor and the order is clear
dat <- dat %>%
  mutate(
    Vitamin_D_cat = factor(Vitamin_D_cat, levels = c("<50", "≥50"))
  )

## 13.2 Full population model: Add interaction items
cox_int_vitd <- coxph(
  Surv(time, status) ~ Activity_Type * Vitamin_D_cat + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + BMI + Diet_quality,
  data = dat
)

summary(cox_int_vitd)

## 13.3 Global interaction Wald test
print("Global interaction verification（Activity_Type × Vitamin_D_cat）:")
anova(cox_int_vitd, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 13.4 Subgroup models: Fit in the <50 and ≥50 groups respectively
cox_sub_under50 <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + BMI + Diet_quality,
  data = dat %>% filter(Vitamin_D_cat == "<50")
)
summary(cox_sub_under50)

cox_sub_over50 <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + BMI + Diet_quality,
  data = dat %>% filter(Vitamin_D_cat == "≥50")
)
summary(cox_sub_over50)

# 14. Activity_Type × BMI_cat multiplicative interaction
## 14.1 Make sure that BMI_cat is a factor and the order is clear
dat <- dat %>%
  mutate(
    BMI_cat = factor(BMI_cat, levels = c("≤25", ">25")) 
  )

## 14.2 Full population model: Add interaction items
cox_int_bmi <- coxph(
  Surv(time, status) ~ Activity_Type * BMI_cat + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + Diet_quality,
  data = dat
)

summary(cox_int_bmi)

## 14.3 Global interaction Wald test
print("Global interaction verification（Activity_Type × BMI_cat）:")
anova(cox_int_bmi, test = "Chisq") %>%
  slice_tail(n = 1) %>%          
  select(-loglik) %>%
  print()

## 14.4 Subgroup models: Fit in the ≤25 and >25 groups respectively
cox_sub_le25 <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + Diet_quality,
  data = dat %>% filter(BMI_cat == "≤25")
)
summary(cox_sub_le25)

cox_sub_gt25 <- coxph(
  Surv(time, status) ~ Activity_Type + Age + Sex + Race +
    Townsend_deprivation_index + Smoking_status + Vitamin_D + Diet_quality,
  data = dat %>% filter(BMI_cat == ">25")
)
summary(cox_sub_gt25)