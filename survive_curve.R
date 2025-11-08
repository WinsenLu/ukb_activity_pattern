# 0. Nagvis -------------------------------------------------------------
library(readr)
library(dplyr)
library(survival)
library(ggplot2)
library(survminer)

# 1. up load ----------------------------------
dat <- read_csv(
  "MS_Cox.csv",
  col_types = cols(
    Date_G35_first_reported = col_character(),
    Date_of_death           = col_character(),
    baseline_date           = col_character()
  )
)

# 2. Manual Date conversion------------------------------------
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
    has_G35    = !is.na(Date_G35_first_reported),   
    time_years = time / 365.25
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

# 6. Drawing (no confidence interval, Y-axis percentage, no grid) --------------------------
surv_fit  <- survfit(Surv(time_years, status) ~ Activity_Type, data = dat)
surv_data <- surv_summary(surv_fit, data = dat) %>% 
  mutate(cuminc = (1 - surv) * 100)   

y_breaks <- seq(0, 0.12, by = 0.03)

ggplot(surv_data, aes(x = time, y = cuminc, colour = strata)) +
  geom_line(size = 1.2) +
  scale_colour_manual(
    name = "Activity Type",
    values = c("Activity_Type=Inactive"       = "#1F77B4",
               "Activity_Type=Active regular" = "#FF7F0E",
               "Activity_Type=Active WW"      = "#2CA02C"),
    labels = c("Inactive", "Active-Regular", "Active-Weekend Warrior")
  ) +
  labs(x = "Years",
       y = "Cumulative Incidence (%)",
       title = "Cumulative Incidence of Atrial Fibrillation by Activity Pattern") +
  scale_y_continuous(
    limits = c(0, 0.15),   
    breaks = y_breaks,
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )