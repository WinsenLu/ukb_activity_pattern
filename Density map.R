library(tidyverse)
library(readr)

## 1. read data------------------------------------------
df <- read_csv("Active_MVPA.csv", locale = locale(encoding = "UTF-8"))

## 2. Variable-length data makes it convenient for ggplot to draw two curves simultaneously --------------------------
df_long <- df %>% 
  pivot_longer(cols = c(Top2_Days_MV_Minutes, Remaining_5_d),
               names_to = "Group",
               values_to = "MVPA_min")

## 3. paint-------------------------
p_final <- ggplot(df_long, aes(x = MVPA_min, fill = Group)) +
  
  geom_density(alpha = .7, colour = "black", na.rm = TRUE) +
  
  #geom_vline(xintercept = 150, linetype = "dashed", colour = "grey30") +
  
  scale_x_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, 250),
                     name = "Daily MVPA (min)") +
  scale_y_continuous(name = "Density") +
  
  scale_fill_manual(values = c("Top2_Days_MV_Minutes" = "#E41A1C",
                               "Remaining_5_d"        = "#377EB8"),
                    labels = c("Top-2 days",
                               "Remaining 5 days")) +
  
  guides(fill = guide_legend(title = NULL)) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid        = element_blank(),      
    panel.border      = element_blank(),      
    axis.line.x       = element_line(colour = "black"), 
    axis.line.y       = element_line(colour = "black"),  
    legend.position   = "top"
  )

print(p_final)