#### REGRESSION (using PA's code)
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

library(readxl)
GSVA_and_MEDIAN <- read_excel("results/GSVA_MEDIAN_comparison_table.xlsx")

GSVA_MED_REGRESSION <- lm(GSVA_MEAN ~ MEDIAN_MEAN, data = GSVA_and_MEDIAN)
summary(GSVA_MED_REGRESSION)

regressor_22 <-
 lm(formula = GSVA_MEAN ~  MEDIAN_MEAN, data = GSVA_and_MEDIAN)

fit_22 <- summary(regressor_22)

p_value_lm_22 <- fit_22$coefficients[2,4]

Rsquared_D2  <- fit_22$r.squared


#PLOT of GSVA and Mread_csv2()#PLOT of GSVA and MEDIAN
GSVA_and_MEDIAN %>%
  group_by(COINCIDENCE) %>%
  ggplot(aes(MEDIAN_MEAN, GSVA_MEAN)) +
  geom_point(colour = "grey") +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth",
              col = "black") +
  stat_regline_equation() +
  xlab("MEDIAN_MEAN") +
  ylab("GSVA_MEAN") + 
  ggtitle(paste0("x vs y, p=", round(p_value_lm_22, digits = 260))) +
  my_theme

##PA's plot:
GSVA_and_MEDIAN %>% ggplot(aes(x =MEDIAN_MEAN,
  y = GSVA_MEAN )) + 
  ggplot2::geom_point(colour = "grey40") +
  geom_smooth(method="lm", col="black") +
  stat_regline_equation() +
  ggtitle(paste0("x vs y, p=", round(p_value_lm_22, digits = 18), 
  "\nR² = ", round(Rsquared_D2, 2))) +
  xlab("MEDIAN_MEAN") +
  ylab("GSVA_MEAN") + 
  my_theme
