# Destructive GRR
# MSA: Destructive GR&R
# Example page 4-30
# Cliff Long
# 2025-01-22




# LOAD PACKAGES ###############################################################

library(here)
library(readxl)
library(janitor)
library(dplyr)
library(ggplot2)
library(broom)
library(pander)


# packages for GRR
library(SixSigma)


# packages for varcomp model
library(lme4)
library(lmerTest)





# LOAD FUNCTIONS ##############################################################

# normal probability plot with confidence bands
res_qqplot_fn <- function(fndat){
  require(qqplotr)
  
  fndat <- as.data.frame(fndat)
  names(fndat) <- "residuals"
  
  ggout <- ggplot(data = fndat, mapping = aes(sample = residuals)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  
  return(ggout)
}




# LOAD DATA ###################################################################

fpath <- "data Perez-Wilson DGRR.xlsx"
fsheet <- "data p4-36"
fsheet <- "data p4-30"

d0 <- read_excel(path = here('book perez wilson example', fpath), 
                 sheet = fsheet)

glimpse(d0)


# variable names
# 'operator' -- if 'appraiser' then change to 'operator'
# 'part'


d0 <- d0 %>% 
  rename(operator = appraiser)


# set specs -------------------------------------------------------------------

lsl <- NA
usl <- NA
tol <- usl - lsl

tol <- 2

study_sigma <- 6


# ANALYSIS ####################################################################

## prep data ==================================================================

glimpse(d0)

d1 <- d0 %>% 
  mutate(sample = factor(sample), 
         trial = factor(trial), 
         operator = factor(operator))

glimpse(d1)

xtabs(~ operator + sample, data = d1)
with(d1, interaction.plot(x.factor = sample, trace.factor = operator, response = reading))



## using aov ==================================================================


### use fixed effects anova to get anova table of SS --------------------------
fit_fixed <- aov(reading ~ operator / sample, data = d1)

# str(fit_fixed)

tidy_out_df <- broom::tidy(fit_fixed)

head(tidy_out_df)


### extract model MS terms and df ---------------------------------------------
# Perez-Wilson p4-24

MS_error <- tidy_out_df |> 
  dplyr::filter(term == 'Residuals') |> 
  select(all_of(c('df', 'meansq')))

print(MS_error)

MS_sample <- tidy_out_df |> 
  dplyr::filter(term == 'operator:sample') |> 
  select(all_of(c('df', 'meansq')))

print(MS_sample)

MS_operator <- tidy_out_df |> 
  dplyr::filter(term == 'operator') |> 
  select(all_of(c('df', 'meansq')))

print(MS_operator)


### get a, b, n by calculating from Df in ANOVA table -------------------------
# NOTE: should match known number of operators, parts, trials

# a = number of operators
# Df_O = a - 1
a <- MS_operator[, 'df'] + 1

print(paste("Operators a = ", a))

# b = number of parts (samples) per operator (not 'trials')
# Df_S = a*(b-1)
b <- ( MS_sample[, 'df'] / a ) + 1

print(paste("Samples b = ", b))

# n = number of trials (per sample)
n <- ( MS_error[, 'df'] / (a*b) ) + 1

print(paste("Trials n = ", n))



### compute GR&R variance components from MS and df terms ---------------------
# Perez-Wilson p4-28

# use brute-force in place of REML to adjust for negative variance component
check_neg_varcomp_fn <- function(fvarcomp){
  if (fvarcomp < 0) {fvarcomp <- 0}
  else {fvarcomp <- fvarcomp}
  
  return(unlist(fvarcomp))
}


varcomp_error <- check_neg_varcomp_fn( MS_error[, 'meansq'] ) 
varcomp_sample <- check_neg_varcomp_fn( (MS_sample[, 'meansq'] - MS_error[, 'meansq']) / n )
varcomp_operator <- check_neg_varcomp_fn( (MS_operator[, 'meansq'] - MS_sample[, 'meansq']) / (n*b) )
varcomp_total <- varcomp_error + varcomp_sample + varcomp_operator

print(paste("variance component ERROR = ", varcomp_error))
print(paste("variance component SAMPLE = ", varcomp_sample))
print(paste("variance component OPERATOR = ", varcomp_operator))
print(paste("variance component TOTAL = ", varcomp_total))


# as percentage of total varcomp
varcomp_error_pct <- varcomp_error / varcomp_total
varcomp_sample_pct <- varcomp_sample / varcomp_total
varcomp_operator_pct <- varcomp_operator / varcomp_total

print(paste("variance component ERROR PCT = ", varcomp_error_pct))
print(paste("variance component SAMPLE PCT = ", varcomp_sample_pct))
print(paste("variance component OPERATOR PCT = ", varcomp_operator_pct))



### compute sigma components from sqrt(varcomp) -------------------------------
# Perez-Wilson p4-28

sigma_error <- sqrt(varcomp_error)
sigma_sample <- sqrt(varcomp_sample)
sigma_operator <- sqrt(varcomp_operator)

print(paste("sigma error = ", sigma_error))
print(paste("sigma sample = ", sigma_sample))
print(paste("sigma operator = ", sigma_operator))


### compute GR&R terms from sigma components ----------------------------------

# see study_sigma above

GV <- sigma_error * study_sigma
PV <- sigma_sample * study_sigma
OV <- sigma_operator * study_sigma

print(paste("study_sigma = ", study_sigma))
print("")
print(paste("GV (gauge, error) = ", GV))
print(paste("PV (sample) = ", PV))
print(paste("OV (operator) = ", OV))




## using lme4 for variance components model ===================================
# https://stat.ethz.ch/~meier/teaching/anova/random-and-mixed-effects-models.html


# without interaction for DGRR - WITH REML
fit_lmer1 <- lmer(reading ~ 1 + (1|operator) + (1|operator:sample), REML = TRUE, data = d1)


# without interaction for DGRR - WITHOUT REML
# fit_lmer2 <- lmer(reading ~ 1+ (1|operator) + (1|operator:sample), REML = FALSE, data = d1)
## turned out to be not very interesting (uses likelihood function)



# assign just one model
fit_lmer <- fit_lmer1
# fit_lmer <- fit_lmer2


# GRR output - model summary
summary(fit_lmer)


# GRR - pseudo anova table
ranova(fit_lmer)

### THESE APPROX MATCH THE AOV OUTPUT ABOVE 

# GRR output - Variance Components
vcor <- data.frame(VarCorr(fit_lmer)) %>% select(-var2)

varcomp_tbl <- vcor %>% 
  select(-sdcor) %>% 
  rename(var_comp = vcov) %>% 
  rename(component = var1) %>% 
  mutate(pct_var_contrib = 100*var_comp/sum(var_comp)) %>% 
  mutate(component = case_when(grp == "sample" ~ "Part", 
                               grp == "operator" ~ "Reproducibility",
                               grp == "Residual" ~ "Error/Repeatability"))

varcomp_tbl %>% pander()


# GRR output - Study Variation

studyvar_tbl <- varcomp_tbl %>% 
  mutate(stddev_comp = sqrt(var_comp), 
         study_var = stddev_comp * study_sigma, 
         pct_study_var = 100*study_var/sum(study_var),
         pct_tolerance = 100*study_var/tol) %>% 
  select(-var_comp) %>% 
  select(-pct_var_contrib)

studyvar_tbl # %>% pander()


# verify model assumptions ####################################################

## pick only one
# fit1 <- fit_aov
# fit1 <- fit_lmer


# verify model assumptions - prep data
d1$fit <- fitted(fit1)
d1$resid <- residuals(fit1)
d1$stdres <- rstudent(fit1)


# plot residuals vs fitted values
d1 %>% 
  ggplot(aes(x = fit, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype = 2, color = 'blue') + 
  ggtitle("Residuals vs Fitted") + 
  NULL



# plot studentized residuals vs fitted values
d1 %>% 
  ggplot(aes(x = fit, y = stdres)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype = 2, color = 'blue') + 
  geom_hline(yintercept = c(-3, 3), linetype = 2, color = 'red') + 
  ggtitle("Studentized Residuals vs Fitted") + 
  NULL



# verify movel assumption - model residual normality
res_qqplot_fn(residuals(fit1))


# =============================================================================
# 7.  TRADITIONAL ANOVA-BASED EMS APPROACH  (VCA package)
# =============================================================================
# VCA fits the nested/crossed structure via ANOVA and directly
# provides EMS-based variance component estimates matching the EMS table.
#
# Formula notation:  y ~ material/(operator)/(sample)
#   material          = fixed
#   operator          = random, crossed with material
#   sample nested in material*operator

cat("\n===== VCA ANOVA-based Variance Components (EMS method) =====\n")

model_vca <- anovaVCA(
  reading ~ operator + operator:sample,
  Data = as.data.frame(d1)
)

print(model_vca)

# output is close to lmer model above


# END CODE ####################################################################
