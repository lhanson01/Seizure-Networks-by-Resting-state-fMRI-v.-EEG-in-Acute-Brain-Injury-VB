#### Cohen's Kappa Function ####
c_kappa <- function(x1, x2){
  N <- length(x1)
  n01 <- N - sum(x1) #amount of negative rsfMRi scans
  n02 <- N - sum(x2) #amount of negative EEG scans
  n11 <- sum(x1) #amount of positive rsfMRI scans
  n12 <- sum(x2) #amount of positive EEG scans
  p_e <- 1/N^2 * (n01 * n02 + n11 * n12)
  p_0 <- sum(x1 == x2) / N
  k <- (p_0 - p_e) / (1 - p_e)
  return(k)
}

#### Estimate Agreement/Kappa, McNemar's Test####

agreement_mcnem <- function(age_group = "All", B){
  if (age_group == "N"){
    index <- as.logical(agr_mcnem_frame$Neonates)
  } else if (age_group == "P") {
    index <- as.logical(agr_mcnem_frame$Children)
  } else if (age_group == "A") {
    index <- as.logical(agr_mcnem_frame$Adults)
  } else {
    index <- as.logical(rep(1, times = nrow(agr_mcnem_frame)))
  }

eeg_pos <- as.numeric(agr_mcnem_frame$Seizure.Positive[index])
mri_pos <- as.numeric(agr_mcnem_frame$Seizure.Network.Positive[index])
#### McNemar's ####
  agree_cont_table <- table(eeg_pos,
                            mri_pos)
  mcnemar_pvalue <- mcnemar.test(agree_cont_table)$p.value

#### Kappa ####
  n <- sum(index)
  percent_agreement <- sum(eeg_pos == mri_pos) / n
  cohens_kappa <- c_kappa(eeg_pos, mri_pos)

  kappa_bs <- vector(mode = "numeric", length = B) # B Bootstrap samples
  for (i in 1:B){
    sample_indeces <- sample(seq_len(n), n, replace = TRUE)
    kappa_bs[i] <- c_kappa(eeg_pos[sample_indeces], mri_pos[sample_indeces])
  }

  kappa_lower_limit <- quantile(kappa_bs, 0.025)
  kappa_upper_limit <- quantile(kappa_bs, 0.975)

  results_frame <- data.frame(age_group = age_group,
                              mcnemar_pvalue = mcnemar_pvalue,
                              percent_agreement = percent_agreement,
                              cohens_kappa = cohens_kappa,
                              kappa_lower_limit = kappa_lower_limit,
                              kappa_upper_limit = kappa_upper_limit)
  return(results_frame)
}


### Set seed for bootstrap in kappa estimate, call function###
set.seed(1234)
estimate_frame <- bind_rows(lapply(c("All", "N", "P", "A"), function(x){
  agreement_mcnem(B = 50000, age_group = x)
}
)
) %>% `row.names<-`(NULL)

#### Agreement Logistic Regression ####
  regr_frame_no_age_pop <- regression_frame %>% select(-Neonates,
                                                       -Children,
                                                       -Adults)
  model <- glm(agree ~ (.-age_at_admit) * age_at_admit,
               data = regr_frame_no_age_pop,
               family = "binomial")
  coefficients <- model$coefficients
  prof_lh_confints <- confint(model)
  p_values <- c(NA,Anova(model, test = "LR")[,3])
  reg_results <- data.frame(coefficients = coefficients,
                            OR = exp(coefficients),
                            OR_lower = exp(prof_lh_confints[,1]),
                            OR_upper = exp(prof_lh_confints[,2]),
                            p_values = p_values)








