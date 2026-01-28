#### Create Demographic Table ####
demo_table_age <- vector(mode="list", length = 0)
for(age_group in c("All", "Neonates", "Children", "Adults")){

  demo_table <- vector(mode = "list",
                       length = 6)

  demo_table <- as.list(as.logical(c(1,0,0,0,1,1)))

  # pick which variables to make demographics table
  names(demo_table) <- c("Sex",
                         "Age at Admission",
                         "Days between EEG and rs-fMRI",
                         "Admission to first EEG",
                         "Seizure Detected on EEG",
                         "SN Detected on rs-fMRI")



  demo_table_data_frame <- agr_mcnem_frame %>%
                          rename(`Age at Admission` = Age.at.Admit,
                                 `Days between EEG and rs-fMRI` = eeg_to_mri,
                                 `Admission to first EEG` = adm_to_eeg,
                                 `Seizure Detected on EEG` = Seizure.Positive,
                                 `SN Detected on rs-fMRI` = Seizure.Network.Positive) %>%
                          mutate(All = 1) %>%
                          filter(.data[[age_group]] == 1)

  demographic_variables <- names(demo_table)
  data_frame_vars <-names(demo_table_data_frame)


  for(i in seq_len(ncol(demo_table_data_frame))){
    if (data_frame_vars[i] %in% demographic_variables){
      if(demo_table[data_frame_vars[i]][[1]]){
        if(data_frame_vars[i] == "Sex"){
          sex <- demo_table_data_frame$Sex
          male_row <- c(sum(sex == 1),
                        sum(sex == 1)/length(sex))
          female_row <- c(sum(sex == 0),
                          sum(sex == 0)/length(sex))
          demo_table[["Sex"]] <- list("Males" = male_row,
                                    "Females" = female_row)
          } else {
            var <- demo_table_data_frame[,i][[1]]
            demo_table[data_frame_vars[i]] <- vector(mode = "list",
                                                                  length = 1)
            demo_table[data_frame_vars[i]][[1]] <- c(sum(var == 1),
                                                     sum(var == 1)/length(var))
          }
        } else {
          var <- as.numeric(demo_table_data_frame[,i][[1]])
          demo_table[data_frame_vars[i]] <- vector(mode = "list",
                                                                length = 1)
          demo_table[data_frame_vars[i]][[1]] <- c(median(var),
                                                            IQR(var))
        }
      }
  }
  demo_table_age[[age_group]] <- list(n = nrow(demo_table_data_frame),
                                      demo_table = demo_table)
}


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

agreement_mcnem <- function(age_group = "All", kappa = TRUE, B = 50000){
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
if (kappa == TRUE){
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
    return(list(results_frame = results_frame,
                agree_cont_table = agree_cont_table))
 }
 return(list(mcnemar_pvalue = mcnemar_pvalue,
             agree_cont_table = agree_cont_table))
}


### Set seed for bootstrap in kappa estimate, call function###
set.seed(1234)
estimate_frame <- bind_rows(lapply(c("All", "N", "P", "A"), function(x){
  agreement_mcnem(B = 50000, age_group = x)$results_frame
}
)
) %>% `row.names<-`(NULL)

contingency_tables <- lapply(lapply(c("All", "N", "P", "A"), agreement_mcnem, kappa = FALSE),
                             function(x){
                               return(x[[2]])
                             }
)

names(contingency_tables) <- c("All", "N", "P", "A")

print(contingency_tables)[[1]]

#### Agreement Logistic Regression ####

regr_frame_no_age_pop <- regression_frame %>% select(-Neonates,
                                                     -Children,
                                                     -Adults)

model <- glm(agree ~. ,
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

#### Agreement Logistic Regression with interaction ####
  interaction_model <- glm(agree ~ (.-age_at_admit) * age_at_admit,
               data = regr_frame_no_age_pop,
               family = "binomial")
  coefficients <- interaction_model$coefficients
  prof_lh_confints <- confint(interaction_model)
  p_values <- c(NA,Anova(interaction_model, test = "LR")[,3])
  interaction_reg_results <- data.frame(coefficients = coefficients,
                            OR = exp(coefficients),
                            OR_lower = exp(prof_lh_confints[,1]),
                            OR_upper = exp(prof_lh_confints[,2]),
                            p_values = p_values)







