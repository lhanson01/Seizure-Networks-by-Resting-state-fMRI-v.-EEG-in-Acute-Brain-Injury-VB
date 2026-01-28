#### Clean Data ####

new_col_names <- make.names(raw_data[2,], unique = TRUE)
renamed_frame <- raw_data
names(renamed_frame) <- new_col_names
choose_columns <- renamed_frame %>% select(ID,
                                           Sex,
                                           Admission.Date,
                                           Population,
                                           Age.at.Admit,
                                           Coded.Diagnosis,
                                           rs.fMRI.Date,
                                           Seizure.Network.Positive,
                                           Seizure.Positive,
                                           Start.Date) %>%
                                    slice(3:n())

agr_mcnem_frame <- choose_columns %>% filter(Seizure.Positive %in% c(0,1),
                                                  Seizure.Network.Positive %in% c(0,1)) %>%
                                      mutate(Neonates = ifelse(as.numeric(Age.at.Admit) < 29, 1, 0),
                                             Children = ifelse(as.numeric(Age.at.Admit) < 6570 & as.numeric(Age.at.Admit) > 28,
                                                               1, 0),
                                             Adults = ifelse(as.numeric(Age.at.Admit) > 6569, 1, 0)) %>%
                                      mutate(adm_to_eeg = as.numeric(Start.Date) - as.numeric(Admission.Date),
                                             eeg_to_mri = as.numeric(rs.fMRI.Date) - as.numeric(Start.Date)) %>%
                                      select(-Start.Date,
                                             -Admission.Date,
                                             -rs.fMRI.Date,
                                             -ID) %>%
                                      select(-Population)

################# Regression Processing #####################

#### Create Diagnosis Coding Dummy Variables ####

coding_vector <- 1:10
coding_variables <- as.data.frame(matrix(data = NA, nrow = nrow(agr_mcnem_frame), ncol = 10))
names(coding_variables) <- paste0("Diagnosis", 1:10)
lapply(1:nrow(agr_mcnem_frame), function(i){
  subj_code <- as.numeric(str_split(agr_mcnem_frame$Coded.Diagnosis[i], ",")[[1]])
  coding_variables[i,] <<- as.numeric(coding_vector %in% subj_code)
}
)

regression_frame <- agr_mcnem_frame %>%
                       mutate(agree = as.numeric(Seizure.Network.Positive == Seizure.Positive)) %>%
                       bind_cols(coding_variables) %>%
                       select(-Coded.Diagnosis)%>%#,
                              #-Diagnosis2,
                              #-Diagnosis7) %>%
                       #filter(abs(adm_to_eeg) < 300 & abs(eeg_to_mri) < 300) %>%
                       mutate(adm_to_eeg = scale(adm_to_eeg),
                              eeg_to_mri = scale(eeg_to_mri),
                              age_at_admit = scale(as.numeric(Age.at.Admit))) %>%
                       select(-Age.at.Admit,
                              -Seizure.Positive,
                              -Seizure.Network.Positive,
                              -Diagnosis8)




