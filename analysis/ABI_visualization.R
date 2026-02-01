
var_names <- c("Stroke and Hemorrhage",
               "Sex",
               "Non-Epileptic Seizure Disorder",
               "Developmental Disorder",
               "Disorders of Consciousness",
               "Age at Admission",
               "Brain abnormalities",
               "EEG - rs-fMRI interval",
               "Hypoxic Brain Injury",
               "Traumatic Brain Injury",
               "Admission-EEG Interval",
               "Developmental Disorders",
               "Miscellaneous Neurological Cond.",
               "Intercept")



reg_results_ORs <- reg_results %>% select(starts_with("OR")) %>% arrange(OR) %>%
                                   mutate(rank = seq_len(n()),
                                          ymin = 0,
                                          ymax = ceiling(max(OR_upper)),
                                          variable = var_names,
                                          variable = factor(variable, levels = variable)) %>%
                                   filter(!(row.names(.) %in% "(Intercept)"))


#### Regression coefficient plots #####

#Odds
ggplot(reg_results_ORs, aes(x=variable, y = OR)) +
  geom_errorbar(aes(ymin = OR_lower, ymax=OR_upper), linewidth = 0.3) +
  geom_point(fill = "lightblue", size = 2, shape = 22, alpha = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() + labs(x = "Covariate", y = "Odds Ratio")

#Log Odds
ggplot(reg_results_ORs, aes(x=variable, y = OR_log)) +
  geom_errorbar(aes(ymin = OR_log_lower, ymax=OR_log_upper), linewidth = 0.3) +
  geom_point(fill = "lightblue", size = 2.5, shape = 21, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + labs(x = "Covariate", y = "Log of Odds Ratio")


# bar graph for population groups

data_w_age_group <- regression_frame %>%
                    mutate("Age Population" = ifelse(Neonates == 1, "Neonates",
                                  ifelse(Children == 1, "Children", "Adults")),
                           )

w_frequency <- data_w_age_group %>% count(`Age Population`, agree) %>%
                group_by(`Age Population`) %>% mutate(rel_freq = n / sum(n))


ggplot(w_frequency, aes(x = factor(`Age Population`, levels = c("Neonates", "Children", "Adults")),
                        y = rel_freq, fill = factor(agree, levels = c(0,1), labels = c("Disagree", "Agree")))) +
  geom_col() +  # Side-by-side proportions
  scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentages
  labs(x = "Age Population", y = "Percent Agreement", fill = "Modality Agreement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))

# bar graphs for TBI, HIE

w_frequency_HIE <- data_w_age_group %>% count(`Age Population`, Diagnosis3) %>%
  group_by(`Age Population`) %>% mutate(rel_freq = n / sum(n))

# bar graph for population groups
ggplot(w_frequency_HIE, aes(x = factor(`Age Population`, levels = c("Neonates", "Children", "Adults")),
                            y = rel_freq,
                            fill = factor(Diagnosis3, levels = c(0,1), labels = c("No HIE", "HIE")))) +
  geom_col() +  # Side-by-side proportions
  scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentages
  labs(x = "Age Population", y = "Percent HIE", fill = "Legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + scale_fill_viridis_d()

#TBI
w_frequency_TBI <- data_w_age_group %>% count(`Age Population`, Diagnosis4) %>%
  group_by(`Age Population`) %>% mutate(rel_freq = n / sum(n))

ggplot(w_frequency_TBI, aes(x = factor(`Age Population`, levels = c("Neonates", "Children", "Adults")),
                            y = rel_freq,
                            fill = factor(Diagnosis4, levels = c(0,1), labels = c("No TBI", "TBI")))) +
  geom_col() +  # Side-by-side proportions
  scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentages
  labs(x = "Age Population", y = "Percent TBI", fill = "Legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + scale_fill_viridis_d()

## Combined Frequency of Detection Bar Graphs ###
vis_frame <- data_w_age_group %>% mutate(EEG_positive = agr_mcnem_frame$Seizure.Positive,
                                         fMRI_positive = agr_mcnem_frame$Seizure.Network.Positive,
                                         EEG_positive = factor(EEG_positive, levels = c(0,1), labels = c("Negative", "Positive")),
                                         fMRI_positive = factor(fMRI_positive, levels = c(0,1), labels = c("Negative", "Positive")),
                                         "Age Population" = factor(`Age Population`, levels = c("Neonates", "Children", "Adults"))) %>%
                                         select(`Age Population`, EEG_positive, fMRI_positive) %>%
              rename("EEG Electrographic Seizure" = EEG_positive, "rs-fMRI SN identification" = fMRI_positive) %>%
              pivot_longer(cols = c(`EEG Electrographic Seizure`, `rs-fMRI SN identification`),
                           names_to = "Modality",
                           values_to = "Result",
                           ) %>%
              count(Modality, Result, `Age Population`) %>%
              mutate(Result = factor(Result, levels = c("Negative", "Positive")),
                     `Age Population` = factor(`Age Population`, levels = c("Neonates", "Children", "Adults")),
                     Group = factor(
                       paste(Result, `Age Population`),
                       levels = c(paste("Negative", c("Neonates", "Children", "Adults")),
                                  paste("Positive", c("Neonates", "Children", "Adults"))
                     )))


text_heights_vector <- numeric(length=nrow(vis_frame))
text_heights_vector[1] <- sum(vis_frame$n[which(vis_frame$Modality == unique(vis_frame$Modality)[1])])
for(i in 2:nrow(vis_frame)){
  if (vis_frame$Modality[i-1] != vis_frame$Modality[i]){
    text_heights_vector[i] <- text_heights_vector[1]
  } else {
    text_heights_vector[i] <- text_heights_vector[i-1] - vis_frame$n[i]
  }
}






ggplot(vis_frame, aes(x=Modality, y = n, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c(
    # Negative (blues)
    "Negative Neonates" = "#bdd7e7",
    "Negative Children" = "#6baed6",
    "Negative Adults"   = "#2171b5",

    # Positive (reds)
    "Positive Neonates" = "#fcae91",
    "Positive Children" = "#fb6a4a",
    "Positive Adults"   = "#cb181d"
  )) + #geom_text(aes(x=Modality, y = n, label = n))
  labs(y = "Number of Subjects") +
  theme(axis.text.x = element_text(size = 14))

