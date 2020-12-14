ukb_df_1000 <- interested_ukb_df %>% head(n=1000)

# Demographics of a UKB sample subset ----

subgroup_of_interest <- (ukb_df_1000$body_mass_index_bmi_f21001_0_0 >= 25)

ukb_context(ukb_df_1000 , subset.var = subgroup_of_interest)

ukb_context(ukb_df_1000, nonmiss.var = "body_mass_index_bmi_f21001_0_0")

# Frequency of an ICD diagnosis by a target variable ----

ukb_icd_freq_by(ukb_df_1000, reference.var = "sex_f31_0_0", freq.plot = FALSE) %>% 
  select(-dx) %>% 
  tidyr::gather(key = "disease", value = "frequency", -categorized_var) %>%
  ggplot2::ggplot(aes(categorized_var, frequency, group = disease,
                      fill = disease)) +
  labs(x = "Reference variable", y = "UKB disease frequency", color = "", fill = "",
       title = "") +
  theme(title = element_text(face = "bold"), panel.grid = element_blank(),
        panel.background = element_rect(color = NULL,
                                        fill = alpha("grey", 0.10)),
        legend.key = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(2))+
  geom_bar(stat = "identity", position = "dodge") +
  guides(fill = guide_legend(ncol = 1), size = FALSE,
         color = FALSE) +
  scale_fill_discrete(labels =  c("coronary artery disease", "cerebrovascular disease",
                                  "lower respiratory tract infection"))



ukb_icd_freq_by(ukb_df_1000, reference.var = "body_mass_index_bmi_f21001_0_0", freq.plot = FALSE) %>% 
  select(-dx) %>% 
  dplyr::mutate(mid = (lower + upper) / 2) %>%
  tidyr::gather(key = "disease", value = "frequency", -categorized_var,
                -lower, -upper, -mid) %>%
  ggplot2::ggplot(aes(mid, frequency, group = disease, color = disease)) +
  labs(x = "Reference variable", y = "UKB disease frequency", color = "", fill = "",
       title = "")  +
  theme(title = element_text(face = "bold"), panel.grid = element_blank(),
        panel.background = element_rect(color = NULL,
                                        fill = alpha("grey", 0.10)),
        legend.key = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(2)) +
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  guides(color = guide_legend(ncol = 1), size = FALSE,
         fill = FALSE) +
  scale_color_discrete(labels = c("coronary artery disease", "cerebrovascular disease",
                                  "lower respiratory tract infection"))

ukb_df_1000 %>% 
#  summarise( Subjects = n() , Male =  sum(sex_f31_0_0 == "Male")/Subjects , 
#             Female = sum(sex_f31_0_0 == "Female")/Subjects) %>%
#  select(-Subjects) %>% 
#  gather(sex, percentage) %>% 
  ggplot(aes(x = sex_f31_0_0,fill=sex_f31_0_0 )) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  width=0.4) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))),
            stat = "count", vjust = -0.25)  +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x = "", y = "", color = "", fill = "",  title = "") +
  theme_pubr()
  theme(title = element_text(face = "bold"), panel.grid = element_blank(),
        panel.background = element_rect(color = NULL,
                                        fill = alpha("grey", 0.10)),
        legend.key = element_blank(), axis.ticks.x = element_blank()) 
  
