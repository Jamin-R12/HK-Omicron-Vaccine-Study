
#------------------------------------------------------------------------------------#
# Prepares figure for fitted parameters, baseline risks, delta amplification risks and VE estimates
# 
# Panel A: Baseline risks by age group and outcome (log scale).
# Panel B: Peak-period relative risks by age and outcome, with a horizontal reference line at 1.
# Panel C: Vaccine effectiveness (VE) by dose, age group, and severity of outcomes, showing confidence intervals.
#------------------------------------------------------------------------------------#

setwd('/Users/byoung/Downloads/HK-Omicron-Vaccine-Study-main')
# READ IN DATA --------------------------------------------------------------------------#######

source("scripts/Model/Model_Inputs.R") # model input
data.input<-readRDS("temp_data/sample_med_init.rds") # Place holder, median  data. 

pacman::p_load('patchwork', 'tidyverse')


#CLEANING and Prepparing dataframes ------------------------------------------------------######

samples_med<-data.frame('indicator'  =names(data.input),
                        'value' =data.input)


p_symp=1-keppa_a1a5
p_infected=c(alpha_uv,alpha_v2, alpha_v3)
DFp_symp = data.frame(age = names(keppa_a1a5), p_symp = (1-keppa_a1a5))


samples_names<-samples_med%>%
  mutate(iter = row_number())%>%
  ungroup %>% 
  mutate(sev = gsub("^p([[:alpha:]]+).*", "\\1", indicator),
                  age = gsub("^.*\\.(.*)", "\\1", indicator),
         sev2= str_extract(indicator,paste0(Sev_indicators, collapse = '|')),  
           )%>%
  mutate(
    sev = factor(sev , levels = c('Hosp', 'SCC', 'Death')),
    eta = str_extract(indicator,"(?!beta)eta"),
    eta = ifelse(is.na(eta),'not',eta ),
    eta = ifelse(grepl('beta',indicator),'not',eta ),
    lrr = ifelse(grepl('lrr',indicator),'lrr', 'not'),
    beta = str_extract(indicator,'^beta'),
    beta = ifelse(grepl('delta',indicator),'beta',beta),
    beta = ifelse(grepl('iota',indicator),'beta',beta),
    beta = ifelse(grepl('seed',indicator),'beta',beta),
    
    beta = ifelse(is.na(beta),'not',beta )
  )

sample_plot2<-samples_names%>%
  filter(!indicator%in%c('beta', 'seed.inf')) %>% 
  mutate(Prob = 'Original') %>%
  bind_rows(samples_names%>%
              filter(!is.na(sev))%>%
              left_join(
                DFp_symp, by =c('age' = 'age')
              ) %>%
              mutate(
                value = ifelse(!is.na(age), value * p_symp, NA ),
                Prob = 'Incl. asymptomatic'
              )
  ) %>%
  mutate(
    age_actual = str_extract(age,'a\\d'),
    dose = str_extract(age,'d\\d'),
    subtitle = ifelse(is.na(sev),'Age-specific Beta', as.character(sev)),
    subtitle = factor(subtitle , levels = c('Hosp', 'SCC', 'Death')),
    age_actual = case_when(
      age_actual == 'a1' ~ '0-17',
      age_actual == 'a2' ~ '18-39',
      age_actual == 'a3' ~ '40-59',
      age_actual == 'a4' ~ '60-79',
      age_actual == 'a5' ~ '80+',
      T~NA))

sample_tab<-sample_plot2%>%
  group_by(indicator,Prob, subtitle, dose, age_actual) %>% 
  summarise(
    meanvalue = mean(value, na.rm= TRUE),
    medvalue = median(value, na.rm= TRUE),
    quant2.5 = quantile(value, 0.025, na.rm= TRUE),
    quant97.5 = quantile( value, 0.975, na.rm= TRUE),
    sdvalue = sd(value, na.rm= TRUE),
    N = n(),
    sev= first(sev),
    age= first(age)
  )


samp_VE_ests<-sample_plot2%>%
  filter(Prob=='Incl. asymptomatic'| grepl('lrr', indicator) ) %>% 
  mutate(sev = ifelse( is.na(sev),as.character(sev2), as.character(sev))) %>% 
  select(-p_symp,-indicator, -Prob,-sev2, -age,-eta,-beta,-Prob,-iter) %>% 
  group_by(sev,age_actual) %>%
  mutate(
    age_actual = ifelse( age_actual=='18-39',  '0-39', age_actual),
    age_actual = factor(age_actual, levels = c('0-39', '40-59', '60-79', '80+')),
    risk =  ifelse(dose=='d0',value, (value[dose=='d0']*value)),
    VE   =  1- (risk/risk[dose=='d0'])) %>%
  select( age_actual,sev, VE, dose) %>% 
  arrange(age_actual,sev, dose) %>% 
  mutate(
    dose = case_when(
      dose == 'd1' ~ 'First dose',
      dose == 'd2' ~ 'Second dose',
      dose == 'd3' ~ 'Booster dose',
      TRUE ~ dose  
    ),
    dose = factor(dose, levels = c('First dose','Second dose', 'Booster dose')),
    sev = factor(sev , levels = c('Hosp', 'SCC', 'Death'))) %>% 
  group_by(age_actual,sev,dose) %>%
  summarise(
    VE_med = median(VE, na.rm= TRUE),
    VE_lower = quantile(VE, 0.025, na.rm= TRUE),
    VE_upper = quantile(VE, 0.975, na.rm= TRUE),
    N = n()
  )


#FIGURES  --------------------------------------------------------------------#######


##Panel A: baseline risks  -----------------------------
BaseRisks_manu<-sample_plot2%>%
  filter(Prob=='Incl. asymptomatic') %>%
  group_by(indicator, subtitle, dose, age_actual) %>%
  summarise(
    medvalue = median(value, na.rm= TRUE),
    quant2.5 = quantile(value, 0.025, na.rm= TRUE),
    quant97.5 = quantile( value, 0.975, na.rm= TRUE),
    N = n()
  ) %>% 
  ungroup() %>% 
  mutate(
    dose = case_when(
      dose == 'd0' ~ 'Unvaccinated',
      dose == 'd2' ~ 'Second dose',
      dose == 'd3' ~ 'Booster dose',
      TRUE ~ dose 
    ),
    dose = factor(dose, levels = c('Unvaccinated', 'Second dose', 'Booster dose'))
  ) %>%
  
  mutate(   subtitle = case_when(
    subtitle == 'Hosp' ~ 'Hospitalisation',
    subtitle == 'SCC' ~ 'SCC',
    subtitle == 'Death' ~ 'Mortality',
    TRUE ~ subtitle  
  ),
  subtitle = factor(subtitle, levels = c('Hospitalisation', 'SCC', 'Mortality')))%>% 
  
  ggplot(
    aes(x=age_actual,y=medvalue, colour = subtitle)
  )+
  geom_point(position=position_dodge2(width = 0.5, preserve = 'single'))+
  labs(x='age group', y='Risk (log-scale)', colour = '')+
  scale_x_discrete(expand=c(0.1, 0))+
  scale_y_log10()+
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#41b6c4")) +
  theme_bw()+
  geom_errorbar(position=position_dodge2(preserve = 'single'), width=0.5, size=.5, aes(x=age_actual, ymin=quant2.5, ymax=quant97.5, colour = subtitle))+
  theme(legend.position = "bottom",
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12))+
  theme(strip.background = element_blank())



## Panel B: delta amplification risks  -----------------------------

delta_tab_manu<-samples_names%>%
  filter(grepl('delta',indicator)) %>%
  mutate(sev_type = str_extract(indicator,'hosp|SCC|Death')) %>%
  mutate(   sev_type = case_when(
    sev_type == 'hosp' ~ 'Hospitalisation',
    sev_type == 'SCC' ~ 'SCC',
    sev_type == 'Death' ~ 'Mortality',
    TRUE ~ sev_type  # If none of the above conditions match, keep the original value
  ),
  sev_type = factor(sev_type, levels = c('Hospitalisation', 'SCC', 'Mortality'))) %>%
  mutate(age_actual = str_extract(indicator,'\\d+'),
         age_actual = case_when(
           age_actual == '1' ~ '0-17',
           age_actual == '2' ~ '18-39',
           age_actual == '3' ~ '40-59',
           age_actual == '4' ~ '60-79',
           age_actual == '5' ~ '80+',
           T~NA)) %>%
  group_by(age_actual, sev_type) %>%
  summarise(medvalue = median(value ),
            quant2.5 = quantile(value, 0.025),
            quant97.5 = quantile(value, 0.975)
  ) %>%
  ggplot(
    aes(x=age_actual,y=medvalue, colour = sev_type)
  )+
  geom_point(position=position_dodge2(width = 0.5, preserve = 'single'))+
  labs(x='age group', y='Peak-period\nrelative risk (log-scale)', colour = "")+
  scale_x_discrete(expand=c(0.1, 0))+
  scale_y_log10()+
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#41b6c4")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbar(position=position_dodge2(preserve = 'single'), width=0.5, size=.5, aes(x=age_actual, ymin=quant2.5, ymax=quant97.5, colour = sev_type))+
  theme_bw()+
  
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12))



## Panel C: VE estimates----------------------------------------#######
samp_VE_gg<-samp_VE_ests %>% 
  filter(!dose=='d0') %>% 
  ggplot(aes(x=age_actual,y=VE_med, colour = sev))+
  geom_point(position=position_dodge2(width = 0.5, preserve = 'single'))+
  labs(x='Age group', y='Vaccine effectiveness', colour = 'Severe outcomes: ')+
  scale_x_discrete(expand=c(0.1, 0))+
  scale_y_continuous(labels = scales::percent)+
  # coord_cartesian(ylim = c(0.01, 1.2))+
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#41b6c4")) +
  theme_bw()+
  geom_errorbar(position=position_dodge2(width = 0.5, preserve = 'single'), width=0.5, size=.5, aes(x=age_actual, ymin=VE_lower, ymax=VE_upper, colour = sev))+
  facet_wrap(~dose,
             ncol = 3)+
  #space out the facets
  
  theme(legend.position = "bottom",
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # remove outside top and right border but keep axis line
        panel.border = element_blank(),
        #add x axis line for upper panel
        # Add lines only to left and bottom by specifying element_line()
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),  
        strip.text = element_text(size = 12, margin = margin(l = 50, r =50)), 
        panel.background = element_blank())



# COMBINE Plots --------------------------------------------------------------#########
#make top plot larger ratio adn remove its xaxis
RiskPlotsCmb<-(BaseRisks_manu+theme( legend.position = "none",
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_text(size=12),
                                     axis.text.x = element_text(size = 9),
                                     axis.text.y = element_text(size = 9),
                                     strip.text = element_text(size = 12)
)+
  delta_tab_manu+ theme(legend.position = "none",
                        axis.title.y = element_text(size=12),
                        axis.title.x = element_blank(),
                        axis.text.x = element_text(size = 9),
                        axis.text.y = element_text(size = 9)))/
  (samp_VE_gg + theme(legend.position = "bottom",
                      axis.title.y = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.x = element_text(size = 9),
                      axis.text.y = element_text(size = 9))
  )+#put guides on the botoom
  # theme(legend.position = "none")+
  plot_annotation(
    tag_levels = 'A',  # Use 'A' for the first plot
    tag_suffix = ")"   # Add a closing parenthesis after the label
  )

RiskPlotsCmb

