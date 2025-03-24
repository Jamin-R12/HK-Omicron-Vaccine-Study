# Model run



#-------------------------------------###
source("scripts/Model/Model_Inputs.R") # model initialisation
source("scripts/Model/Model_TransObs_ftns.R") # seir and observational model functions.
data.input<-readRDS("temp_data/sample_med_init.rds") #median calibrated data. 






##------------------------------------##
#RUN - EPI Func###############

# run SEIR and severe health outcome observational proccess
DaysIndicators_l<-TransRep(
    beta_age = data.input[1],
    iseed = data.input[2], 
    pHosp_base =data.input[3:7],
    pSCC_base =data.input[8:12],
    pDeath_base=data.input[13:17],
    lrrHosp = data.input[18:29], lrrSCC = data.input[30:41], lrrDeath = data.input[42:53],
    deltas_l = data.input[105:119],
    store_output=TRUE
  )




##-----------------------------------##
#PLOT model output



# Severe health outcomes by date of event.  
dt.list <- lapply(1:3, function(i) {
  dt <- as.data.table(cbind(day = 1:nrow(DaysIndicators_l[[1]][[1]]), DaysIndicators_l[[1]][[i]]))
  melt(dt, id.vars = 'day', variable.name = "name", value.name = "incidence")
})
# Combine all melted data tables into one
dt.lg <- rbindlist(dt.list)

SHO_incidence<-dt.lg%>%
  mutate(
    dose = str_extract(name, "d\\d"),
    age = str_extract(name, "a\\d"),
    indi = str_extract(name, "hosp|scc|death"),
    name2 = paste0(age,dose))




#Infection incidence from SEIR output
dt<-as.data.table(cbind(day = 1:nrow(seir.L.out[[2]]),seir.L.out[[2]]))
dt.inf <- melt(dt, id.vars = 'day', variable.name = "name", value.name = "value")


Inf_incidence<-dt.inf %>%
  mutate(
    dose = str_extract(name, "d\\d"),
    age = str_extract(name, "\\d+_\\d+|\\d+\\+"), # Extract age ranges like "0_5" or "80+"
    indi= 'infections',
    name2 = paste0(age,dose)
  )  %>% group_by(day, indi,name2)%>%
  summarise(value = sum(value))%>%
  group_by(name2) %>%
  arrange(day) %>%
  mutate(incidence = ifelse(row_number()==1,value, value - lag(value)))





#combing and plot. 
Model_inci<-SHO_incidence %>% 
  bind_rows(Inf_incidence) %>% 
  mutate(   
    indi = factor(indi, levels = c("infections", "hosp", "scc", "death"))
) %>% 
  group_by(day, indi)%>%
  summarise(incidence = sum(incidence))%>%
  ggplot(aes(x=day, y=incidence, colour = indi))+
  geom_point()+
  facet_wrap(~indi,
             scales = 'free',
             ncol = 1
             )+
  #set colour
  scale_colour_manual(values = c("#1f78b4", "#33a02c", "#41b6c4", "#ff7f00")) %>% 

  #make publishable
  theme(legend.position = 'none')+
  #strip plot backgrounds & grid
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # remove outside top and right border but keep axis line
        panel.border = element_blank(),
        # Add lines only to left and bottom by specifying element_line()
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),        
        panel.background = element_blank())

Model_inci



