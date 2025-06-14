# INTRODUCTION ------------------------------------------------------------

#Modified from Corbella et al 2017 script. 

# EPIDEMIC FUNCTION--
#       simulator of the daily severe health outcomes
#       (i.e. simulator of the discrete time 
#       deterministic SEEIR model and the observational
#       process)

'Function:
  1) seiir.model
  2) TransRep = 1) seiir.model + Observational process ( lines 613- 660)'

#PACKAGES##########
pacman::p_load(msm,car,tmvtnorm,deSolve,dplyr, pscl,surveillance,  readr, progress, tidyverse,matrixStats, data.table, coda,socialmixr)


###1) TRANSMISSION COMPONENT---------------------------------------------
# function to solve the SEEIIR deterministic continuous time system of 
# equation at specific time points.
#
# Args:
#   t      : vector of time at which to solve the system
#   x      : vector containing the starting state (S0, E10, Ia0, Is0, R0)
#   params : list of the parameters of the system 
#            containing transition rates (beta, sigma, gamma),
#             susceptibility, (alpha) and proportion asymptomatic (keppa)
#           
#
#

seiir.model <- function (t, x, params) {  
  Suv<- x[1:10]
  Sv1<- x[11:20]
  Sv2<- x[21:30]
  Sv3<- x[31:40]
  
  Euv<-x[41:50]
  Ev1<-x[51:60]
  Ev2<-x[61:70]
  Ev3<-x[71:80]
  
  Iauv<-x[81:90]
  Iav1<-x[91:100]
  Iav2<-x[101:110]
  Iav3<-x[111:120]
  
  Isuv<-x[121:130]
  Isv1<-x[131:140]
  Isv2<-x[141:150]
  Isv3<-x[151:160]
  
  Ruv<-x[161:170]
  Rv1<-x[171:180]
  Rv2<-x[181:190]
  Rv3<-x[191:200]
  
  with(                             
    as.list(params),                
    {
      
      t_r<-floor(t)
      
      # set hosprate for peak/non-peak periods
      if (t_r > 24 & t_r < 68 ) {
        hosprate_uv<-hosprate_uv[11:20]
        hosprate_v2<-hosprate_v2[11:20]
        hosprate_v3<-hosprate_v3[11:20]
        
      } else {
        hosprate_uv<-hosprate_uv[1:10]
        hosprate_v2<-hosprate_v2[1:10]
        hosprate_v3<-hosprate_v3[1:10]      
        
      }
      
      
      #set vaccination rates for age groups.
      # Need to scale vaccines, to include those that get vaccinated after infection, 25% case ascertainment, CANT get vac. some will still despite symp infections. assumes 60% of symptomatic cases get vaccinated, so need to scale. 
      
      ratuv_S  <- Suv/(Suv+Euv+Iauv+(Ruv*0.1)) #scale by 0.75 to account for 25% case ascertainment
      ratv1_S  <- Sv1/(Sv1+Ev1+Iav1+(Rv1*0.1))
      ratv2_S  <- Sv2/(Sv2+Ev2+Iav2+(Rv2*0.1))
      
      #if any are NaN, replace with 0
      ratuv_S[is.nan(ratuv_S)] <- 0
      ratv1_S[is.nan(ratv1_S)] <- 0
      ratv2_S[is.nan(ratv2_S)] <- 0
      
      #set vaccination rates variation
      pi_v1_mod <- pi_v1$pi_v1_bef * ratuv_S
      pi_v2_mod <- pi_v2$pi_v2_bef * ratv1_S
      pi_v3_mod <- pi_v3$pi_v3_bef * ratv2_S
      
   
      #children
      if(t_r>=Day_incrVR_d1_u5){ # day 37
        pi_v1_mod[1:2]<- pi_v1$pi_v1_aft[1:2]* ratuv_S[1:2] #u5 & 6-11
        
        if(t_r>=Day_incrVR_d2_u5){# 65
          pi_v2_mod[1:2]<-  pi_v2$pi_v2_aft[1:2]* ratv1_S[1:2] #u5 & 6-11
          
          if(t_r>=Day_incrVR_d3_u17){ # day 62
            pi_v3_mod[1:3]<- pi_v3$pi_v3_aft[1:3]* ratv2_S[1:3] #u5 & 6-11 & 12-17
            
            if(t_r>=Day_incrVR_d2_12_17){# 65
              pi_v1_mod[3]<- pi_v1$pi_v1_aft[3]* ratuv_S[3] #12-17
              pi_v2_mod[3]<- pi_v2$pi_v2_aft[3]* ratv1_S[3] #12-17
            }
          }
        }
      }
      
      #adults
      if(t_r>=Day_incrVR_d3_18_59){ # day 26  day 27
        pi_v3_mod[4:7]<- pi_v3$pi_v3_aft[4:7]* ratv2_S[4:7] #18 & 60-69
        pi_v2_mod[4:7]<- pi_v2$pi_v2_aft[4:7]* ratv1_S[4:7] #18 & 60+ actually Day_incrVR_d2_18_59 27.. red
        if(t_r>=Day_incrVR_d1_18_59){
          pi_v1_mod[4:7]<- pi_v1$pi_v1_aft[4:7]* ratuv_S[4:7] #18+ also day 52
        }
      }
      
      #seniors
      if(t_r>=Day_incrVR_d3_60_69){ # day 26 & 27
        pi_v3_mod[8]   <- pi_v3$pi_v3_aft[8]* ratv2_S[8] #60-69
        pi_v2_mod[8:10]<- pi_v2$pi_v2_aft[8:10]* ratv1_S[8:10]
        
        if(t_r>=Day_incrVR_d3_70pls){ # day 51 & 52
          
          pi_v1_mod[8:10]<- pi_v1$pi_v1_aft[8:10]* ratuv_S[8:10] #60 
          pi_v3_mod[9:10]<- pi_v3$pi_v3_aft[9:10]* ratv2_S[9:10] #70+
        }
      }
    
    
      
      
      
      #get time varying iota for respective contact matric
      iota<-time_varying_iota[t_r,] 
      
      #multiply each contact matrix by their respective iota element
      contact_school[prim_school,]    <- contact_school[prim_school,]*iota[1]
      contact_school[,prim_school]    <- contact_school[,prim_school]*iota[1]
      
      contact_school[sec_school,]     <- contact_school[sec_school,]*iota[2]
      contact_school[,sec_school]     <- contact_school[,sec_school]*iota[2]
      contact_school[teacher,teacher] <- contact_school[teacher,teacher]*mean(iota[1:2])
      
      contact_other                   <- contact_other*iota[3]
      contact_work                    <- contact_work*iota[4]
      contact_home                    <- contact_home*iota[5]
      contact_transport               <- contact_transport*iota[6]
      
      C_mat_iota<-contact_work+contact_home+contact_transport+contact_other+contact_school
      
      
      
      
      # FI X Cmat#########
      
      
      #multiply the vector across each row of the contact matrix
      FI_uv_uv   <- rowSums(C_mat_iota * matrix(Isuv * beta_d0, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v1   <- rowSums(C_mat_iota * matrix(Isv1 * beta_d0, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v2   <- rowSums(C_mat_iota * matrix(Isv2 * beta_d2, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v3   <- rowSums(C_mat_iota * matrix(Isv3 * beta_d3, nrow = 10,ncol = 10, byrow = TRUE))

      
      FI_uv_uvasp<- rowSums(C_mat_iota * matrix(Iauv * beta_d0asyp, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v1asp<- rowSums(C_mat_iota * matrix(Iav1 * beta_d0asyp, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v2asp<- rowSums(C_mat_iota * matrix(Iav2 * beta_d2asyp, nrow = 10,ncol = 10, byrow = TRUE))
      FI_uv_v3asp<- rowSums(C_mat_iota * matrix(Iav3 * beta_d3asyp, nrow = 10,ncol = 10, byrow = TRUE))
      
      #Force of infection
      # Multiply FOI by alpha (suscept of infected)
      FI_v1_uv <- FI_uv_uv
      FI_v1_v1 <- FI_uv_v1
      FI_v1_v2 <- FI_uv_v2
      FI_v1_v3 <- FI_uv_v3
      
      FI_v2_uv <- alpha_v2*FI_uv_uv
      FI_v2_v1 <- alpha_v2*FI_uv_v1
      FI_v2_v2 <- alpha_v2*FI_uv_v2
      FI_v2_v3 <- alpha_v2*FI_uv_v3
      
      FI_v3_uv <- alpha_v3*FI_uv_uv
      FI_v3_v1 <- alpha_v3*FI_uv_v1
      FI_v3_v2 <- alpha_v3*FI_uv_v2
      FI_v3_v3 <- alpha_v3*FI_uv_v3
      
      
      FI_v1_uvasp <- FI_uv_uvasp
      FI_v1_v1asp <- FI_uv_v1asp
      FI_v1_v2asp <- FI_uv_v2asp
      FI_v1_v3asp <- FI_uv_v3asp
      
      FI_v2_uvasp <- alpha_v2*FI_uv_uvasp
      FI_v2_v1asp <- alpha_v2*FI_uv_v1asp
      FI_v2_v2asp <- alpha_v2*FI_uv_v2asp
      FI_v2_v3asp <- alpha_v2*FI_uv_v3asp
      
      FI_v3_uvasp <- alpha_v3*FI_uv_uvasp
      FI_v3_v1asp <- alpha_v3*FI_uv_v1asp
      FI_v3_v2asp <- alpha_v3*FI_uv_v2asp
      FI_v3_v3asp <- alpha_v3*FI_uv_v3asp
      
      FI_uv       <- FI_uv_uv+FI_uv_uvasp+FI_uv_v1+FI_uv_v1asp+FI_uv_v2+FI_uv_v2asp+FI_uv_v3+FI_uv_v3asp
      FI_v1       <- FI_v1_uv+FI_v1_uvasp+FI_v1_v1+FI_v1_v1asp+FI_v1_v2+FI_v1_v2asp+FI_v1_v3+FI_v1_v3asp
      FI_v2       <- FI_v2_uv+FI_v2_uvasp+FI_v2_v1+FI_v2_v1asp+FI_v2_v2+FI_v2_v2asp+FI_v2_v3+FI_v2_v3asp
      FI_v3       <- FI_v3_uv+FI_v3_uvasp+FI_v3_v1+FI_v3_v1asp+FI_v3_v2+FI_v3_v2asp+FI_v3_v3+FI_v3_v3asp
      
      dSuv <- -(Suv*FI_uv) - (pi_v1_mod)
      dSv1 <- -(Sv1*FI_v1) + (pi_v1_mod) - (pi_v2_mod)
      dSv2 <- -(Sv2*FI_v2) + (pi_v2_mod) - (pi_v3_mod)
      dSv3 <- -(Sv3*FI_v3) + (pi_v3_mod)
      
      dEuv <- +(Suv*FI_uv) - (sigma*Euv)
      dEv1 <- +(Sv1*FI_v1) - (sigma*Ev1)
      dEv2 <- +(Sv2*FI_v2) - (sigma*Ev2)
      dEv3 <- +(Sv3*FI_v3) - (sigma*Ev3)
      
      
      dIauv <- (sigma*Euv*keppa_uv) - (gamma*Iauv)
      dIav1 <- (sigma*Ev1*keppa_uv) - (gamma*Iav1)
      dIav2 <- (sigma*Ev2*keppa_v2) - (gamma*Iav2)
      dIav3 <- (sigma*Ev3*keppa_v3) - (gamma*Iav3)
      
      #only remove hospital individuals from system, assume all hosp are also SCC & Death
      dIsuv <- (sigma*Euv*(1-keppa_uv)) - (Isuv*pHosp_uv*hosprate_uv)  - (Isuv*(1-pHosp_uv)* gamma)
      dIsv1 <- (sigma*Ev1*(1-keppa_uv)) - (Isv1*pHosp_v1*hosprate_uv)  - (Isv1*(1-pHosp_v1)* gamma)
      dIsv2 <- (sigma*Ev2*(1-keppa_v2)) - (Isv2*pHosp_v2*hosprate_v2)  - (Isv2*(1-pHosp_v2)* gamma) 
      dIsv3 <- (sigma*Ev3*(1-keppa_v3)) - (Isv3*pHosp_v3*hosprate_v3)  - (Isv3*(1-pHosp_v3)* gamma)
      
      dinfect.d.uv=(Suv*FI_uv)
      
      dinfect.d.v1=(Sv1*FI_v1)
      
      dinfect.d.v2=(Sv2*FI_v2)
      
      dinfect.d.v3=(Sv3*FI_v3)
      
      #hospital
      dHuv <- +(Isuv*pHosp_uv*hosprate_uv)
      dHv1 <- +(Isv1*pHosp_v1*hosprate_uv) #>> Use unvaccinated rate. 
      dHv2 <- +(Isv2*pHosp_v2*hosprate_v2)
      dHv3 <- +(Isv3*pHosp_v3*hosprate_v3)
      
      dRuv <- (gamma*Iauv) + (gamma*Isuv*(1-pHosp_uv))
      dRv1 <- (gamma*Iav1) + (gamma*Isv1*(1-pHosp_v1))
      dRv2 <- (gamma*Iav2) + (gamma*Isv2*(1-pHosp_v2))
      dRv3 <- (gamma*Iav3) + (gamma*Isv3*(1-pHosp_v3))
      
    
      dx <- c(dSuv,dSv1,dSv2,dSv3,
              dEuv,dEv1,dEv2,dEv3,
              dIauv,dIav1,dIav2,dIav3,
              dIsuv,dIsv1,dIsv2,dIsv3,
              dRuv,dRv1,dRv2,dRv3,
              dHuv,dHv1, dHv2, dHv3,
              dinfect.d.uv,dinfect.d.v1, dinfect.d.v2, dinfect.d.v3
      )
      
      list(dx)                    
    })
  
}















# Funciton to map the amplification factor delta for a set of amplification over time. 
Age_Spec_Slope<-function(DELTA, LENGTH){
  sequences_list <- map(DELTA, ~seq(1, .x, length.out = LENGTH))
  slope <- as.matrix(do.call(cbind, sequences_list))
  
  return(slope)
}

#Function to reduce the 10 bands to 5 age groups
Reduce10bands_a5bands<-function(mat10bands){
  colnames(mat10bands)<-gsub("\\.0_5","\\.a1",colnames(mat10bands))
  colnames(mat10bands)<-gsub("6_11","a1",colnames(mat10bands))
  colnames(mat10bands)<-gsub("12_17","a1",colnames(mat10bands))
  colnames(mat10bands)<-gsub("18_29","a2",colnames(mat10bands))
  colnames(mat10bands)<-gsub("30_39","a2",colnames(mat10bands))
  colnames(mat10bands)<-gsub("40_49","a3",colnames(mat10bands))
  colnames(mat10bands)<-gsub("50_59","a3",colnames(mat10bands))
  colnames(mat10bands)<-gsub("60_69","a4",colnames(mat10bands))
  colnames(mat10bands)<-gsub("70_79","a4",colnames(mat10bands))
  colnames(mat10bands)<-gsub("80\\+","a5",colnames(mat10bands))
  
  # reduce and sum where columns match a5bds_10bds
  mat5bands<-matrix(0,nrow=nrow(mat10bands),ncol=length(id_v))
  colnames(mat5bands)<-c(id_v[uvid_v],id_v[v1id_v], id_v[v2id_v],id_v[v3id_v])
  
  for(x_id in  id_v){
    mat5bands[,x_id]<-rowSums(matrix(mat10bands[,grepl(x_id,colnames(mat10bands))], nrow=nrow(mat10bands)))
  }
  return(mat5bands)
}






# EPIDMIC MODEL ------------------------------------------------------------- ############
# Combines SEIIR model (severe disease process) and the observationa process 
# Function to organise inputs, solve SEIIR (transmission component) where the Infected symptomatic or used to estiamt the severe disease incidence   
# Outputs the estimated Daily, Hosp, SCC, Deaths incidence
# epidemic model
#
#Parameters
# Args:
#   reff    : effective reprod. number 
#   ages    : age classes (child, adult & senior)
#   doses   : vac dose classes (0/uv, 2, 3)
#   ps      : dist of age and dose amonsgt pop.
#   alpha   : age and dose susceptibility
#   keppa   : age and dose proportion asymptomatic
#   rNaught : Basic reproduction number
#   end     : end of the epidemic
#   C_mat   : age-specific contact matrix
#   iseed   : Infections seeding
#   sigma   : Latent period
#   gamma   : infectious (symptomatic/asymp) period
#   PHSM    : Public health and safety measures, age-dose specfic reduction in contacts
#   iota    : magnifying factor of PHSMs
#   pi      : age specific vaccination rates, v2 getting 2nd dose, v3 getting boosted. 
#   phi_v2_reduc : reduction in onwards transmissible for 2 dose breaktrough cases
#   phi_v3_reduc : reduction in onwards transmissible for 3 dose breaktrough cases  
#   iota_times    :
#   iota_v        :



#SEIIR + OBS Model via convolution 

# # rm(N)
#   pHosp_base=pHosp_base0; pSCC_base=pSCC_base0; pDeath_base = pDeath_base0;beta_age = sample_med_init[1:5];phi_v2_reduc = sample_med_init[52];phi_v3_reduc = sample_med_init[53];iseed = sample_med_init[6];deltas_l = sample_med_init[99:113];iota_values =iota_values0;store_output = TRUE;pi_v1 = pi_v1_0;pi_v2 = pi_v2_0;pi_v3 = pi_v3_0;time_varying_iota = time_varying_iota_0;baseline = TRUE;base_2nd_cov=NULL;base_3rd_cov=NULL;  lrrHosp = lrrHosp0; lrrSCC = lrrSCC0; lrrDeath = lrrDeath0


TransRep <- function(
    pHosp_base, pSCC_base, pDeath_base,
    lrrHosp, lrrSCC, lrrDeath,
    beta_age,
    iseed ,
    deltas_l,
    store_output = FALSE,
    pi_v1 = pi_v1_0,
    pi_v2 = pi_v2_0,
    pi_v3 = pi_v3_0,
    time_varying_iota = time_varying_iota_0,
    baseline = TRUE,
    base_1st_cov=NULL,
    base_2nd_cov=NULL,
    base_3rd_cov=NULL)
{
  
  
  
  
  if(!is.null(base_2nd_cov)){
    
    #create new staring baseline population.
    
    start_d0_cov        <- (1-(base_1st_cov))
    start_d1_cov        <- (base_1st_cov-base_2nd_cov)
    start_d2_cov        <- (base_2nd_cov - base_3rd_cov)
    names(start_d0_cov) <- paste0(ages_10bds,'d0')
    # order vector by names with id_v10bds
    N.vacp.temp         <- c(start_d0_cov,start_d1_cov,start_d2_cov,base_3rd_cov)[match(id_v10bds,names(c(start_d0_cov,start_d1_cov,start_d2_cov,base_3rd_cov)))]
    #new starting pops. 
    N.pop<-rep(N,each = n.doses) *N.vacp.temp
    names(N.pop)<-id_v10bds
    
  }
  
  
  pop_dist   <- N.pop/sum(N.pop)
  
  
  #latent period
  sigma      <- 1/latent 
  gamma      <- 1/c(5.2,4.2,6.1)[IE]    # 6 days infectious period Yu 2023 10.1016/j.ijid.2023.02.011
  
  
  
  #ve effects##########
  
  #hosp------
  pHosp_l    = list(
    pHosp_base,
    pHosp_base*lrrHosp[v1lrr],
    pHosp_base*lrrHosp[v2lrr],
    pHosp_base*lrrHosp[v3lrr])
  
  #change to 10 elements for SEIR model. 
  pHosp_l_10 = list(pHosp_l[[1]][a5bds_10bds],
                    pHosp_l[[2]][a5bds_10bds],
                    pHosp_l[[3]][a5bds_10bds],
                    pHosp_l[[4]][a5bds_10bds])
  #scc------
  pSCC_l     = list(
    pSCC_base,
    pSCC_base*lrrSCC[v1lrr],
    pSCC_base*lrrSCC[v2lrr],
    pSCC_base*lrrSCC[v3lrr])
  
  
  #death------
  pDeath_l   = list(
    pDeath_base, #unvaccinated
    pDeath_base*lrrDeath[v1lrr],#v1
    pDeath_base*lrrDeath[v2lrr],#v2
    pDeath_base*lrrDeath[v3lrr])#v3
  
  
  
  
  ######deltas#########
  #hosp delta
  slope_hosp1<-Age_Spec_Slope(deltas_l[1:5], lngthup)
  slope_hosp2<-Age_Spec_Slope(deltas_l[1:5], lngthdn)
  #SCC delta
  slope_SCC1<-Age_Spec_Slope(deltas_l[6:10], lngthup)
  slope_SCC2<-Age_Spec_Slope(deltas_l[6:10], lngthdn)
  #Death delta
  slope_Death1<-Age_Spec_Slope(deltas_l[11:15], lngthup)
  slope_Death2<-Age_Spec_Slope(deltas_l[11:15], lngthdn)
  
  
  slopeUp_l<-list(slope_hosp1[,rep(1:5,times = 4)],
                  slope_SCC1[,rep(1:5,times = 4)],
                  slope_Death1[,rep(1:5,times = 4)])
  
  slopeDn_l<-list(slope_hosp2[,rep(1:5,times = 4)],
                  slope_SCC2[,rep(1:5,times = 4)],
                  slope_Death2[,rep(1:5,times = 4)])
  
  
  ##Transmission rates. 
  #beta calcs#########
  
  # beta = beta_age[1]
  beta   = beta_age[1] 
  
  beta_d0     <- beta/N #normalisig over population aspect.
  
  beta_d2     <- beta_d0 * 0.8 # 20% reduction infectivity
  beta_d3     <- beta_d0 * 0.7 # 30% reduction infectivity
  
  beta_d0asyp <- beta_d0 * phi_asymp_v
  beta_d2asyp <- beta_d2 * phi_asymp_v
  beta_d3asyp <- beta_d3 * phi_asymp_v
  
  
  
  #initial states###########
  I.0       <- round(pop_dist*iseed*0.5 , 0)
  
  Ia.0      <- round(I.0 * keppa, 0)
  Is.0      <- round(I.0 * (1-keppa),0) 
  
  E.0       <- round(pop_dist*iseed *0.5,0)
  R.0       <- 0*N.pop
  S.0       <- N.pop - E.0 - Ia.0 - Is.0 - R.0 
  
  
  times     <- seq(1, end , by=1)
  
  
  
  ###params#######
  params       <- list(sigma=sigma,
                       
                       gamma=gamma,
                       
                       beta_d0=beta_d0,
                       beta_d2=beta_d2,
                       beta_d3=beta_d3,
                       
                       beta_d0asyp=beta_d0asyp,
                       beta_d2asyp=beta_d2asyp,
                       beta_d3asyp=beta_d3asyp,
                       
                       keppa_uv=keppa_uv,
                       keppa_v2=keppa_v2,
                       keppa_v3=keppa_v3,
                       
                       alpha_v2=alpha_v2,
                       alpha_v3=alpha_v3,
                       
                       pi_v1 = pi_v1,
                       pi_v2 = pi_v2,
                       pi_v3 = pi_v3,
                       Day_incrVR_d1_u5 = Day_incrVR_d1_u5,
                       
                       Day_incrVR_d2_6_11 = Day_incrVR_d2_6_11,
                       Day_incrVR_d2_12_17 = Day_incrVR_d2_12_17,
                       Day_incrVR_d2_18_59 = Day_incrVR_d2_18_59,
                       
                       Day_incrVR_d3_u17 = Day_incrVR_d3_u17,
                       Day_incrVR_d3_18_59 = Day_incrVR_d3_18_59,
                       Day_incrVR_d3_70pls = Day_incrVR_d3_70pls,
                       
                       
                       
                       #peak and non-peak rates, age dose specific rates. 
                       #hosp delay d0
                       hosprate_uv = c(rep(med_hosp_days_peak$non_pk[1],7),med_hosp_days_peak$non_pk[2],med_hosp_days_peak$non_pk[2],med_hosp_days_peak$non_pk[3],
                                       rep(med_hosp_days_peak$peak[1],7),med_hosp_days_peak$peak[2],med_hosp_days_peak$peak[2],med_hosp_days_peak$peak[3]),
                       
                       #hosp delay d2
                       hosprate_v2 = c(rep(med_hosp_days_peak$non_pk[4],7),med_hosp_days_peak$non_pk[5],med_hosp_days_peak$non_pk[5],med_hosp_days_peak$non_pk[6],
                                       rep(med_hosp_days_peak$peak[4],7),med_hosp_days_peak$peak[5],med_hosp_days_peak$peak[5],med_hosp_days_peak$peak[6]),
                       
                       #hosp delay d3 
                       hosprate_v3 = c(rep(med_hosp_days_peak$non_pk[4],7),med_hosp_days_peak$non_pk[5],med_hosp_days_peak$non_pk[5],med_hosp_days_peak$non_pk[6],
                                       rep(med_hosp_days_peak$peak[4],7),med_hosp_days_peak$peak[5],med_hosp_days_peak$peak[5],med_hosp_days_peak$peak[6]),
                       
                       pHosp_uv    = pHosp_l_10[[1]],
                       pHosp_v1    = pHosp_l_10[[2]],
                       pHosp_v2    = pHosp_l_10[[3]],
                       pHosp_v3    = pHosp_l_10[[4]],
                       
                       
                       
                       time_varying_iota   = time_varying_iota, 
                       contact_home = contact_home,
                       contact_work = contact_work,
                       contact_other = contact_other,
                       contact_school = contact_school,
                       contact_transport = contact_transport,
                       prim_school = prim_school,
                       sec_school = sec_school
  )
  
  
  xstart       <- c(Suv=S.0[uvid_v10bds],Sv1 = S.0[v1id_v10bds],Sv2 = S.0[v2id_v10bds],Sv3 = S.0[v3id_v10bds],
                    Euv=E.0[uvid_v10bds], Ev1=E.0[v1id_v10bds], Ev2=E.0[v2id_v10bds], Ev3=E.0[v3id_v10bds],
                    Iauv=Ia.0[uvid_v10bds], Iav1=Ia.0[v1id_v10bds], Iav2=Ia.0[v2id_v10bds], Iav3=Ia.0[v3id_v10bds], 
                    Isuv=Is.0[uvid_v10bds], Isv1=Is.0[v1id_v10bds], Isv2=Is.0[v2id_v10bds], Isv3=Is.0[v3id_v10bds],
                    Ruv=R.0[uvid_v10bds], Rv1=R.0[v1id_v10bds], Rv2=R.0[v2id_v10bds], Rv3=R.0[v3id_v10bds],
                    Huv=R.0[uvid_v10bds], Hv1=R.0[v1id_v10bds], Hv2=R.0[v2id_v10bds], Hv3=R.0[v3id_v10bds],
                    
                    infect.d.uv =Is.0[uvid_v10bds]+Ia.0[uvid_v10bds],
                    infect.d.v1 =Is.0[v1id_v10bds]+Ia.0[v1id_v10bds],
                    infect.d.v2 =Is.0[v2id_v10bds]+Ia.0[v2id_v10bds],
                    infect.d.v3 =Is.0[v3id_v10bds]+Ia.0[v3id_v10bds]
  )
  
  
  
  
  
  ##desolve - ode rk4######
  out <- as.matrix(ode(xstart,times,seiir.model,params, method ='rk4' )) #desolve, solving equation.
  # out <- as.matrix(lsoda(xstart,times,seiir.model,params)) #desolve, solving equation.
  
  #extract symptomatice and asymptomatic infections
  D_Totinf_10 <-out[,grepl("^infect.d.*",colnames(out))]
  
  #get new daily increase in infections include 0 as starting point
  D_Totinf_10<-rbind(rep(0,ncol(D_Totinf_10)),D_Totinf_10)
  D_Totinf_10<-D_Totinf_10[-1,] - D_Totinf_10[-nrow(D_Totinf_10),]
  
  D_Totinf<-Reduce10bands_a5bands(D_Totinf_10)
  
  
  #for analysis script
  if(store_output==TRUE){
    #output for anaylsis.
    seir.L.out<<-list(out[,grepl("^S[a-z].*$",colnames(out))],
                      # out[,grepl("^infect.d..*",colnames(out))],
                      out[,grepl("^infect.d.*",colnames(out))] ,
                      out[,grepl("^R.*$",colnames(out))],
                      out[,grepl("^E.*$",colnames(out))],
                      out[,grepl("^H.*$",colnames(out))],
                      out[,grepl("^Ia.*$",colnames(out))]+out[,grepl("^Is.*$",colnames(out))])
    
    
  } else{
    #output for iterations
    
    # d_Sympinf<-d_Sympinf[rat_frtday:length(d_Sympinf)]
    # d_Asympinf<-d_Asympinf[rat_frtday:length(d_Asympinf)]
  }
  
  
 ##OBSERVATIONAL COMPONENT-------------------------------------------------------------------------------------#
   ##Severe Health - convolution ---output daily counts----
  
  DaysIndy_0<-matrix(0,nrow = nrow(D_Totinf), ncol = 20 )
  pIndy_l<-list(pHosp_l,pSCC_l,pDeath_l)
  # ind = Sev_indicators[1]
  
  DaysIndicators_l<-lapply(Sev_indicators, function(ind) {
    
    #reset output.
    DaysIndy<-DaysIndy_0
    DaysIndy<-cbind(1:nrow(DaysIndy), DaysIndy)
    
    #call appropriate indicator
    ind_index<-which(Sev_indicators == ind)
    NonPeakprobwtdays <- as.matrix(NonPeak_ProwtNDays_mat_list[[ind_index]])
    Peakprobwtdays <- as.matrix(Peak_ProwtNDays_mat_list[[ind_index]])
    
    #call appropriate prob/slope modifier.
    pIndy<-unlist(pIndy_l[[ind_index]])
    SlopeUp<-unlist(slopeUp_l[[ind_index]])
    SlopeDn<-unlist(slopeDn_l[[ind_index]])
    # row<-DaysIndy[3,]
    DaysIndy<-apply(DaysIndy, 1, function(row){
      row_result <- row[2:21]
      
      for(jIndy in 0:(nrow(NonPeakprobwtdays)-1)) {
        if(jIndy<row[1]){
          
          # during peak -------------------------------------------------------------
          if (row[1]>(peak_tstart) & row[1]<=(peak_all)){
            
            t_1  <- row[1]-(peak_tstart)
            # multiply element of vector by number
            pIndy_mod = pIndy * 0.96
            pIndy_mod = pIndy * (SlopeUp[t_1,])
            
          } else if (row[1]>peak_all & row[1]<=(peak_tend)){
            t_2<-(peak_tend+1)-row[1]
            pIndy_mod = pIndy * (SlopeDn[t_2,])
            
          }else{
            
            pIndy_mod = pIndy
            
          }
          
          if ( between(row[1], peak_tstart, peak_tend) ){
            row_result <- row_result + as.matrix(pIndy_mod*as.matrix(D_Totinf[row[1]-jIndy,])*Peakprobwtdays[(1+jIndy),])
          } else {
            row_result <- row_result + as.matrix(pIndy_mod*as.matrix(D_Totinf[row[1]-jIndy,])*NonPeakprobwtdays[(1+jIndy),])
          }
          
          
          
        }
      }
      return(row_result)
    })
    
    DaysIndy <- t(DaysIndy)
    colnames(DaysIndy)<-str_extract(colnames(D_Totinf), 'a\\dd\\d$')
    
    DaysIndy<-DaysIndy[,match(id_v,colnames(DaysIndy))]
    #rowsum a1 &a2 for each d1,d2,d3
    DaysIndy[,"a2d1"]<-rowSums(DaysIndy[,c("a2d1","a1d1")])
    DaysIndy[,"a2d2"]<-rowSums(DaysIndy[,c("a2d2","a1d2")])
    DaysIndy[,"a2d3"]<-rowSums(DaysIndy[,c("a2d3","a1d3")])
    DaysIndy<-DaysIndy[,!(colnames(DaysIndy) %in% c("a1d1","a1d2","a1d3"))]
    
    colnames(DaysIndy)<- paste0(tolower(ind), colnames(DaysIndy)) #leave out indi as to match easier later.
    
    
    return(DaysIndy)
  })
  #SH outcomes and Prevalence output.
  return(list(DaysIndicators_l))
  
}


##----------------------------------------------------------------------------##







