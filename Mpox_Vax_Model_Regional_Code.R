library(GillespieSSA2)
library(ggplot2)
library(tidyverse)
library(gdata)
library(reshape2)
library(tools)
# This code re-codes model into GillespieSSA2 package, which speeds up processing significantly.

# Regions by clade
# 
list_clade1b <- c("SudKivu","NordKivu","Tanganyika","HautKatanga",
                  "Lomami","Western","Eastern","Northern","Central",
                  "WesternSL","EasternSL","NorthernSL","SouthernSL","Buhumuza","Bujumbura",    
                  "Burunga","Butanyerera","Gitega")

list_clade1a1b <- c("KongoCentral","Kinshasa","Kasai","MaiNdombe","Tshopo")

list_clade1a <-  c("BasUele","Equateur","HautLomami","HautUele","Ituri",      
                   "KasaiCentral","KasaiOriental" ,"Kwango",       
                   "Kwilu", "Lualaba", "Maniema","Mongala",     
                   "NordUbangi","Sankuru" ,  "SudUbangi" ,"Tshuapa")

list_clade1bcountries <-  c("Western","Eastern","Northern","Central",
                            "Buhumuza","Bujumbura",    
                            "Burunga","Butanyerera","Gitega")

list_burundi <- c("Buhumuza","Bujumbura",    
                  "Burunga","Butanyerera","Gitega")

list_uganda <-  c("Western","Eastern","Northern","Central")


params <- c(
  rho = 1/21,      # Recovery rate
  beta_sex1a = r0_sex1a*(recoveryrate),
  beta_sex1b = r0_sex1b*(recoveryrate),
  mu_a_U5 = mortality_a_u5/21, #death rate, under 5
  mu_a_U15 = mortality_a_u15/21, #death rate, under 15
  mu_a_O15LR = mortality_a_o15/21, #death rate, over 15
  mu_a_O15HR = mortality_a_o15/21, #death rate, over 15
  mu_b_U5 = mortality_b_u5/21, #death rate, under 5
  mu_b_U15 = mortality_b_u15/21, #death rate, under 15
  mu_b_O15LR = mortality_b_o15/21, #death rate, over 15
  mu_b_O15HR = mortality_b_o15/21, #death rate, over 15
  omega_U5= 0, #vaccination rate, under 5
  omega_U15= 0, #vaccination rate, under 15
  omega_O15LR= 0, #vaccination rate, over 15 low risk
  omega_O15HR= 0,  #vaccination rate, over 15 high risk 
  exog_shock_multi_1a = exog_shock_multi_1a,
  r0_multi_clade1a=1
)


locations <- demog$location

age_groups <- c("U5", "U15", "O15HR", "O15LR")


y_init <- c()
for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age <- age_groups[b]
    location <- locations[c]
    
    clade = init_cases[init_cases$province==location,2]

    y_init[paste0("S_", age, "_", location)] <- demog[which(demog$location==location), age] -
      #Subtract infected clade 1a
      ifelse(clade=="clade1a" & age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")]),
             ifelse(clade=="clade1a1b"& age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")])*0.5,
                    ifelse(clade=="clade1a" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")]),
                           ifelse(clade=="clade1a1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")])*0.5,
                                  ifelse(clade=="clade1a" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")]),
                                         ifelse(clade=="clade1a1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")])*0.5,0))))))-
      #Subtract infected clade 1b
      ifelse(clade=="clade1b" & age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")]),
             ifelse(clade=="clade1a1b"& age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")])*0.5,
                    ifelse(clade=="clade1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")]),
                           ifelse(clade=="clade1a1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")])*0.5,
                                  ifelse(clade=="clade1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")]),
                                         ifelse(clade=="clade1a1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")])*0.5,0)))))) -
      #Subtract vaccinated
      ifelse((Vaccination_scenario==1) & age=="U5", 
             Vax_percent*Vax_effect*demog[which(demog$location==location), age],
             ifelse((Vaccination_scenario==2) & age=="O15HR",
                    Vax_percent*Adult_vax_prec*Vax_effect*demog[which(demog$location==location), "O15HR"],
                    ifelse((Vaccination_scenario==2) & age=="O15LR",
                           Vax_percent*(1-Adult_vax_prec)*Vax_effect*demog[which(demog$location==location), "O15HR"], 
                           ifelse(Vaccination_scenario==3 & location %in% c(list_clade1a,list_clade1a1b) & age=="U5",
                                  Vax_percent*Vax_effect*demog[which(demog$location==location), age],
                                  ifelse(Vaccination_scenario==3 & location %in% c(list_clade1b,list_clade1a1b) & age=="O15HR",
                                         Vax_percent*Adult_vax_prec*Vax_effect*demog[which(demog$location==location), "O15HR"],
                                         ifelse(Vaccination_scenario==3 & location %in% c(list_clade1b,list_clade1a1b) & age=="O15LR",
                                                Vax_percent*(1-Adult_vax_prec)*Vax_effect*demog[which(demog$location==location), "O15HR"],
                                                ifelse(Vaccination_scenario==4, Vax_percent*Vax_effect * demog[which(demog$location==location), age],
                                                       0))))))) -
      #Subtract recovered
      ifelse(age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_R")]),
             ifelse(age=="O15HR", round(as.integer(init_cases[which(init_cases$province==location),paste0("O15","_R")])*0.8,0),round(init_cases[which(init_cases$province==location),paste0("O15","_R")]*0.2,0))) -
      #Subtract dead
      ifelse(age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_D")]),
             ifelse(age=="O15HR", round(as.integer(init_cases[which(init_cases$province==location),paste0("O15","_D")])*0.8,0),round(init_cases[which(init_cases$province==location),paste0("O15","_D")]*0.2,0)))
    
    
    
    y_init[paste0("Ia_", age, "_", location)] <- ifelse(clade=="clade1a" & age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")]),
                                                        ifelse(clade=="clade1a1b"& age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")])*0.5,
                                                               ifelse(clade=="clade1a" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")]),
                                                                      ifelse(clade=="clade1a1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")])*0.5,
                                                                             ifelse(clade=="clade1a" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")]),
                                                                                    ifelse(clade=="clade1a1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")])*0.5,0))))))
    
    y_init[paste0("Ib_", age, "_", location)] <- ifelse(clade=="clade1b" & age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")]),
                                                        ifelse(clade=="clade1a1b"& age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_I")])*0.5,
                                                               ifelse(clade=="clade1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")]),
                                                                      ifelse(clade=="clade1a1b" & age %in% c("O15LR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15LR","_I")])*0.5,
                                                                             ifelse(clade=="clade1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")]),
                                                                                    ifelse(clade=="clade1a1b" & age %in% c("O15HR"),as.integer(init_cases[which(init_cases$province==location),paste0("O15HR","_I")])*0.5,0))))))
    
    
    y_init[paste0("V_", age, "_", location)] <- ifelse((Vaccination_scenario==1) & age=="U5", 
                                                       Vax_percent*Vax_effect*demog[which(demog$location==location), age],
                                                       ifelse((Vaccination_scenario==2) & age=="O15HR",
                                                              Vax_percent*Adult_vax_prec*Vax_effect*demog[which(demog$location==location), "O15HR"],
                                                              ifelse((Vaccination_scenario==2) & age=="O15LR",
                                                                     Vax_percent*(1-Adult_vax_prec)*Vax_effect*demog[which(demog$location==location), "O15HR"], 
                                                                     ifelse(Vaccination_scenario==3 & location %in% c(list_clade1a,list_clade1a1b) & age=="U5",
                                                                            Vax_percent*Vax_effect*demog[which(demog$location==location), age],
                                                                            ifelse(Vaccination_scenario==3 & location %in% c(list_clade1b,list_clade1a1b) & age=="O15HR",
                                                                                   Vax_percent*Adult_vax_prec*Vax_effect*demog[which(demog$location==location), "O15HR"],
                                                                                   ifelse(Vaccination_scenario==3 & location %in% c(list_clade1b,list_clade1a1b) & age=="O15LR",
                                                                                          Vax_percent*(1-Adult_vax_prec)*Vax_effect*demog[which(demog$location==location), "O15HR"],
                                                                                          ifelse(Vaccination_scenario==4, Vax_percent*Vax_effect * demog[which(demog$location==location), age],
                                                                                                 0)))))))
    
    y_init[paste0("R_", age, "_", location)] <- ifelse(age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_R")]),
                                                       ifelse(age=="O15HR", round(as.integer(init_cases[which(init_cases$province==location),paste0("O15","_R")])*0.8,0),round(init_cases[which(init_cases$province==location),paste0("O15","_R")]*0.2,0)))
    
    y_init[paste0("D_", age, "_", location)] <- ifelse(age %in% c("U5","U15"), as.integer(init_cases[which(init_cases$province==location),paste0(age,"_D")]),
                                                       ifelse(age=="O15HR", round(as.integer(init_cases[which(init_cases$province==location),paste0("O15","_D")])*0.8,0),round(init_cases[which(init_cases$province==location),paste0("O15","_D")]*0.2,0)))
  }
}

y_init <- round(y_init,0)


# Create reactions and effects datasets

reactions= list()
nu <- matrix(0, nrow = 0, ncol = length(y_init))
colnames(nu) <- names(y_init)


# Create all reactions and effects that happen within the model

for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age<- age_groups[b]
    location <- locations[c]
    
    
    # With R0s 
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(",location,"_R0_",age, "* S_", age, "_", location, " * Ia_", age_groups[1], "_", location, 
             " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+(",location,"_R0_",age,  "* S_", age, "_", location, " * Ia_", age_groups[2], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (",location,"_R0_",age,  "* S_", age, "_", location, " * Ia_", age_groups[3], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (",location,"_R0_",age,  "* S_", age, "_", location, " * Ia_", age_groups[4], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    # Human-to-human transmission: clade 1b
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(r0_multi_clade1a*",location,"_R0_",age, "* S_", age, "_", location, " * Ib_", age_groups[1], "_", location, 
             " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+(r0_multi_clade1a*",location,"_R0_",age,  "* S_", age, "_", location, " * Ib_", age_groups[2], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (r0_multi_clade1a*",location,"_R0_",age, "* S_", age, "_", location, " * Ib_", age_groups[3], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (r0_multi_clade1a*",location,"_R0_",age, "* S_", age, "_", location, " * Ib_", age_groups[4], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    # 
   
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(0.95*",location, "_exogshock_", age,"* S_", age, "_", location,") - (0.95*",location, "_exogshock_", age,"* V_", age, "_", location,")"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
    
    if (b>2) {
      
      for (i in 1:length(locations)) {
        
        other_country <- locations[i]
        
        if (location %in% c("BasUele","Equateur","HautLomami","HautUele","Ituri",
                            "Kasai","KasaiCentral","KasaiOriental","Kinshasa","KongoCentral",
                            "Kwango","Kwilu","Lualaba","MaiNdombe","Maniema",
                            "Mongala" ,"NordUbangi","Sankuru","SudUbangi","Tshopo",
                            "Tshuapa"))
          
          
          
        {
          
          reactions <- append(reactions, list(c(
            #paste0("S_", age, "_", location, " -> I_", age, "_", location),
            paste0(location,"_", other_country, "*0.2* S_",age,"_", location, "*(Ia_", age_groups[1], "_", other_country,
                   "+Ia_", age_groups[2], "_", other_country,"+ Ia_", age_groups[3], "_", other_country,"+ Ia_", age_groups[4], "_", other_country,')'))))
          
          
          new_row <- rep(0, length(y_init))
          names(new_row) <- names(y_init)
          new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
          new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
          nu <- rbind(nu, new_row)
          
        }
        
        if (location %in% locations_1b_2)
          # if (locations %in% c("Burundi","CAR","Kenya","Rwanda","Uganda","HautKatanga",
          #                     "Kasai","Kinshasa","KongoCentral","Lomami","MaiNdombe","NordKivu",
          #                     "SudKivu","Tanganyika","Tshopo"))
          
          
          
          
        {
          
          reactions <- append(reactions, list(c(
            #paste0("S_", age, "_", location, " -> I_", age, "_", location),
            paste0(location,"_", other_country, "*0.2* S_",age,"_", location, "*(Ib_", age_groups[1], "_", other_country,
                   "+Ib_", age_groups[2], "_", other_country,"+ Ib_", age_groups[3], "_", other_country,"+ Ib_", age_groups[4], "_", other_country,')'))))
          
          
          new_row <- rep(0, length(y_init))
          names(new_row) <- names(y_init)
          new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
          new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
          nu <- rbind(nu, new_row)
          
        }
      }
      
    }
    
    
    # Sexual transmission among high-risk adults
    
    if (b == 3 && any(!(location %in% c(list_uganda)))) {
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("beta_sex1a* S_", age, "_", location, "*Ia_",age,"_",location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
               "+Ib_", age, "_", location,"+R_", age, "_", location,
               "+V_", age, "_", location,")"))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("beta_sex1b* S_", age, "_", location, "*Ib_",age,"_",location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
               "+Ib_", age, "_", location,"+R_", age, "_", location,
               "+V_", age, "_", location,")"))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    
    if (b == 3 && any((location %in% c(list_uganda)))) {
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("beta_sex1a* S_", age, "_", location, "*Ia_",age,"_",location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
               "+Ib_", age, "_", location,"+R_", age, "_", location,
               "+V_", age, "_", location,")"))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("beta_sex1b* S_", age, "_", location, "*Ib_",age,"_",location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
               "+Ib_", age, "_", location,"+R_", age, "_", location,
               "+V_", age, "_", location,")"))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    
    # Vaccination
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> V", age, "_", location),
      paste0("omega * S_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("V_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
    # Recovery, clade 1a
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
      paste0("rho * Ia_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
    # Recovery, clade 1b
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
      paste0("rho * Ib_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
    # Death, clade 1a
    
    if (any(!(location %in% c("Sankuru","SudUbangi")))) {
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> D", age, "_", location),
        paste0("mu_a_",age ,"* Ia_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    if (any((location %in% c("Sankuru","SudUbangi")))) {
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> D", age, "_", location),
        paste0("0.3*mu_a_",age ,"* Ia_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    # Death, clade 1b
    
    if (any((location %in% c(list_uganda,"Tshopo")))) {
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> D", age, "_", location),
        paste0("10*mu_b_", age, "* Ib_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    if (any(!(location %in% c(list_uganda, "Tshopo","SudKivu")))) {
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> D", age, "_", location),
        paste0("mu_b_", age, "* Ib_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    if (any((location %in% c("SudKivu")))) {
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> D", age, "_", location),
        paste0("0.4*mu_b_", age, "* Ib_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
  }
}

reactions2 <- as.character(reactions)
nu2 <- t(nu)
#
out <-
  GillespieSSA2::ssa(
    initial_state = y_init,
    reactions = port_reactions(x0= y_init, a = reactions2, nu = nu2),
    params = params7,
    method = ssa_exact(),
    final_time = 20, # running for 273 days: July 28, 2024 (seeded case) to April 27, 2025 (DRC total case data)
    census_interval = 0.001,#.001,
    verbose = FALSE
  )

max(out$time)

data <- out$state



