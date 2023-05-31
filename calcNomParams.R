calcNomParams_fit = function(beta) {
  
  ####################################### Known Parameters ###############################################
  
  # ------ bigET
  #https://doi.org/10.1016/j.lfs.2012.08.008
  #Big ET-1    4 pg/ml in healthy subjects, 17 pg/ml in hemodialysis patients
  mw_BigET = 4283.0 #g/mol
  BigET0 = 4*1000/mw_BigET #pg/ml converted to pmol/l  https://doi.org/10.1016/j.lfs.2012.08.008  # value: 0.9339248
  
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1218999/
  #ECE - Big ET1 Kcat/Km = 4.6 - 6 /s/uM
  
  Km = 0.75*1e6 #pM
  Kcat = 3.3*60 #(/min)
  E_ECE_BigET = Kcat/Km  #L/pmol/min
  
  
  
  # ------ ET-1
  #https://doi.org/10.1016/j.lfs.2012.08.008
  #Plasma ET-1  1.5 +- 0.1 fmol/mL or pmol/L
  mw_ET1 = 2491.9 #g/mol
  ET1_cent0 = beta["ET1_cent0"] #3.2 #1.2 #pmol/L
  
  
  #Bacon and Davenport 1996
  #ET-1 is equipotent for ETA and ETB and will therefore
  #activate both
  #Kd = 0.4 #nM
  Kd_ET1 = 400 #pmol/L
  
  
  BigET_infusion_rate =0 #pmol/min
  
  ####################################### Parameters to be estimated ##########################################
  
  #Fraction of internalization that occurs through ETB vs ETA (1 - all ETB, 0 - all ETA)
  ETB_internalization_fraction = beta["ETB_internalization_fraction"] #0.8
  ETB_internalization_fraction_cent = beta["ETB_internalization_fraction_cent"] #0.8 #ETB_internalization_fraction #0.75
  
  V_cent = beta["V_cent"]  #L
  V_peri = beta["V_peri"]   #L
  V_bigET = beta["V_peri"]  #L  needed for BigET dosing in Hunter et al. 
  
  
  #10.1081/PRG-120024025
  #ECE activity 172 +/- 8 pmol ET/ml/h - human
  ECE_conc =beta["ECE_conc"] #172*1000/60 # 172/60 #pmol ET/L/min   # value: 3.666667
  
  
  K_ET1_pc = beta["K_ET1_pc"]
  K_ET1_cp = beta["K_ET1_cp"]
  
  Q_ET1_pc = K_ET1_pc*V_peri  #L/min
  Q_ET1_cp = K_ET1_cp*V_cent  #to be estimated
  
  
  Rtot_cent = beta["Rtot_cent"]
  
  
  Kint = beta["Kint"] #/min internalization rate
  
  Qint_cent= Kint*Rtot_cent*V_cent #pmol/min --> Kint*V_cent*Rtot_cent
  
  
  
  ####################################### Calculated Parameters ###############################################
  
  
  ECE_activity = ECE_conc*(E_ECE_BigET*BigET0) #pmol/L value: 13086.94   (pmol/L/min)/(L/min/pmol * pmol/L)
  BigET_prod_rate = ECE_activity*V_bigET  #pmol/min , from steady state condition:
  
  #Can solve for Qint_peri and ET1_peri0 from steady state constraints:
  #Let
  
  ET1_peri0 = (Qint_cent*(ET1_cent0/(Kd_ET1 + ET1_cent0)) + Q_ET1_cp*ET1_cent0) / Q_ET1_pc
  
  
  Qint_peri = (BigET_prod_rate - Q_ET1_pc*ET1_peri0 + Q_ET1_cp*ET1_cent0) /
    (ET1_peri0/(Kd_ET1 + ET1_peri0))
  
  Rtot_peri = Qint_peri / (Kint*V_peri)
  
  
  #ETB receptor concentrations, calculated from estimated internalization rates and volumes
  ETB_total_peri0 = Rtot_peri*ETB_internalization_fraction  #pmol/L total receptor concenntration
  ETB_total_cent0 = Rtot_cent*ETB_internalization_fraction_cent
  
  ETA_total_peri0 = Rtot_peri*(1-ETB_internalization_fraction)   #pmol/L total receptor concenntration
  ETA_total_cent0 = Rtot_cent*(1-ETB_internalization_fraction_cent)
  
  #Calculate Ligend-receptor complex
  ET1_ETA_cent0 = ETA_total_cent0*ET1_cent0/(Kd_ET1 + ET1_cent0) #
  ET1_ETB_cent0 = ETB_total_cent0*ET1_cent0/(Kd_ET1 + ET1_cent0) #pmol/L
  
  ET1_total_cent0 = ET1_cent0 + ET1_ETA_cent0 +  ET1_ETB_cent0 #pmol/L(total = free + bound)
  
  #pmol/L
  ET1_ETA_peri0 = ETA_total_peri0*ET1_peri0/(Kd_ET1 + ET1_peri0) #
  ET1_ETB_peri0 = ETB_total_peri0*ET1_peri0/(Kd_ET1 + ET1_peri0) #pmol/L
  ET1_total_peri0 = ET1_peri0 + ET1_ETA_peri0 +  ET1_ETB_peri0 #(total = free + bound)
  
  
  #ETA and ETA inhibitor parameters
  ETA_inhibition = 0
  ETB_inhibition = 0
  ETA_inhib_slope = 30 #min
  ETB_inhib_slope = 15 #min
  ETA_ramp = 0
  ETB_ramp = 0
  ETA_halflife = 60  #12 min
  ETB_halflife = 60  #60 min
  C_Ib_peri = 1 #delay constant
  
  
  Kia = 1#I/I_Kia
  Kib = 1#I/I_Kib
  
  I = 0 #inhibitor concentration
  
  t=sort(ls())
  param=sapply(t,names)
  for (i in 1:length(t)){
    param[i]=get(t[i])
  }
  param$param=NULL
  param = data.frame(param)
  return(param)
}







