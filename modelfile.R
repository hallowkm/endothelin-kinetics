
#Define Model
ode <- "

#Concentration = Amount / Volume,  pmol/L
  BigET = BigET_amt/V_bigET
  ET1_total_peri = ET1_total_peri_amt/V_peri
  ET1_total_cent = ET1_total_cent_amt/V_cent
  ET1_total_peri_labeled =ET1_total_peri_labeled_amt/V_peri
  ET1_total_cent_labeled = ET1_total_cent_labeled_amt/V_cent


#ETA and ETB Inhibitor Effects
if (ETA_ramp == 1) {
  ETA_inhibitor_effect = min(ETA_inhibition, ETA_inhibition*time/ETA_inhib_slope)
} else {
  if (ETA_ramp ==-1 ) { #turn off drug
      ETA_inhibitor_effect = ETA_inhibition*exp(-time*log(2)/ETA_halflife)
  }   else {
      ETA_inhibitor_effect = ETA_inhibition
  }
}


if (ETB_ramp == 1) {
  ETB_inhibitor_effect = min(ETB_inhibition, ETB_inhibition*time/ETB_inhib_slope)
} else {
  if ( ETB_ramp == -1) {  #turn off drug over time
      ETB_inhibitor_effect = ETB_inhibition*exp(-time*log(2)/ETB_halflife)
 }   else {
      ETB_inhibitor_effect = ETB_inhibition
 }
}

I = ETB_inhibitor_effect + ETA_inhibitor_effect

I_Kia = I/Kia
I_Kib = I/Kib
I_Kib_peri = I_peri/Kib



#Assume constant receptor concentrations, pmol/L
  ETA_total_cent = ETA_total_cent0
  ETA_total_peri = ETA_total_peri0
  ETB_total_cent = ETB_total_cent0
  ETB_total_peri = ETB_total_peri0

#Receptor internalization rates
  Qint_cent= Kint*Rtot_cent*V_cent #pmol/min 
  Qint_peri= Kint*Rtot_peri*V_peri #pmol/min 
  
  QintB_peri = Qint_peri*(ETB_internalization_fraction)
  QintA_peri = Qint_peri*(1-ETB_internalization_fraction)
  QintB_cent = Qint_cent*(ETB_internalization_fraction_cent)
  QintA_cent = Qint_cent*(1-ETB_internalization_fraction_cent)


#Track radiolabeled and non-labeled entities, pmol/L
  ET1_total_cent_combined = ET1_total_cent + ET1_total_cent_labeled
  ET1_total_peri_combined = ET1_total_peri + ET1_total_peri_labeled
  ET1_cent_fraction_labeled = ET1_total_cent_labeled/ET1_total_cent_combined
  ET1_peri_fraction_labeled = ET1_total_peri_labeled/ET1_total_peri_combined


#Calculate Free ET1 from total ET1 and total receptor concentrations, under quasiequilibrium assumption (Mager 2005)
#pmol/L
  # ET1_peri_combined =(1/2)*(ET1_total_peri_combined - (ETA_total_peri + ETB_total_peri) - Kd_ET1*(1+I_Kia + I_Kib_peri) +
  #         sqrt( (ET1_total_peri_combined - (ETA_total_peri + ETB_total_peri) - Kd_ET1*(1+I_Kia + I_Kib_peri)  )^2 + 4*ET1_total_peri_combined*Kd_ET1*(1+I_Kia + I_Kib_peri)))
  # 
  # ET1_cent_combined =(1/2)*(ET1_total_cent_combined - (ETA_total_cent + ETB_total_cent) - Kd_ET1*(1+I_Kia + I_Kib) +
  #         sqrt( (ET1_total_cent_combined - (ETA_total_cent + ETB_total_cent) - Kd_ET1*(1+I_Kia + I_Kib)  )^2 + 4*ET1_total_cent_combined*Kd_ET1*(1+I_Kia + I_Kib)))


############# Solve for free ET-1 in peripheral compartment ################
#Equations result in a cubic polynomail (see manuscript)
#Use cubic formula to solve
#Cubic formula results in complex numbers in some steps, which cancel later
#However, since RxODE doesn't recognize complex numbers, it is necessary to convert to polar coordinates while solving

#Cubic polynomial coefficient
a1 = -1/Kd_ET1^2 
b1 = ET1_total_peri/Kd_ET1^2 - 2/Kd_ET1 - I_peri/(Kd_ET1*Kib) - I_peri/(Kd_ET1*Kia) - (1/Kd_ET1^2)*(ETA_total_peri + ETB_total_peri)
cc1 = 2*ET1_total_peri/Kd_ET1 + ET1_total_peri*I_peri/(Kd_ET1*Kia) + ET1_total_peri*I_peri/(Kd_ET1*Kib) - (1+I_peri/Kia+I_peri/Kib+I_peri^2/(Kia*Kib)) - (1/Kd_ET1)*(ETA_total_peri*(1+I_peri/Kib) + ETB_total_peri*(1+I_peri/Kia))
d1 = ET1_total_peri*(1+I_peri/Kia + I_peri /Kib + I_peri^2/(Kia*Kib))

#Make leading coefficient one
a = a1/a1
b = b1/a1
cc = cc1/a1
d = d1/a1

#Terms in cubic formula
term1 = (-b^3)/27 + b*cc/6 - d/2
term2 = cc/3 - (b^2)/9

#convert to polar form
x= term1
y = (sqrt(abs(term1^2 + term2^3)))
r = sqrt(x^2 + y^2)


theta1 = atan(y/x) + pi

if (x > 0) {
  theta1 = atan(y/x) 
}


#Calculate solution
ET1_peri_combined = 2*(r^(1/3))*cos(theta1/3) - b/3

############# Solve for free ET-1 in central compartment ################
#Equations result in a cubic polynomail (see manuscript)
#Use cubic formula to solve
#Cubic formula results in complex numbers in some steps, which cancel later
#However, since RxODE doesn't recognize complex numbers, it is necessary to convert to polar coordinates while solving

#Cubic polynomial coefficient
a2 = -1/Kd_ET1^2 
b2 = ET1_total_cent/Kd_ET1^2 - 2/Kd_ET1 - I/(Kd_ET1*Kib) - I/(Kd_ET1*Kia) - (1/Kd_ET1^2)*(ETA_total_cent + ETB_total_cent)
cc2 = 2*ET1_total_cent/Kd_ET1 + ET1_total_cent*I/(Kd_ET1*Kia) + ET1_total_cent*I/(Kd_ET1*Kib) - (1+I/Kia+I/Kib+I^2/(Kia*Kib)) - (1/Kd_ET1)*(ETA_total_cent*(1+I/Kib) + ETB_total_cent*(1+I/Kia))
d2 = ET1_total_cent*(1+I/Kia + I /Kib + I^2/(Kia*Kib))

#Make leading coefficient one
aa = a2/a2
bb = b2/a2
ccc = cc2/a2
dd = d2/a2

#Terms in cubic formula
term1c = (-bb^3)/27 + bb*ccc/6 - dd/2
term2c = ccc/3 - (bb^2)/9

#convert to polar form
xc= term1c
yc = (sqrt(abs(term1c^2 + term2c^3)))
rc = sqrt(xc^2 + yc^2)

#Handle discontinuities in atan when xc = 0 and changes signs
if (xc == 0) {
  xc = xc-0.0001
}

theta1c = atan(yc/xc) + pi

if (xc > 0) {
  theta1c = atan(yc/xc)   
}

# 
# if (term1c < 0) {
#   theta1c = theta1c+pi
# } else {
#   theta1c = theta1c - pi
# }

#Calculate solution
ET1_cent_combined = 2*(rc^(1/3))*cos(theta1c/3) - bb/3



#Separate unlabeled and radiolabeled portions
  ET1_cent = ET1_cent_combined*(1-ET1_cent_fraction_labeled)
  ET1_cent_labeled = ET1_cent_combined*ET1_cent_fraction_labeled
  ET1_peri = ET1_peri_combined*(1-ET1_peri_fraction_labeled)
  ET1_peri_labeled = ET1_peri_combined*ET1_peri_fraction_labeled


#Receptor-Bound Fraction
  ETA_cent_bound_fraction_combined =   ET1_cent_combined/(Kd_ET1*(1+I_Kia) + ET1_cent_combined)
  ETA_peri_bound_fraction_combined =   ET1_peri_combined/(Kd_ET1*(1+I_Kia) + ET1_peri_combined)
  ETB_cent_bound_fraction_combined =   ET1_cent_combined/(Kd_ET1*(1+I_Kib) + ET1_cent_combined)
  ETB_peri_bound_fraction_combined =   ET1_peri_combined/(Kd_ET1*(1+I_Kib_peri) + ET1_peri_combined)

#Ligand-receptor complex, pmol/L
  ET1_ETA_peri_combined = ETA_total_peri*ETA_peri_bound_fraction_combined
  ET1_ETB_peri_combined = ETB_total_peri*ETB_peri_bound_fraction_combined
  ET1_ETA_cent_combined = ETA_total_cent*ETA_cent_bound_fraction_combined
  ET1_ETB_cent_combined = ETB_total_cent*ETB_cent_bound_fraction_combined

#Separate unlabeled and radiolabeled portions
  ET1_ETA_cent = ET1_ETA_cent_combined*(1-ET1_cent_fraction_labeled)  
  ET1_ETA_cent_labeled = ET1_ETA_cent_combined*ET1_cent_fraction_labeled
  ET1_ETB_cent = ET1_ETB_cent_combined*(1-ET1_cent_fraction_labeled)
  ET1_ETB_cent_labeled = ET1_ETB_cent_combined*ET1_cent_fraction_labeled

  ET1_ETA_peri = ET1_ETA_peri_combined*(1-ET1_peri_fraction_labeled)
  ET1_ETA_peri_labeled = ET1_ETA_peri_combined*ET1_peri_fraction_labeled
  ET1_ETB_peri = ET1_ETB_peri_combined*(1-ET1_peri_fraction_labeled)
  ET1_ETB_peri_labeled = ET1_ETB_peri_combined*ET1_peri_fraction_labeled


#------------------------- ET1 production from Big-ET
#pmol/min
  ET1_production_from_BIGET =  E_ECE_BigET * BigET * ECE_conc * V_bigET

#------------------------- ET1 distribution between compartments 
#pmol/min

  Q_ET1_pc = K_ET1_pc*V_peri  #L/min
  Q_ET1_cp = K_ET1_cp*V_cent  #L/min
  
  ET1_distribution = - Q_ET1_pc*ET1_peri + Q_ET1_cp*ET1_cent      
  ET1_distribution_labeled =   - Q_ET1_pc*ET1_peri_labeled + Q_ET1_cp*ET1_cent_labeled

#------------------------- ET1 clearance by receptor internalization
#pmol/min
  ET1_ETB_internalization_peri = QintB_peri*ETB_peri_bound_fraction_combined*(1-ET1_peri_fraction_labeled)
  ET1_ETB_internalization_peri_labeled = QintB_peri*ETB_peri_bound_fraction_combined*ET1_peri_fraction_labeled
  ET1_ETA_internalization_peri = QintA_peri*ETA_peri_bound_fraction_combined*(1-ET1_peri_fraction_labeled)
  ET1_ETA_internalization_peri_labeled = QintA_peri*ETA_peri_bound_fraction_combined*ET1_peri_fraction_labeled

  
  ET1_ETB_internalization_cent = QintB_cent*ETB_cent_bound_fraction_combined*(1-ET1_cent_fraction_labeled)
  ET1_ETB_internalization_cent_labeled = QintB_cent*ETB_cent_bound_fraction_combined*(ET1_cent_fraction_labeled)
  ET1_ETA_internalization_cent = QintA_cent*ETA_cent_bound_fraction_combined*(1-ET1_cent_fraction_labeled)
  ET1_ETA_internalization_cent_labeled = QintA_cent*ETA_cent_bound_fraction_combined*(ET1_cent_fraction_labeled)


#------------------------------------ State Variables ------------------------------------------------------#


d/dt(BigET_amt) =  
   
    #Production of bigET
    BigET_prod_rate -  
    
    #Conversion to ET1
    ET1_production_from_BIGET; 
                  

d/dt(ET1_total_peri_amt)=   
  
    # production from bigET:
    ET1_production_from_BIGET +  
    
    # ET1 transfer between cent and peri:
    ET1_distribution -  

    #Internalization and clearance through receptor binding
    ET1_ETB_internalization_peri - ET1_ETA_internalization_peri
  
d/dt(ET1_total_cent_amt)=   

   # ET1 transfer between cent and peri:
    - ET1_distribution -
    
      #Internalization and clearance through receptor binding
    ET1_ETB_internalization_cent - ET1_ETA_internalization_cent

#Radiolabeled ET-1
d/dt(ET1_total_cent_labeled_amt) =
    
    # ET1 transfer between cent and peri:
    Q_ET1_pc*ET1_peri_labeled - Q_ET1_cp*ET1_cent_labeled - 
    
    #Internalization and clearance through receptor binding
    ET1_ETB_internalization_cent_labeled - ET1_ETA_internalization_cent_labeled
  
d/dt(ET1_total_peri_labeled_amt) =  

    # ET1 transfer between cent and peri:
    - Q_ET1_pc*ET1_peri_labeled + Q_ET1_cp*ET1_cent_labeled -

    #Internalization and clearance through receptor binding
    ET1_ETB_internalization_peri_labeled - ET1_ETA_internalization_peri_labeled
  
d/dt(I_peri) = C_Ib_peri*(I - I_peri)

"

#Compile model
mod = RxODE(model = ode)


