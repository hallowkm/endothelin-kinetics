############################# Endothelin-1 Kinetics Model #########################################

#Authors: KM Hallow, University of Georgia
#June 22, 2022


#This file simulates and plots the following study:

# 
# Hunter, R.W., et al., First-in-Man Demonstration of Direct Endothelin-Mediated Natriuresis and Diuresis. 
# Hypertension, 2017. 70(1): p. 192-200.
# 
# Parker, J.D., et al., Human endothelin-1 clearance kinetics revealed by a radiotracer technique. 
# J Pharmacol Exp Ther, 1999. 289(1): p. 261-5.
# 
# Kaasjager, K.A., et al., Role of endothelin receptor subtypes in the systemic and renal responses to endothelin-1 in humans. 
# J Am Soc Nephrol, 1997. 8(1): p. 32-9.


###################################################################################################



library(RxODE)
library(tidyverse)
library(gridExtra)
library(ggpubr)

#Load Experimental Data
ETdat = read.csv("ExpDataET.csv")  
hunter = filter(ETdat, AuthorYear == "Hunter2017" & Variable == "Plasma ET1" )
kaas = filter(ETdat, AuthorYear == "Kaasjaver 1997" & Variable == "Plasma ET1")
parker = filter(ETdat, AuthorYear == "Parker1999" & Variable == "Labeled ET1" & Time <= 100 )

#Load model
source("modelfile.R")
mod <- RxODE(model = ode)
#Load fitted parameters
load("betafit.rds")

#Load other parameters
source("calcNomParams.R")






############------------ Simulate Parker 1999 ----------------################



#Use parker-specific ECE value
eta = beta
eta["ET1_cent0"] = 2.5
eta["ECE_conc"] = beta["ECE_conc"]*beta["ECE_adj_parker"]
theta = calcNomParams_fit(eta)

#Define initial conditions
inits = c(
  BigET = theta$BigET0,
  ET1_total_peri_amt = theta$ET1_total_peri0*theta$V_peri,
  ET1_total_cent_amt = theta$ET1_total_cent0*theta$V_cent,
  ET1_total_cent_labeled_amt = 0,
  ET1_total_peri_labeled_amt = 0,
  I_peri = 0
)

ev = eventTable()
ev = ev$add.sampling(parker$Time)

#Add dosing of radiolabeled ET-1, dose is arbitrary
dose = 0.0001
ev$add.dosing(dose = dose, dosing.to = 4, rate = dose/3, stop.time = 5)

#Simulate
x = data.frame(mod$run(theta, ev, inits))

#Plot
x$parker = parker$Value
maxval = max(x$parker)
h1 = ggplot(x) + 
  geom_point(aes(x=time, y = parker/maxval, color = "Parker1999")) + 
  geom_path(aes(x=time, y = ET1_cent_labeled/max(ET1_cent_labeled), color = "Model"),size= 1) +
  scale_y_continuous(trans='log', breaks = c(0.001, 0.01, 0.1, 1,10), limits = c(0.001, 10) ) +
  ylab("% of Peak Plasma ET-1") + 
  xlab("Time (min)") + 
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.7,0.8), legend.background = element_rect(fill='transparent')) + 
  annotate("segment", x = 0, xend = 0, y = 2e-3, yend = 6e-3,
           arrow = arrow(length = unit(0.3, "cm")), size = 0.75) +
  annotate("text", x= 0, y = 1.25e-3, label = "Radiolabeled\nET-1", fontface = "bold", size = 3, hjust = 0) + 
  labs(tag = "A")


############------------ Simulate Hunter -------------#############

eta = beta
eta["ET1_cent0"] = 0.5
theta = calcNomParams_fit(eta)

inits = c(
  BigET = theta$BigET0,
  ET1_total_peri_amt = theta$ET1_total_peri0*theta$V_peri,
  ET1_total_cent_amt = theta$ET1_total_cent0*theta$V_cent,
  ET1_total_cent_labeled_amt = 0,
  ET1_total_peri_labeled_amt = 0,
  I_peri = 0
)

evh = eventTable()
evh$add.sampling(seq(0, 180, by=1))
evh$add.dosing(dose = .75*30, dosing.to = 1, rate = .75)
evh$add.dosing(dose = 15*30, dosing.to = 1, rate = 15, start.time = 30)
evh$add.dosing(dose = 300*30, dosing.to = 1, rate = 300, start.time= 60)


xh = data.frame(mod$run(theta, evh, inits))

h2 = ggplot() + geom_path(data = xh, mapping = aes(x=time, y = ET1_cent - xh$ET1_cent[1], color = "Model"), size=1) +
  geom_point(data = hunter, mapping = aes(x=Time, y = Value*transformation, color = "Hunter 2017")) + 
  geom_errorbar(data = hunter, mapping = aes(x=Time, ymin = (Value - SD)*transformation, ymax =(Value + SD)*transformation), width = 10) +
  ylab("% of Peak Plasma ET-1") + 
  xlab("Time (min)") + 
  ylim(c(-0.8, 1)) + 
  scale_x_continuous( breaks = c(0, 30,60,90,120,150)) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8),legend.background = element_rect(fill='transparent')) + 
  annotate("rect", xmin = 30, xmax = 60, ymin = -0.5, ymax = -0.25,
           alpha = .3)+
  annotate("rect", xmin = 60, xmax = 90, ymin = -0.5, ymax = -0.25,
           alpha = .4)+
  annotate("rect", xmin = 0, xmax = 30, ymin = -0.5, ymax = -0.25,
           alpha = .1)+
  annotate("text",x=15, y = -0.375,label= "0.75" )+
  annotate("text",x=45, y = -0.375,label= "15" )+
  annotate("text",x=75, y = -0.375,label= "300" )+
  annotate("text", x= 45, y = -0.65, label = "Big ET-1 (pmol/min)", fontface = "bold", size = 3) +
  labs(tag = "B")


############------------ Simulate Kaasjager 1997  -------------#############

#Define parameters, use study-specific ECE
eta = beta
eta["ET1_cent0"] = 3.2
eta["ECE_conc"] = beta["ECE_conc"]*beta["ECE_adj_kaas"]
theta = calcNomParams_fit(eta)

#Define initial conditions
inits = c(
  BigET = theta$BigET0,
  ET1_total_peri_amt = theta$ET1_total_peri0*theta$V_peri,
  ET1_total_cent_amt = theta$ET1_total_cent0*theta$V_cent,
  ET1_total_cent_labeled_amt = 0,
  ET1_total_peri_labeled_amt = 0,
  I_peri = 0
)

evK = eventTable()
evK$add.sampling(seq(0,180+60))

#Specify ET-1 infusion
ET1_infusion_rate_cent = (0.5/theta$mw_ET1)*77*1000 #ng/kg/min * mol/g *kg *(1000 pmol/nmol) /L  = pmol/min/L  70 kg, 5L blood
evK$add.dosing(dose = ET1_infusion_rate_cent*60, rate = ET1_infusion_rate_cent, dosing.to = 3)
evK$add.dosing(dose = ET1_infusion_rate_cent*2*60, rate = ET1_infusion_rate_cent*2, dosing.to = 3, start.time = 60)
evK$add.dosing(dose = ET1_infusion_rate_cent*4*60, rate = ET1_infusion_rate_cent*4, dosing.to = 3, start.time = 120)

#Simulate
xk = data.frame(mod$run(theta, evK, inits))

#plot
h3 = ggplot() + geom_path(data = xk, mapping = aes(x=time, y = ET1_cent - xk$ET1_cent[1], color = "Model"), size = 1) +
  geom_point(data = kaas, mapping = aes(x=Time, y = (Value - kaas$Value[1])*transformation, color = "Kaasjager 1997"))+
  geom_errorbar(data = kaas, mapping = aes(x=Time, ymin = (Value - kaas$Value[1] - SD)*transformation, ymax =(Value - kaas$Value[1] + SD)*transformation), width = 20) +
  ylab("Plasma ET-1 (pmol/L)") + 
  xlab("Time (min)") + 
  ylim(c(-2.5,8))+
  scale_x_continuous( breaks = c(0, 30,60,90,120,150, 180, 210,240)) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8), legend.background = element_rect(fill='transparent'))+#, color = NA))+
  annotate("rect", xmin = 0, xmax = 60, ymin = -1.5, ymax = -0.5,
           alpha = .1)+
  annotate("rect", xmin = 60, xmax = 120, ymin = -1.5, ymax = -0.5,
           alpha = .2)+
  annotate("rect", xmin = 120, xmax = 180,  ymin = -1.5, ymax = -0.5,
           alpha = .4)+
  annotate("text",x=30, y = -1,label= "0.5" )+
  annotate("text",x=90, y = -1,label= "1.0" )+
  annotate("text",x=150, y = -1,label= "2.0" )+
  annotate("text", x= 90, y = -2, label = "ET-1 (ng/kg/min)", fontface = "bold", size = 3) + 
  labs(tag = "C")
h3

dev.new()
grid.arrange(h1, h2,h3)

ggarrange(h1,h2,h3,nrow = 3) %>% ggexport(filename = "calibrationFig.pdf", width = 5, height = 10)





