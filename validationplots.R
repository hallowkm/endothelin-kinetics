############################# Endothelin-1 Kinetics Model #########################################

#Authors: KM Hallow, University of Georgia
#June 22, 2022


#This file simulates and plots the following study:

# 
# Bohm, F., et al., ETA receptors mediate vasoconstriction, whereas ETB receptors clear endothelin-1 in the splanchnic and renal circulation of healthy men. 
# Clin Sci (Lond), 2003. 104(2): p. 143-51.


###################################################################################################


setwd("C:/models/renalmodel/Endothelin/Endothelin Kinetics - Copy/PublicationScripts")
#Load experimental data
ETdat = read.csv("ExpDataET.csv")  
bohm = filter(ETdat, AuthorYear == "Bohm2003" & Variable == "ET1")


#Load model
source("modelfile.R")

#Load fitted parameters
load("betafit.rds")


#Load other parameters
source("calcNomParams.R")


eta = beta

eta["ET1_cent0"] = 3.2
eta["ECE_conc"] = beta["ECE_conc"]*0.4 
eta["ETB_internalization_fraction_cent"] = 0.8
theta = calcNomParams_fit(eta)

theta$ETA_inhib_slope = 30
theta$ETB_inhib_slope = 15
theta$ETA_halflife = 60
theta$C_Ib_peri = 1


inits = c(
  BigET_amt = theta$BigET0*theta$V_bigET,
  ET1_total_peri_amt = theta$ET1_total_peri0*theta$V_peri,
  ET1_total_cent_amt = theta$ET1_total_cent0*theta$V_cent,
  ET1_total_cent_labeled_amt = 0,
  ET1_total_peri_labeled_amt = 0,
  I_peri = 0
)

evf1 = eventTable()
evf1$add.sampling(seq(0,30,by=0.01))

evf2 = eventTable()
evf2$add.sampling(seq(0,15,by=0.1))


#Set up event table
evf = eventTable()
evf$add.sampling(seq(0,20,by=0.01))
evf$add.dosing(dose = 4*77*30, rate = 4*77, dosing.to = "ET1_total_cent_amt")

xf = NULL

#Simulate 30 min no drug
xf1 = data.frame(mod$run(theta, evf1, inits))
inits1 = as.list(xf1[dim(xf1)[1], names(xf1) %in% names(inits)])

#Simulate ET1 infusion
xf = data.frame(mod$run(theta, evf, inits1))
xf$time = xf$time + 30
xf = rbind(xf1, xf)
xf$TREAT = "Placebo"

#Simulate BQ124 - ETA antagonist arm
theta$ETA_inhibition = 90
theta$ETB_inhibition = 0
theta$ETA_ramp = 1
theta$ETB_ramp = 1
theta$Kia = 0.73 #https://pubmed.ncbi.nlm.nih.gov/7768260/, https://bpspubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1476-5381.1996.tb15212.x
theta$Kib = 24300 #theta$Kia*4e4
xf1 = data.frame(mod$run(theta, evf1, inits))
inits1 = as.list(xf1[dim(xf1)[1], names(xf1) %in% names(inits)])

theta$ETA_ramp = -1
xfA = data.frame(mod$run(theta, evf, inits1))
xfA$time = xfA$time + 30
xfA = rbind(xf1, xfA)
xfA$TREAT = "BQ123"

#Simulate BQ788 - ETB antagonist arm
theta$ETA_inhibition = 0
theta$ETB_inhibition = 90
theta$Kib = 9.8   #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1915765/
theta$Kia = 1000#theta$Kib*1e3
theta$ETA_ramp = 1
theta$ETB_ramp = 1
xf1 = data.frame(mod$run(theta, evf2, inits))
inits1 = as.list(xf1[dim(xf1)[1], names(xf1) %in% names(inits)])

#Simulate BQ788 - ETB antagonist arm
theta$ETA_ramp = 0
theta$ETB_ramp = 0
xf2 = data.frame(mod$run(theta, evf2, inits1))
inits1 = as.list(xf2[dim(xf2)[1], names(xf2) %in% names(inits)])
xf2$time = xf2$time + 15
xf1 = rbind(xf1, xf2)



theta$ETA_ramp = 0
theta$ETB_ramp = -1
xfB = data.frame(mod$run(theta, evf, inits1))
xfB$time = xfB$time + 30
xfB = rbind(xf1, xfB)
xfB$TREAT = "BQ788"





#Simulate BQ788 - ETB antagonist arm - slow
theta$C_Ib_peri = 0.05
theta$ETA_inhibition = 0
theta$ETB_inhibition = 90
theta$Kib = 9.8   #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1915765/
theta$Kia = 1000#theta$Kib*1e3
theta$ETA_ramp = 1
theta$ETB_ramp = 1
xf1 = data.frame(mod$run(theta, evf2, inits))
inits1 = as.list(xf1[dim(xf1)[1], names(xf1) %in% names(inits)])

#Simulate BQ788 - ETB antagonist arm
theta$ETA_ramp = 0
theta$ETB_ramp = 0
xf2 = data.frame(mod$run(theta, evf2, inits1))
inits1 = as.list(xf2[dim(xf2)[1], names(xf2) %in% names(inits)])
xf2$time = xf2$time + 15
xf1 = rbind(xf1, xf2)



theta$ETA_ramp = 0
theta$ETB_ramp = -1
xfBs = data.frame(mod$run(theta, evf, inits1))
xfBs$time = xfBs$time + 30
xfBs = rbind(xf1, xfBs)
xfBs$TREAT = "BQ788 - slow onset"


xf=rbind(xf, xfA, xfB, xfBs)

xf$TREAT = factor(xf$TREAT, levels = c("Placebo", "BQ123","BQ788", "BQ788 - slow onset"))
bohm$Group = factor(bohm$Group, levels = c("Placebo", "BQ123","BQ788"))
Valcol<-rep(c("#FFAB00", "#830051",  "#edbbdb",  "#3F4445", "#C4D600","#3C1053"),20)#"#003865",



h4 = ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_cent - xf$ET1_cent[1], color = TREAT), size= 1) +
  geom_point(data = bohm, mapping = aes(x=Time, y = Value*transformation - bohm$Value[1], color = Group)) + 
  geom_errorbar(data = bohm, mapping = aes(x=Time, ymin = Value - bohm$Value[1] - SD, ymax = Value- bohm$Value[1] + SD, color = Group), width = 5) + 
  ylab("Change in Plasma ET-1 (pmol/L)") + 
  xlab("Time (min)") + 
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(4,1,2,3)]))+
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8), legend.background = element_rect(fill='transparent')) +
  ylim(c(-15, 40)) + 
  annotate("rect", xmin = 30, xmax = 50, ymin = -14, ymax = -2,
           alpha = .3)+
  annotate("text",x=40, y = -8,label= "ET-1 (4 pmol/kg/min)" ) +
  annotate("rect", xmin = 0, xmax = 30, ymin = -6, ymax = -2,
           alpha = .15, fill = Valcol[4])+
  annotate("text",x=15, y = -4,label= "Placebo: Saline" ) + 
  annotate("rect", xmin = 0, xmax = 15, ymin = -14, ymax = -10,
           alpha = .2, fill = Valcol[2])+
  annotate("text",x=15, y = -8,label= "ETAi: BQ123" )+    
  annotate("rect", xmin = 0, xmax = 30, ymin = -10, ymax = -6,
           alpha = .2, fill = Valcol[1])+
  annotate("text",x=7.5, y = -12,label= "ETBi: BQ788" )# + 
#geom_path(data = xfB, mapping = aes(x=time, y = ET1_cent - xf$ET1_cent[1], color = "BQ788 delayed"))

h4


xf$drugeffect = xf$ETA_inhibitor_effect
xf$drugeffect[xf$TREAT == "BQ788"] = xf$ETB_inhibitor_effect[xf$TREAT == "BQ788"]
xf$drugeffect[xf$TREAT == "Placebo"] = NA

ggplot(xf[xf$TREAT != "Placebo",]) + geom_path(aes(x=time, y = drugeffect, color =  TREAT), size = 1) +
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(1,2)])) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8), legend.background = element_rect(fill='transparent')) +
  xlab("Time (min)")  


h5 = ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_ETA_peri, color = TREAT), size= 1) +
  ggtitle(expression("[ET1-ET"[A]*"]"[tissue]))+
  ylab("(pmol/L)") +
  xlab("Time (min)") +
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(4,1,2,3)])) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8),
                        legend.background = element_rect(fill='transparent'),
                        plot.title = element_text(hjust=0.5))


h6 = ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_ETB_peri, color = TREAT), size= 1) +
  ggtitle(expression("[ET1-ET"[B]*"]"[tissue]))+
  ylab("(pmol/L)") +
  xlab("Time (min)") +
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(4,1,2,3)])) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8),
                        legend.background = element_rect(fill='transparent'),
                        plot.title = element_text(hjust=0.5))

h7 = ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_ETA_cent, color = TREAT), size= 1) +
  ggtitle(expression("[ET1-ET"[A]*"]"[plasma]))+
  ylab("(pmol/L)") +
  xlab("Time (min)") +
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(4,1,2,3)])) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8),
                        legend.background = element_rect(fill='transparent'),
                        plot.title = element_text(hjust=0.5))


h8 = ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_ETB_cent, color = TREAT), size= 1) +
  ggtitle(expression("[ET1-ET"[B]*"]"[plasma]))+
  ylab("(pmol/L)") +
  xlab("Time (min)") +
  scale_x_continuous( breaks = c(0, 30,60)) +
  scale_color_manual(values = c(Valcol[c(4,1,2,3)])) +
  theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8),
                        legend.background = element_rect(fill='transparent'),
                        plot.title = element_text(hjust=0.5))

h9 = ggarrange(h7,h8,h5,h6, common.legend = TRUE, nrow = 2, ncol = 2)

ggarrange(h4, h9, nrow  = 2, labels = c("A", "B"),common.legend = FALSE) # %>% ggexport(filename = "bohmfig.pdf", height = 9, width = 6)
# ggplot() + geom_path(data = xf, mapping = aes(x=time, y = ET1_ETA_internalization_cent, color = TREAT), size= 1) +
#   ylab("Change in Plasma ET-1 (pmol/L)") + 
#   xlab("Time (min)") + 
#   scale_x_continuous( breaks = c(0, 30,60)) +
#   theme_light() + theme(legend.title = element_blank(), legend.position = c(0.25,0.8), legend.background = element_rect(fill='transparent')) 
# 


# ggplot() + geom_path(data = xf, aes(x=time, y = ETB_inhibitor_effect, color = TREAT))
# ggplot() + geom_path(data = xf, aes(x=time, y = ETA_inhibitor_effect, color = TREAT))
#See plot from Peter - drugs with different affinities for ETA vs ETB. Does it matter?
#How much ETB needed? etc.

#Link Atrasentan PK?


#ggplot(xf[xf$TREAT == "BQ788",]) + ggplot() + geom_path(data = xf, aes(x=time, y = I_Kib))






