# Analysis of testing data
# Kucharski et al (2024)
# https://github.com/adamkucharski/pl-testing

# Load libraries -------------------------------------------------------------------

library(dplyr)
library(readxl)
library(incidence2)
library(lubridate)
library(readr)
library(epitools)

# Load data ---------------------------------------------------------------

# Read in incidence data
incidence_data <- read_csv("data_out/incidence_data.csv")
incidence_data_lft <- read_csv("data_out/incidence_data_lft.csv")

# Load Ct data
ct_data <- read_csv("data_out/ct_values.csv")

# Load variant-specific Ct distributions
first_test_wt <- read_csv("data_out/first_test_wt.csv")
first_test_alpha <- read_csv("data_out/first_test_alpha.csv")

# Load PCR vs LFT comparison data
proportion_data_pcr <- read_csv("data_out/proportion_data_pcr.csv")
proportion_data_lft <- read_csv("data_out/proportion_data_lft.csv")

# ONS community infection survey estimated prevalence by age
# Source: Abbott & Funk, 2022: https://www.medrxiv.org/content/10.1101/2022.03.29.22273101v1
data_cis_est <- read_csv("https://raw.githubusercontent.com/epiforecasts/inc2prev/main/outputs/estimates_age.csv")

data_cis_est_16 <- data_cis_est |> filter(name=="est_prev" & variable=="16-24")
data_cis_est_25 <- data_cis_est |> filter(name=="est_prev" & variable=="25-34")


# Plot prevalence data -------------------------------------------------------------------

# Helper function: plot data and binom CI
plot_CI <- function(dates,xx,nn,colA="black",cex_in=1) {
  
  for(ii in 1:length(nn)){
    test_r <- binom.test(xx[ii],nn[ii])
    CI1 <- as.numeric(test_r$conf.int)[1]
    CI2 <- as.numeric(test_r$conf.int)[2]
    points(dates[ii],xx[ii]/nn[ii],col=colA,pch=19,cex=cex_in); lines(c(dates[ii],dates[ii]),c(CI1,CI2),col=colA)
  }
  
}

# Plot positivity across 3 timeseries
par(mfcol=c(2,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
xrange <- as.Date(c("2020-05-01","2022-01-01")) #"2022-08-15"))
yrange <- c(0,0.1)

# PL data
plot((incidence_data$date_index),0*incidence_data$pos_ratio,main="SARS-CoV-2 prevalence",
     col="white",xlab="",ylab="proportion positive",type="l",xlim=xrange,ylim=yrange,yaxt="n",yaxs="i")
grid(ny=NULL,nx=NA,col="lightgray")
title(LETTERS[1],adj=0)

# Select LFT period
pick_lft <- (incidence_data$pos_ratio<0.05 | incidence_data$date_index>as.Date("2021-12-01"))

# Plot PCR data only
incidence_data_plot <- incidence_data[pick_lft,]

polygon(c(incidence_data[!pick_lft,]$date_index,rev(incidence_data[!pick_lft,]$date_index)),
        c(rep(0,sum(!pick_lft)),rep(1,sum(!pick_lft))),lty=0,col=rgb(0,0,0,0.15))

plot_CI((incidence_data_plot$date_index),incidence_data_plot$n_pos,incidence_data_plot$count,cex_in=0.8)


# CIS data
polygon(c(data_cis_est_25$date,rev(data_cis_est_25$date)),c(data_cis_est_25$q5,rev(data_cis_est_25$q95))/100,col=rgb(0,0,1,0.5),lty=0)
lines(data_cis_est_25$date,data_cis_est_25$q50/100,col="blue")

polygon(c(data_cis_est_16$date,rev(data_cis_est_16$date)),c(data_cis_est_16$q5,rev(data_cis_est_16$q95))/100,col=rgb(1,1,0,0.5),lty=0)
lines(data_cis_est_16$date,data_cis_est_16$q50/100,col="orange")

axis(2, at=seq(0,1,0.02),labels = sapply(seq(0,100,2),function(x){paste(x,"%",sep="")}),col = "black") 
title(ylab="Percent positive", line=3, cex.lab=1)

graphics::text(x=xrange[1]+220,y=0.065,labels="- ONS, age 25-29",col="blue")
graphics::text(x=xrange[1]+220,y=0.072,labels="- ONS, age 16-24",col="orange")
graphics::text(x=as.Date("2021-06-21")+10,y=0.092,labels="LFT testing",col=rgb(0,0,0,0.7),adj=0)


# Plot individual trajectories -------------------------------------------------------------------

# Plot on normal scale
plot(ct_data$diff_date,ct_data$CT_value,pch=19,col=rgb(0,0,0,1),cex=0.8,
     ylim=c(50,20),main="Individual Ct values",ylab="Ct value",xlab="days since first positive test")
grid(ny=NULL,nx=NA,col="lightgray")

symbols(110, 34, circles=c(21), inches=FALSE, add=TRUE, bg=NA, fg="red", lty=2, lwd=2, xlab="", ylab="")
symbols(255, 27, circles=c(25), inches=FALSE, add=TRUE, bg=NA, fg="purple", lty=2, lwd=2, xlab="", ylab="")

graphics::text(x=90,y=28,labels="One individual reinfected\n in Alpha wave (3 samples)",col="red",adj=0,cex=0.8)
graphics::text(x=275,y=34,labels="Three individuals reinfected\n in Delta wave",col="purple",adj=1,cex=0.8)

title(LETTERS[2],adj=0)

dev.copy(png,paste0("outputs/ct_dynamics.png"),units="cm",width=15,height=20,res=200)
dev.off()
  
# Get variant Ct distributions -------------------------------------------------------------------

# Calculate mean Ct values
c(mean(first_test_wt$CT_value),sd(first_test_wt$CT_value),length(first_test_wt$CT_value))
c(mean(first_test_alpha$CT_value),sd(first_test_alpha$CT_value),length(first_test_alpha$CT_value))

# Plot distributions
break_n <- seq(15,50,5)
par(mfcol=c(2,1),mar=c(4,3,1,1),mgp=c(2,0.6,0),las=1)

hist(first_test_wt$CT_value,xlab="Ct value",main="Wildtype",prob=F,breaks=break_n)
title(LETTERS[1],adj=0)

hist(first_test_alpha$CT_value,xlab="Ct value",main="Alpha",prob=F,breaks=break_n)
title(LETTERS[2],adj=0)
dev.copy(png,paste0("outputs/Ct_period.png"),units="cm",width=12,height=10,res=200)
dev.off()


# Reinfection analysis ----------------------------------------------------

# Define 2x2 matrix
matrix_1 <- matrix(c(1,263,164,3694),nrow=2)

# Calculate odds ratio
oddsratio(matrix_1)

# Compare PCR vs antigen ----------------------------------------------------------

par(mfcol=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
plot(proportion_data_pcr$Days_Since_First_Positive,proportion_data_pcr$Proportion_Positive,
     pch=19,col="dark blue",xlab="days since first positive PCR",ylab="proportion positive",xlim=c(-5,20),ylim=c(0,1),yaxs="i")
lines(proportion_data_pcr$Days_Since_First_Positive,proportion_data_pcr$Proportion_Positive,col="dark blue",lwd=0.5)
points(proportion_data_lft$Days_Since_First_Positive,proportion_data_lft$Proportion_Positive,
       pch=19,col="dark orange")
lines(proportion_data_lft$Days_Since_First_Positive,proportion_data_lft$Proportion_Positive,col="dark orange",lwd=0.5)

dev.copy(png,paste0("outputs/PCR_vs_LFT.png"),units="cm",width=12,height=10,res=200)
dev.off()





