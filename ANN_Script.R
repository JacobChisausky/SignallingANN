#The point of this code is to read ANN stats from the cpp program and allow us to visualize the phenotypes
#set.seed(1234567) 

library(ggplot2)
library(plotly)
library(reshape2)

rlu <- function(x) {
  if (x < 0) {
    return(0)
  } else if (x > 1) {
    return(1)
  } else {
    return(x)
  }
} #Return 0 if < 0 and 1 if > 1

#Sender Functions ####
senderANN_output <- function(s,q,annS){
  n1 <- rlu(annS[1+1] + s) # the '+1's are to account for the 0 index in cpp
  n2 <- rlu(annS[2+1] + q)
  
  n3 <- rlu(annS[3+1] + n1 * annS[11+1] + n2 * annS[15+1])
  n4 <- rlu(annS[4+1] + n1 * annS[12+1] + n2 * annS[16+1])
  n5 <- rlu(annS[5+1] + n1 * annS[13+1] + n2 * annS[17+1])
  n6 <- rlu(annS[6+1] + n1 * annS[14+1] + n2 * annS[18+1])
  
  n7 <- rlu(annS[7+1] + n3 * annS[19+1] + n4 * annS[23+1] + n5 * annS[27+1] + n6 * annS[31+1])
  n8 <- rlu(annS[8+1] + n3 * annS[20+1] + n4 * annS[24+1] + n5 * annS[28+1] + n6 * annS[32+1])
  n9 <- rlu(annS[9+1] + n3 * annS[21+1] + n4 * annS[25+1] + n5 * annS[29+1] + n6 * annS[33+1])
  n10 <- rlu(annS[10+1] + n3 * annS[22+1] + n4 * annS[26+1] + n5 * annS[30+1] + n6 * annS[34+1])
  
  output <- rlu(annS[0+1] + n7 * annS[35+1] + n8 * annS[36+1] + n9 * annS[37+1] + n10 *annS[38+1])
  return(output)
} #Gives output of a sender ANN
senderPhenotype <- function(s_num, q_num, ANN_S){
  #s_num = number of strength values to plot. More = finer plot
  #q_num = number of quality values to calculate
  
  phenotype_S <- data.frame()
  
  for (q_cur in 0:q_num){
    q <- q_cur/q_num
    for (s_cur in 0:s_num){
      s <- s_cur/s_num
      out <- senderANN_output(s,q,ANN_S)
      row <- c(q,s,out)
      phenotype_S <- rbind(phenotype_S,row)
    }
  }
  colnames(phenotype_S) <- c("q","s","output")
  return(phenotype_S)
} #Returns a data frame for s_num s values and q_num q. For plotting
generate_ANNs_random <- function(lim){
  return(runif(39, min=-abs(lim), max=abs(lim)))
} #Returns sender ANN with random vars from -lim to lim
generate_ANNs_randomComplex <- function(lim){
  a<-generate_ANNs_random(lim)
  phenotype <- senderPhenotype(10,10,a)
  while (sum(unique(phenotype$output)) < 10){
    a<-generate_ANNs_random(1)
    phenotype <- senderPhenotype(10,10,a)
    #I could do this in cpp at initialization to make more interesting starting phenotypes?
  }
  return(a)
} #Returns sender which is complex, i.e, not flat
generate_ANNs_null <-function(){
  return(c(rep(0,39)))
} #Returns sender ANN with all 0s
senderANN_mutate <- function(gen, mut_rate, mut_step, ANN_S){
  ann<-ANN_S
  for (g in i:gen){
    for (i in 1:39){ #For each variable
      if (runif(1,0,1) < mut_rate){ #check if mut happens
        mutSize <- rcauchy(1, 0, mut_step)  #generate mut size
        ann[i] <- ann[i] + mutSize   #change variable
      }
    }
  }
  return(ann)
}#mutates an ANN gen number of times
senderANN_printSurface <- function(resolution, ANN_S){
  phenotype <- senderPhenotype(resolution,resolution,ANN_S)
  data <- dcast(phenotype, q ~ s, value.var = "output")
  rownames(data) <- data[,1]
  data<-data[,-1]
  colnames(data) <- as.numeric(colnames(data))
  rownames(data) <- as.numeric(rownames(data))
  slabs<-colnames(data)
  qlabs<-rownames(data)
  matrix <- as.matrix(data)
  
  p <- plot_ly(z = ~matrix, type = "surface")
  
  
  p <- layout(p, scene = list(xaxis = list(title = "s", 
                                           ticketmode = 'array',
                                           ticktext = round(as.numeric(slabs),2),
                                           tickvals = (as.numeric(slabs)*(ncol(matrix)-1)),
                                           tickformat = '.2f',
                                           autorange="reversed"),
                              yaxis = list(title = "q", 
                                           ticketmode = 'array',
                                           ticktext = round(as.numeric(qlabs),2),
                                           tickvals = (as.numeric(qlabs)*(nrow(matrix)-1)),
                                           tickformat = '.2f'),
                              zaxis = list(title = "Pr")))
  
  print(p)
} 

#Receiver Functions ####
receiverANN_output <- function(s,annR){
  n1 = rlu(annR[1+1] + s);
  
  n2 = rlu(annR[2+1] + n1 * annR[10+1]);
  n3 = rlu(annR[3+1] + n1 * annR[11+1]);
  n4 = rlu(annR[4+1] + n1 * annR[12+1]);
  n5 = rlu(annR[5+1] + n1 * annR[13+1]);
  
  n6 = rlu(annR[6+1] + n2 * annR[14+1] + n3 * annR[18+1] + n4 * annR[22+1] + n5 * annR[26+1]);
  n7 = rlu(annR[7+1] + n2 * annR[15+1] + n3 + annR[19+1] + n4 * annR[23+1] + n5 * annR[27+1]);
  n8 = rlu(annR[8+1] + n2 * annR[16+1] + n3 * annR[20+1] + n4 * annR[24+1] + n5 * annR[28+1]);
  n9 = rlu(annR[9+1] + n2 * annR[17+1] + n3 * annR[21+1] + n4 * annR[25+1] + n5 * annR[29+1]);
  
  output = rlu(annR[0+1] + n6 * annR[30+1] + n7 * annR[31+1] + n8 * annR[32+1] + n9 * annR[33+1]);
  return(output)
} #Gives output of a receiver ANN
receiverPhenotype <- function(s_num, ANN_R){
  #s_num = number of strength values to plot. More = finer plot
  #q_num = number of quality values to calculate
  
  phenotype_R <- data.frame()
  
    for (s_cur in 0:s_num){
      s <- s_cur/s_num
      out <- receiverANN_output(s,ANN_R)
      row <- c(s,out)
      phenotype_R <- rbind(phenotype_R,row)
    }
  colnames(phenotype_R) <- c("s","output")
  return(phenotype_R)
} #Returns a data frame for s_num s values. For plotting
generate_ANNr_random <- function(lim){
  return(runif(34, min=-abs(lim), max=abs(lim)))
}
generate_ANNR_randomComplex <- function(lim){
  a<-generate_ANNr_random(lim)
  phenotype <- receiverPhenotype(20,a)
  while (sum(unique(phenotype$output)) < 2){
    a<-generate_ANNr_random(1)
    phenotype <- receiverPhenotype(20,a)
    #I could do this in cpp at initialization to make more interesting starting phenotypes?
  }
  return(a)
} #Returns receiver which is complex, i.e, not flat
generate_ANNr_null <-function(){
  return(c(rep(0,34)))
} #Returns receiver ANN with all 0s
receiverANN_mutate <- function(gen, mut_rate, mut_step, ANN_R){
  ann<-ANN_R
  for (g in i:gen){
    for (i in 1:34){ #For each variable
      if (runif(1,0,1) < mut_rate){ #check if mut happens
        mutSize <- rcauchy(1, 0, mut_step)  #generate mut size
        ann[i] <- ann[i] + mutSize   #change variable
      }
    }
  }
  return(ann)
}#mutates an ANN gen number of times


#Senders tests ####

#testing the effect of mutations on random ANNs
a<-generate_ANNs_randomComplex(1)
senderANN_printSurface(10,a)

a<-senderANN_mutate(10,0.01,0.05,a)
senderANN_printSurface(10,a)

# Receivers tests ####

# Investigating mutation sizes on receiver phenotypes
for (j in 1:10){
a<-generate_ANNR_randomComplex(1)
d<-receiverPhenotype(100,a)
d$gen<-0
d_mut <- d
gen <- 0
for (i in c(1,10,100,200,300,500,1000,2000,3000,5000)){
  a<-receiverANN_mutate(i-gen,0.01,0.05,a)  
  d<-receiverPhenotype(100,a)
  gen<-i  
  d$gen <- gen
  d_mut<-rbind(d_mut,d)
}

p<-ggplot(d_mut,aes(x=s,y=output)) + 
  geom_line(aes(color=log10(gen),group=gen),linewidth=1.5,alpha=.6) +
  ylim(0,1) + xlim(0,1) +
  theme_bw() + scale_color_viridis_c(option="viridis")
print(p)
}


##Read data ####
directory <- "C:/Users/owner/Documents/S4/Simulation_ANN"
file <- "117_14_05_29_annVars_forR.csv"  #null receivers - good
file <- "117_18_33_44_annVars_null_Senders.CSV" #null senders - good
#file <- "120_11_06_06_annVars_first_real_test.csv" #Pooling equilibrium - all signal, all respond
file <- "120_17_55_35_annVars_first_real_test.csv"  #They are all sending signal 0 since it is tried first
    #So - all signals are s = 0. This causes no decrease in survival.
    #Why do receivers respond? Not responding should be the same as responding: payoff of q or 1-q will be the same
    #Same problem as in other sim: If senders don't send, no inventive to respond
        #If responders don't respond, no incentive to send.
file <- "121_11_00_26_annVars_test_honest_start_smallN.csv" #Honest at 100000. By 120000 when coevolution start, a few receivers respond to very strong signals. We reach pooling equilibrium soon.
file <- "121_13_38_28_annVars_test_honest_start_smallN.csv" #By g=220000 when honest until 200000: receivers A1 always.
    #How can I increase the incentive to A1 only high q and A2 only low q?
    #Implement k coefficient to receiver payoffs
file <- "honestANNs/123_19_26_08_ANNvars_get_honest_ANNs_9.csv"

annFile <- read.csv(paste0(directory,"/",file))
params <- read.csv(paste0(directory,"/",pFile))
sFile <- subset(annFile,indType == "Sender")
rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
sEnd <- subset(sFile,gen==max(sFile$gen))
rEnd <- subset(rFile,gen==max(rFile$gen))
##

unique(sFile$gen)
unique(annFile$rep)
# Senders ####

n<-4
n_annS<-as.numeric(sEnd[n,7:45])
senderANN_printSurface(30,n_annS)


#Null Receiver analytical predictuon
#prediction for null receivers
data<-senderPhenotype(100,100,n_annS)
pred<-function(x){
  return( (1/(2+10*x))^(1/(1+10*x)) )
}
ggplot(data) +
  geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2) + 
  scale_color_viridis_c() +
  theme_bw() +
  geom_function(fun = pred, colour = "red", linewidth=3)


#Real data - receivers
unique(rFile$gen)
rEnd <- subset(rFile,gen==220000)

n<-1
n_annR<-as.numeric(rEnd[n,-(1:6)])
data<-receiverPhenotype(100,n_annR)

ggplot(data,aes(x=s,y=output)) + geom_line(color="red",size=1.4) +
  theme_bw() + ylim(0,1) + labs(y="Pr")

dataAll<-data.frame()
for (n in 1:100){
n_annR<-as.numeric(rEnd[n,-(1:6)])
data<-receiverPhenotype(100,n_annR)
data$n<-n
dataAll<-rbind(dataAll,data)
}
ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
  theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals")








#### TestParams ####
directory <- "C:/Users/owner/Documents/S4/Simulation_ANN"
#f 1 had more honest receivers at 100 000 and some at 125 000
#f 2 had lots more variability at 100 000
#f 3 looked okay at 100 000
#f 4 was really noisy 

fNum <- 1
file <- paste0("testParams/123_00_05_42_annVars_testParams_",fNum,".csv")
pFile <- paste0("testParams/123_00_05_42_params_testParams_",fNum,".csv")
params <- read.csv(paste0(directory,"/",pFile))

params1 <- read.csv(paste0(directory,"/",paste0("testParams/123_00_05_42_params_testParams_",1,".csv")))
params2 <- read.csv(paste0(directory,"/",paste0("testParams/123_00_05_42_params_testParams_",2,".csv")))
params3 <- read.csv(paste0(directory,"/",paste0("testParams/123_00_05_42_params_testParams_",3,".csv")))
params4 <- read.csv(paste0(directory,"/",paste0("testParams/123_00_05_42_params_testParams_",4,".csv")))

#Clean
params1
params3
#Messy
params2
params4
df <- rbind(params1,params3,params2,params4)
df

unique(annFile$gen)
annFile <- read.csv(paste0(directory,"/",file))
sFile <- subset(annFile,indType == "Sender")
rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
sEnd <- subset(sFile,gen==max(sFile$gen))
rEnd <- subset(rFile,gen==max(rFile$gen))
##

unique(sFile$gen)
unique(annFile$rep)
# Senders ####
params


n<-1
n_annS<-as.numeric(sEnd[n,7:45])
senderANN_printSurface(30,n_annS)


#Null Receiver analytical predictuon
#prediction for null receivers
data<-senderPhenotype(100,100,n_annS)
pred<-function(x){
  return( (1/(2+10*x))^(1/(1+10*x)) )
}
ggplot(data) +
  geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2) + 
  scale_color_viridis_c() +
  theme_bw() +
  geom_function(fun = pred, colour = "red", linewidth=3)


#Real data - receivers
unique(rFile$gen)
rEnd <- subset(rFile,gen==100000)

n<-1
n_annR<-as.numeric(rEnd[n,-(1:6)])
data<-receiverPhenotype(100,n_annR)

ggplot(data,aes(x=s,y=output)) + geom_line(color="red",size=1.4) +
  theme_bw() + ylim(0,1) + labs(y="Pr")

dataAll<-data.frame()
for (n in 1:100){
  n_annR<-as.numeric(rEnd[n,-(1:6)])
  data<-receiverPhenotype(100,n_annR)
  data$n<-n
  dataAll<-rbind(dataAll,data)
}
ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
  theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals")

