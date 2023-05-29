#The point of this code is to read ANN stats from the cpp program and allow us to visualize the phenotypes
#set.seed(1234567) 

#Exper 1 is looking at multiple s vales
#Exper 2 is investigating what mut rates and steps promote hybrid equilibria


library(ggplot2)
library(plotly)
library(reshape2)

#General Functions ####
pred<-function(x){
  return( (1/(2+10*x))^(1/(1+10*x)) )
}

rowProduct <- function(row){
  total <- row[1]
  for (i in 2:length(row)){
    total <- total*row[i]
  }
  return(total)
}

rlu <- function(x) {
  if (x < 0) {
    return(0)
  } else if (x > 1) {
    return(1)
  } else {
    return(x)
  }
} #Return 0 if < 0 and 1 if > 1

get_levels <- function(levels){
  increment = 1.0/(levels-1)
  vals <- c()
  for (i in 0:(levels-1)){
    vals<-c(vals,increment*i)
  }
  return(vals);
}

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

senderPhenotypeDiscrete <- function(s_levels, q_levels, ANN_S){
  s_vals<-get_levels(s_levels)
  q_vals<-get_levels(q_levels)
  
  phenotype_S <- data.frame()
  
  for (q in q_vals){
    for (s in s_vals){
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
receiverANN_output <- function(s, r,annR){
  n1 <- rlu(annR[1+1] + r) # the '+1's are to account for the 0 index in cpp
  n2 <- rlu(annR[2+1] + s)
  
  n3 <- rlu(annR[3+1] + n1 * annR[11+1] + n2 * annR[15+1])
  n4 <- rlu(annR[4+1] + n1 * annR[12+1] + n2 * annR[16+1])
  n5 <- rlu(annR[5+1] + n1 * annR[13+1] + n2 * annR[17+1])
  n6 <- rlu(annR[6+1] + n1 * annR[14+1] + n2 * annR[18+1])
  
  n7 <- rlu(annR[7+1] + n3 * annR[19+1] + n4 * annR[23+1] + n5 * annR[27+1] + n6 * annR[31+1])
  n8 <- rlu(annR[8+1] + n3 * annR[20+1] + n4 * annR[24+1] + n5 * annR[28+1] + n6 * annR[32+1])
  n9 <- rlu(annR[9+1] + n3 * annR[21+1] + n4 * annR[25+1] + n5 * annR[29+1] + n6 * annR[33+1])
  n10 <- rlu(annR[10+1] + n3 * annR[22+1] + n4 * annR[26+1] + n5 * annR[30+1] + n6 * annR[34+1])
  
  output <- rlu(annR[0+1] + n7 * annR[35+1] + n8 * annR[36+1] + n9 * annR[37+1] + n10 *annR[38+1])
  return(output)
} #Gives output of a sender ANN

receiverPhenotype <- function(s_num, r_num, ANN_R){
  #s_num = number of strength values to plot. More = finer plot
  #q_num = number of quality values to calculate
  
  phenotype_R <- data.frame()
  
  for (r_cur in 0:r_num){
    r <- r_cur/r_num
    for (s_cur in 0:s_num){
      s <- s_cur/s_num
      out <- receiverANN_output(s,r,ANN_R)
      row <- c(r,s,out)
      phenotype_R <- rbind(phenotype_R,row)
    }
  }
  colnames(phenotype_R) <- c("r","s","output")
  return(phenotype_R)
} #Returns a data frame for s_num s values and q_num q. For plotting

receiverPhenotypeDiscrete <- function(s_levels, r_levels, ANN_R){
  s_vals<-get_levels(s_levels)
  r_vals<-get_levels(r_levels)
  
  phenotype_R <- data.frame()
  
  for (s in s_vals){
    for (r in r_vals){
      out <- receiverANN_output(s,r,ANN_R)
      row <- c(s,r,out)
      phenotype_R <- rbind(phenotype_R,row)
    }
  }
  colnames(phenotype_R) <- c("s","r","output")
  return(phenotype_R)
}

generate_ANNr_random <- function(lim){
  return(runif(39, min=-abs(lim), max=abs(lim)))
} #Returns ANN with random vars from -lim to lim
generate_ANNr_randomComplex <- function(lim){
  a<-generate_ANNr_random(lim)
  phenotype <- receiverPhenotype(10,10,a)
  while (sum(unique(phenotype$output)) < 10){
    a<-generate_ANNr_random(1)
    phenotype <- receiverPhenotype(10,10,a)
    #I could do this in cpp at initialization to make more interesting starting phenotypes?
  }
  return(a)
} #Returns ann which is complex, i.e, not flat
generate_ANNr_null <-function(){
  return(c(rep(0,39)))
} #Returns ANN with all 0s
receiverANN_mutate <- function(gen, mut_rate, mut_step, ANN_R){
  ann<-ANN_R
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
receiverANN_printSurface <- function(resolution, ANN_R){
  phenotype <- receiverPhenotype(resolution,resolution,ANN_R)
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

#Functions for Prob Distributions ####

calcProbs<-function(phenotype,tries_max){
  receiver<-FALSE
  if (colnames(phenotype)[1]=="s"){
    receiver<-TRUE
    colnames(phenotype) <- c("q","s","output")
  }
  
  dataAll<-data.frame()
  for (q_cur in unique(phenotype$q)){
    data<-subset(phenotype,q==q_cur)
    data$probRaw <- data$output/sum(data$output)
    if (sum(data$output)==0){data$probRaw<-0}
    
    #Chance of selecting each by normal method (in first ten rounds)
    data$probSelectionNormal <- (1-(1-data$output)^tries_max)
    
    #chance neither selected by normal method = rowProduct(1-data$probSelectionNormal)
    #ADD this to s=0
    
    #chance one is selected by normal method = 1-rowProduct(1-data$probSelectionNormal)
    
    #chance each one is selected by normal method
    #Calculated as chance either will be selected * chance each is selected vs the other per round
    data$probNormalEach <- data$probRaw*(1-rowProduct(1-data$probSelectionNormal))
    
    #prob final is corrected for no event occuring
    data$probFinal <- ifelse(data$s==0,rowProduct(1-data$probSelectionNormal) + data$probNormalEach,data$probNormalEach)
    dataAll<-rbind(dataAll,data)
  }
  
  output<-dataAll[,-c(4:6)]
  if (receiver==TRUE){
    colnames(output) <- c("s","r","output","probFinal")
  }
  return(output)
}




#Test 1: Initial ####
directory_test1 <- "C:/Users/owner/eclipse-workspace/SignallingANN/data"
directory_test1 <- "C:/Users/owner/Documents/DONE/hybrid_mutRates"
annFiles <- list.files(directory_test1,"*annVars*")
paramFiles <- list.files(directory_test1,"*params_t*")
annFiles
annFile <- read.csv(paste0(directory_test1,"/",annFiles[2]))

paramsAll<-data.frame()
for (i in paramFiles){
  paramFile <- read.csv(paste0(directory_test1,"/",i))
  paramsAll<-rbind(paramsAll,paramFile)
}
paramsAll

unique(annFile$gen)


#for (G in unique(annFile$gen)){
#for (num in 1:length(annFiles)){
annFiles
num<-1
annFile <- read.csv(paste0(directory_test1,"/",annFiles[num]))
paramFile <- read.csv(paste0(directory_test1,"/",paramFiles[num]))
s_vals<-get_levels(paramFile$s_levels)
q_vals<-get_levels(paramFile$q_levels)
r_vals<-get_levels(paramFile$r_levels)

sFile <- subset(annFile,indType == "Sender")

unique(annFile$gen)
G<-max(annFile$gen)

sGen <- subset(sFile,gen==G)
rGen <- subset(rFile,gen==G)

#Put individuals together
dataMult<-data.frame()
for (n in 1:50){
  n_annS<-as.numeric(sGen[n,7:45])
  data<-senderPhenotypeDiscrete(2,2,n_annS)
  data$n <- n  
  dataMult<-rbind(dataMult,data)
}
p <- ggplot(dataMult) +
  geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=4) + 
  scale_color_viridis_c() +
  theme_bw() +
  facet_wrap(~n) +
  labs(title=paste0("Senders: ",num)) +
  labs(subtitle=paste0("gen = ",G,
                       "\nk = ",paramFile$k,
                       "\ncMin = ",paramFile$cMin,
                       "\ncMax = ",paramFile$cMax
  ))
p

paramFile


ggsave(plot=p,paste0(num,"_",G/10000,"_S.png"),
       device="png",path=directory_test1,height=8,width=10,unit="in")


if (1==2){
  n<-1
  n_annS<-as.numeric(sGen[n,7:45])
  
  #senderANN_printSurface(30,n_annS)
  data<-senderPhenotype(70,70,n_annS)
  ggplot(data) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw()
}


#Receivers
paramFile$cMax
paramFile$cMin

unique(annFile$gen)
G<-2500

rFile <- subset(annFile,indType == "Receiver")
rGen <- subset(rFile,gen==G)
rEnd <- subset(rFile,gen==max(rFile$gen))

dataMult<-data.frame()
for (n in 1:40){
  n_annR<-as.numeric(rGen[n,7:45])
  data<-receiverPhenotypeDiscrete(2,2,n_annR)
  data$n <- n  
  dataMult<-rbind(dataMult,data)
}
p <- ggplot(dataMult) +
  geom_point(aes(x=s,y=r,color=output),alpha=0.8,size=4) + 
  scale_color_viridis_c() +
  theme_bw() +
  facet_wrap(~n) +
  labs(title=paste0("Receivers: ",num)) +
  labs(subtitle=paste0("gen = ",G,
                       "\nk = ",paramFile$k,
                       "\ncMin = ",paramFile$cMin,
                       "\ncMax = ",paramFile$cMax
  ))
p

ggsave(plot=p,paste0(num,"_",G/10000,"_R.png"),
       device="png",path=directory_test1,height=8,width=10,unit="in")


#### Discrete phenotypes - probabilities ####
unique(sFile$gen)
G<-20000
sGen <- subset(sFile,gen==G)
rGen <- subset(rFile,gen==G)

probsAll<-data.frame()
for (n in 1:nrow(sGen)){
  n_annS<-as.numeric(sGen[n,7:45])
  d<-senderPhenotypeDiscrete(2,2,n_annS) 
  probs<-calcProbs(d,10)
  colnames(probs)<-c("state","strength","output","probFinal")
  probs$n <- n
  probs$indType = "Sender"
  probsAll<-rbind(probsAll,probs)
}
for (n in 1:nrow(rGen)){
  n_annR<-as.numeric(rGen[n,7:45])
  d<-receiverPhenotypeDiscrete(2,2,n_annR) 
  probs<-calcProbs(d,10)
  colnames(probs)<-c("state","strength","output","probFinal")
  probs$n <- n
  probs$indType = "Receiver"
  probs
  probsAll<-rbind(probsAll,probs)
}

# Plot for all levels = 2
ggplot(subset(probsAll,strength==1)) + 
  geom_jitter(aes(x=state,y=probFinal),height=0,width=0.05,alpha=.6) +
  facet_grid(indType~.) +
  theme_bw() +
  stat_summary(
    fun = "mean",
    geom = "point",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red",
    aes(x=state,y=probFinal) ) + 
  labs(subtitle="Received signal strength (receivers)\nQuality (senders)") +
  labs(y="Probability") +
  labs(x="Quality (senders) or \nReceived Signal Strength (receivers)") +
  ylim(0,1) +
  #Expected beta
  geom_point(aes(y=ifelse(indType=="Receiver",ifelse(paramFile$cMax<1,paramFile$cMax,1),-5)),x=1,size=5,color="blue") +
  #Expected alpha
  geom_point(aes(y=ifelse(indType=="Sender",ifelse(paramFile$cMax<1,(paramFile$m/(1-paramFile$m)),0),-5)),x=0,size=5,color="blue")

#### Discrete phenotypes - collate files ####
#for (G in unique(annFile$gen)){
unique(annFile$gen)
G<-20000
probsMaster<-data.frame()
for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directory_test1,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directory_test1,"/",paramFiles[num]))
  s_vals<-get_levels(paramFile$s_levels)
  q_vals<-get_levels(paramFile$q_levels)
  r_vals<-get_levels(paramFile$r_levels)
  sFile <- subset(annFile,indType == "Sender")
  rFile <- subset(annFile,indType == "Receiver")
  sGen <- subset(sFile,gen==G)
  rGen <- subset(rFile,gen==G)
  
  probsAll<-data.frame()
  for (n in 1:nrow(sGen)){
    n_annS<-as.numeric(sGen[n,7:45])
    d<-senderPhenotypeDiscrete(2,2,n_annS) 
    probs<-calcProbs(d,10)
    colnames(probs)<-c("state","strength","output","probFinal")
    probs$n <- n
    probs$indType = "Sender"
    probsAll<-rbind(probsAll,probs)
  }
  for (n in 1:nrow(rGen)){
    n_annR<-as.numeric(rGen[n,7:45])
    d<-receiverPhenotypeDiscrete(2,2,n_annR) 
    probs<-calcProbs(d,10)
    colnames(probs)<-c("state","strength","output","probFinal")
    probs$n <- n
    probs$indType = "Receiver"
    probs
    probsAll<-rbind(probsAll,probs)
  }
  probsAll$num<-num
  probsAll$mut_rate_ann <- paramFile$mut_rate_ann_S
  probsAll$mut_step_ann <- paramFile$mut_step_ann_S
  
  probsMaster<-rbind(probsMaster,probsAll)
}

paramFile

ggplot(subset(probsMaster,strength==1)) + 
  geom_jitter(aes(x=state,y=probFinal),height=0,width=0.15,alpha=.4) +
  theme_bw() +
  stat_summary(
    fun = "mean",
    geom = "point",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red",
    aes(x=state,y=probFinal) ) + 
  labs(subtitle="Received signal strength (receivers)\nQuality (senders)") +
  labs(y="Probability") +
  labs(x="Quality (senders) or \nReceived Signal Strength (receivers)") +
  ylim(0,1) +
  #Expected beta
  geom_point(aes(y=ifelse(indType=="Receiver",ifelse(paramFile$cMax<1,paramFile$cMax,1),-5)),x=1,size=5,color="blue") +
  #Expected alpha
  geom_point(aes(y=ifelse(indType=="Sender",ifelse(paramFile$cMax<1,(paramFile$m/(1-paramFile$m)),0),-5)),x=0,size=5,color="blue") +
  facet_grid(mut_rate_ann+indType~mut_step_ann)

#Exper 1: More s levels ####
#May have too low of mut rates or steps for this to work...
directory_exper1 <- "C:/Users/owner/Documents/DONE/exper1b_many_s"
annFiles <- list.files(directory_exper1,"*annVars*")
paramFiles <- list.files(directory_exper1,pattern = ".*_params_.*\\.csv$")

annFile <- read.csv(paste0(directory_exper1,"/",annFiles[2]))
paramsAll<-data.frame()
for (i in paramFiles){
  paramFile <- read.csv(paste0(directory_exper1,"/",i))
  paramsAll<-rbind(paramsAll,paramFile)
}
paramsAll

unique(annFile$gen)


#for (G in unique(annFile$gen)){
#for (num in 1:length(annFiles)){
annFiles
num<-1
annFile <- read.csv(paste0(directory_exper1,"/",annFiles[num]))
paramFile <- read.csv(paste0(directory_exper1,"/",paramFiles[num]))
s_vals<-get_levels(paramFile$s_levels)
q_vals<-get_levels(paramFile$q_levels)
r_vals<-get_levels(paramFile$r_levels)


sFile <- subset(annFile,indType == "Sender")

unique(annFile$gen)
G<-max(annFile$gen)

sGen <- subset(sFile,gen==G)
rGen <- subset(rFile,gen==G)

#Put individuals together
dataMult<-data.frame()
for (n in 1:50){
  n_annS<-as.numeric(sGen[n,7:45])
  data<-senderPhenotypeDiscrete(2,2,n_annS)
  data$n <- n  
  dataMult<-rbind(dataMult,data)
}
p <- ggplot(dataMult) +
  geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=4) + 
  scale_color_viridis_c() +
  theme_bw() +
  facet_wrap(~n) +
  labs(title=paste0("Senders: ",num)) +
  labs(subtitle=paste0("gen = ",G,
                       "\nk = ",paramFile$k,
                       "\ncMin = ",paramFile$cMin,
                       "\ncMax = ",paramFile$cMax
  ))
p

paramFile


ggsave(plot=p,paste0(num,"_",G/10000,"_S.png"),
       device="png",path=directory_exper1,height=8,width=10,unit="in")


if (1==2){
  n<-1
  n_annS<-as.numeric(sGen[n,7:45])
  
  #senderANN_printSurface(30,n_annS)
  data<-senderPhenotype(70,70,n_annS)
  ggplot(data) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw()
}


#Receivers
paramFile$cMax
paramFile$cMin

unique(annFile$gen)
G<-2500

rFile <- subset(annFile,indType == "Receiver")
rGen <- subset(rFile,gen==G)
rEnd <- subset(rFile,gen==max(rFile$gen))

dataMult<-data.frame()
for (n in 1:40){
  n_annR<-as.numeric(rGen[n,7:45])
  data<-receiverPhenotypeDiscrete(2,2,n_annR)
  data$n <- n  
  dataMult<-rbind(dataMult,data)
}
p <- ggplot(dataMult) +
  geom_point(aes(x=s,y=r,color=output),alpha=0.8,size=4) + 
  scale_color_viridis_c() +
  theme_bw() +
  facet_wrap(~n) +
  labs(title=paste0("Receivers: ",num)) +
  labs(subtitle=paste0("gen = ",G,
                       "\nk = ",paramFile$k,
                       "\ncMin = ",paramFile$cMin,
                       "\ncMax = ",paramFile$cMax
  ))
p

ggsave(plot=p,paste0(num,"_",G/10000,"_R.png"),
       device="png",path=directory_exper1,height=8,width=10,unit="in")


#### Discrete phenotypes - probabilities ####
unique(sFile$gen)
G<-20000
sGen <- subset(sFile,gen==G)
rGen <- subset(rFile,gen==G)

probsAll<-data.frame()
for (n in 1:nrow(sGen)){
  n_annS<-as.numeric(sGen[n,7:45])
  d<-senderPhenotypeDiscrete(2,2,n_annS) 
  probs<-calcProbs(d,10)
  colnames(probs)<-c("state","strength","output","probFinal")
  probs$n <- n
  probs$indType = "Sender"
  probsAll<-rbind(probsAll,probs)
}
for (n in 1:nrow(rGen)){
  n_annR<-as.numeric(rGen[n,7:45])
  d<-receiverPhenotypeDiscrete(2,2,n_annR) 
  probs<-calcProbs(d,10)
  colnames(probs)<-c("state","strength","output","probFinal")
  probs$n <- n
  probs$indType = "Receiver"
  probs
  probsAll<-rbind(probsAll,probs)
}

# Plot for all levels = 2
ggplot(subset(probsAll,strength==1)) + 
  geom_jitter(aes(x=state,y=probFinal),height=0,width=0.05,alpha=.6) +
  facet_grid(indType~.) +
  theme_bw() +
  stat_summary(
    fun = "mean",
    geom = "point",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red",
    aes(x=state,y=probFinal) ) + 
  labs(subtitle="Received signal strength (receivers)\nQuality (senders)") +
  labs(y="Probability") +
  labs(x="Quality (senders) or \nReceived Signal Strength (receivers)") +
  ylim(0,1) +
  #Expected beta
  geom_point(aes(y=ifelse(indType=="Receiver",ifelse(paramFile$cMax<1,paramFile$cMax,1),-5)),x=1,size=5,color="blue") +
  #Expected alpha
  geom_point(aes(y=ifelse(indType=="Sender",ifelse(paramFile$cMax<1,(paramFile$m/(1-paramFile$m)),0),-5)),x=0,size=5,color="blue")

#### Discrete phenotypes - collate files ####
#for (G in unique(annFile$gen)){
unique(annFile$gen)
G<-20000
probsMaster<-data.frame()
for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directory_exper1,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directory_exper1,"/",paramFiles[num]))
  s_vals<-get_levels(paramFile$s_levels)
  q_vals<-get_levels(paramFile$q_levels)
  r_vals<-get_levels(paramFile$r_levels)
  sFile <- subset(annFile,indType == "Sender")
  rFile <- subset(annFile,indType == "Receiver")
  sGen <- subset(sFile,gen==G)
  rGen <- subset(rFile,gen==G)
  
  probsAll<-data.frame()
  for (n in 1:nrow(sGen)){
    n_annS<-as.numeric(sGen[n,7:45])
    d<-senderPhenotypeDiscrete(2,2,n_annS) 
    probs<-calcProbs(d,10)
    colnames(probs)<-c("state","strength","output","probFinal")
    probs$n <- n
    probs$indType = "Sender"
    probsAll<-rbind(probsAll,probs)
  }
  for (n in 1:nrow(rGen)){
    n_annR<-as.numeric(rGen[n,7:45])
    d<-receiverPhenotypeDiscrete(2,2,n_annR) 
    probs<-calcProbs(d,10)
    colnames(probs)<-c("state","strength","output","probFinal")
    probs$n <- n
    probs$indType = "Receiver"
    probs
    probsAll<-rbind(probsAll,probs)
  }
  probsAll$num<-num
  probsAll$mut_rate_ann <- paramFile$mut_rate_ann_S
  probsAll$mut_step_ann <- paramFile$mut_step_ann_S
  
  probsMaster<-rbind(probsMaster,probsAll)
}

paramFile

ggplot(subset(probsMaster,strength==1)) + 
  geom_jitter(aes(x=state,y=probFinal),height=0,width=0.15,alpha=.4) +
  theme_bw() +
  stat_summary(
    fun = "mean",
    geom = "point",
    col = "black",
    size = 3,
    shape = 24,
    fill = "red",
    aes(x=state,y=probFinal) ) + 
  labs(subtitle="Received signal strength (receivers)\nQuality (senders)") +
  labs(y="Probability") +
  labs(x="Quality (senders) or \nReceived Signal Strength (receivers)") +
  ylim(0,1) +
  #Expected beta
  geom_point(aes(y=ifelse(indType=="Receiver",ifelse(paramFile$cMax<1,paramFile$cMax,1),-5)),x=1,size=5,color="blue") +
  #Expected alpha
  geom_point(aes(y=ifelse(indType=="Sender",ifelse(paramFile$cMax<1,(paramFile$m/(1-paramFile$m)),0),-5)),x=0,size=5,color="blue") +
  facet_grid(mut_rate_ann+indType~mut_step_ann)


