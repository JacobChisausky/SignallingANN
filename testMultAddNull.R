#The point of this code is to read ANN stats from the cpp program and allow us to visualize the phenotypes
#set.seed(1234567) 

library(ggplot2)
library(plotly)
library(reshape2)

#General Functions ####
pred<-function(x){
  return( (1/(2+10*x))^(1/(1+10*x)) )
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

##Read data - before code changes####
directory <- "C:/Users/owner/Documents/S4/Simulation_ANN/Test_MultAdd"
annFiles <- list.files(directory,"*annVars*")
paramFiles <- list.files(directory,"*params_t*")

#1 ff:0 k:2 bad. Signal only when q>.9,s>.9
#2 ff:1 k:2 Good - though s is a bit high for larger qs.
#3 ff:0 k:10 bad: signal only when high q. A bit diagonal (meaning s and q corellate positively) but only at very high q
#4 ff:1 k:10 bad. Threshold - always signal strongly
#5 ff:0 k:2 bad. Only signal at very high q. Blob in top right corner
#6 ff:1 k:2 bad: diagonal but send 0 first always means always 0
#7 ff:0 k:10 bad: top right corner
#8 ff:1 k:10 bad - same as 5
#9 ff:0 k:2 top right corner
#10 ff:1 k:2 good - diagonal
#11 ff:1 k:2 good -diagonal
#12 ff:0 k:10 spot in top right corner. A little diagonal

#To test - try more honest generations
#To test - just null senders and receivers
#Problem - c was 10 for these

num <- 12
annFile <- read.csv(paste0(directory,"/",annFiles[num]))
paramFile <- read.csv(paste0(directory,"/",paramFiles[num]))
paramFile
paramFile$fitnessFunction
paramFile$nullHonestBeginG

sFile <- subset(annFile,indType == "Sender")
sPreNull <- subset(sFile,gen==150000)
sNull <- subset(sFile,gen==paramFile$nullHonestBeginG)
sEnd <- subset(sFile,gen==max(sFile$gen))

n<-2
n_annS<-as.numeric(sPreNull[n,7:45])
senderANN_printSurface(30,n_annS)

data<-senderPhenotype(70,70,n_annS)
ggplot(data) +
  geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
  scale_color_viridis_c() +
  theme_bw() +
  geom_function(fun = pred, colour = "red", linewidth=3)


#Real data - receivers
directory <- "C:/Users/owner/Documents/S4/Simulation_ANN/Test_MultAdd"
annFiles <- list.files(directory,"*annVars*")
paramFiles <- list.files(directory,"*params_t*")

num <- 12
annFile <- read.csv(paste0(directory,"/",annFiles[num]))
paramFile <- read.csv(paste0(directory,"/",paramFiles[num]))
paramFile$fitnessFunction
paramFile$nullHonestBeginG
paramFile
#1 ff 0. k 2: good. Low s cut off
#2 ff 1. k 2: good. Low s cut off
#3 ff 0. k 10: very good. Less cut off.
#6 ff 1. k 2. good. cutoff at s=.15
#10 ff 1. k 2. pretty good. cutoff .15
#11 ff 1. k 2. good. 

#4 ff 1. k 10. decent. Pr < s, should be =
#5 ff 0. k 2. decent. Large cutoff. good above s=.375
#7 ff 0. k 10. poor. Cutoff at s=.5. some threshold, some fit expected
#8 ff 1. k 10. decent. cutoff at .375
#9 ff 0. k 2. decent. cutoff .25
#12 ff 0. k 10. decent. cutoff .3



rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
rNull <- subset(rFile,gen==paramFile$nullHonestBeginG)
rEnd <- subset(rFile,gen==max(rFile$gen))

#n<-1
#n_annR<-as.numeric(rNull[n,-(1:6)])
#data<-receiverPhenotype(100,n_annR)
#ggplot(data,aes(x=s,y=output)) + geom_line(color="red",size=1.4) +
#  theme_bw() + ylim(0,1) + labs(y="Pr")

dataAll<-data.frame()
for (n in 1:100){
  n_annR<-as.numeric(rNull[n,-(1:6)])
  data<-receiverPhenotype(100,n_annR)
  data$n<-n
  dataAll<-rbind(dataAll,data)
}
ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
  theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals")

#Test 1 - Mid. LongHonest ####
#Note on these results:
#ff=1 and k=2 is best. k=5 or ff=0 just gives flat 'always send strong signal'.
#Why? I don't know. I can redo analytics with additive fitness in desmos
directoryT1_LongHonest <- "D:/StAndrews/SignallingANN/code_MultAdd/codeMultAdd_longHonest"
annFiles <- list.files(directoryT1_LongHonest,"*annVars*")
paramFiles <- list.files(directoryT1_LongHonest,"*params_t*")

#1,2,3 look bad
#but 4 looks good!
for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directoryT1_LongHonest,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directoryT1_LongHonest,"/",paramFiles[num]))
  paramFile$fitnessFunction
  
  paramFile$nullHonestBeginG
  
  sFile <- subset(annFile,indType == "Sender")
  sNull <- subset(sFile,gen==390000)
  sEnd <- subset(sFile,gen==max(sFile$gen))
  
  if (1==2){
    n<-4
    n_annS<-as.numeric(sNull[n,7:45])
    senderANN_printSurface(30,n_annS)
    
    data<-senderPhenotype(70,70,n_annS)
    ggplot(data) +
      geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
      scale_color_viridis_c() +
      theme_bw() +
      geom_function(fun = pred, colour = "red", linewidth=3)
  }
  
  #Put 30 individuals together
  dataMult<-data.frame()
  for (n in 1:20){
    n_annS<-as.numeric(sNull[n,7:45])
    data<-senderPhenotype(30,30,n_annS)
    data$n <- n  
    dataMult<-rbind(dataMult,data)
  }
  p <- ggplot(dataMult) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw() +
    facet_wrap(~n) +
    geom_function(fun = pred, colour = "red", linewidth=2) +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders
    ))
  
  ggsave(plot=p,paste0(num,"_S.png"),
         device="png",path=directoryT1_LongHonest,height=8,width=10,unit="in")
  
  rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
  rNull <- subset(rFile,gen==390000)
  rEnd <- subset(rFile,gen==max(rFile$gen))
  
  dataAll<-data.frame()
  for (n in 1:100){
    n_annR<-as.numeric(rNull[n,-(1:6)])
    data<-receiverPhenotype(100,n_annR)
    data$n<-n
    dataAll<-rbind(dataAll,data)
  }
  pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
    theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders))
  
  ggsave(plot=pR,paste0(num,"_R.png"),
         device="png",path=directoryT1_LongHonest,height=8,width=8,unit="in")
  
}

#Test 1 - End. LongHonestEnd ####
#Note on these results:
#ff=1 and k=2 is best. k=5 or ff=0 just gives flat 'always send strong signal'.
#Why? I don't know. I can redo analytics with additive fitness in desmos
directoryT1_LongHonest <- "D:/StAndrews/SignallingANN/code_MultAdd/codeMultAdd_longHonest"
annFiles <- list.files(directoryT1_LongHonest,"*annVars*")
paramFiles <- list.files(directoryT1_LongHonest,"*params_t*")

for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directoryT1_LongHonest,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directoryT1_LongHonest,"/",paramFiles[num]))
  paramFile$fitnessFunction
  
  paramFile$nullHonestBeginG
  
  sFile <- subset(annFile,indType == "Sender")
  sNull <- subset(sFile,gen==390000)
  sEnd <- subset(sFile,gen==max(sFile$gen))
  
  #Put 30 individuals together
  dataMult<-data.frame()
  for (n in 1:20){
    n_annS<-as.numeric(sEnd[n,7:45])
    data<-senderPhenotype(30,30,n_annS)
    data$n <- n  
    dataMult<-rbind(dataMult,data)
  }
  p <- ggplot(dataMult) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw() +
    facet_wrap(~n) +
    geom_function(fun = pred, colour = "red", linewidth=2) +
    labs(title=paste0("End: ",num)) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders
    ))
  
  ggsave(plot=p,paste0(num,"_EndS.png"),
         device="png",path=directoryT1_LongHonest,height=8,width=10,unit="in")
  
  rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
  rNull <- subset(rFile,gen==390000)
  rEnd <- subset(rFile,gen==max(rFile$gen))
  
  dataAll<-data.frame()
  for (n in 1:100){
    n_annR<-as.numeric(rEnd[n,-(1:6)])
    data<-receiverPhenotype(100,n_annR)
    data$n<-n
    dataAll<-rbind(dataAll,data)
  }
  pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
    theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
    labs(title=paste0("End: ",num)) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders))
  
  ggsave(plot=pR,paste0(num,"_EndR.png"),
         device="png",path=directoryT1_LongHonest,height=8,width=8,unit="in")
  
}


#Test 1 - null senders ####
directoryT1_NullSenders <- "D:/StAndrews/SignallingANN/code_MultAdd/codeMultAdd_nullSenders"
annFiles <- list.files(directoryT1_NullSenders,"*annVars*")
paramFiles <- list.files(directoryT1_NullSenders,"*params_t*")

for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directoryT1_NullSenders,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directoryT1_NullSenders,"/",paramFiles[num]))
  paramFile$fitnessFunction
  
  paramFile$nullHonestBeginG
  
  sFile <- subset(annFile,indType == "Sender")
  #sNull <- subset(sFile,gen==390000)
  sEnd <- subset(sFile,gen==max(sFile$gen))
  
  #Put 30 individuals together
  dataMult<-data.frame()
  for (n in 1:20){
    n_annS<-as.numeric(sEnd[n,7:45])
    data<-senderPhenotype(30,30,n_annS)
    data$n <- n  
    dataMult<-rbind(dataMult,data)
  }
  p <- ggplot(dataMult) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw() +
    facet_wrap(~n) +
    geom_function(fun = pred, colour = "red", linewidth=2) +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders
    ))
  
  ggsave(plot=p,paste0(num,"_S.png"),
         device="png",path=directoryT1_NullSenders,height=8,width=10,unit="in")
  
  rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
  #rNull <- subset(rFile,gen==390000)
  rEnd <- subset(rFile,gen==max(rFile$gen))
  
  dataAll<-data.frame()
  for (n in 1:100){
    n_annR<-as.numeric(rEnd[n,-(1:6)])
    data<-receiverPhenotype(100,n_annR)
    data$n<-n
    dataAll<-rbind(dataAll,data)
  }
  pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
    theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders))
  
  ggsave(plot=pR,paste0(num,"_R.png"),
         device="png",path=directoryT1_NullSenders,height=8,width=8,unit="in")
  
}

#Test 1 - null receivers ####
#Try with null receivers respons with Pr = 0.5s.
#This will result in senders signalling if q>0.5 only
#Then have receivers evolve to this: senders signal at s = q

directoryT1_NullReceivers <- "D:/StAndrews/SignallingANN/code_MultAdd/codeMultAdd_nullReceivers"
annFiles <- list.files(directoryT1_NullReceivers,"*annVars*")
paramFiles <- list.files(directoryT1_NullReceivers,"*params_t*")

for (num in 1:length(annFiles)){
  annFile <- read.csv(paste0(directoryT1_NullReceivers,"/",annFiles[num]))
  paramFile <- read.csv(paste0(directoryT1_NullReceivers,"/",paramFiles[num]))
  paramFile$fitnessFunction
  
  sFile <- subset(annFile,indType == "Sender")
  
  unique(sFile$gen)
  paramFile$nullHonestBeginG
  
  #sNull <- subset(sFile,gen==390000)
  sEnd <- subset(sFile,gen==max(sFile$gen))
  
  if (1==2){
    n<-4
    n_annS<-as.numeric(sEnd[n,7:45])
    senderANN_printSurface(30,n_annS)
    
    data<-senderPhenotype(70,70,n_annS)
    ggplot(data) +
      geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
      scale_color_viridis_c() +
      theme_bw() +
      geom_function(fun = pred, colour = "red", linewidth=3)
  }
  
  #Put 30 individuals together
  dataMult<-data.frame()
  for (n in 1:20){
    n_annS<-as.numeric(sEnd[n,7:45])
    data<-senderPhenotype(30,30,n_annS)
    data$n <- n  
    dataMult<-rbind(dataMult,data)
  }
  p <- ggplot(dataMult) +
    geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
    scale_color_viridis_c() +
    theme_bw() +
    facet_wrap(~n) +
    geom_function(fun = pred, colour = "red", linewidth=2) +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders
    ))
  
  ggsave(plot=p,paste0(num,"_S.png"),
         device="png",path=directoryT1_NullReceivers,height=8,width=10,unit="in")
  
  rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
  #rNull <- subset(rFile,gen==390000)
  rEnd <- subset(rFile,gen==max(rFile$gen))
  
  dataAll<-data.frame()
  for (n in 1:100){
    n_annR<-as.numeric(rEnd[n,-(1:6)])
    data<-receiverPhenotype(100,n_annR)
    data$n<-n
    dataAll<-rbind(dataAll,data)
  }
  pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
    theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
    labs(title=num) +
    labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                         "\nk = ",paramFile$k,
                         "\nnullReceivers = ", paramFile$nullReceivers,                
                         "\nnullSenders = ", paramFile$nullSenders))
  
  ggsave(plot=pR,paste0(num,"_R.png"),
         device="png",path=directoryT1_NullReceivers,height=8,width=8,unit="in")
  
}



#Test 2 LongHonest, after code changes ####
directoryT2_longHonest <- "D:/StAndrews/SignallingANN/code_MA_tests2/test_2_multFitness_doubles"
annFiles <- list.files(directoryT2_longHonest,"*annVars*")
paramFiles <- list.files(directoryT2_longHonest,"*params_t*")

for (G in c(400000,450000,500000,550000,600000)){
  for (num in 1:length(annFiles)){
    annFile <- read.csv(paste0(directoryT2_longHonest,"/",annFiles[num]))
    paramFile <- read.csv(paste0(directoryT2_longHonest,"/",paramFiles[num]))
    paramFile$fitnessFunction
    
    paramFile$nullHonestBeginG
    unique(annFile$gen)
    
    sFile <- subset(annFile,indType == "Sender")
    sNull <- subset(sFile,gen==G)
    sEnd <- subset(sFile,gen==max(sFile$gen))
    
    if (1==2){
      n<-1
      n_annS<-as.numeric(sNull[n,7:45])
      
      senderANN_printSurface(30,n_annS)
      
      data<-senderPhenotype(70,70,n_annS)
      ggplot(data) +
        geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
        scale_color_viridis_c() +
        theme_bw() +
        geom_function(fun = pred, colour = "red", linewidth=3)
    }
    
    #Put 30 individuals together
    dataMult<-data.frame()
    for (n in 1:20){
      n_annS<-as.numeric(sNull[n,7:45])
      data<-senderPhenotype(30,30,n_annS)
      data$n <- n  
      dataMult<-rbind(dataMult,data)
    }
    p <- ggplot(dataMult) +
      geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
      scale_color_viridis_c() +
      theme_bw() +
      facet_wrap(~n) +
      geom_function(fun = pred, colour = "red", linewidth=2) +
      labs(title=num) +
      labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                           "\ngen = ",G,
                           "\nk = ",paramFile$k,
                           "\nnullReceivers = ", paramFile$nullReceivers,                
                           "\nnullSenders = ", paramFile$nullSenders
      ))
    
    ggsave(plot=p,paste0(num,"_",G,"_S.png"),
           device="png",path=directoryT2_longHonest,height=8,width=10,unit="in")
    
    rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
    rNull <- subset(rFile,gen==G)
    rEnd <- subset(rFile,gen==max(rFile$gen))
    
    dataAll<-data.frame()
    for (n in 1:100){
      n_annR<-as.numeric(rNull[n,-(1:6)])
      data<-receiverPhenotype(100,n_annR)
      data$n<-n
      dataAll<-rbind(dataAll,data)
    }
    pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
      theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
      labs(title=num) +
      labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                           "\ngen = ",G,
                           "\nk = ",paramFile$k,
                           "\nnullReceivers = ", paramFile$nullReceivers,                
                           "\nnullSenders = ", paramFile$nullSenders))
    
    ggsave(plot=pR,paste0(num,"_",G,"_R.png"),
           device="png",path=directoryT2_longHonest,height=8,width=8,unit="in")
    
  }
}

#Test 2 additive ####
directoryT2_add <- "D:/StAndrews/SignallingANN/code_MA_tests2/test_addFit_p"
annFiles <- list.files(directoryT2_add,"*annVars*")
paramFiles <- list.files(directoryT2_add,"*params_t*")

paramsAll<-data.frame()
for (num in 1:length(annFiles)){
  paramFile <- read.csv(paste0(directoryT2_add,"/",paramFiles[num]))
  paramFile$num <- num
  paramsAll<-rbind(paramFile,paramsAll)
}
paramsAll
#Parmas - mut steps, mut rates, k
#Data at 50000 100000 150000 200000 250000 300000 3500000 400000
#Null until 200000

for (G in c(200000, 250000, 300000, 350000, 400000)){
  for (num in 1:length(annFiles)){
    annFile <- read.csv(paste0(directoryT2_add,"/",annFiles[num]))
    paramFile <- read.csv(paste0(directoryT2_add,"/",paramFiles[num]))
    paramFile$fitnessFunction
    paramFile$nullHonestBeginG
    unique(annFile$gen)
    
    sFile <- subset(annFile,indType == "Sender")
    sNull <- subset(sFile,gen==G)
    sEnd <- subset(sFile,gen==max(sFile$gen))
    
    if (1==2){
      n<-1
      n_annS<-as.numeric(sNull[n,7:45])
      
      senderANN_printSurface(30,n_annS)
      
      data<-senderPhenotype(70,70,n_annS)
      ggplot(data) +
        geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
        scale_color_viridis_c() +
        theme_bw() +
        geom_function(fun = pred, colour = "red", linewidth=3)
    }
    
    #Put 30 individuals together
    dataMult<-data.frame()
    for (n in 1:20){
      n_annS<-as.numeric(sNull[n,7:45])
      data<-senderPhenotype(30,30,n_annS)
      data$n <- n  
      dataMult<-rbind(dataMult,data)
    }
    p <- ggplot(dataMult) +
      geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
      scale_color_viridis_c() +
      theme_bw() +
      facet_wrap(~n) +
      geom_function(fun = pred, colour = "red", linewidth=2) +
      labs(title=num) +
      labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                           "\ngen = ",G,
                           "\nmutRates = ", paramFile$mut_rate_ann_R,
                           "\nmutSteps = ", paramFile$mut_step_ann_R,
                           "\nk = ",paramFile$k
                            ))
    
    ggsave(plot=p,paste0(num,"_",G/10000,"_S.png"),
           device="png",path=directoryT2_add,height=8,width=10,unit="in")
    
    rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
    rNull <- subset(rFile,gen==G)
    rEnd <- subset(rFile,gen==max(rFile$gen))
    
    dataAll<-data.frame()
    for (n in 1:100){
      n_annR<-as.numeric(rNull[n,-(1:6)])
      data<-receiverPhenotype(100,n_annR)
      data$n<-n
      dataAll<-rbind(dataAll,data)
    }
    pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
      theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
      labs(title=num) +
      labs(subtitle=paste0("ff = ",paramFile$fitnessFunction,
                           "\ngen = ",G,
                           "\nmutRates = ", paramFile$mut_rate_ann_R,
                           "\nmutSteps = ", paramFile$mut_step_ann_R,
                           "\nk = ",paramFile$k))
    ggsave(plot=pR,paste0(num,"_",G/10000,"_R.png"),
           device="png",path=directoryT2_add,height=8,width=8,unit="in")
    
  }
}
