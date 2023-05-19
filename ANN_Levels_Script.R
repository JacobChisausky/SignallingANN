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

#Test 1: Initial ####
directory_test1 <- "D:/StAndrews/ANN_Levels/tests1/testDiscrete_2"
annFiles <- list.files(directory_test1,"*annVars*")
paramFiles <- list.files(directory_test1,"*params_t*")

annFile <- read.csv(paste0(directory_test1,"/",annFiles[1]))

paramsAll<-data.frame()
for (i in paramFiles){
  paramFile <- read.csv(paste0(directory_test1,"/",i))
  paramsAll<-rbind(paramsAll,paramFile)
}
paramsAll

unique(annFile$gen)


#for (G in unique(annFile$gen)){
  for (num in 1:length(annFiles)){
    annFile <- read.csv(paste0(directory_test1,"/",annFiles[num]))
    paramFile <- read.csv(paste0(directory_test1,"/",paramFiles[num]))
    
    unique(annFile$gen)
    G <- max(sFile$gen)
    
    
    sFile <- subset(annFile,indType == "Sender")
    sGen <- subset(sFile,gen==G)
    sEnd <- subset(sFile,gen==max(sFile$gen))
    
    if (1==2){
      n<-1
      n_annS<-as.numeric(sGen[n,7:45])
      
      senderANN_printSurface(30,n_annS)
      
      data<-senderPhenotype(70,70,n_annS)
      ggplot(data) +
        geom_point(aes(x=q,y=s,color=output),alpha=0.8,size=2.7) + 
        scale_color_viridis_c() +
        theme_bw() +
        geom_function(fun = pred, colour = "red", linewidth=3)
    }
    
    #Put 20 individuals together
    dataMult<-data.frame()
    for (n in 1:20){
      n_annS<-as.numeric(sGen[n,7:45])
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
      labs(subtitle=paste0("gen = ",G,
                           "\nk = ",paramFile$k,
                           "\nc0 = ",paramFile$c0,
                           "\nc1 = ",paramFile$c1
                          ))
    
    ggsave(plot=p,paste0(num,"_",G/10000,"_S.png"),
           device="png",path=directory_test1,height=8,width=10,unit="in")
    
    
    rFile <- subset(annFile,indType == "Receiver")[,-(41:45)]
    rGen <- subset(rFile,gen==G)
    rEnd <- subset(rFile,gen==max(rFile$gen))
    
    dataAll<-data.frame()
    for (n in 1:100){
      n_annR<-as.numeric(rGen[n,-(1:6)])
      data<-receiverPhenotype(100,n_annR)
      data$n<-n
      dataAll<-rbind(dataAll,data)
    }
    pR<-ggplot(dataAll,aes(x=s,y=output)) + geom_path(aes(group=n),alpha=0.1,size=.8) +
      theme_bw() + ylim(0,1) + labs(y="Pr", title="100 individuals") +
      labs(title=num) +
      labs(subtitle=paste0("gen = ",G,
                           "\nk = ",paramFile$k,
                           "\nc0 = ",paramFile$c0,
                           "\nc1 = ",paramFile$c1
                           ))
    pR
    ggsave(plot=pR,paste0(num,"_",G/10000,"_R.png"),
           device="png",path=directory_test1,height=8,width=8,unit="in")
    
  }
#}

