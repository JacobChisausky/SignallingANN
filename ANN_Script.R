#The point of this code is to read ANN stats from the cpp program and allow us to visualize the phenotypes
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
}
senderANN <- function(s,q,annS){
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
}
randANN_S <- function(lim){
  return(runif(39, min=-abs(lim), max=abs(lim)))
}
senderPhenotype <- function(s_num, q_num, ANN_S){
  #s_num = number of strength values to plot. More = finer plot
  #q_num = number of quality values to calculate
  
  phenotype_S <- data.frame()
  
  for (q_cur in 0:q_num){
    q <- q_cur/q_num
    for (s_cur in 0:s_num){
      s <- s_cur/s_num
      out <- senderANN(s,q,ANN_S)
      row <- c(q,s,out)
      phenotype_S <- rbind(phenotype_S,row)
    }
  }
  colnames(phenotype_S) <- c("q","s","output")
  return(phenotype_S)
} #Returns a data fram for s_num s values and q_num q vales








#Senders ####

phenotype <- senderPhenotype(10,10,randANN_S(1))
while (sum(unique(phenotype$output)) < 10){
  phenotype <- senderPhenotype(10,10,randANN_S(1))
  #I could do this in cpp at initialization to make more interesting starting phenotypes?
}

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
                                         tickformat = '.2f'),
                            yaxis = list(title = "q", 
                                         ticketmode = 'array',
                                         ticktext = round(as.numeric(qlabs),2),
                                         tickvals = (as.numeric(qlabs)*(nrow(matrix)-1)),
                                         tickformat = '.2f'),
                            zaxis = list(title = "Pr")))

p






##Reading ANNs from file and displaying
directory <- "C:/Users/owner/Documents/S4/Simulation_ANN"
file <- "116_17_38_46_annVars_forR.csv"

annFile <- read.csv(paste0(directory,"/",file))
sFile <- subset(annFile,indType == "Sender")
rFile <- subset(annFile,indType == "Receiver")[,-(((ncol(rFile))-4):ncol(rFile))]
sEnd <- subset(sFile,gen==max(sFile$gen))
rEnd <- subset(rFile,gen==max(rFile$gen))

# Senders ####
n<-2
n_annS<-as.numeric(sFile[n,7:45])

phenotype <- senderPhenotype(10,10,n_annS)

{
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
                                           tickformat = '.2f'),
                              yaxis = list(title = "q", 
                                           ticketmode = 'array',
                                           ticktext = round(as.numeric(qlabs),2),
                                           tickvals = (as.numeric(qlabs)*(nrow(matrix)-1)),
                                           tickformat = '.2f'),
                              zaxis = list(title = "Pr")))
}  
  p






