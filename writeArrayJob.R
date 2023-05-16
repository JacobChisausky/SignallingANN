library("lubridate")

sh_BinaryName      <- "Program"  #compiled program name

hpc_time <- "0-12:00:00"
hpc_mem <- 8000

#auto_time <- TRUE #If true, ignore hpc_time and automatically calculate expected time

jobFolder <- "testNewNullBehavior"  #The folder path when you want the JOBS to be written to. To be created in this directiory

programPath <- "~/../../scratch/s5429412/SignallingANN/src/Program"  #Where is the compiled program located? Include program name

parallelReplicates <- 5   #This is different from 'replicates' below. In the below replicates, all jobs are set up to run sequentially
#In contrast, this will make replicates to be run in parallel. So, it will copy the entire job (including the number of replicates as given below)
# and write them each to a different job

replicates <- c(1)
k <- c(2)
seed <- c(0)  #If 0, then give each job a different random seed
N <- c(1000.0)
G <- c(400000)
c <- c(1.0)
p <- c(0.5)
d <- c(5.0)
init_ann_range <- c(1.0)
mut_rate_ann_S <- c(0.01)
mut_rate_ann_R <- c(0.01)
mut_step_ann_S <- c(0.01)
mut_step_ann_R <- c(0.01)
send_0_first <- c("true")  #as string
s_max <- c(10)  
interactionPartners <- c(10) 
fitnessFunction <- c(1)
complexInit <- c("true")
nullReceivers <- c("false")
nullSenders <- c("false")
nullHonestBeginG <- c(200000)
Report_annVar <- c(8)
Report_annVar_N <- c(100)
Report_annInit <- c("false")   #As string
dataFileName <- "testNewNull"   #as string
dataFileFolder <- "."   #as string


sameMutRates <- TRUE  #If true, mut_rate_ann_R will be overrided and replaced by mut_rate_ann_S 
sameMutSteps <- TRUE  #If true, mut_step_ann_R will be overrided and replaced by mut_step_ann_S

##-------------------------------------

#dataFileFolder<-paste0("\"",dataFileFolder,"\"")
#dataFileName <- paste0("\"",dataFileName,"\"")
fileName<-jobFolder

#Time allocation estimate####
# if (auto_time == TRUE){
#   minuteEst <- (0.000000005) * max(G) * max(N) * max(interactionPartners) * (1.5 / (max(interactionPartners)+5)) * max(replicates)
#   minuteEst <- max(minuteEst,2)
#   
#   #Give twice as much time as expected to be safe
#   minuteAlloc <- minuteEst * 2
#   dAlloc <- floor(minuteAlloc/1440)
#   rem <- minuteAlloc - dAlloc*1440
#   hAlloc <- floor(rem/60)
#   rem <- rem - hAlloc*60
#   if (hAlloc < 10){
#     hAlloc <- paste0("0",hAlloc)
#   }
#   if (rem < 10){
#     rem <- paste0("0",rem)
#   }
#   hpc_time <- paste0(dAlloc,"-",hAlloc,":",rem,":00")
# } ####

hpc_time

#Add time to folder path in same way as cpp - so that alphabetically, the newest folder is always last
t_hr <- format(Sys.time(), "%H")
t_min <- format(Sys.time(), "%M")
t_sec <- format(Sys.time(), "%S")
t_yday <- yday(Sys.Date())-1
date <- paste(t_yday,t_hr,t_min,t_sec,sep="_")

#folderPath<-paste0("./",date,"_",jobFolder)
folderPath<-paste0("./",jobFolder)

dir.create(folderPath)
fileNameList<-c()
iterator<-1;  #add this to each job title in case two jobs get submitted with same name at same second

dataFileName <- fileName
#dataFileFolder <- folderPath

#Move compiled program into new folder
file.copy(programPath,folderPath)

#Check to see if any files in folder share name and if so, increase iterator so as to not overwrite them

params<-list()
params<-list(replicates,
             k,
             seed,
             N,
             G,
             c,
             p,
             d,
             init_ann_range,
             mut_rate_ann_S,
             mut_rate_ann_R,
             mut_step_ann_S,
             mut_step_ann_R,
             send_0_first,
             s_max,
             interactionPartners,
             fitnessFunction,
             complexInit,
             nullReceivers,
             nullSenders,
             nullHonestBeginG,
             Report_annVar,
             Report_annVar_N,
             Report_annInit,
             dataFileName,
             dataFileFolder
)
paramsNames <- c("replicates",
                 "k",
                 "seed",
                 "N",
                 "G",
                 "c",
                 "p",
                 "d",
                 "init_ann_range",
                 "mut_rate_ann_S",
                 "mut_rate_ann_R",
                 "mut_step_ann_S",
                 "mut_step_ann_R",
                 "send_0_first",
                 "s_max",
                 "interactionPartners",
                 "fitnessFunction",
                 "complexInit",
                 "nullReceivers",
                 "nullSenders",
                 "nullHonestBeginG",
                 "Report_annVar",
                 "Report_annVar_N",
                 "Report_annInit",
                 "dataFileName",
                 "dataFileFolder"
)

#names<-c()  #This chunk is for setting the .json file names, based on params_in_file_name
#i<-1
#for (i in 1:length(params_in_file_name)){
#  names<-append(names,match(params_in_file_name[i],paramsNames))
#}
#names<-sort(names)

for (pR in 1:parallelReplicates){
  
  for (x1 in replicates){
    for (x2 in k){
      for (x3 in seed){
        for (x4 in N){
          for (x5 in G){
            for (x6 in c){
              for (x6a in p){
                for (x6b in d){
                for (x7 in init_ann_range){
                  for (x8 in mut_rate_ann_S){
                    if (sameMutRates == TRUE){
                      mut_rate_ann_R <- x8
                    }
                    for (x9 in mut_rate_ann_R){
                      for (x10 in mut_step_ann_S){
                        if (sameMutSteps == TRUE){
                          mut_step_ann_R <- x10
                        }
                        for (x11 in mut_step_ann_R){
                          for (x12 in send_0_first){
                            for (x13 in s_max){
                              for (x14 in interactionPartners){
                                for (x14a in fitnessFunction){
                                  for (x15 in complexInit){
                                    for (x16 in nullReceivers){
                                      for (x17 in nullSenders){
                                        for (x18 in nullHonestBeginG){
                                          for (x19 in Report_annVar){
                                            for (x20 in Report_annVar_N){
                                              for (x21 in Report_annInit){
                                                
                                                if (x3 == 0){
                                                  x3 <- sample(1:99999999, 1)
                                                }
                                                
                                                Xparams <- c(x1,x2,x3,x4,x5,x6,x6a,x6b,x7,x8,x9,x10,
                                                             x11,x12,x13,x14,x14a,x15,x16,x17,x18,x19,x20,
                                                             x21,paste0("\"",dataFileName,"_",iterator,"\""),paste0("\"",dataFileFolder,"/","\""))
                                                
                                                writeTo <- paste0(folderPath,"/",fileName,"_",iterator)
                                                txt <- "{\n\t"
                                                i <- 1
                                                for (i in 1:(length(paramsNames)-1)){
                                                  txt <- paste(txt,"\"",paramsNames[i],"\": ",Xparams[i],",\n\t",sep="")
                                                }
                                                txt <- paste(txt,"\"",paramsNames[i+1],"\": ",Xparams[i+1],"\n",
                                                             "}",sep="")
                                                
                                                fileConn<-file(paste0(writeTo,"_params.json"))
                                                writeLines(txt,fileConn)
                                                close(fileConn)
                                                
                                                fileNameList<-append(fileNameList,paste0(fileName,"_",iterator))
                                                
                                                iterator <- iterator + 1
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } 
}

txt<-paste0("#!/bin/bash\n#SBATCH --job-name=",fileName,
            "\n#SBATCH --mail-user=jacob.chisausky@evobio.eu\n#SBATCH --mail-type=ALL\n#SBATCH --time=",
            hpc_time,"\n#SBATCH --mem=",hpc_mem,"\n#SBATCH --array=1-",iterator-1,
            "\n\nmodule load foss/2022a\n\necho \"SLURM_JOBID: \" $SLURM_JOBID\necho \"SLURM_ARRAY_TASK_ID: \" $SLURM_ARRAY_TASK_ID\necho \"SLURM_ARRAY_JOB_ID: \" $SLURM_ARRAY_JOB_ID\n\n./",sh_BinaryName," `ls -d *.json | awk NR==$SLURM_ARRAY_TASK_ID`")

fileConn<-file(paste(folderPath,"/SignallingModel",".sh",sep=""))
writeLines(txt,fileConn)
close(fileConn)
