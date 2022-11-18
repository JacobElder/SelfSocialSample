##Variables of interest######################################################################################################################################################################

#working directory, use slashes and not backslashes
setwd("C:\\Users\\Alex Koch\\Desktop\\Q-SpAM\\Processing and visualizing similarity data collected with Q-SpAM")
inputFile<-"Q-SpAM (mapping the domain of influential leaders).xml"
numberOfLabels<-62
numberOfConditions<-0

#vector containing desired dimensions for MDS
dimensions<-c(2:5)

##Importing##########################################################################################################################################################################################

#Check if required packages are installed
if ("ggplot2" %in% installed.packages() == FALSE) {install.packages("ggplot2")}
#if ("car" %in% installed.packages() == FALSE) {install.packages("car")}
if ("smacof" %in% installed.packages() == FALSE) {install.packages("smacof")}
if ("scatterplot3d" %in% installed.packages() == FALSE) {install.packages("scatterplot3d")}
if ("rio" %in% installed.packages() == FALSE) {install.packages("rio")}
if ("ggrepel" %in% installed.packages() == FALSE) {install.packages("ggrepel")}

#Load libraries
library(ggplot2)
#library(car)
library(smacof)
library(scatterplot3d)
library(rio)
library(ggrepel)

#import data
spam<-import(inputFile)

#Define missing values
spam[spam==-9999]<-NA

#reclassing
spam$xwindowsize<-as.numeric(spam$xwindowsize)
spam$ywindowsize<-as.numeric(spam$ywindowsize)

#Extract label names
selection<-NULL
for (i in 1:numberOfLabels){
  selection<-c(selection, paste0("label",i,"name"))
}
labelNames<-vector("list", numberOfLabels)

labelNameSelection<-as.data.frame(spam[selection])

for (j in 1:numberOfLabels){
  for (i in 1:nrow(spam))
    if (labelNameSelection[i,j]!="leer"){
      labelNames[[j]]<-as.character(labelNameSelection[i,j])
      break
    }
}

#Flip y-Axis for coordinates
selection<-NULL
for (i in 1:numberOfLabels){
  selection<-c(selection, paste0("label",i,"y"))
}
for (i in 1:nrow(spam)){
  spam[i,selection]<-sapply(as.numeric(spam[i,selection]), function(x){x*-1+spam$ywindowsize[i]})
}

#Extract coordinate matrix
selection<-NULL
for (i in 1:numberOfLabels){
  selection<-c(selection, paste0("label",i,"x"))
  selection<-c(selection, paste0("label",i,"y"))
}

#No conditions
coordinates<-as.data.frame(spam[selection])
names(coordinates)<-names(spam[selection])
#testing
coordinates$id<-1:nrow(coordinates)

#with conditions
if (numberOfConditions!=0){
  coordinatesByCondition<-vector("list", numberOfConditions)
  for (i in 1:numberOfConditions){
    assign(paste0("coordinates",i),coordinates[which(spam$condition==i),])
    coordinatesByCondition[[i]]<-coordinates[which(spam$condition==i),]
  }
}

#Extract label input
selection<-NULL
for (i in 1:numberOfLabels){
  selection<-c(selection, paste0("label",i,"input"))
}
labelInputs<-as.data.frame(spam[selection])

#create folders
dir.create("participant coordinates", showWarnings=FALSE)
dir.create("participant graphs", showWarnings=FALSE)
dir.create("participant distances", showWarnings=FALSE)
dir.create("MDS graphs", showWarnings=FALSE)
dir.create("MDS solutions", showWarnings=FALSE)
if (numberOfConditions!=0){
  for (i in 1:numberOfConditions){
    dir.create(paste0("participant coordinates/condition",i), showWarnings=FALSE)
    dir.create(paste0("participant graphs/condition",i), showWarnings=FALSE)
    dir.create(paste0("participant distances/condition",i), showWarnings=FALSE)
  }
}

##Calculating distances###########################################################################################################################################################################################

#No conditions
if (numberOfConditions==0){
  #Calculate coordinate and distance matrices and save them in a list
  allCoordinates<-vector("list", nrow(spam))
  allDistances<-vector("list", nrow(spam))
  
  for (i in 1:nrow(coordinates)){
    tempCoordinates<-matrix(as.numeric(coordinates[i,1:(numberOfLabels*2)]), ncol=2, byrow=TRUE, dimnames=list(labelNames, c("x","y")))
    tempParticipantNumber<-coordinates[i,(numberOfLabels*2+1)]
    #Export graph
    try(png(filename=paste0("participant graphs/",tempParticipantNumber,".png")))
    print(ggplot(as.data.frame(tempCoordinates), aes(x=x, y=y)) +
            geom_point(shape=1) +
            xlim(0, spam$xwindowsize[tempParticipantNumber]) +
            ylim(0, spam$ywindowsize[tempParticipantNumber]) +
            #scale_x_continuous(limits=c(0, spam$xwindowsize[tempParticipantNumber]), breaks=seq(0, spam$xwindowsize[tempParticipantNumber], 100)) +
            #scale_y_continuous(limits=c(0, spam$ywindowsize[tempParticipantNumber]), breaks=seq(0, spam$ywindowsize[tempParticipantNumber], 100)) +
            #coord_fixed(ratio=spam$xwindowsize[tempParticipantNumber]/spam$ywindowsize[tempParticipantNumber], xlim=c(0, spam$xwindowsize[tempParticipantNumber]) , ylim=c(0, spam$ywindowsize[tempParticipantNumber])) +
            geom_text_repel(aes(x=x, y=y, label=rownames(tempCoordinates))))
    dev.off()
    #Save and export coordinates that are relative to screen size
    tempCoordinates<-apply(tempCoordinates, 2, function(x) x/sqrt(spam$xwindowsize[tempParticipantNumber]^2+spam$ywindowsize[tempParticipantNumber]^2))
    assign(paste0("ppn",i,"coordinates"),tempCoordinates)
    write.table(tempCoordinates,file = paste("participant coordinates/ppn",tempParticipantNumber,".csv", sep=""),sep=",",row.names=FALSE,col.names=TRUE)
    allCoordinates[[tempParticipantNumber]]<-(tempCoordinates)
    
    #Save and export distances
    tempDistances<-as.matrix(dist(tempCoordinates))
    assign(paste0("ppn",tempParticipantNumber,"distances"),tempDistances)
    write.table(as.matrix(tempDistances),file = paste("participant distances/ppn",tempParticipantNumber,".csv", sep=""),sep=",",row.names=FALSE,col.names=TRUE)
    allDistances[[tempParticipantNumber]]<-tempDistances
  }
} else {
#With conditions  
  #Calculate coordinate and distance matrices and save them in a list
  allCoordinates<-vector("list", numberOfConditions)
  allDistances<-vector("list", numberOfConditions)
  
  allCoordinatesByCondition<-vector("list", numberOfConditions)
  allDistancesByCondition<-vector("list", numberOfConditions)
  for (j in 1:numberOfConditions){
    tempCoordinatesList<-vector("list", nrow(coordinatesByCondition[[j]]))
    tempDistancesList<-vector("list", nrow(coordinatesByCondition[[j]]))
    for (i in 1:nrow(coordinatesByCondition[[j]])){
      tempCoordinates<-matrix(as.numeric(coordinatesByCondition[[j]][i,1:(numberOfLabels*2)]), ncol=2, byrow=TRUE, dimnames=list(labelNames, c("x","y")))
      tempParticipantNumber<-coordinatesByCondition[[j]][i,(numberOfLabels*2+1)]
      #Export graph
      try(png(filename=paste0("participant graphs/condition", j, "/", tempParticipantNumber,".png")))
      
      print(ggplot(as.data.frame(tempCoordinates), aes(x=x, y=y)) +
              geom_point(shape=1) +
              xlim(0, spam$xwindowsize[tempParticipantNumber]) +
              ylim(0, spam$ywindowsize[tempParticipantNumber]) +
              #scale_x_continuous(limits=c(0, spam$xwindowsize[tempParticipantNumber]), breaks=seq(0, spam$xwindowsize[tempParticipantNumber], 100)) +
              #scale_y_continuous(limits=c(0, spam$ywindowsize[tempParticipantNumber]), breaks=seq(0, spam$ywindowsize[tempParticipantNumber], 100)) +
              #coord_fixed(ratio=spam$xwindowsize[tempParticipantNumber]/spam$ywindowsize[tempParticipantNumber], xlim=c(0, spam$xwindowsize[tempParticipantNumber]) , ylim=c(0, spam$ywindowsize[tempParticipantNumber])) +
              geom_text_repel(aes(x=x, y=y, label=rownames(tempCoordinates))))
      dev.off()
      #ggsave(path=paste0("participant graphs/condition", j, "/"), filename=paste0(tempParticipantNumber, ".jpg"), device="jpg", dpi=50)
      #Save and export coordinates that are relative to screen size
      tempCoordinates<-apply(tempCoordinates, 2, function(x) x/sqrt(spam$xwindowsize[tempParticipantNumber]^2+spam$ywindowsize[tempParticipantNumber]^2))
      assign(paste0("ppn",tempParticipantNumber,"coordinates"),tempCoordinates)
      write.table(tempCoordinates,file = paste("participant coordinates/condition",j,"/ppn",rownames(coordinatesByCondition[[j]])[i],".csv", sep=""),sep=",",row.names=FALSE,col.names=TRUE)
      tempCoordinatesList[[tempParticipantNumber]]<-(tempCoordinates)
      
      #Save and export distances
      tempDistances<-as.matrix(dist(tempCoordinates))
      tempDistancesList[[tempParticipantNumber]]<-as.matrix(dist(tempCoordinates))
      assign(paste0("ppn",tempParticipantNumber,"distances"),tempDistances)
      write.table(as.matrix(tempDistances),file = paste("participant distances/condition",j,"/ppn",rownames(coordinatesByCondition[[j]])[tempParticipantNumber],".csv", sep=""),sep=",",row.names=FALSE,col.names=TRUE)
    }
    assign(paste0("Coordinates",j),tempCoordinatesList)
    assign(paste0("Distances",j),tempDistancesList)
    allCoordinatesByCondition[[j]]<-tempCoordinatesList
    allDistancesByCondition[[j]]<-tempDistancesList
    
  }
}

##Averaging################################################################################################################################################################################################

#No conditions
if (numberOfConditions==0){
  averageCoordinates<-apply(array(unlist(allCoordinates), dim=c(numberOfLabels, numberOfLabels, nrow(spam))), 1:2, mean, na.rm=TRUE)
  write.table(averageCoordinates,file = "averageCoordinates.csv",sep=",",row.names=FALSE,col.names=TRUE)
  
  averageDistances<-apply(array(unlist(allDistances), dim=c(numberOfLabels, numberOfLabels, nrow(spam))), 1:2, mean, na.rm=TRUE)
  write.table(averageDistances,file = "averageDistances.csv",sep=",",row.names=FALSE,col.names=TRUE)
  
  #Export graph
  #try(png(filename="averageCoordinates.png"))
  #try(scatterplot(x=as.numeric(averageCoordinates[,1]), y=as.numeric(averageCoordinates[,2]), cex=1, cex.sub=4, smoother=FALSE, xlim=c(0,1), ylim=c(0,1), id.n=numberOfLabels, labels=labelNames, reg.line=FALSE, boxplots=FALSE))
  #dev.off()
  
} else {
  #With conditions
  averageCoordinatesByCondition<-vector("list", numberOfConditions)
  averageDistancesByCondition<-vector("list", numberOfConditions)
  
  for (i in 1:numberOfConditions){
    tempAverageCoordinates<-apply(array(unlist(allCoordinatesByCondition), dim=c(numberOfLabels, numberOfLabels, nrow(spam))), 1:2, mean, na.rm=TRUE)
    assign(paste0("averageCoordinates",i), tempAverageCoordinates)
    averageCoordinatesByCondition[[i]]<-tempAverageCoordinates
    write.table(tempAverageCoordinates,file = paste0("averageCoordinatesCond",i,".csv"),sep=",",row.names=FALSE,col.names=TRUE)
    
    tempAverageDistances<-apply(array(unlist(allDistancesByCondition[[i]]), dim=c(numberOfLabels, numberOfLabels, nrow(spam))), 1:2, mean, na.rm=TRUE)
    assign(paste0("averageDistances",i), tempAverageDistances)
    averageDistancesByCondition[[i]]<-tempAverageDistances
    write.table(tempAverageDistances,file = paste0("averageDistancesCond",i,".csv"),sep=",",row.names=FALSE,col.names=TRUE)
  }
  #Export graph
  #try(png(filename=paste0("participant graphs/condition",j,"/",rownames(coordinatesByCondition[[j]])[i],".png")))
  #try(scatterplot(x=as.numeric(tempCoordinates[,1]), y=as.numeric(tempCoordinates[,2]), cex=1, cex.sub=4, smoother=FALSE, xlim=c(0,spam$xwindowsize[i]), ylim=c(0,spam$ywindowsize[i]), id.n=numberOfLabels, labels=labelNames, reg.line=FALSE, boxplots=FALSE))
  #dev.off()
  
  
}

##MDS#####################################################################################################################################################################################

#No conditions
if (numberOfConditions==0){
  stress<-matrix(ncol=1, nrow=length(dimensions), dimnames=list(dimensions, "stress"))
  
  #Calculate and save MDS, stress, and MDS graphs for all dimensions
  k<-0
  for (i in dimensions){
    k<-k+1
    tempMDS<-smacofSym(averageDistances, ndim=i, type="interval")
    assign(paste0("dim",i,"MDS"), tempMDS)
    write.csv(tempMDS$conf, paste0("MDS solutions/","dim",i,".csv"))
    stress[k,1]<-tempMDS$stress
    if (i==3){
      jpeg(filename = paste0("MDS graphs/","dim",i,".png"), quality=100, width=1200, height=1200)
      tempMDS2d<-scatterplot3d(x=tempMDS$conf[,1], y=tempMDS$conf[,2], z=tempMDS$conf[,3], type = "h", cex.axis=2, cex.lab=2)$xyz.convert(tempMDS$conf[,1],tempMDS$conf[,2],tempMDS$conf[,3])
      text(tempMDS2d, labels = rownames(tempMDS$conf), cex=1.5, pos=3)
      dev.off()
    } else if (i==2){
      jpeg(filename = paste0("MDS graphs/","dim",i,".png"), quality=100, width=1200, height=1200)
      print(ggplot(as.data.frame(tempMDS$conf), aes(x=D1, y=D2)) +
              geom_point(shape=1) +
              geom_text_repel(aes(x=D1, y=D2, label=rownames(tempCoordinates))))
      dev.off()
    }
  }
  write.table(stress, "stress.txt")
} else {
  
  #With conditions
  for (i in 1:numberOfConditions){
    tempStress<-matrix(ncol=1, nrow=length(dimensions), dimnames=list(dimensions, "stress"))
    
    #Calculate and save MDS, stress, and MDS graphs for all dimensions in all conditions
    k<-0
    for (j in dimensions){
      k<-k+1
      tempMDS<-smacofSym(averageDistancesByCondition[[i]], ndim=j, type="interval")
      assign(paste0("cond",i,"dim",j,"MDS"), tempMDS)
      write.csv(tempMDS$conf, paste0("MDS solutions/","cond",i,"dim",j,".csv"))
      tempStress[k,1]<-tempMDS$stress
      if (j==3){
        jpeg(filename = paste0("MDS graphs/","cond",i,"dim",j,".png"), quality=100, width=1200, height=1200)
        tempMDS2d<-scatterplot3d(x=tempMDS$conf[,1], y=tempMDS$conf[,2], z=tempMDS$conf[,3], type = "h", cex.axis=2, cex.lab=2)$xyz.convert(tempMDS$conf[,1],tempMDS$conf[,2],tempMDS$conf[,3])
        text(tempMDS2d, labels = rownames(tempMDS$conf), cex=1.5, pos=3)
        dev.off()
      } else if (j==2){
        jpeg(filename = paste0("MDS graphs/","cond",i,"dim",j,".png"), quality=100, width=1200, height=1200)
        #scatterplot(x=tempMDS$conf[,1], y=tempMDS$conf[,2], id.n=numberOfLabels, labels=rownames(tempMDS$conf), smoother=FALSE, cex.lab=2, id.cex=2, cex.axis=2, cex=1.5, reg.line=FALSE, boxplots=FALSE)
        
        print(ggplot(as.data.frame(tempMDS$conf), aes(x=D1, y=D2)) +
                geom_point(shape=1) +
                geom_text_repel(aes(x=D1, y=D2, label=rownames(tempCoordinates))))
        
        dev.off()
      }
    }
    assign(paste0("stress",i), tempStress)
    write.table(tempStress, paste0("stress_cond",i,".txt"))
  }  
}