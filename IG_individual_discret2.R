# source("IG.R") - changes in line 47, 56, 53 to to modify output file names

rm(list=ls())

source("/Users/Raghu/Documents/asthma/modified_code/individual_dir_info.R")
source(libfile)
res.dir="/Users/Raghu/Documents/asthma/modified_results/discret2/cluster%"

## get parameters ##

# read clinical data #
clinic.file=file.path(data.dir, "original", "dataclusterkm_IGsort_samplescl6var10_NOalvmacro5.txt")
clinic.data=read.delim(clinic.file)
#data=apply(clinic.data[, 4:ncol(clinic.data)], 2, as.numeric)
data=as.matrix(clinic.data[, 4:ncol(clinic.data)])
col.data=colnames(data) #names of columns in data matrix
overallres.dir = res.dir #main results folder with "%" to be gsubed for current cluster analysis
# get phenotype data #
cls=as.numeric(clinic.data[, "Cluster"]) #get all the numbers in the cluster column
originalphenotypes=cls #the column of orignial cluster determinations for each patient
allclusters = unique(originalphenotypes) #unique clusters in initial data = c(1,2,3,4,5,6)

for (cluster.var in allclusters) { #does each cluster analysis one at a time
  phenotypes = originalphenotypes #original clusters that is going to be modified 
  cat(cluster.var,"\n") #output which cluster is being worked on
  change = which(phenotypes > cluster.var) 
  phenotypes[change] = 0
  change = which(phenotypes < cluster.var) #anything that is not equal to the cluster we are interested in is going to be 0
  phenotypes[change] = 0
  change = which(phenotypes == cluster.var) #
  phenotypes[change] = cluster.var #change this equal to cluster.var when doing second case where discret breaks into 6 clusters/ change to be equal to 1 in first case of discret2
  clusterlabels=unique(originalphenotypes) #change to unique(originalphenotypes) for discret6 #change to unique(phenotypes) for discret2
  n.class.pheno=length(clusterlabels) #counts the number of unique cluster numbers (discret6 = 6 & discret2 = 2)

    infogain.data=discret.data=su.data=NULL #set all of these objects to NULL
    for (i.var in 1:ncol(data)) {
    # i.var=1

      #cat(i.var, "\n") # write i.var to file which is the object that is being iterated....reason we see 1-120 in output

      data1=NULL
      cl.data1=length(unique(data[,i.var])) #number of unique elements in a given column
      if (cl.data1>n.class.pheno) { #if the number of unique elements is more than the 6 clusters we have set up
        data1=discret.equalwidth(data[,i.var], n.class.pheno) #clusters the varibles we see into n.class.pheno(number of clusters) -> bins
      } else {
        data1=data[,i.var] #if not greater keep the column the same since it has <= the number of clusters we want
      }

      infogain.data=c(infogain.data, infogain(phenotypes, data1)) #adds the infogain value into matrix for each variable columns
      su.data=c(su.data, infogain.normalized(phenotypes, data1))  #have the normalized info gain value as well
      discret.data=rbind(discret.data, data1) #converts linear matrix into each has an indivudal column
    }
  
    # at end of loop discret.data has each column represent all the variables for one patient
    res.dir=overallres.dir
    res.dir = gsub("%",cluster.var,res.dir) #change results to put into c
    idx.varig=1:ncol(data) #matrix of the sequence from 1 to the number of columns
    idx.varig.reord=idx.varig[order(infogain.data[idx.varig], decreasing=T)] #sorts the varibles from largest to smallest infogain #the largest are the variables that best edict cluster labels
    infogain.data.all=cbind(col.data[idx.varig.reord], infogain.data[idx.varig.reord]) # has the variables and infogain values from largest to smallest (Best to worst variables for estimating clusters)
    colnames(infogain.data.all)=c("Variable", "InfoGain") #set column headers
    filenameprefix = "Infogain_"
    filenamespecific = gsub("_",cluster.var,filenameprefix)
    infogain.file=gsub("datacluster", filenamespecific, basename(clinic.file)) #change begining of name of infile into InfogainEW
    infogain.file.fp=file.path(res.dir, infogain.file) #put results file in results directory
    write.table(infogain.data.all, infogain.file.fp, row.names=F,col.names=T, quote=F, sep="\t") #generate table in file
  
    #do the same thing for the normalized inforgain data 
    idx.varsu=1:ncol(data)
    idx.varsu.reord=idx.varsu[order(su.data[idx.varsu], decreasing=T)]
    su.data.all=cbind(col.data[idx.varsu.reord], su.data[idx.varsu.reord])
    colnames(su.data.all)=c("Variable", "SU")
    filenameprefix = "SuEW_Individual_"
    filenamespecific = gsub("_",cluster.var,filenameprefix)
    su.file=gsub("datacluster", filenamespecific, basename(clinic.file))
    su.file.fp=file.path(res.dir, su.file)
    write.table(su.data.all, su.file.fp, row.names=F,col.names=T, quote=F, sep="\t")

    #make file for the matrix of varibales for each patient
    discret.data.all=cbind(col.data[idx.varig.reord], discret.data[idx.varig.reord, ])
    colnames(discret.data.all)=c("Variable", phenotypes)
    filenameprefix = "discretEW_Individual_"
    filenamespecific = gsub("_",cluster.var,filenameprefix)
    discret.data.file=gsub("datacluster", filenamespecific, basename(clinic.file))
    discret.data.file.fp=file.path(res.dir, discret.data.file)
    write.table(discret.data.all, discret.data.file.fp, row.names=F,col.names=T, quote=F, sep="\t")
}




