# source("FCBFig005_6clkmdata_NOalvmacro5.R")

#input-> file produced by IG.R where columns are patients and the rows are variables and the discretized clusters
#identifies redundant variables for each cluster type - modified from "FCBFig005_6clkmdata_NOalvmacro5.R"

rm(list=ls())

source("/Users/Raghu/Documents/asthma/modified_code/individual_dir_info.R")
source(libfile) #infogain_function.R
#source(libfile2)

for (discretization.type in c(2,6)){
discret.folder = gsub("%",discretization.type,"discret%")

  for (cluster.type in c(1:6)){

    ## read discrete data ##
    cluster.folder = gsub("%",cluster.type,"cluster%")
    inputfile = gsub("%",cluster.type,"C%_discret_indiviudal_d!.txt")
    inputfile = gsub("!",discretization.type,inputfile)
    #overall.inputfile = "discretEW%Individual%km_IGsort_samplescl6var10_NOalvmacro5.txt"
    #cluster.inputfile = gsub("%",cluster.type,overall.inputfile)


    discret.datafile.fp=file.path(res.dir,discret.folder,cluster.folder,inputfile) #text file of data 
    output.datafile.path = file.path(res.dir,discret.folder,cluster.folder)
    #discret.datafile="discretEWkm_6clpats_NOalvmacro5neword.txt"
    #discret.datafile.fp=file.path(kmeansclass.newdir, discret.datafile)


      discret.data=var.name=data=idxrm.varname=data=n.sample=n.var=cls=phenotypes=NULL
      su.ic.var=idx.ord=S.list=S.list.var=idx.Fj=Fj=idx.Fi=Fi=su1.ij=su1.ic=idx.Fi.rm=NULL
      var.reord=S.list.var=NULL

      #discret.data=readdata(discret.datafile.fp, cat="text")
      discret.data=as.matrix(read.delim(discret.datafile.fp)) #reads in data from file as matrix
      #file is results from discret.data; rows are variables ; columns are patients and the cluster
      #they have been assigned to; if less than 6 unique values of variable, it remains original 
      #cluster demarcation

      var.name=discret.data[,"Variable"]
      data=apply(discret.data[, -1],2,as.numeric)#converts the cluster numbers which are strings initially into numbers
      rownames(data)=var.name #sets rownames of data to be the variables names

      # remove "Phenotype_Group", and "subjclass_screen"   -> # why do we remove these??? 
      idxrm.varname=c(grep("Phenotype_Group", var.name), grep("subjclass_screen", var.name),
      		  grep("dxasthma", var.name), grep("Total_iCS", var.name),
    		  grep("Total_BetaAgonists", var.name))  #, grep("icort", var.name)
      if (any(idxrm.varname)) {
        data=data[-idxrm.varname, ]
      }
      n.sample=ncol(data)
      n.var=nrow(data)

      #removed X in column names so we can convert into numbers ex. X1.2 -> 1.2 which represents cluster for patient
      columns = colnames(discret.data)[-1]
      for (i in 1:length(columns)) {
        columns[i] = gsub("X","",columns[i])
      }
      
      cls=as.integer(columns) #converts the string of cluster number into integers
      phenotypes=cls

      su.ic=NULL # calculate infogain comparing each variable to the phenotypes columns values
      for (i.var in 1:n.var) {
      # i.var =1
        cat(i.var, "\n") #displays which variable it caclulates it for to the console
        su.ic=c(su.ic, infogain(phenotypes, data[i.var, ]))
      }

    ## order su.data in descending order ## why do we do this if the file is already ordered from best variable to worst based on infogain variable
     
    su.ic.var=rownames(data) #variable names 
    idx.ord=order(su.ic, decreasing=T) # indicies of order based on decreasing infogain value
    S.list=su.ic[idx.ord] #order infograin values
    S.list.var=su.ic.var[idx.ord] #order variables based on infogain
    names(S.list)=S.list.var

    S.list=S.list[S.list>0.05] #remove variables that have infogain value less .05  -> why? they don't matter anyway?
    S.list.var=names(S.list) 

    ##
    idx.Fj=1
    Fj=S.list.var[idx.Fj] # start with highest rated variable 
    var.reord=NULL
    while (!is.na(Fj)) { #iterate through all the variables in S.list.var until NA element is selected(at end of list)

     cat(Fj, "\n") #print out line number to console
     idx.Fi=idx.Fj+1 #increase index to look at next variable
     Fi=S.list.var[idx.Fi] #get next variable 
     idx.Fi.rm=NULL
     while (!is.na(Fi)) { #compare variable j with all other possible combinations of varibles i
      su1.ij=infogain(data[Fi, ], data[Fj, ]) #calculate the similarity between the two varibales Fi and Fj
      #su1.ic=S.list[idx.Fi]
      su1.ic=infogain(data[Fi, ], phenotypes) #compare that similarity to the phenotypes 
      if (su1.ij>=su1.ic) {  #if closer infogain match to variable than to cluster 
       idx.Fi.rm=c(idx.Fi.rm, idx.Fi) #adds the index of the variable that is redundant to the variable at hand 
      }
      idx.Fi=idx.Fi+1 # look a the next variable
      Fi=S.list.var[idx.Fi] # new variable to compare 
     } # while loop ends when compared to all other variables 
     if (!is.null(idx.Fi.rm)) { # if some redundant variable(s) was identified 
      var.reord=c(var.reord, paste("(", Fj, ")", sep=""), S.list.var[idx.Fi.rm]) 
      S.list.var=S.list.var[-idx.Fi.rm] #remove varible from list of total variables 
     } else {
      var.reord=c(var.reord, paste("(", Fj, ")", sep=""))
     }

     idx.Fj=idx.Fj+1
     Fj=S.list.var[idx.Fj]
    }

    # the infogain value for each variable when compared to the phenotypes cluster determination for all the variables that are not redundant 
    S.list.f=cbind(S.list.var, formatC(S.list[S.list.var]))
    colnames(S.list.f)=c("Variable", "Su")
    #S.list.ffile=gsub("discretE", "FCBFig005nonredundvar", discret.datafile.fp)
    S.list.filename = gsub("%",cluster.type,"C%_FCBnonredund_d!.txt")
    S.list.filename = gsub("!",discretization.type,S.list.filename)
    S.list.ffile = file.path(output.datafile.path,S.list.filename)
    write.table(S.list.f, S.list.ffile, col.names=T, row.names=F, quote=F, sep="\t")

    #writees the 
    #var.reord.file=gsub("discretE", "varreord.FCBFig005nonredund", discret.datafile.fp)
    var.reord.filename = gsub("%",cluster.type,"C%_varreord.FCBnonredund_d!.txt")
    var.reord.filename = gsub("!",discretization.type,var.reord.filename)
    var.reord.file = file.path(output.datafile.path,var.reord.filename)
    write.table(var.reord, var.reord.file, col.names=F, row.names=F, quote=F, sep="\t")
  }
}
