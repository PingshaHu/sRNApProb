library(edgeR,lib.loc='~/R/x86_64-unknown-linux-gnu-library/3.0')
#library(edgeR)
args<- commandArgs(TRUE)
# 7-28-2015 fixed null locus estimation bug by adding 1. 

### Read the table ###
infile=args[1] # infile="smRNA_mergecounts.txt"
samples=sort(unlist(strsplit(args[2],","))) # samples=C24,C24,C24,Ler,Ler,Ler,flc,flc,flc,fcl,fcl
#nlines=as.numeric(unlist(strsplit(system(paste("wc -l",infile,sep=" "),intern=T)," "))[1]) # intern to capture the return value
headline=readLines(infile,n=1)
headcols=unlist(strsplit(headline,"\\t"))
countcols=grep("count",headcols) # I need sample name in the Total, sample-rep#.count
x=read.delim(infile,header=TRUE)

### random sampling form locus and null 5000 samples ###  inputs, the count data per sample, directly from loci file. 
x.null <- x[x$Content=="null",]
x.smRNA <- x[x$Content=="smRNA",]
#N_sample=x.null[sample(nrow(x.null),5000,replace = FALSE),]  
N_sample=x.null[sample(nrow(x.null),500,replace = TRUE),]  
#S_sample=x.smRNA[sample(nrow(x.smRNA),5000,replace = FALSE),]
S_sample=x.smRNA[sample(nrow(x.smRNA),500,replace = TRUE),]

# get lib size
y <- DGEList(counts=x[,countcols]) # y <- DGEList(counts=x[,-c(1,2,3,4)],group=group)
y <- calcNormFactors(y,method="TMM") # y[1]] is/ the sample number list, y[2]] is the lib.size and normfactors
lib.size=y$samples$lib.size

# get the count data devideded by length
RPKM <-function(x){
  y=x[countcols]
  pos=unlist(strsplit(as.character(x[1]),"\\."))
  Length=as.numeric(pos[3])-as.numeric(pos[2])+1
  y.RPK=as.numeric(y)/(Length/1000)
  y.RPKM=y.RPK/(lib.size/1000000)
  return(y.RPKM)	
}	
x.RPKM=data.frame(t(apply(x,1,RPKM)))
colnames(x.RPKM)=headcols[countcols]
x.RPKM=round(x.RPKM*100,0)

N_sample.RPKM=data.frame(t(apply(N_sample,1,RPKM)))
colnames(N_sample.RPKM)=headcols[countcols]
N_sample.RPKM=round(N_sample.RPKM*100,0)

S_sample.RPKM=data.frame(t(apply(S_sample,1,RPKM)))
colnames(S_sample.RPKM)=headcols[countcols]
S_sample.RPKM=round(S_sample.RPKM*100,0)


### decide the experiment design ###
# using headline info to decide the experiment design
group <- factor(samples)
design <- model.matrix(~group)
grouplevel=levels(group)

### estimate phi and u 
# phi from edge
#library(edgeR)
phi_est <-function(x){
y <- DGEList(counts=x,group=group) # y <- DGEList(counts=x[,-c(1,2,3,4)],group=group)
y <- calcNormFactors(y,method="TMM") # y[1]] is the sample number list, y[2]] is the lib.size and normfactors
y <- estimateCommonDisp(y)
#y <- estimateGLMCommonDisp(y,design)
phi=y$common.dispersion
return(phi)
}



# phi for null and small RNA
#phi_n=phi_est(N_sample.RPKM)
phi_n=phi_est(N_sample.RPKM+1) # fixed null bugs by add 1. 
phi_s=phi_est(S_sample.RPKM)
r_n=1/phi_n
r_s=1/phi_s

# u_i for each group i by MLE
#Probcount <- function(k_i,r,u) {
#	P=gamma(k_i + r)/((gamma(r))*gamma(k_i+1))*(1/(1+u*(1/r)))^r*(u/(r+u))^k_i
#}
Probcount_1 <- function(k,r_s,u_s){
p=pnbinom(k,size=r_s,mu=u_s)
}
Probcount_0 <- function(k,r_n,u_n){
p=1-pnbinom(k,size=r_n,mu=u_n)
}

mu_est <- function(x,r){
# x is N_sample or S_sample
# r is 1/phi
k=p=u=NULL
for ( i in 1:length(grouplevel)) {
	patn=paste(grouplevel[i],".+counts",sep="")
	#x.group=x[,grep(patn,headcols)] # this is for not normalized data with whole x as input
	x.group=x[,grep(patn,colnames(x))] # this is for normalized RPKM data with x.RPKM as input
	k[[i]]=apply(x.group,1,mean)
	N=length(k[[i]])
	p[i]=sum(k[[i]]/N)/(r+sum(k[[i]]/N))
	u[i]=p[i]*r/(1-p[i])	
}	
	return(list(k=k,u=u))
}
S_est=mu_est(S_sample.RPKM,r_s)
N_est=mu_est(N_sample.RPKM,r_s)

Prob_est <- function(x,r_s,u_s,r_n,u_n){
# x is x.RPKM sample
# u_s from S_est$u and u_n from N_est$u
k=P_1=P_0=NULL
for ( i in 1:length(grouplevel)) {
	patn=paste(grouplevel[i],".+counts",sep="")
	x.group=x[,grep(patn,colnames(x))] # this is for normalized RPKM data with x.RPKM as input
	k[[i]]=apply(x.group,1,mean)
	P_1[[i]]=sapply(k[[i]],Probcount_1,r_s=r_s,u_s=u_s[i])	
	P_0[[i]]=sapply(k[[i]],Probcount_0,r_n=r_n,u_n=u_n[i])	
}	
	return(list(P_1=P_1,P_0=P_0,k=k))
}
Probest=Prob_est(x.RPKM,r_s=r_s,u_s=S_est$u,u_n=N_est$u,r_n=r_n)


# get Problength Probability
Problength <- function(x){
prob=NULL
for ( i in 1:length(grouplevel)) {
	patn=paste(grouplevel[i],".+_uniformdis",sep="")
	x.group=x[,grep(patn,headcols)]
	prob[[i]]=apply(x.group,1,mean)
  }	
	#return(list(problen=prob))
	return(prob)
}
# get the Problength_1 for small RNA locus
p=Problength(x)
Problength_1=lapply(p,function(x){1-x})
# got Problenth_0 for null locus
Problength_0=Problength(x)


### get the posterior Prob ###
# 1. Prob(data|I=1)=Probcount_1*Problength_1
# S_est$P_1[[i]]*S_est$Problength_1[[i]]

# 2. Prob(data|I=0)=Probcount_0*Problength_0
# S_est$P_0[[i]]*S_est$Problength_0[[i]]

# 3. the posterior Prob
#Prob(I=1|data)=p*Prob(data|I=1)/(p*Prob(data|I=1)+(1-p)*Prob(data|I=0))
#pProb=0.5*Probcount_1*Problength_1/(0.5*Probcount_1*Problength_1+0.5*Probcount_0*Problength_0)

pProb=NULL
for ( i in 1:length(grouplevel)){
	pProb[[i]]=0.5*Probest$P_1[[i]]*Problength_1[[i]]/(0.5*Probest$P_1[[i]]*Problength_1[[i]]+0.5*Probest$P_0[[i]]*Problength_0[[i]])
}	

#p1=Probest$P_1[[i]]
#plen1=Problength_1[[i]]  # problem
#p0=Probest$P_0[[i]]
#plen0=Problength_0[[i]]  # problem


posteriorProb=matrix(unlist(pProb),ncol=length(grouplevel))
colnames(posteriorProb)=paste(grouplevel,"posteriorProb",sep="_")
mydata=data.frame(x,posteriorProb)
write.table(mydata,"smRNA_posteriorProb.txt",sep="\t",quote=F)
# reformat smRNA_posteriprProb.txt as a new file with count and posteriorProb 
pProb_head=colnames(mydata)
tmp=grep("uniformdis",pProb_head)
mydata2=mydata[,-tmp]
newhead=colnames(mydata2)
write.table(mydata2,"sRNA_loci_wprob.txt",sep="\t",quote=F,row.names=F)

# P1 is the probability to be as small RNA locus. 
P1=matrix(unlist(Probest$P_1),ncol=length(grouplevel),byrow=F)
colnames(P1)=paste(grouplevel,"Probcount1",sep="_")
rownames(P1)=paste(x[,1],x[,2],sep=":")
write.table(P1,"smRNA_Probcount_1.txt",row.names=F,sep="\t",quote=F)

# P0 is the probability to be as null locus. 
P0=matrix(unlist(Probest$P_0),ncol=length(grouplevel),byrow=F)
colnames(P0)=paste(grouplevel,"Probcount0",sep="_")
rownames(P0)=paste(x[,1],x[,2],sep=":")
write.table(P0,"smRNA_Probcount_0.txt",row.names=F,sep="\t",quote=F)




