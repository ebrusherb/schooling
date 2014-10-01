library(igraph)
library(R.matlab)

lap<-function(G){
	n=length(V(G))
	M=array(0,c(n,n))
	for(i in 1:n){
		M[i,neighbors(G,i,mode="in")]=-1/degree(G,i,mode="in")
	}
	diag(M)<-1+diag(M)
	M
}

id<-function(n){
	m=array(0,c(n,n))
	diag(m)=1
	return(m)
}

N=50

its=1000

vals=array(,its)
vals2=array(,its)
sums=array(,its)

j=1
while(j<=its){

positions<-array(,c(N,2))
positions[,1]=runif(N,0,1)
positions[,2]=runif(N,0,1)
d=as.matrix(dist(positions))

strategy=sample(1:(N-1),N,replace=TRUE)

M<-array(0,c(N,N))
for(i in 1:N){
	neighbors=order(d[i,])[2:(strategy[i]+1)]
	M[i,neighbors]=1/strategy[i]
}
L=id(N)-M

# G=erdos.renyi.game(N,p=.4,type="gnp",directed=TRUE)
# L=lap(G)
i=which.max(abs(Re(eigen(t(L))$vectors[,N])))
x=Re(eigen(L)$values)[N-1]
x=round(x,5)
y=array(0,N)
for(k in 1:N){
y[k]=Re(eigen(L[-k,-k])$values)[N-1]
}
sums[j]=sum(y)
y=round(y,5)
z=x-y[i]
vals[j]=x
vals2[j]=y[i]
j=j+1
print(c(x,z))
if(z<0){
	print("Negative")
	j=its+1
	}
# if(z==0){
	# print("interesting")
	# j=its+1
# }
}

# par(oma=c(0,4,0,0))
plot(vals,vals2,bty="n",cex.lab=1.5,xlab=expression(bar(lambda)),ylab=expression(lambda^i),mgp=c(2.5,1,0))
abline(0,1)

claim<-function(L,v){
	lambda=eigen(L)$values
	us=eigen(L)$vectors
	i=which.max(abs(v%*%us/norm(v,type='2')))
	lambda=lambda[-i]
	gap=round(min(diff(rev(lambda))),5)
	p=perp(v)
	lambdatilde=eigen(p%*%L%*%t(p))$values
	diffs=abs(lambda-lambdatilde)
}

# H=readMat('H2vals.mat')
# Hvals=H$H2vals
# Hmod=H$H2valsmod
# Hmod=apply(Hmod,2,mean)
# Hmod=Hmod[-c(140,439)]
# Hvals=Hvals[-c(140,439)]

# quartz()
# pdf(file='deviations.pdf')

# par(oma=c(0,2,0,0))
# plot(Hvals,Hmod,xlab="Longterm deviation from consensus",ylab="Longterm deviation from signal",bty='n',cex.lab=1.5)

# graphics.off()