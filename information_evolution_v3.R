library(igraph)
library(gplots)
library(fields)
library(network)

#####FUNCTIONS 
dynamics<-function(x,G){
	n<-length(x)
	new<-x
	for(i in 1:n){
		new[i]=x[i]+sum((x[neighbors(G,i,mode="in")]-x[i])/degree(G,i,mode="in"))
	}
	new
}

lap<-function(G){
	n=length(V(G))
	M=array(0,c(n,n))
	for(i in 1:n){
		M[i,neighbors(G,i,mode="in")]=-1/degree(G,i,mode="in")
	}
	diag(M)<-1+diag(M)
	M
}

Q<-function(n){
	a=rep(1,n)
	a=a/norm(a,type="2")
	m=array(0,c(n-1,n))
	m[1,1:2]=c(1/sqrt(2),-1/sqrt(2))
	for(j in 2:(n-1)){
		c=runif(n,-1,1)
		c=c-sum(c)/n
		for(k in 1:(j-1)){
			c<-c-(c%*%m[k,])*m[k,]
		}
		c<-c/norm(c,type="2")
		m[j,]=c
	}
	m
}

perp<-function(v){
	n=length(v)
	a=v/norm(v,type="2")
	m=array(0,c(n-1,n))
	c=runif(n,-1,1)
	c=c-(c%*%a)*a
	c=c/norm(c,type="2")
	m[1,]=c
	for(j in 2:(n-1)){
		c=runif(n,-1,1)
		c=c-(c%*%a)*a
		for(k in 1:(j-1)){
			c<-c-(c%*%m[k,])*m[k,]
		}
		c<-c/norm(c,type="2")
		m[j,]=c
	}
	m
}

H2<-function(G){
	vals<-eigen(lap(G))$values
	vals<-Re(vals)
	vals<-vals[-length(vals)]
	h=sqrt(1/2*sum(1/vals))
	h
}

nextgen<-function(strategy,fitness){
	l<-floor(length(fitness)/2)
	replace=order(fitness)[(length(fitness)-l+1):length(fitness)]
	v=strategy
	v[replace]=strategy[order(fitness)[1:l]]
	v
}

id<-function(n){
	m=array(0,c(n,n))
	diag(m)=1
	return(m)
}



################# evolution
N<-50

rainb<-rainbow(N)

positions<-array(,c(N,2))
positions[,1]=runif(N,0,1)
positions[,2]=runif(N,0,1)
d=as.matrix(dist(positions))

generations<-1

dist.gens<-array(,c(generations,N,N))
strategy<-array(,c(generations+1,N))
strategy[1,]<-sample(1:N,N,replace=TRUE)
H2vals<-array(,generations)

for(k in 1:generations){

M<-array(0,c(N,N))
for(i in 1:N){
	neighbors=order(d[i,])[2:(strategy[k,i]+1)]
	M[i,neighbors]=1
}
G<-graph.adjacency(t(M),mode="directed")
H2vals[k]=H2(G)

signal=runif(1,0,1)
runs=100

	for(ind in 1:N){
	chosen=ind
	x=array(,c(N,runs))

	x[,1]=runif(N,0,1)
	x[chosen,1]=signal


	for(i in 2:runs){
		x[,i]=dynamics(x[,i-1],G)
		x[chosen,i]=signal
		# x=x+rnorm(N,0,.01)
		}

	# plot(x[1,],type="l",ylim=c(0,1),col=rainb[1])
	# for(i in 2:N){lines(x[i,],col=rainb[i])}

	for(i in 1:N){
		dist.gens[k,ind,i]=sqrt(sum((x[i,]-x[chosen,])^2))
		}

	mean.dists<-array(,N)
	for(i in 1:N){
		mean.dists[i]=mean(dist.gens[k,-i,i])
	}
	strategy[k+1,]=nextgen(strategy[k,],mean.dists)
	
	}
}

######## analytical speed of learning
convergence<-function(D,strategy){
	N=length(strategy)
	M<-array(0,c(N,N))
	for(i in 1:N){
		neighbors=order(D[i,])[2:(strategy[i]+1)]
		M[i,neighbors]=1
	}
	# weighted<-M
	# for(i in 1:N){
		# weighted[,i]=M[,i]/sum(M[,i])
	# }
	G=graph.adjacency(t(M),mode='directed')
	L=lap(G)
	delta=array(0,c(N,N))
	for(i in 1:N){
		path_dists=shortest.paths(G,v=i,mode='out')
		not_children=which(path_dists==Inf)
		children=setdiff(which(path_dists!=Inf),i)
		delta[not_children,i]=0
		if(length(children)>0){
			L2=L[children,children]
			vals=Re(eigen((L2))$values)
			vecs=Re(eigen((L2))$vectors)
			delta[children,i]=min(vals)*abs(vecs[,which.min(vals)])
		}
	}
	means=apply(delta,1,mean)
	return(means)
}

convergence2<-function(M){
	N=dim(M)[2]
	G=graph.adjacency(t(M),mode='directed')
	L=lap(G)
	delta=array(0,c(N,N))
	for(i in 1:N){
		path_dists=shortest.paths(G,v=i,mode='out')
		not_children=which(path_dists==Inf)
		children=setdiff(which(path_dists!=Inf),i)
		delta[not_children,i]=0
		if(length(children)>0){
			L2=L[children,children]
			vals=Re(eigen((L2))$values)
			vecs=Re(eigen((L2))$vectors)
			delta[children,i]=min(vals)*abs(vecs[,which.min(vals)])
		}
	}
	means=apply(delta,1,mean)
	return(means)
}

convergence3<-function(M){
	N=dim(M)[2]
	G=graph.adjacency(t(M),mode='directed')
	L=lap(G)
	delta=array(0,c(N,N))
	for(i in 1:N){
		path_dists=shortest.paths(G,v=i,mode='out')
		not_children=which(path_dists==Inf)
		children=setdiff(which(path_dists!=Inf),i)
		delta[not_children,i]=0
		if(length(children)>0){
			L2=L[children,children]
			vals=Re(eigen((L2))$values)
			vecs=Re(eigen((L2))$vectors)
			delta[children,i]=(exp(-min(vals)))*abs(vecs[,which.min(vals)])
		}
	}
	means=apply(delta,1,mean)
	return(means)
}

########ESS for random networks

N=50

Native=seq(1,N-1,by=5)
Invader=seq(1,N-1,by=5)
count1=length(Native)
count2=length(Invader)

advantage=array(,c(count1,count2))
colnames(advantage)=paste("Invader",Invader)
rownames(advantage)=paste("Native",Native)
for(i in 1:count1){
	for(j in 1:count2){
		play1=Native[i]
		play2=Invader[j]
		
		strategy=array(play1,N)
		strategy[1]=play2
		
		its=10
	
		rates=array(,c(its,N))
		
		for(k in 1:its){
		positions<-array(,c(N,2))
		positions[,1]=runif(N,0,1)
		positions[,2]=runif(N,0,1)
		d=as.matrix(dist(positions))	
		rates[k,]=convergence(d,strategy)
		}
		native_avg=mean(rates[,2:N])
		invader_avg=mean(rates[,1])
		advantage[i,j]=invader_avg-native_avg
	}
}

image.plot(native,invader,advantage,xlab="Native",ylab="Invader")

# rainb=rainbow(count1)

# plot(invader,-ESS[1,],type="l",ylim=range(-ESS),xlab="Invader Strategy",ylab="Benefit to Invader",col=rainb[1])
# for(i in 2:count1){
	# lines(-ESS[i,],col=rainb[i])
# }


#### find an ESS 

findESS<-function(U){
	n=dim(U)[2]
	yes=array(0,n)
	for(i in 1:n){
		if(U[i,i]>max(U[-i,i])){yes[i]=1}
		for(j in 1:n){
			if(U[i,i]==U[j,i]&&U[i,j]>U[j,j]){yes[i]=1}
		}
	}
	return(yes)
}

cost<-function(j){
	j
}


costmat=array(,c(N,N))
for(i in 1:N){
	costmat[,i]=cost(i)
}

Utility=advantage-.5*costmat
findESS(Utility)


#############################
#### ESS for circular networks 
circle<-function(n,k){
	m=array(0,c(n,n))
	l=ceiling(k/2)
	indices=rep(1:n,3)
	if(k==1){
		for(i in 1:(n-1)){
			m[i,i+1]=1
		}
		m[n,1]=1
	} else
	for(i in 1:n){
		m[i,indices[(i+1:l)+n]]=1
		m[i,indices[(i-1:(k-l))+n]]=1
	}
	return(m)
}

numcount=10
count1=10
count2=10

advantage=array(,c(numcount,count1,count2))

N=floor(seq(20,150,length.out=numcount))

for(k in 1:numcount){

Native=floor(seq(2,N[k]-1,length.out=count1))
Invader=floor(seq(2,N[k]-1,length.out=count2))


indices=rep(1:N[k],3)

cheater=floor(N[k]/2)

for(i in 1:count1){
	for(j in 1:count2){
M=circle(N[k],Native[i])
a=Invader[j]
l=ceiling(a/2)
M[cheater,]=0
if(a>1){
M[cheater,indices[(cheater+1:l)+N[k]]]=1
M[cheater,indices[(cheater-1:(a-l))+N[k]]]=1
}
if(a==1){
	M[cheater,cheater+1]=1
}
v=convergence3(M)
advantage[k,i,j]=(-v[cheater]+mean(v[-cheater]))/mean(v[-cheater])
	}
}
image.plot(Native,Invader,advantage[k,,])
}

setwd('/Users/eleanorbrush/Desktop')

######plot the advantage the invader has over natives for various network sizes
quartz()
pdf(file='invader_advantage.pdf')
# postscript("invader_advantage.eps",bg="white",horizontal=FALSE)

# par(oma=c(0,0,0,0))
layout(matrix(1:2,nrow=2,byrow=TRUE))
for(i in c(1,6)){
	Native=floor(seq(2,N[i]-1,length.out=count1))
	Invader=floor(seq(2,N[i]-1,length.out=count2))
	image.plot(Native,Invader,advantage[i,,],xlab="Native Strategy",ylab="Invader Strategy",main=paste("N = ",N[i],sep=""),col=tim.colors(),cex.lab=1.5,cex.main=1.5,cex.axis=1,mgp=c(2.5,1,0))
	abline(0,1)
}
# par(oma=c(0,0,0,4))# reset margin to be much smaller.
# image.plot( legend.only=TRUE, zlim=range(advantage[c(1,5,10),,]))

graphics.off()

#### plot ESS strategy as a function of network size


ESS=array(,numcount)
for(i in 1:numcount){
	Native=floor(seq(2,N[i]-1,length.out=count1))
	Invader=floor(seq(2,N[i]-1,length.out=count2))
	ESS[i]=Native[min(which(Invader[apply(advantage[i,,],1,which.max)]-Native<=0))]
}

quartz()
pdf(file='ESS_strategy.pdf')

layout(1)
par(oma=c(0,2,0,0))
plot(N,ESS,xlab="Number of Nodes",ylab="Optimal Strategy",cex=.9,pch=1,cex.lab=1.5,bty="n")
f<-lm(ESS~N)
y<-predict(f,newdata=data.frame(x=N))
lines(N,y)

graphics.off()

###### plot convergence scores for one network 
n=20	
native=4
invader=18

indices=rep(1:n,3)

cheater=floor(n/2)

M=circle(n,native)
a=invader
l=ceiling(a/2)
M[cheater,]=0
if(a>1){
M[cheater,indices[(cheater+1:l)+n]]=1
M[cheater,indices[(cheater-1:(a-l))+n]]=1
}
if(a==1){
	M[cheater,cheater+1]=1
}


quartz()

pdf(file='convergence_scores_one_network.pdf')
par(oma=c(0,4,0,0))
v=convergence3(M)
plot((1:n)-cheater,v,xlab="Distance from Invader",ylab="Convergence Score",main=paste("N = ",n,", native strategy = ",native,", invader strategy = ",invader,sep=""),cex.lab=1.5,mgp=c(2.5,1,0),cex.main=1.5,bty="n",axes=FALSE)
axis(1,at=seq(1-cheater,n-cheater,by=2))
axis(2,at=round(seq(min(v),max(v),length.out=4),4))
graphics.off()

#########network plots
n=7 

positions<-array(,c(n,2))
positions[,1]=runif(n,0,1)
positions[,2]=runif(n,0,1)
d=as.matrix(dist(positions))

strategy=sample(1:(n-1),n,replace=TRUE)
M<-array(0,c(n,n))
for(i in 1:n){
	neighbors=order(d[i,])[2:(strategy[i]+1)]
	M[i,neighbors]=1
}
net=as.network(M)
plot(net,vertex.cex=3,vertex.col='#3962F7',label=strategy,label.pos=5,coord=positions)

n=20	
native=4
invader=18

pis=seq(0,2*pi,length.out=n+1)[-(n+1)]
xcoords=cos(pis);ycoords=sin(pis)
coords=cbind(xcoords,ycoords)

indices=rep(1:n,3)

cheater=1

M=circle(n,native)
a=invader
l=ceiling(a/2)
M[cheater,]=0
if(a>1){
M[cheater,indices[(cheater+1:l)+n]]=1
M[cheater,indices[(cheater-1:(a-l))+n]]=1
}
if(a==1){
	M[cheater,cheater+1]=1
}

net=as.network(M)
plot(net,vertex.cex=2,vertex.col='#3962F7',coord=coords)
