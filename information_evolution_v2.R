library(igraph)
library(gplots)

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
		M[neighbors(G,i,mode="in"),i]=-1/degree(G,i,mode="in")
	}
	diag(M)<-1
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

N<-50

rainb<-rainbow(N)

positions<-array(,c(N,2))
positions[,1]=runif(N,0,1)
positions[,2]=runif(N,0,1)
d=as.matrix(dist(positions))

its<-1

dist.its<-array(,c(its,N,N))
strategy<-array(,c(its+1,N))
strategy[1,]<-sample(1:N,N,replace=TRUE)
H2vals<-array(,its)

for(k in 1:its){

M<-array(0,c(N,N))
for(i in 1:N){
	neighbors=order(d[i,])[2:(strategy[k,i]+1)]
	M[neighbors,i]=1
}
G<-graph.adjacency(M,mode="directed")
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
		dist.its[k,ind,i]=sqrt(sum((x[i,]-x[chosen,])^2))
		}

	mean.dists<-array(,N)
	for(i in 1:N){
		mean.dists[i]=mean(dist.its[k,-i,i])
	}
	strategy[k+1,]=nextgen(strategy[k,],mean.dists)
	
	}
}

l=length(unique(strategy))
means=array(,l)
cis=array(,l)
for(i in 1:l){
	w=which(strategy[-1]==i)
	means[i]=mean(dist.its[-1,][w,])
	cis[i]=qt(0.975,length(w)*its)*apply(matrix(dist.its[-1,][w,],nrow=1),1,sd)/sqrt(length(w)*its)
}

plotCI(x=unique(strategy),y=means,uiw=cis)

convergence<-function(D,strategy){
	N=length(strategy)
	M<-array(0,c(N,N))
	for(i in 1:N){
		neighbors=order(D[i,])[2:(strategy[i]+1)]
		M[neighbors,i]=1
	}
	weighted<-M
	for(i in 1:N){
		weighted[,i]=M[,i]/sum(M[,i])
	}
	delta=array(0,c(N,N))
	for(i in 1:N){
		W=weighted[-i,-i]
		delta[-i,i]=Re(eigen(t(W))$values[1])*abs(Re(eigen(t(W))$vectors[,1]))
	}
	means=apply(delta,1,mean)
	return(means)
}