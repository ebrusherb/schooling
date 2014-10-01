library(igraph)

f<-function(k,i,m,t){
	x=0
	if(is.element(i,c(1:(t+1),m-t-1,m-1))){
		x=tan(pi/(m-2)*k*(i-1))+1
	}
	if(is.element(i,(t+2):(m-t-2))){
		x=tan(pi/(m-2)*k*(i-m/2))+1
	}
	x
}

lambda<-function(k,m,t){
	vec=matrix(1:t,nrow=1)
	vec=apply(vec,2,f,k=k,m=m,t=t)
	vec=cumprod(vec)
	x=1-1/(2*t)*sum(vec)
	return(x)
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


#####
####
m=21
t=2
G=graph.adjacency(circle(n=m,k=2*t))
L=lap(G)
Lhat=L
Lhat[1,]=0

V=array(0,c(m,m))
V[,1]=1
for(i in 1:((m-1)/2)){
	V[,i+1]=cos(2*pi*i/m*(0:(m-1)))
	V[,i+(m-1)/2+1]=sin(2*pi*i/m*(0:(m-1)))
}

Lambda=(L%*%V/V)[m,]

Vhat=eigen(Lhat)$vectors
Vhat=cbind(Vhat[,m],Vhat)
Vhat=Vhat[,-(m+1)]

Lambdahat=(Lhat%*%Vhat/Vhat)[m,]

M=V%*%diag(1-Lambda)

# A=array(,c(m,(m-1)/2))
# B=array(,c(m,(m-1)/2))
# for(i in 1:((m-1)/2)){
	# A[,i]=V[1+i,]+V[m-i+1,]
	# B[,i]=V[1+i,]-V[m-i+1,]
# }


# C=array(0,c(m,m))
# for(i in 1:m){
	# for(j in 1:m){
		# C[i,j]=sum(Vhat[,j]*V[,i])/(sum(V[,i]^2))
	# }
# }


# C=array(0,c(m,m))()
# C[,1]=c(1,rep(0,m-1))
# for(j in 1:((m-1)/2)){
	# C[((m-1)/2+1):(m-1)+1,j+1]=2*sin(2*pi*j*(1:((m-1)/2))/m)
	# C[1:((m-1)/2)+1,j+((m-1)/2)+1]=2*cos(2*pi*j*(1:((m-1)/2))/m)
# }

# Vhat2=V%*%(C)

A=array(,c(m,m))
A[1,]=-V[1,]
for(i in 2:m){
	for(j in 1:((m-1)/2+1)){
		A[i,j]=-1/(2*t)*sum(cos(2*pi/m*(j-1)*((i-1+1):(i-1+t))))-1/(2*t)*sum(cos(2*pi/m*(j-1)*((i-1-1):(i-1-t))))
	}
	for(j in ((m-1)/2+2):(m)){
		A[i,j]=-1/(2*t)*sum(sin(2*pi/m*(j-(m-1)/2-1)*((i-1+1):(i-1+t))))-1/(2*t)*sum(sin(2*pi/m*(j-(m-1)/2-1)*((i-1-1):(i-1-t))))
	}
}

Vnorm=V
Vnorm[,1]=V[,1]/sqrt(m)
Vnorm[,2:m]=V[,2:m]/sqrt(m/2)

M2=t(V)%*%A
M=solve(V)%*%A
C=eigen(M)$vectors
Vhat=V%*%C

F<-function(t,m){
	v=array(,(m-1)/2)
	for(i in 1:((m-1)/2)){
		v[i]=1/t*sum(cos(2*pi*i*(1:t)/m))
	}
	return(v)
}

VA=-matrix(1,ncol=1,nrow=floor((m-1)/2))%*%matrix(1,nrow=1,ncol=floor((m-1)/2))+matrix(1,ncol=1,nrow=floor((m-1)/2))%*%matrix(F(t,m),nrow=1)-m/2*diag(F(t,m))

C=eigen(M)$vectors
Vhat=V[,c(1,((m-1)/2+2):m)]
Vhat=cbind(Vhat,V[,2:((m-1)/2+1)]%*%t(C))


######
mult<-function(i,m){
	v=array(,m-1)
	v[1]=1
	mu=-2*cos(pi*i/m)
	# mu=runif(1)
	for(k in 2:(m-1)){
		c=p=array(,ceiling(k/2))
		for(j in 1:ceiling(k/2)){
 				c[j]=choose(k-j,j-1)*(-1)^(j-1)
 				p[j]=(k-2*j+1)}
		v[k]=sum(c*mu^p)
	}
	v=c(0,v)
	return(v)
}

V=matrix(1:(m-1),nrow=1)
V=apply(V,2,mult,m=m)

mult2<-function(i,m){
	v=array(,m-1)
	v[1]=1
	mu=-2*cos(pi*i/m)
	# mu=runif(1)
	for(k in 2:(m-1)){
		c=p=array(,floor(k/2)+1)
		for(r in 0:floor(k/2)){
 				c[r+1]=choose(k-1-r,r)*(-1)^(r)
 				p[r+1]=(k-1-2*r)}
		v[k]=sum(c*mu^p)
	}
	v=c(0,v)
	return(v)
}

U<-function(i,m){
	v=array(,m-1)
	v[1]=1
	mu=-2*cos(pi*i/m)
	for(k in 2:(m-1)){
		c=p=array(,floor(k/2)+1)
		for(r in 0:floor(k/2)){
 				c[r+1]=choose(k-1-r,r)*(-1)^(r)
 				p[r+1]=(k-1-2*r)}
		v[k]=sum(c*mu^p)
	}
	v=c(0,v)
	return(v)
}

vecs<-function(m){
	v=array(0,c(m,m))
	v[,m]=rep(1,m)
	for(i in 1:(m-1)){
		v[2:m,i]=sin(pi*(1:(m-1))*(m-i)/m)/sin(pi*(m-i)/m)
	}
	return(v)
}


