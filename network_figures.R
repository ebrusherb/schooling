library(network)

n=40
k=18

c=circle(n,k)
net=network(c)

x=cos(2*pi*(1:n)/n)
y=sin(2*pi*(1:n)/n)
coord=cbind(x,y)

plot(net,coord=coord,usearrows=FALSE,col=rgb(113,117,239))