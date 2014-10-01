from __future__ import division
import numpy as np
import scipy
import scipy.spatial
import math, random

def dynamics(x,M):
	n = len(x)
	new = list()
	for i in range(n):
		neighbors = 0
		summed_diffs = 0
		for j in range(n):
			if M[j,i]>0:
				neighbors += 1
				summed_diffs += M[j,i]*(x[j]-x[i])
		to_add = x[i] + summed_diffs / neighbors
		new.append(to_add)
	return new
	
def lap(M):
	n = M.shape[1]
	m = np.zeros_like(M)
	m = m.astype(np.float)
	for i in range(n):
		neighbors = 0
		for j in range(n):
			if M[j,i]>0:
				neighbors +=1
		for j in range(n):
			if M[j,i]>0:
				m[j,i]=-1/neighbors
		m[i,i]=1
	return m
	
def H2(M):
	vals,vec = np.linalg.eig(lap(M))
	vals = vals.real
	vals = sorted(vals)[1:]
	vals = np.array(vals)
	h = math.sqrt(.5*sum(1/vals))
	return h
	
def nextgen(strategy,fitness):
	l = len(fitness)//2	
	to_return = list(strategy)
	temp = np.array(fitness).argsort()
	ranks = temp.argsort()
	ranks = list(ranks)
	w = 0
	for j in range(len(fitness)):
		if ranks[j]>=(len(fitness)-l):
			to_return[j]=strategy[ranks.index(w)]
			w += 1
	return to_return

def dist(N):
	positions = np.empty(shape=(N,2))
	for i in range(N):
		positions[i,0]=random.uniform(0,1)
		positions[i,1]=random.uniform(0,1)
	D = scipy.spatial.distance.cdist(positions,positions)
	return D

N = 20
its = 10
D = dist(N)
runs = 500

dists = np.empty(shape=(its,N,N))
strategy = np.empty(shape=(its+1,N))
for j in range(N):
	strategy[0,j] = random.randint(1,N-1)

H2vals = np.empty(shape=its)

for k in range(its):
	strategy2 = strategy[k,:]
	
	#construct communication matrix
	M = np.zeros(shape=(N,N))
	for i in range(N):
		dist_rank = list(D[:,i].argsort().argsort())
		v = range(N)
		v.remove(i)
		for s in v:
			if dist_rank[s] <= strategy2[i]:
				M[s,i]=1
	H2vals[k]=H2(M)
	
	signal = random.uniform(0,1)
	
	#update opinions according to dynamical system, where signaler receiver maintain opinion about signal
	for ind in range(N):
		chosen = ind
		x = np.empty(shape=(N,runs))
		for l in range(N):
			x[l,0]=random.uniform(0,1)
		x[chosen,0]=signal
		
		for i in range(1,runs):
			x[:,i]=dynamics(x[:,i-1],M)
			
		for i in range(N):
			dists[k,ind,i] = math.sqrt(sum((x[i,:]-x[chosen,:])**2))
			
	meandists = np.empty(N)
	v = list(dists[k,:,0])
	v = np.array(v[1:N])
	meandists[0]=np.mean(v)
	for i in range(1,N):
		v = list(dists[k,:,i])
		v = v[:i-1]+v[i:N]
		v = np.array(v)
		meandists[i]=np.mean(v)
			
	strategy[k+1,:]=nextgen(strategy[k,:],meandists)
	


	
