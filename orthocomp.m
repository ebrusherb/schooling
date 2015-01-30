% function C = orthocomp(reduced)

N=size(reduced,1);
empty=N-size(reduced,2);
C=[zeros(N,empty) reduced];

for i=empty:-1:1
   vec=rand(N,1);
   for j=N:-1:(i+1)
      vec=vec-(reshape(vec,1,[])*C(:,j))*C(:,j); 
   end
   vec=vec/norm(vec);
   C(:,i)=vec;
end