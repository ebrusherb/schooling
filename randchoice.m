function choice=randchoice(options,probabilities)
c=cumsum(probabilities);
x=rand(1);
v=x-c;
v=sum(v>=0)+1;
choice=options(v);

