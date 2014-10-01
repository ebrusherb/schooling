figure
hold on
n=size(ESSseaten,2);

for i=1:n
    l=size(ESSseaten{i},2);
    plot(radvals(i)*ones(l,1),ESSseaten{i},'o');
end
