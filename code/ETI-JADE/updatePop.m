function pop2 = updatePop(n, pop_si, popold, s1_Ind, DM)

sizep=size(pop_si,1);
for i=1 : sizep
    d_max=randperm(n, DM(i));
    popold(s1_Ind(i),d_max)=pop_si(i,d_max);
end
pop2=popold;

 