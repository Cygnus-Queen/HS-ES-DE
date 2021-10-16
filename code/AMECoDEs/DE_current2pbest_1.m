function [p2,fit2,Fm2,CRm2,counter] = DE_current2pbest_1(p,fit,Fm2,CRm2,A,counter,popsize,q,n,c,problem,Xmax,Xmin)
%DE/rand-to-pbest/1
[~,sortid]=sort(fit,'ascend');
T=ceil(q*popsize);
pbestid=sortid(1:T);
[F,CR] = randFCR(popsize, CRm2, 0.1, Fm2,  0.1);
U=[p;A];
for i=1:popsize
    Rindex=randi([1,T],1);
    pbest=pbestid(Rindex);
%     while pbestid==i
%         Rindex=randi([1,T],1);
%         pbest=pbestid(Rindex);
%     end
    r1=randi([1,popsize],1);
    while  r1==i 
        r1=randi([1,popsize],1);
    end
    
    r2=randi([1,popsize],1);
    while  r2==i || r2==r1 
        r2=randi([1,popsize],1);
    end    
    
    k=randi([1,popsize+size(A,1)],1);
    while  k==r2 || k==r1  || k==i
        k=randi([1,popsize+size(A,1)],1);
    end  
    V=p(i,:)+F(i)*(p(pbest,:)-p(i,:))+F(i)*(p(r2,:)-U(k,:));
    %crossover
    j_rand = floor(rand * n) + 1;
    t = rand(1, n) < CR(i);
    t(1, j_rand) = 1;
    index=find(t==1);
    CR(i)=length(index)/n;
    t_ = 1 - t;
    V = t .* V + t_ .* p(i, :);
    indexup=find(V>Xmax(problem));
    V(indexup)=max(Xmin(problem), 2*Xmax(problem)-V(indexup));   
    indexlow=find(V<Xmin(problem));
    V(indexlow)=min(Xmax(problem), 2*Xmin(problem)-V(indexlow)); 
    uSet(i,:)=V;
end
fitSet = cec14_func(uSet',problem);%评价新种群
SF=[];SCR=[];dfit=[];
p2=uSet;
fit2=fitSet;
for i=1:popsize
    if fitSet(i)<=fit(i)
%         fit1(i)=fitSet(i);
%         p1(i,:)=uSet(i,:);
        SF=[SF,F(i)]; SCR=[SCR,CR(i)];
        dfit=[dfit, abs(fitSet(i)-fit(i))];
        counter(i)=0;
    else
        counter(i)=counter(i)+1;
    end
end
if ~isempty(SF)
    w=dfit/sum(dfit);
    Fm2=(1-c)*Fm2+c*sum(w.*(SF.^2))/sum(w.*SF);
    CRm2=(1-c)*CRm2+c*sum(w.*SCR);  
else
    Fm2=(1-c)*Fm2+c*rand;
    CRm2=(1-c)*CRm2+c*rand;
end
end

