function [p1,fit1,Fm1,CRm1,counter] = DE_current2nbest_1__(p,fit,Fm1,CRm1,A,counter, popsize,q,n,c,problem,Xmax,Xmin)
%DE_current2nbest_1
T=ceil(q*popsize);
if mod(T,2)~=0
   T=T+1;
end
[F,CR]=randFCR(popsize, CRm1, 0.1, Fm1,  0.1);
%U=[p;A];
A=[];
U=[];
for i=1:popsize
    rn=mod([i:i+T/2-1],popsize)+1;
    ln=mod([i-T/2-1+popsize:i-2+popsize],popsize)+1;
    nid=[ln,i,rn];
    nfit=fit(nid);
    [~,bestnid]=min(nfit);
    nbest=nid(bestnid);      
    r1=randi([1,popsize],1);
    while  r1==i 
           r1=randi([1,popsize],1);
    end 
    r2=randi([1,popsize],1);
    while  r2==i || r2==r1
           r2=randi([1,popsize],1);
    end     
% % % % % % % % % %     k=randi([1,popsize+size(A,1)],1);
% % % % % % % % % %     while  k==r1 || k==i
% % % % % % % % % %         k=randi([1,popsize+size(A,1)],1);
% % % % % % % % % %     end    
    %v=p(i,:)+F(i)*(p(nbest,:)-p(i,:))+F(i)*(p(r1,:)-U(k,:));      
    v=p(i,:)+F(i)*(p(r1,:)-p(r2,:));   
    %crossover
    j_rand = floor(rand * n) + 1;
    t = rand(1, n) < CR(i);
    t(1, j_rand) = 1;
    index=find(t==1);
    CR(i)=length(index)/n;
    t_ = 1 - t;
    V = t .* v + t_ .* p(i, :);   
    indexup=find(V>Xmax(problem));
    V(indexup)=max(Xmin(problem), 2*Xmax(problem)-V(indexup));   
    indexlow=find(V<Xmin(problem));
    V(indexlow)=min(Xmax(problem), 2*Xmin(problem)-V(indexlow)); 
    uSet(i,:)=V;
end
fitSet = cec14_func(uSet',problem);%评价新种群
%选择
SF=[];SCR=[]; dfit=[];
p1=uSet;
fit1=fitSet;
for i=1:popsize
    counter(i)=0;
% % % % % % % % % %     if fitSet(i)<=fit(i)
% % % % % % % % % % %         fit1(i)=fitSet(i);
    p1(i,:)=uSet(i,:);
% % % % % % % % % %         SF=[SF,F(i)];
% % % % % % % % % %         SCR=[SCR,CR(i)];
% % % % % % % % % %         dfit=[dfit, abs(fitSet(i)-fit(i))];   
% % % % % % % % % %         counter(i)=0;
% % % % % % % % % %     else
% % % % % % % % % %         counter(i)=counter(i)+1;
% % % % % % % % % %     end
end
% % % % % % % % % % if ~isempty(SF)
% % % % % % % % % %     w=dfit/sum(dfit);
% % % % % % % % % %     Fm1=(1-c)*Fm1+c*sum(w.*(SF.^2))/sum(w.*SF);
% % % % % % % % % %     CRm1=(1-c)*CRm1+c*sum(w.*SCR);
% % % % % % % % % % else
% % % % % % % % % %     Fm1=(1-c)*Fm1+c*rand;
% % % % % % % % % %     CRm1=(1-c)*CRm1+c*rand;
% % % % % % % % % % end

end

