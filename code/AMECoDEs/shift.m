function [p,fit,FES,counter] = shift(p,fit,A,popsize,counter,FES,q,n,problem,Xmax,Xmin)
 mp=mean(p);
 dp=1/popsize*sqrt(sum(pdist2(p,mp)));  
 [~,bestid]=min(fit);
if dp<0.001
    for i=1:popsize
        if counter(i)>=n && i~=bestid;           
            nj=randi([1,n],1);
            Rnd=randperm(n);
            j=Rnd(1:nj);
            p(i,j)=Xmin(problem)+(Xmax(problem)-Xmin(problem)).*rand(1,nj);
            fit(i)=cec14_func(p(i,:)',problem);
            counter(i)=0;
            FES=FES+1;
        end
    end
else
    U=[p;A];
    [~,sortid]=sort(fit,'ascend');
    pbestid=sortid(1:q*popsize);
    for i=1:popsize
        if counter(i)>n
            rnd=randi([1,floor(q*popsize)],1);
            pbest=pbestid(rnd);
            
            r1=randi([1,popsize],1);
            while  r1==i || r1==pbest
                r1=randi([1,popsize],1);
            end
            r2=randi([1,popsize],1);
            while  r1==i || r2==r1 || r2==pbest
                r2=randi([1,popsize],1);
            end          
            V=p(pbest,:);            
            Rnd=randperm(n);
            rj=randi([1,n],1);
            j=Rnd(1:rj);
            V(j)=p(pbest,j)+rand*(p(r1,j)-p(r2,j))+0.1*randn(1,rj);
             
            indexup=find(V>Xmax(problem));
            V(indexup)=max(Xmin(problem), 2*Xmax(problem)-V(indexup));
            indexlow=find(V<Xmin(problem));
            V(indexlow)=min(Xmax(problem), 2*Xmin(problem)-V(indexlow));
            fitV=cec14_func(V',problem);
            if fitV<=fit(i)
                p(i,:)=V;
                fit(i)=fitV;
                counter(i)=0;
            else
                counter(i)=counter(i)+1;
            end
            FES=FES+1;
        end
    end
end
end





