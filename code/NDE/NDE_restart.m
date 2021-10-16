%%..............................%%
clc;
clear all;
tic;

format long;
format compact;

%'NDE'
%'CEC2014'

D=30;
times = 10;
NP_ini=10*D;
NP_end=5;
gm=10;
FES_max=D*10000*times;
c=0.1;
Fm=0.5;
CRm=0.5;

funcset=[1:30];
lu=[-100*ones(1,D);100*ones(1,D)];

for funcindex=[17:30]
    
    F_index=funcset(funcindex);
    filename = strcat(strcat('f',num2str(funcindex)),'_restart.txt');
    fp = fopen(filename,'a+');
    time=1;
    totaltime=30;
    
    outcome=[];
    while time<=totaltime
        
        NP=NP_ini;
        FES=0;
        %Initialize population
        P=repmat(lu(1,:),NP,1)+rand(NP,D).*repmat((lu(2,:)-lu(1,:)),NP,1);
        fitness=cec14_func(P',F_index)-100*F_index;
        FES=NP;
        reFES = NP;
        
        N_rsize=ones(1,NP);
        Numg=zeros(1,NP); % '0' denotes the guider is good and '1' denotes the guider is worst
        Nums=zeros(1,NP); % '0' denotes the individual is good and '1' denotes the individual is worst
        
        curdiv = 1;
        everbestfitinrun = realmax();
        curbestfit = fitness(1);
        prebestfit = realmax();
        sign = 0;
        everbestfitall = min(fitness);
        xxx = 0;
        curdiv = 0.0;
        for x =1 : D
            midpoint(x) = median(P(:,x));
        end
        distobest = 1 : NP;
        for x = 1: NP
            distobest (x)= 0;
            for y = 1 : D
                distobest(x) = distobest(x) + abs((P(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
            end
            distobest (x) = distobest (x) / D;
            curdiv = curdiv + distobest (x);
        end
        curdiv = curdiv / NP;
        everbestfitall = min(fitness);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fitness), min(fitness),everbestfitall);% %g %g %d %d              %,everbestfitall,everbestfitinrun,mycount,ourOpt);
        
        %% Main Loop
        while FES<=FES_max
            
            v=[];
            u=[];
            
            F=[];
            CR=[];
            %   F=randcauchy(1,NP,Fm,0.1);
            
            
            F=Fm+0.1*tan((pi*rand(1,NP)-0.5));
            
            CR=normrnd(CRm,0.1,1,NP);
            
            F=min(1,F);
            pos=[];
            pos=find(F<0);
            while ~isempty(pos)
                F(pos)=Fm+0.1*tan((pi*rand(1,length(pos))-0.5));
                F=min(1,F);
                pos=find(F<0);
            end
            CR=min(1,max(0,CR));
            %%%%%%%%%%
            prebestfit = curbestfit;
            curbestfit = min(fitness);
            %curbestchrom = popold(ind, :);
            if prebestfit < everbestfitinrun
                everbestfitinrun = prebestfit;
            end
            
            if prebestfit < everbestfitall
                everbestfitall = prebestfit;
            end
            if sign == 0
                if (curbestfit >= everbestfitinrun) ||  (curbestfit < everbestfitinrun && (everbestfitinrun - curbestfit) / everbestfitinrun < 1e-5)
                    mycount = mycount +1;
                    %fprintf('%e %e\n',curbestfit,everbestfitinrun);
                else
                    mycount = 0;
                end
                upcount =500;
                if abs(everbestfitinrun - everbestfitall) < 1e-10
                    upcount = 1000;
                end
                if mycount >= upcount
                    sign = 1;%表示下一轮重启
                    mycount = 0;
                end
            else
                sign = 0;%表示下一轮恢复正常
            end
            %%%%%%%%%%
            for i=1:NP
                % construct the neighborhood for each individual
                neighbindex=[];
                neighb=[];
                fitness_neighb=[];
                
                neighbindex=mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1;
                neighb=P(neighbindex,:);
                fitness_neighb=fitness(neighbindex);
                
                [fnbest,nindex]=min(fitness_neighb);
                nbest=neighb(nindex,:);
                
                indexset=[1:NP];
                indexset(i)=[];
                temp=floor(rand*(NP-1))+1;
                index1=indexset(temp);
                indexset(temp)=[];
                temp=floor(rand*(NP-2))+1;
                index2=indexset(temp);
                
                lindexset=[1:2*N_rsize(i)+1];
                lindexset(N_rsize(i)+1)=[];
                ltemp=floor(rand*2*N_rsize(i))+1;
                lindex1=lindexset(ltemp);
                lindexset(ltemp)=[];
                ltemp=floor(rand*(2*N_rsize(i)-1))+1;
                lindex2=lindexset(ltemp);
                
                fnworst=max(fitness_neighb);
                fnmean=mean((fnworst-fitness_neighb)/(fnworst-fnbest));
                
                % Mutation
                if rand<(1+exp(20*((fnworst-fitness(i))/(fnworst-fnbest)-fnmean)))^-1
                    v(i,:)= neighb(lindex1,:) +F(i)*(P(index1,:)- P(index2, :));
                else
                    v(i,:)=P(i,:)+F(i)*(nbest-P(i,:)) +F(i)*(neighb(lindex1,:)-neighb(lindex2,:)) +F(i)*(P(index1,:)- P(index2, :));
                end
                
                % Bound constration handling
                vioLow = find(v(i,:) < lu(1, :));
                if ~isempty(vioLow)
                    v(i, vioLow) = 2 .* lu(1, vioLow) - v(i, vioLow);
                    vioLowUpper = find(v(i, vioLow) > lu(2, vioLow));
                    if ~isempty(vioLowUpper)
                        v(i, vioLow(vioLowUpper)) = lu(2, vioLow(vioLowUpper));
                    end
                end
                
                vioUpper = find(v(i,:)> lu(2, :));
                if ~isempty(vioUpper)
                    v(i, vioUpper) = 2 .* lu(2, vioUpper) - v(i, vioUpper);
                    vioUpperLow = find(v(i, vioUpper) < lu(1, vioUpper));
                    if ~isempty(vioUpperLow)
                        v(i, vioUpper(vioUpperLow)) = lu(1, vioUpper(vioUpperLow));
                    end
                end
                
                % Crossover
                j_rand = floor(rand * D) + 1;
                t = rand(1, D) < CR(i);
                t(1, j_rand) = 1;
                t_ = 1 - t;
                u(i, :) = t .* v(i,:) + t_ .* P(i, :);
            end
            
            fitness_u=[];
            fitness_u=cec14_func(u',F_index)-100*F_index;
            
            FES=FES+NP;
            reFES = reFES + NP;
            % Selection
            if sign == 0
                Tempfitness=[];
                TempP=[];
                Tempfitness=fitness;
                TempP=P;
                
                suc_F=[];
                suc_CR=[];
                increasement=[];
                for i=1:NP
                    if fitness_u(i)<Tempfitness(i)
                        increasement=[increasement Tempfitness(i)-fitness_u(i)];
                        TempP(i,:)=u(i,:);
                        Tempfitness(i)=fitness_u(i);
                        suc_F=[suc_F F(i)];
                        suc_CR=[suc_CR CR(i)];
                    else
                        if  fitness_u(i)==Tempfitness(i)
                            TempP(i,:)=u(i,:);
                            Tempfitness(i)=fitness_u(i);
                        end
                    end
                end
                
                % Updata control parameters
                if ~isempty(suc_F)
                    weight=increasement/sum(increasement);
                    meanCR=sum(weight.*suc_CR);
                    meanF=sum(weight.*suc_F.^2)/sum(weight.*suc_F);
                    Fm=(1-c)*Fm+c*meanF;
                    CRm=(1-c)*CRm+c*meanCR;
                end
                var=[];
                for i=1:NP
                    var=[var std(Tempfitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1))];
                    if min(Tempfitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1))<min(fitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1))
                        Numg(i)=0;
                        Nums(i)=0;
                    else
                        Numg(i)=Numg(i)+1;
                        if  mean(Tempfitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1))>=mean(fitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1)) && Tempfitness(i)~=min(Tempfitness)
                            Nums(i)=Nums(i)+1;
                        end
                    end
                end
                
                % avelliate the evolutionay dilemmas
                for i=1:NP
                    if Numg(i)==gm
                        if rand > Nums(i)/gm
                            N_rsize(i)=N_rsize(i)+1;
                            N_rsize(i)=min(N_rsize(i),floor(0.5*(NP-1)));
                        else
                            t = rand(1, D) <1-min((reFES/FES_max),(max(Tempfitness)-Tempfitness(i))/(max(Tempfitness)-min(Tempfitness)));
                            t_ = 1 - t;
                            if var(i)<mean(var)
                                I=lu(1,:)+rand(1,D).*(lu(2,:)-lu(1,:));
                                TempP(i,:) = t.*I + t_.*TempP(i,:);
                            else
                                [~,bindex]=min(Tempfitness(mod([NP+i-1-N_rsize(i):NP+i-1+N_rsize(i)],NP)+1));
                                gbest=TempP(mod(NP+i-1-N_rsize(i)+bindex-1,NP)+1,:);
                                TempP(i,:) = t.*gbest + t_.*TempP(i,:);
                            end
                            Tempfitness(i)=cec14_func(TempP(i,:)',F_index)-100*F_index;
                            FES=FES+1;
                            reFES = reFES + 1;
                        end
                        Numg(i)=0;
                        Nums(i)=0;
                    end
                end
                
                %Adaptively adjust the population size by the reduction method
                NPnew=round(NP_ini+(NP_end-NP_ini)*reFES/(0.1 * FES_max));
                if NPnew < NP_end
                    NPnew = NP_end;
                end
                if NPnew<NP
                    Divnum=NP-NPnew;
                    [~,sindex]=sort(Tempfitness,'descend');
                    TempP(sindex(1:Divnum),:)=[];
                    Tempfitness(sindex(1:Divnum))=[];
                    Numg(sindex(1:Divnum))=[];
                    Nums(sindex(1:Divnum))=[];
                    N_rsize(sindex(1:Divnum))=[];
                    var(sindex(1:Divnum))=[];
                    NP=NPnew;
                end
                %% Execute the neighborhood-based adaptive evolution mechanism
                %Track and record the performance and diversity of neighborhood of each individual
                
                
                P=[];
                fitness=[];
                P=TempP;
                fitness=Tempfitness;
            else
                %disp('duck');
                NP=NP_ini;
                
                %Initialize population
                P=repmat(lu(1,:),NP,1)+rand(NP,D).*repmat((lu(2,:)-lu(1,:)),NP,1);
                fitness=cec14_func(P',F_index)-100*F_index;
                
                
                N_rsize=ones(1,NP);
                Numg=zeros(1,NP); % '0' denotes the guider is good and '1' denotes the guider is worst
                Nums=zeros(1,NP); % '0' denotes the individual is good and '1' denotes the individual is worst
                sign = 0;
                reFES = 0;
            end          
            if (FES/(D * 10000 * times) >= (xxx / 100))
                if xxx > 0
                    curdiv = 0.0;
                    for x =1 : D
                        midpoint(x) = median(P(:,x));
                    end
                    distobest = 1 : NP;
                    for x = 1: NP
                        distobest (x)= 0;
                        for y = 1 : D
                            distobest(x) = distobest(x) + abs((P(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                        end
                        distobest (x) = distobest (x) / D;
                        curdiv = curdiv + distobest (x);
                    end
                    curdiv = curdiv / NP;
                    if min(fitness) < everbestfitall
                        everbestfitall = min(fitness);
                    end
                    fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fitness), min(fitness),everbestfitall);% %g %g %d %d              %,everbestfitall,everbestfitinrun,mycount,ourOpt);
                end
                xxx = xxx + 1;
            end
            %disp(NP);
            
        end
        outcome=[outcome everbestfitall];
        time=time+1;
    end
    %sort(outcome)
    %mean(outcome)
    %std(outcome)
    fprintf(fp,'%e (%e): \n',mean(outcome),std(outcome));
    for x = 1 : totaltime
        fprintf(fp,'%s ', num2str(outcome(x)));
    end
    fclose(fp);
end
toc;