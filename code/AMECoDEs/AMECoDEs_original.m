% clear all
clear;
clc;
format short e
format compact
warning off
%disp('AMECoDEs')
problemset = [1:30];% The index of test functions
fbias = [1:30]*100; %The actual optimal objctive function value of each functions
n = 30; %dimension
Xmin = -100*ones(1,30); % lower boundary
Xmax = 100*ones(1,30); %up boundary
popsize =200; % population size
totaltime =30; % the total number of runs
maxFES=n*10000 * 1; % th1e maximum number of function evaluation
c=0.1;
q=0.1;
for fun_num=[1:30] % fun_num is the index of the functions
    problem=problemset(fun_num);
    filename = strcat(strcat('f',num2str(problem)),'_original.txt');
    fp = fopen(filename,'a+');
     filename2 = strcat(strcat('results','_2017x1.txt'));
    fp2 = fopen(filename2,'a+');
    outcome = [];
    for time = 1:totaltime
        xxx = 0;
        curbestfit = realmax();%%%%%%%%%%%%%%%%%%our parameters
        prebestfit = realmax();
        curdiv = 1.0;
        prediv = 1.0;
        flag = 0;
        count = 0;

        sign = 0;
        sign2 = 0;
        everbestfitall = realmax();
        everbestfitinrun = realmax();%%%%%%%%%%%%%%%%%%our parameters
        rand('seed', sum(100*clock));  % Random seed
        p= Xmin(problem)+(Xmax(problem)-Xmin(problem))*rand(popsize,n); % the initial population
        fit = cec14_func(p',problem); % evaluate the population
        FES =popsize; % the number of function evaluation
        A=[];
        g = 1;   % the generation count
        fbest(g,time)=min(fit)-fbias(problem);% the value of function error in each generation
        Fm1=0.5;CRm1=0.5;Fm2=0.5;CRm2=0.5;Fm3=0.5;CRm3=0.5; % control parameters of each strategies  
        everbestfitall = min(fit) - fbias(problem);
        %FES = NP;
        curdiv = 0.0;
        for x =1 : n
            midpoint(x) = median(p(:,x));
        end
        distobest = 1 : popsize;
        for x = 1: popsize
            distobest (x)= 0;
            for y = 1 : n
                distobest(x) = distobest(x) + abs((p(x,y) - midpoint(y))/(Xmax(1,y) - Xmin(1,y)));
            end
            distobest (x) = distobest (x) / n;
            curdiv = curdiv + distobest (x);
        end
        curdiv = curdiv / popsize;
        %disp(curdiv);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fit)-fbias(problem), min(fit)-fbias(problem), everbestfitall);                        
        
        counter=zeros(1,popsize);       
        while FES < maxFES                  
             [p1,fit1,Fm1,CRm1,counter]=DE_current2nbest_1(p,fit,Fm1,CRm1,A,counter, popsize,q,n,c,problem,Xmax,Xmin); % current-to-nbest
             [p2,fit2,Fm2,CRm2,counter]=DE_current2pbest_1(p,fit,Fm2,CRm2,A,counter,popsize,q,n,c,problem,Xmax,Xmin); % current-to-pbest                   
             FES=FES+2*popsize;         
             A0=[];
             
             %selection************************* 
             for i=1:popsize
                UU=[p1(i,:);p2(i,:);]; Ufit=[fit1(i);fit2(i)];
                [minfit,minid]=min(Ufit); 
                                                                
                if minfit<=fit(i);
                    A0=[A0;p(i,:)];
                    fit(i)=minfit;
                    p(i,:)=UU(minid,:);
                end
             end   
            
            %shift mechanism
            if FES<maxFES              
               [p,fit,FES,counter] = shift(p,fit,A,popsize,counter,FES,q,n,problem,Xmax,Xmin);                      
            end          
            %*****************************
            
            %update archive
            A=[A;A0];
            LA=size(A,1);
            if LA>popsize
                L=LA-popsize;
                Rnd=randperm(LA);
                A(Rnd(1:L),:)=[];
            end
            
            g=g+1;
            fd=min(fit)-fbias(problem);
            if fd<fbest(g-1,time)
                fbest(g,time)=min(fit)-fbias(problem);
            else
                fbest(g,time)=fbest(g-1,time);
            end
            if fd < everbestfitall
                everbestfitall = fd;
            end
            
            if (FES / maxFES) >= (xxx / 100)
                if xxx > 0
                    curdiv = 0.0;
                    for x =1 : n
                        midpoint(x) = median(p(:,x));
                    end
                    distobest = 1 : popsize;
                    for x = 1: popsize
                        distobest (x)= 0;
                        for y = 1 : n
                            distobest(x) = distobest(x) + abs((p(x,y) - midpoint(y))/(Xmax(y) - Xmin(y)));
                        end
                        distobest (x) = distobest (x) / n;
                        curdiv = curdiv + distobest (x);
                    end
                    curdiv = curdiv / popsize;
                    fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fit)-fbias(problem), min(fit)-fbias(problem), everbestfitall);
                end
                xxx = xxx + 1;
            end             
            
            
            
        end
         outcome = [outcome, everbestfitall];
    end
    %Fbest{problem}=fbest; %save the best result at each generation   
    %meanER = mean(DEoutcome)
    %stdER=std(DEoutcome)   
    fprintf(fp,'%e (%e): ', mean(outcome), std(outcome));
     fprintf(fp2,'%e (%e): ', mean(outcome), std(outcome));
    for x = 1 : 30
        fprintf(fp,'%s ', num2str(outcome(x)));
         fprintf(fp2,'%s ', num2str(outcome(x)));
    end
     fprintf(fp2,'\n');
end


