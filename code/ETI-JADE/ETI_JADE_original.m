%**************************************************************************************************
% Author: Wei Du
% Last Edited: Jan 25, 2016
% Email: duwei0203@gmail.com
% Reference: Differential Evolution With Event-Triggered Impulsive Control
%                  IEEE Transactions on Cybernetics, doi: 10.1109/TCYB.2015.2512942
%**************************************************************************************************

% The MATLAB source codes of JADE was provided by Dr. J. Zhang.
% References: J. Zhang and A. C. Sanderson, "JADE: adaptive differential evolution with optional external archive," IEEE Trans. Evolut. Comput., vol. 13, no. 5, pp. 945-958, 2009.

clc;
clear all;
tic;

format long;
format compact;
warning off

'ETI-JADE'

n=30;

popsize=100;

lu = [-100 * ones(1, n); 100 * ones(1, n)];

fhd=str2func('cec14_func');

problemSet = 1:30;

for problemIndex = [1:30]

    problem = problemSet(problemIndex);
    filename = strcat(strcat('f',num2str(problem)),'_original.txt');   
    fp = fopen(filename,'a+');    
    func_num = problem;

    % Record the best results
    outcome = []; 
    beishu = 1;
    MaxGen = n * 10000 * beishu / popsize;
    curdiv = 0.0;
    time = 1;
    % The total number of runs
    totalTime = 30;

    while time <= totalTime

        rand('state', sum(100 * clock));

        % Initialize the main population
        popold = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

        valParents = feval(fhd, popold', func_num) -  problemIndex * 100;
		%disp(valParents);
        valParents = valParents';
        
        gbest=min(valParents);

        c = 1/10;
        p = 0.05;

        CRm = 0.5;
        goodCR = [];
        Fm = 0.5;

        Afactor = 1;
        
        srpop=0;
        
        pbest_STG=zeros(popsize,1);
        
        LN=1;
        UN=popsize;
        num=LN;
        
        archive.NP = Afactor * popsize; % the maximum size of the archive
        archive.pop = zeros(0, n); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        %% the values and indices of the best solutions
        [valBest, indBest] = sort(valParents, 'ascend');
        gen = 1;
        FES = popsize;
        xxx = 0;
        curdiv = 0.0;   
        for x =1 : n
           midpoint(x) = median(popold(:,x));
        end
        distobest = 1 : popsize;
        for x = 1: popsize
            distobest (x)= 0;
            for y = 1 : n
                distobest(x) = distobest(x) + abs((popold(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
            end
            distobest (x) = distobest (x) / n;
            curdiv = curdiv + distobest (x);
        end                   
        curdiv = curdiv / popsize;
        everbestfitall = min(valParents);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(valParents), min(valParents),everbestfitall);% %g %g %d %d              %,everbestfitall,everbestfitinrun,mycount,ourOpt);
        ourOpt = 0; 
        everbestfitall = realmax();
        everbestfitinrun = realmax();
        noprogressgen = 0;
        curbestfit = valParents(1);
        prebestfit = realmax();	
        curbestchrom = popold(1, :);
        preind = 0;
        while FES < n * 10000 * beishu
            
            pop = popold; % the old population becomes the current population
            
            valParents_temp=valParents;

            if FES > 1 && ~isempty(goodCR) && sum(goodF) > 0 % If goodF and goodCR are empty, pause the update
                CRm = (1 - c) * CRm + c * mean(goodCR);
                Fm = (1 - c) * Fm + c * sum(goodF .^ 2) / sum(goodF); % Lehmer mean 
            end

            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F, CR] = randFCR(popsize, CRm, 0.1, Fm, 0.1);

            r0 = [1 : popsize];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);

            % Find the p-best solutions
            pNP = max(round(p * popsize), 2); % choose at least two best solutions
            randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(indBest(randindex), :); % randomly choose one of the top 100p% solutions

            prebestfit = curbestfit;
            prebestchrom = curbestchrom;
            %preind = ind;
            [curbestfit,ind] = min(valParents);
            curbestchrom = popold(ind, :);
%           [curbestfit, index] = min(mixVal);
           %noprogressFES = noprogressFES + mixPopSizexx;
           if prebestfit < everbestfitinrun  
               everbestfitinrun = prebestfit;
               everbestchrominrun = prebestchrom;
           end
 
           if prebestfit < everbestfitall
               everbestfitall = prebestfit;
            end                        
            
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            vi = pop + F(:, ones(1, n)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));

            vi = boundConstraint(vi, pop, lu);

            % == == == == = Crossover == == == == =
            mask = rand(popsize, n) > CR(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize n], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);

            valOffspring = feval(fhd, ui', func_num)  -  problemIndex * 100;
            valOffspring = valOffspring';

            FES = FES + popsize;
            
            gbest_temp=gbest;
            srpop_temp=srpop;

            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents, I] = min([valParents, valOffspring], [], 2);
            popold = pop;

            archive = updateArchive(archive, popold(I == 2, :), valParents(I == 2));

            popold(I == 2, :) = ui(I == 2, :);

            goodCR = CR(I == 2);
            goodF = F(I == 2);
            
            gbest=min(valParents);
            
            % ranking assignment
            [~,indFit] = sort(valParents);
            KP=indFit(1);
            gbest_Pop=popold(KP,:);
            [~, indFitTrans]=sort(indFit);
            pop_R1=indFitTrans; %the ranking of the fitness value
            
            pbest_STG(valParents==valParents_temp) = pbest_STG(valParents==valParents_temp)+1;
            pbest_STG(valParents~=valParents_temp) = 0;
            [~, indSTG]=sort(pbest_STG);
            [~, indSTGTrans]=sort(indSTG);
            pop_R2=indSTGTrans; %the ranking of the number of consecutive stagnation generation
           
            pop_R=pop_R1+pop_R2;
            
            srpop=sum(I==2)/popsize; %the update rate (UR) of the population at each generation
           
            if gbest ~= gbest_temp 
                num=randi(num,1,1);  %adaptive mechanism to determine the number of individuals taking impulsive control
            end
            
            if srpop==0 %UR is 0
                num=min(num, UN); %adaptive mechanism to determine the number of individuals taking impulsive control
                [~, rank_ind]=sort(pop_R, 'descend');
                s1_Ind=rank_ind(1:num);
                s1_Ind=setdiff(s1_Ind, KP);
                if ~isempty(s1_Ind)
                    s_size=size(s1_Ind,1);
                    randmx=rand(s_size,1);
                    pr=0.2;
                    while sum(randmx<=pr)==0 %random selection of individuals
                        pr=pr+0.2;
                        pr=min(pr,1.0);
                    end
                    ds_size=sum(randmx<=pr);
                    lu2(1,:)=min(popold);
                    lu2(2,:)=max(popold);
                    popold(s1_Ind(randmx<=pr), :)=repmat(lu2(1, :), ds_size, 1) + rand(ds_size, n) .* (repmat(lu2(2, :)-lu2(1, :), ds_size, 1)); %injecting destabilizing impulses
                    valParents3=feval(fhd, popold(s1_Ind(randmx<=pr), :)', func_num)  -  problemIndex * 100;
                    valParents(s1_Ind(randmx<=pr))=valParents3';
                    pbest_STG(s1_Ind(randmx<=pr))=0;
                    FES=FES+ds_size;
                end
            elseif srpop~=0 && srpop<srpop_temp %UR is not 0, but UR reduces
                num=min(num, UN); %adaptive mechanism to determine the number of individuals taking impulsive control
                [~, rank_ind]=sort(pop_R, 'descend'); 
                s1_Ind=rank_ind(1:num);
                s_size=size(s1_Ind,1);
                
                %stabilizing impulsive control
                DM=randi(n,s_size,1);
                rd_ind=randi([1 popsize], s_size, 1);
                [~, IK] = min([valParents(s1_Ind), valParents(rd_ind)], [], 2); % IK == 1: the stagnant one is better; IK == 2: the random one is better

                Rk=rand(s_size,n); 

                pop_si=Rk.*repmat(gbest_Pop, s_size, 1)+(1-Rk).*popold(s1_Ind,:); 
                pop_si(IK==2, :)=popold(rd_ind(IK==2), :); 
                pop2=updatePop(n, pop_si, popold, s1_Ind, DM); 
                valParents2=feval(fhd, pop2', func_num) -  problemIndex * 100;
                valParents2=valParents2';

                %check if the fitness of the new one is better
                [~,IR] = min([valParents(s1_Ind), valParents2(s1_Ind)], [], 2); % IR == 2: the new one is better
                popold(s1_Ind(IR==2), :)=pop2(s1_Ind(IR==2), :);
                valParents(s1_Ind(IR==2))=valParents2(s1_Ind(IR==2));

                FES=FES+s_size;
                
                pbest_STG(s1_Ind(IR==2))=0;
                
                if sum(IR==2)==0 
                    num=num+1; %adaptive mechanism to determine the number of individuals taking impulsive control
                    s2_Ind=s1_Ind(IR==1);
                    s2_Ind=setdiff(s2_Ind, KP);
                    if ~isempty(s2_Ind) 
                        randmx=rand(length(s2_Ind),1); 
                        pr=0.2;
                        while sum(randmx<=pr)==0 %random selection of individuals
                            pr=pr+0.2;
                            pr=min(pr,1.0);
                        end
                        ds_size=sum(randmx<=pr);
                        lu2(1,:)=min(popold);
                        lu2(2,:)=max(popold);
                        popold(s2_Ind(randmx<=pr), :)=repmat(lu2(1, :), ds_size, 1) + rand(ds_size, n) .* (repmat(lu2(2, :)-lu2(1, :), ds_size, 1)); %injecting destabilizing impulses
                        valParents3=feval(fhd, popold(s2_Ind(randmx<=pr), :)', func_num) -  problemIndex * 100;
                        valParents(s2_Ind(randmx<=pr))=valParents3';
                        pbest_STG(s2_Ind(randmx<=pr))=0;
                        FES=FES+ds_size;
                    end
                end

                if min(valParents) < gbest
                    KP_temp=find(valParents==min(valParents));
                    KP=KP_temp(1); 
                    gbest=min(valParents);
                    num=randi(num,1,1); %adaptive mechanism to determine the number of individuals taking impulsive control 
                end
            end
            [valBest indBest] = sort(valParents, 'ascend');
            gen = gen + 1;
             if (FES/(n * 10000 * beishu)) >= (xxx / 100)
                if xxx > 0
                curdiv = 0.0;   
                for x =1 : n
                   midpoint(x) = median(pop(:,x));
                end
                distobest = 1 : popsize;
                for x = 1: popsize
                    distobest (x)= 0;
                    for y = 1 : n
                        distobest(x) = distobest(x) + abs((pop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                    end
                    distobest (x) = distobest (x) / n;
                    curdiv = curdiv + distobest (x);
                end                   
                curdiv = curdiv / popsize;
                fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(valParents), min(valParents),everbestfitall);% %g %g %d %d              %,everbestfitall,everbestfitinrun,mycount,ourOpt);
                end
                xxx = xxx + 1;
            end
        end
        %disp (valParents);
        outcome = [outcome min(valParents)];

        time = time + 1;

    end

    sort(outcome)
    mean(outcome)
    std(outcome)
    fprintf(fp,'%e %e: ',mean(outcome),std(outcome));
    for x = 1 : totalTime
        fprintf(fp,'%s ', num2str(outcome(x)));
    end
    fclose(fp);    
end
toc;
