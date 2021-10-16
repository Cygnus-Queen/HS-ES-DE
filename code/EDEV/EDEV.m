% EDEV
% Developed by Guohua Wu, College of Information Systems and Management,
% National University of Defense Technology, Changsha, 410073,  guohuawu.nudt@gmail.com

clc;
clear all;
tic;
warning off
format shortG;
format compact;

'EDEV'
mixPopSizexx = 100
genForChangexx = 20
leastSelectionProxx = 0.1

D = 30;
n = 30;
beishu = 1;
lu = [-80 * ones(1, n); 80 * ones(1, n)];
for fun = [1:30]
    filename = strcat(strcat('f',num2str(fun)),'_30D.txt');  
    fp = fopen(filename,'a+'); 
     filename2 = strcat(strcat('results','_2014x1.txt'));
    fp2 = fopen(filename2,'a+');
    Solve = [];
    for run=1:30
        % Define the dimension of the problem
        rand('seed', sum(100*clock));
        %% general initial
        FES = 0;
        xxx = 0;
       
        
        everbestfit = realmax();
        curbestfit = realmax();
        curdiv = 1;
        
        
        leastSelectionPro =leastSelectionProxx;
        arrayGbestChange = [1,1,1];
        arrayGbestChangeRate = [0,0,0];
        arrayCondidateLN = [1,5,20];
        genForChange = genForChangexx;
        MaxFES = D*10000*beishu;%%%%%%%%%%%%%%%%
        mixPopSize = mixPopSizexx;
        MaxGen = MaxFES/mixPopSize;
        indexBestLN = 1;
        numViaLN = [0,0,0];
        rateViaLN = [];

        
        mixPop = repmat(lu(1, :), mixPopSize, 1) + rand(mixPopSize, D) .* (repmat(lu(2, :) - lu(1, :), mixPopSize, 1));
        temp = cec14_func(mixPop', fun);
        [mmm, nnn] = size(temp);       
        mixVal = (temp - ones(mmm, nnn) * fun *100)';%benchmark_func(mixPop, fun, o, A, M, a, alpha, b);
        %disp(mixVal);
        %disp('-----------');
%         valOffspringxx = cec14_func(ui', fun);
        overallBestVal = min(mixVal);
        popsize=mixPopSize;
        permutation = randperm(mixPopSize);
        arrayThird= permutation(1:leastSelectionPro*mixPopSize);  % for EPSDE
        arraySecond = permutation(leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize); % for CoDE
        arrayFirst = permutation(2*leastSelectionPro*mixPopSize+1:end);  % for JADE
     
        %% Initialize JADE related
        popold = mixPop(arrayFirst,:) ;
        valParents =mixVal(arrayFirst);
        c = 1/10;
        pj = 0.05;
        CRm = 0.5;
        Fm = 0.5;
        Afactor = 1;
        archive.NP = length(arrayFirst); % the maximum size of the archive
        archive.pop = zeros(0, D); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        %the values and indices of the best solutions
        [valBest, indBest] = sort(valParents, 'ascend');
        
     
        
        %% Initialize CoDE related
        popCODE = mixPop(arraySecond,:) ;
        valCODE = mixVal(arraySecond);
        
        %% Initialize EPSDE related
        I_D = D;
        FM_pop = mixPop(arrayThird,:) ;
        FM_popold  = zeros(size(FM_pop)); % toggle population
        I_NP = length(arrayThird);
        val   = zeros(I_NP, 1);    %create and reset the "cost array"
        FVr_bestmem = zeros(1, I_D); % best population member ever
        FVr_bestmemit = zeros(1, I_D); % best population member in iteration
        I_nfeval  = 0;      % number of function evaluations
        
        FVr_minbound = lu(1, :);
        FVr_maxbound = lu(2, :);
        Lbound = lu(1, :);
        Ubound = lu(2, :);
        %------Evaluate the best member after initialization----------------------
        
        I_best_index = 1;     % start with first population member
        temp = cec14_func(FM_pop(I_best_index, :)',fun);
        [mmm, nnn] = size(temp);       
        valxx  = temp - ones(mmm, nnn) * fun *100;%benchmark_func(FM_pop(I_best_index, :), fun, o, A, M, a, alpha, b);    
        val(1) =valxx';
        F_bestval = val(1);     % best objective function value so far
        FES = FES + 1;
        for k = 2 : I_NP       % check the remaining members
            valxx = cec14_func(FM_pop(k, :)',fun)  - fun * 100;%benchmark_func(FM_pop(k, :), fun, o, A, M, a, alpha, b);
            val(k) = valxx';
            FES = FES + 1;
            if (val(k) < F_bestval)
                
                I_best_index = k;    % save its location
                F_bestval  = val(k);
            end
        end
        FVr_bestmem = FM_pop(I_best_index, :); % best member of current iteration
        
        %------DE-Minimization---------------------------------------------
        %------FM_popold is the population which has to compete. It is--------
        %------static through one iteration. FM_pop is the newly--------------
        %------emerging population.----------------------------------------
        I_iter = 1;
        FF = [0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
        CR = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
        spara = randperm(3)';
        Para = zeros(I_NP, 3);
        inde = 1 : I_NP;
        inde = inde';
        RR = zeros(mixPopSize, 1);
        PPara = [];
        gen = 1;
        FESj = 0;
        cc = 1;
        curbestchrom = mixPop(1, :);
        everbestfitinrun = realmax();
        everbestfitall = realmax();
        
         for x =1 : n
           midpoint(x) = median(mixPop(:,x));
        end
        distobest = 1 : popsize;
        for x = 1: popsize
            distobest (x)= 0;
            for y = 1 : n
                distobest(x) = distobest(x) + abs((mixPop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
            end
            distobest (x) = distobest (x) / n;
            curdiv = curdiv + distobest (x);
        end                   
        curdiv = curdiv / popsize;
        everbestfitall = min(valParents);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(valParents), min(valParents),everbestfitall);%
        
        
        while FES<=MaxFES%gen<=MaxGen && 
            
            %% update the index for adjust different algorithms
            %         [gen,overallBestVal]
            if mod(gen,genForChange) == 0
                arrayGbestChangeRate(1) = arrayGbestChange(1)/length(arrayFirst);
                arrayGbestChangeRate(2) = arrayGbestChange(2)/(length(arraySecond)*3);%(length(arraySecond)*3);
                arrayGbestChangeRate(3) = arrayGbestChange(3)/length(arrayThird);
                [~,indexBestLN]=max(arrayGbestChangeRate);
%                 if sum(arrayGbestChangeRate == arrayGbestChangeRate(1)) == 3
%                     indexBestLN = randi([1,3],1);
%                 end              
                arrayGbestChange = [1,1,1];
                arrayGbestChangeRate =  [0,0,0];
            end
            permutation = randperm(mixPopSize);
            if indexBestLN == 1
                arrayThird= permutation(1:leastSelectionPro*mixPopSize);
                arraySecond = permutation(leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize);
                arrayFirst = permutation(2*leastSelectionPro*mixPopSize+1:end);
                numViaLN(1) = numViaLN(1) + 1;
            elseif indexBestLN == 2
                arrayThird = permutation(1:leastSelectionPro*mixPopSize);
                arrayFirst = permutation(leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize);
                arraySecond  = permutation(2*leastSelectionPro*mixPopSize+1:end);
                numViaLN(2) = numViaLN(2) + 1;
            elseif indexBestLN == 3
                arrayFirst = permutation(1:leastSelectionPro*mixPopSize);
                arraySecond = permutation(leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize);
                arrayThird  = permutation(2*leastSelectionPro*mixPopSize+1:end);
                numViaLN(3) = numViaLN(3) + 1;
            end
            rateViaLN(gen,:) = numViaLN/sum(numViaLN);
            %prebestfit = curbestfit;
            %curbestfit = min(mixVal);   
            prebestfit = curbestfit;
            prebestchrom = curbestchrom;
            %preind = ind;
            [curbestfit,ind] = min(mixVal);
            curbestchrom = mixPop(ind, :);
%           [curbestfit, index] = min(mixVal);
           %noprogressFES = noprogressFES + mixPopSizexx;
           if prebestfit < everbestfitinrun  
               everbestfitinrun = prebestfit;
               everbestchrominrun = prebestchrom;
           end
 
           if prebestfit < everbestfitall
               everbestfitall = prebestfit;
            end                        
            %% forllowing is JADE***************************************************************************
            %*************************************************************************************************************
            pop = mixPop(arrayFirst,:); % the old population becomes the current population
            valParents = mixVal(arrayFirst);
            popsize = length(arrayFirst);        
            if FESj > 1 && ~isempty(goodCR) && sum(goodF) > 0 % If goodF and goodCR are empty, pause the update
                CRm = (1 - c) * CRm + c * mean(goodCR);
                Fm = (1 - c) * Fm + c * sum(goodF .^ 2) / sum(goodF); % Lehmer mean
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [Fj, CRj] = randFCR(popsize, CRm, 0.1, Fm, 0.1);
            r0 = [1 : popsize];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents, 'ascend');
            pNP = max(round(pj * popsize), 2); % choose at least two best solutions
            randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            vi = pop + Fj(:, ones(1, n)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            if sum(fun ==[7,25]) == 0
                vi = boundConstraint(vi, pop, lu);
            end
            % == == == == = Crossover == == == == =
            mask = rand(popsize, n) > CRj(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize n], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            
%             valOffspring = benchmark_func(ui, fun, o, A, M, a, alpha, b);
            [mmm, nnn] = size(cec14_func(ui', fun));
            valOffspringxx = cec14_func(ui', fun) - ones(mmm, nnn) * fun *100;%benchmark_func(ui, fun, o, A, M, a, alpha, b);
            valOffspring = valOffspringxx';
            FESj = FESj + popsize;
            FES = FES + popsize;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents, I] = min([valParents, valOffspring], [], 2);
            popold = pop;
            archive = updateArchive(archive, popold(I == 2, :), valParents(I == 2));
            popold(I == 2, :) = ui(I == 2, :);
            goodCR = CRj(I == 2);
            goodF = Fj(I == 2);
            
            if min(valParents)< overallBestVal
                overallBestVal = min(valParents);
            end
            arrayGbestChange(1) = arrayGbestChange(1) + sum(mixVal(arrayFirst)- valParents);
            mixPop(arrayFirst,:) = popold;
            mixVal(arrayFirst) = valParents;
            
            
            %% following is for CODE****************************************************************************
            %*************************************************************************************************************
            popCODE = mixPop(arraySecond,:) ;
            valCODE = mixVal(arraySecond);
            popsizeC= length(arraySecond);              
            pTemp = popCODE;
            fitTemp = valCODE;
            % uSet: the set of trial vectors
            uSet = zeros(3 * popsizeC, n);
            for i = 1 : popsizeC
                % the three control parameter settings
                F    = [1.0 1.0 0.8];
                CRC = [0.1 0.9 0.2];
                % Uniformly and randomly select one of the control
                % parameter settings for each trial vector generation strategy
                paraIndex = floor(rand(1, 3) * length(F)) + 1;
                % Generate the trail vectors
                u = generator(popCODE, lu, i, F, CRC, popsizeC, n, paraIndex,fun);
                uSet(i * 3 - 2 : 3 * i, :) = u;
                FES = FES + 3;
            end
            % Evaluate the trial vectors
%             fitSet = benchmark_func(uSet, fun, o, A, M, a, alpha, b);
            %disp (cec14_func(uSet', fun));
            [mmm, nnn] = size(cec14_func(uSet', fun));
            fitSetxx = cec14_func(uSet', fun) - ones(mmm, nnn) * fun *100; %benchmark_func(uSet, fun, o, A, M, a, alpha, b);
            fitSet = fitSetxx';
            for i = 1 : popsizeC
                % Choose the best trial vector from the three trial vectors
                [minVal, minID] = min(fitSet(3 * i - 2 : 3 * i, :));
                bestInd = uSet(3 * (i - 1) + minID, :);
                bestIndFit = fitSet(3 * (i - 1) + minID, :);
                % Choose the better one between the trial vector and the
                % target vector
                if valCODE(i) >= bestIndFit
                    pTemp(i, :) = bestInd;
                    fitTemp(i, :) = bestIndFit;
                end
            end
            if  min(fitTemp) < overallBestVal
                overallBestVal = min(fitTemp) ;
            end
            popCODE = pTemp;
            valCODE = fitTemp;
            arrayGbestChange(2) = arrayGbestChange(2) + sum(mixVal(arraySecond)- valCODE);
            mixPop(arraySecond,:) = popCODE;
            mixVal(arraySecond) = valCODE;
            
            %% following is for EPSDE******************************************************************************
            %*************************************************************************************************************
            
            FM_pop = mixPop(arrayThird,:) ;
            val = mixVal(arrayThird);
            I_NP = length(arrayThird);
                    
            inde = 1 : I_NP;
            inde = inde';
            I_iter =gen;
            if(I_iter == 1 || length(Para)<I_NP)
                Para(inde, :) = [spara(randi(size(spara, 1),  [1, size(inde, 1)])), CR(randi(size(CR, 1), [1,size(inde, 1)])), FF(randi(size(FF, 1),  [1, size(inde, 1)]))];
            else
                for k = 1 : size(inde, 1)
                    if((rand <= RATE && ~isempty(PPara))|| length(RR) < I_NP)
                        RR(k) = randi(size(PPara, 1),  [1,1]);
                        Para(inde(k), :) = PPara(randi(size(PPara, 1),  [1, 1]), :);
                    else
                        RR(k) = 0;
                        Para(inde(k), :) = [spara(randi(size(spara, 1), [1, 1])), CR(randi(size(CR, 1),  [1, 1])), FF(randi(size(FF, 1), [1, 1]))];
                    end
                end
            end
            RRR = [];
            count = 0;
            FM_popold = FM_pop;
            
            for i = 1 : I_NP
                FM_mui = rand(1, I_D) < Para(i, 2);
                dd = find(FM_mui == 1);
                if isempty(dd)
                    ddd = ceil(rand * I_D);
                    FM_mui(ddd) = 1;
                end
                FM_mpo = FM_mui < 0.5;
                FM_bm = FVr_bestmem;
                para(i, :) = normrnd(Para(i, 3), 0.001, 1, I_D);
                if(Para(i, 1) == 1)
                    %DE/best/2/bin
                    ind = randperm(I_NP);
                    FM_pm3 = FM_popold(ind(1), :);
                    FM_pm4 = FM_popold(ind(2), :);
                    FM_pm5 = FM_popold(ind(3), :);
                    FM_pm6 = FM_popold(ind(4), :);
                    FM_ui(i, :) = FM_bm +(FM_pm3 - FM_pm4 + FM_pm5 - FM_pm6) .* para(i, :);
                    FM_ui(i, :) = FM_popold(i, :) .* FM_mpo + FM_ui(i, :) .* FM_mui;
                end
                if(Para(i, 1) == 2)
                    %DE/rand/1/bin
                    ind = randperm(I_NP);
                    FM_pm7 = FM_popold(ind(1), :);
                    FM_pm8 = FM_popold(ind(2), :);
                    FM_pm9 = FM_popold(ind(3), :);
                    FM_ui(i, :) = FM_pm7 +para(i, :) .* (FM_pm8 - FM_pm9);
                    FM_ui(i, :) = FM_popold(i, :) .* FM_mpo+ FM_ui(i, :) .* FM_mui;  % crossover
                end
                if(Para(i, 1) == 3)
                    % DE/current-to-rand/1/bin/
                    ind = randperm(I_NP);
                    FM_pm21 = FM_popold(ind(1), :);
                    FM_pm22 = FM_popold(ind(2), :);
                    FM_pm23 = FM_popold(ind(3), :);
                    FM_ui(i, :) = FM_popold(i, :) + rand(1, I_D) .* (FM_pm21-FM_popold(i, :)) + para(i, :) .* (FM_pm22 - FM_pm23);  % differential variation
                end
                if sum(fun ==[7,25]) == 0
                    if(FM_ui(i, :) < Lbound)
                        FM_ui(i, :) = FVr_minbound + (FVr_maxbound-FVr_minbound) .* rand(1, I_D);
                    end
                    if(FM_ui(i, :) > Ubound)
                        FM_ui(i, :) = FVr_minbound + (FVr_maxbound-FVr_minbound) .* rand(1, I_D);
                    end
                end
            end
%             tempval = benchmark_func(FM_ui, fun, o, A, M, a, alpha, b);
            %disp(cec14_func(FM_ui', fun));
            [mmm, nnn] = size(cec14_func(FM_ui', fun));
            tempvalxx = cec14_func(FM_ui', fun) - ones(mmm, nnn) * fun *100;
            tempval = tempvalxx';
            FES = FES + length(arrayThird);
            for i = 1 : I_NP
                I_nfeval = I_nfeval + 1;
                if (tempval(i) < val(i))
                    FM_pop(i, :) = FM_ui(i, :);
                    val(i) = tempval(i);
                    PPara = [Para(i, :); PPara];
                    if(RR(i) ~= 0)
                        RRR = [RRR; RR(i)];
                    end
                    if (tempval(i) < F_bestval)
                        F_bestval = tempval(i);
                        FVr_bestmem = FM_ui(i, :);
                        I_best_index = i;
                    end
                else
                    count = count + 1;
                end
            end
            if min(val)< overallBestVal
                overallBestVal = min(val);
            end
            arrayGbestChange(3) = arrayGbestChange(3) + sum(mixVal(arrayThird)- val);
            mixPop(arrayThird,:) = FM_pop;
            mixVal(arrayThird) = val;
            
            PPara(RRR, :) = [];
            rate(I_iter, 1) = count / I_NP;
            if(I_iter > 10)
                RATE = mean(rate((I_iter-10) : I_iter), 1);
            else
                RATE = mean(rate, 1);
            end
            I_iter = I_iter + 1;
            gen = gen+1;
            %         [gen,overallBestVal]
           if everbestfit > curbestfit
               everbestfit = curbestfit;
           end
           %cc = cc + 1;
           %disp(curbestfit);
           %if FES >= cc * (MaxFES / 100 )
           %    fprintf(fp,'%d %g %g %g\r\n', gen, curbestfit, curdiv, everbestfit);
           %end
 %           disp(gen);
 %%%%%%%
%             if rem(gen, MaxGen / 800) == 0 
%                 curdiv = 0.0;   
%                 for x =1 : n
%                    midpoint(x) = median(mixPop(:,x));
%                 end
%                 distobest = 1 : mixPopSizexx;
%                 for x = 1: mixPopSizexx
%                     distobest (x)= 0;
%                     for y = 1 : n
%                         distobest(x) = distobest(x) + abs((mixPop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
%                     end
%                     distobest (x) = distobest (x) / n;
%                     curdiv = curdiv + distobest (x);
%                 end                   
%                 curdiv = curdiv / mixPopSizexx;
%                 %disp(curdiv);
%            end  
           % 
           if everbestfit > curbestfit
               everbestfit = curbestfit;
           end
           %cc = cc + 1;
           %disp(curbestfit);
           %if FES >= cc * (MaxFES / 100 )
           %    fprintf(fp,'%d %g %g %g\r\n', gen, curbestfit, curdiv, everbestfit);
           %end
 %           disp(gen);
            
         if (FES/ MaxFES) >= (xxx / 100)
              if xxx > 0
                curdiv = 0.0;   
                for x =1 : n
                   midpoint(x) = median(mixPop(:,x));
                end
                distobest = 1 : popsize;
                for x = 1: popsize
                    distobest (x)= 0;
                    for y = 1 : n
                        distobest(x) = distobest(x) + abs((mixPop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                    end
                    distobest (x) = distobest (x) / n;
                    curdiv = curdiv + distobest (x);
                end                   
                curdiv = curdiv / popsize;
                fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(valParents), min(valParents),everbestfitall);% %g %g %d %d              %,everbestfitall,everbestfitinrun,mycount,ourOpt);
             end
                xxx = xxx + 1;
         end
             
             
 %%%%%%%%
        end
        disp('one time.');
        disp(everbestfit);
        %overallBestVal
        overallBestValSet(fun,run) = overallBestVal;
        Solve = [Solve everbestfit];
    end
    fprintf(fp,'%s %s: ', num2str(mean(Solve)), num2str(std(Solve)));
     fprintf(fp2,'%s %s: ', num2str(mean(Solve)), num2str(std(Solve)));
    for x = 1 : run
        fprintf(fp,'%s ', num2str(Solve(x)));
           fprintf(fp2,'%s ', num2str(Solve(x)));
    end
    fprintf(fp2,'\n');
    fclose(fp);        
end
%genForChangexx