%**************************************************************************************************
% Developed by Guohua Wu, College of Information Systems and Management, National University of Defense Technology,
% Changsha, China, 410073
% guohuawu@nudt.edu.cn

% Reference: 
% G. Wu, R. Mallipeddi, P. N. Suganthan, R. Wang, and H. Chen, "Differential evolution with multi-population based ensemble of mutation strategies," Information Sciences, vol. 329, pp. 329-345, 2016.
%**************************************************************************************************

clc;
clear all;
tic;

warning off
format shortG;
format compact;

strParameterDescription = 'xxpopsize = 250 '
% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
% problemSet = 18:20;
xxpopsize = 250;
totalTime = 30;
xxleastSelectionPro = 0.2;

for problemIndex = [1:30]
    
    problem = problemIndex; %problemSet(problemIndex);
    filename = strcat(strcat('f',num2str(problem)),'_original.txt');    
    fp = fopen(filename,'a+');   
         filename2 = strcat(strcat('results','x1.txt'));
    fp2 = fopen(filename2,'a+');
    Solve = [];
    % Define the dimension of the problem
    n = 30;
    
    lu = [-100 * ones(1, n); 100 * ones(1, n)];
    
    time = 1;
    times = 1;
    D = n;
    % The total number of runs
    while time <= totalTime
        
        rand('seed', sum(100 * clock));    
        CRm2 = 0.5;
        Fm2 = 0.5;
        CRm3 = 0.5;
        Fm3 = 0.5;       
        %% the values and indices of the best solutions
%         [valBest, indBest] = sort(valParents, 'ascend');       
        FES = 0;
        leastSelectionPro = xxleastSelectionPro;
        arrayGbestChange = [1,1,1];
        arrayGbestChangeRate = [0,0,0];
        genForChange = 20;
        MaxFES = n * 10000 * times;
        mixPopSize = xxpopsize;
        MaxGen = MaxFES/mixPopSize;
        indexBestLN = 1;
        numViaLN = [0,0,0];
        rateViaLN = zeros(MaxGen,3);
        
        mixPop = repmat(lu(1, :), mixPopSize, 1) + rand(mixPopSize, D) .* (repmat(lu(2, :) - lu(1, :), mixPopSize, 1));
        %mixVal = benchmark_func(mixPop, problem, o, A, M, a, alpha, b);
         mixVal = cec14_func(mixPop',problem)'-100*problem;
        overallBestVal = min(mixVal);
 
        FES = xxpopsize;
        everbestfitall = min(mixVal);
        curdiv = 0.0;
        for x =1 : n
            midpoint(x) = median(mixPop(:,x));
        end
        distobest = 1 : xxpopsize;
        for x = 1: xxpopsize
            distobest (x)= 0;
            for y = 1 : n
                distobest(x) = distobest(x) + abs((mixPop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
            end
            distobest (x) = distobest (x) / n;
            curdiv = curdiv + distobest (x);
        end
        curdiv = curdiv / xxpopsize;
        %fprintf(fp,'%d %e %e %e\r\n', FES, min(mixVal),curdiv,everbestfitall);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(mixVal), min(mixVal), everbestfitall);                
        
        permutation = randperm(mixPopSize);
        arraySecond = permutation(1:leastSelectionPro*mixPopSize);  % for EPSDE
        arrayFirst = permutation(leastSelectionPro*mixPopSize+1: end); %  for JADE
        % arrayFirst = permutation(2*leastSelectionPro*mixPopSize+1:end);  % for JADE
        
        popold = mixPop(arrayFirst,:) ;
        valParents =mixVal(arrayFirst);
        c = 1/10;
        pj = 0.04;
        CRm = 0.5;
        Fm = 0.5;
        Afactor = 1;
        archive.NP = mixPopSize; % the maximum size of the archive
        archive.pop = zeros(0, D); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        
        gen = 0;
        FESj = 0;
        initialPopSize = mixPopSize;
        previousPopSize = mixPopSize;
        consumedFES = [1,1,1];
        everbestfit = realmax();
        curbestfit = realmax();
        curdiv = 1.0;
      %%
        while  FES < MaxFES
            gen = gen +1;         
%             if mod(gen,1) == 0
% %                 mixPopSize =  round(((4-initialPopSize)*(FES/(n * 10000))+ initialPopSize));
% %                 alpha = 0.5-(0.5-lowvalue)*((i/MaxGen)).^8;
%                 mixPopSize =  initialPopSize-(initialPopSize-4)* (FES/(n * 10000))^2;%      round(((4-initialPopSize)*(FES/(n * 10000))+ initialPopSize));
%                 mixPopSize = round(mixPopSize);
%                 [~,I]=sort(mixVal, 'descend');              
%                 mixVal(I(1:previousPopSize - mixPopSize))= [];
%                 mixPop(I(1:previousPopSize - mixPopSize),:) = [];
%                 previousPopSize = mixPopSize;
%                 [row,colum] = size(archive.pop);
%                 if ~(archive.NP >  row)
%                 archive.NP = round(0.6*mixPopSize);
%                 tempNP = archive.NP;
%                 perm = randperm(tempNP);
%                 archive.pop = archive.pop(perm,:);
%                 archive.funvalues = archive.funvalues(perm);
%                 end
%             end            
            if mod(gen,genForChange) == 0
                arrayGbestChangeRate(1) = arrayGbestChange(1)/consumedFES(1);
                arrayGbestChangeRate(2) = arrayGbestChange(2)/consumedFES(2);
                arrayGbestChangeRate(3) = arrayGbestChange(3)/consumedFES(3);
                [~,indexBestLN]=max(arrayGbestChangeRate);
                if sum(arrayGbestChangeRate == arrayGbestChangeRate(1)) == 3
                    indexBestLN = randi([1,3],1);
                end
                arrayGbestChange = [0.1,0.1,0.1];
                arrayGbestChangeRate =  [0,0,0];
                consumedFES = [1,1,1];
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
            consumedFES = consumedFES + [length(arrayFirst),length(arraySecond),length(arrayThird)];
            if mixPopSize<20
               arrayFirst =  permutation;    
               arraySecond = [];
               arrayThird = [];
            end      
%             pop = popold; % the old population becomes the current population
%%  %% ===========================mutation 1=====================================%%%%
            prebestfit = curbestfit;
            curbestfit = min(mixVal);
           if ~isempty(arrayFirst)               
            pop1 = mixPop(arrayFirst,:); % the old population becomes the current population
            valParents1 = mixVal(arrayFirst);
            popsize = length(arrayFirst);
            [~,I1]=sort(mixVal, 'ascend');
            [~,I2]=sort(valParents1, 'descend');
            pop1(I2(1),:) = mixPop(I1(1),:);
            valParents1(I2(1)) = mixVal(I1(1));
            prevalParents1 = valParents1;
            
            if FESj > 1 && ~isempty(goodCR) && sum(goodF) > 0 % If goodF and goodCR are empty, pause the update
                CRm = (1 - c) * CRm + c * mean(goodCR);
                Fm = (1 - c) * Fm + c * sum(goodF .^ 2) / sum(goodF); % Lehmer mean
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [Fj, CRj] = randFCR(popsize, CRm, 0.1, Fm, 0.1);
            r0 = [1 : popsize];
            popAll = [pop1; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents1, 'ascend');
            pNP = max(round(pj * popsize), 2); % choose at least two best solutions
            randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop1(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            vi = pop1 + Fj(:, ones(1, n)) .* (pbest - pop1 + pop1(r1, :) - popAll(r2, :));
%             if sum(problem ==[7,25]) == 0
%                 vi = boundConstraint(vi, pop1, lu);
%             end
            % == == == == = Crossover == == == == =
            mask = rand(popsize, n) > CRj(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize n], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop1(mask);
            
            %valOffspring1 = benchmark_func(ui, problem, o, A, M, a, alpha, b);
            valOffspring1 = cec14_func(ui',problem)'-100*problem;
            FESj = FESj + popsize;
            FES = FES + popsize;
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents1, I] = min([valParents1, valOffspring1], [], 2);
            popold1 = pop1;
            archive = updateArchive(archive, popold1(I == 2, :), valParents1(I == 2));
            popold1(I == 2, :) = ui(I == 2, :);
            goodCR = CRj(I == 2);
            goodF = Fj(I == 2);  
            if min(valParents1)< overallBestVal
                overallBestVal = min(valParents1);
            end
            arrayGbestChange(1) = arrayGbestChange(1) + sum(prevalParents1- valParents1);
            if prevalParents1(I2(1)) == valParents1(I2(1))
               popold1(I2(1),:) =  mixPop(arrayFirst(I2(1)),:);
               valParents1(I2(1)) = mixVal(arrayFirst(I2(1)));
            end
            mixPop(arrayFirst,:) = popold1;
            mixVal(arrayFirst) = valParents1;
           end
%             [valBest indBest] = sort(valParents1, 'ascend');
            
            %% ===========================mutation 2=====================================%%%%
           if ~isempty(arraySecond)
            pop2 = mixPop(arraySecond,:); % the old population becomes the current population 
            valParents2 = mixVal(arraySecond);
            popsize2 = length(arraySecond);
            [~,I1]=sort(mixVal, 'ascend');
            [~,I2]=sort(valParents2, 'descend');
            pop2(I2(1),:) = mixPop(I1(1),:);
            valParents2(I2(1)) = mixVal(I1(1));          
            prevalParents2 = valParents2;
            
            if gen > 1 && ~isempty(goodCR2) && sum(goodF2) > 0 % If goodF and goodCR are empty, pause the update
                CRm2 = (1 - c) * CRm2 + c * mean(goodCR2);
                Fm2 = (1 - c) * Fm2 + c * sum(goodF2 .^ 2) / sum(goodF2); % Lehmer mean
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F2, CR2] = randFCR(popsize2, CRm2, 0.1, Fm2, 0.1);
            rot = (0:1:popsize2-1);
            ind = randperm(2);
            a1  = randperm(popsize2);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize2);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize2);
            a3  = a2(rt+1);
            pm1 = pop2(a1,:);             % shuffled population 1
            pm2 = pop2(a2,:);             % shuffled population 2
            pm3 = pop2(a3,:);             % shuffled population 3
            vi =pop2 + repmat(rand(popsize2,1),1,D) .* (pm1 - pop2) + F2(:, ones(1, D)) .* (pm2 - pm3);
%             if sum(problem==[7,25]) ==0
%                 vi = boundConstraint(vi, pop2, lu);
%             end
            ui = vi;
            %valOffspring2 = benchmark_func(ui, problem, o, A, M, a, alpha, b);
            valOffspring2 = cec14_func(ui',problem)'-100*problem;
            FES = FES + popsize2;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents2, I] = min([valParents2, valOffspring2], [], 2);
            popold2 = pop2;
%             archive = updateArchive(archive, popold(I == 2, :), valParents2(I == 2));
            popold2(I == 2, :) = ui(I == 2, :);
            goodCR2 = CR2(I == 2);
            goodF2 = F2(I == 2);
            
            arrayGbestChange(2) = arrayGbestChange(2) + sum(prevalParents2- valParents2);
            if prevalParents2(I2(1)) == valParents2(I2(1))
                popold2(I2(1),:) =  mixPop(arraySecond(I2(1)),:);
                valParents2(I2(1)) = mixVal(arraySecond(I2(1)));
            end
            mixPop(arraySecond,:) = popold2;
            mixVal(arraySecond) = valParents2;   
           end
            
%             [valBest indBest] = sort(valParents2, 'ascend');
           %% ===========================mutation 3 =====================================%%%%
           if ~isempty(arrayThird)
            pop3 = mixPop(arrayThird,:); % the old population becomes the current population 
            valParents3 = mixVal(arrayThird);
            popsize3 = length(arrayThird);
            [~,I1]=sort(mixVal, 'ascend');
            [~,I2]=sort(valParents3, 'descend');
            pop3(I2(1),:) = mixPop(I1(1),:);
            valParents3(I2(1)) = mixVal(I1(1));                    
            prevalParents3 = valParents3;
            if gen > 1 && ~isempty(goodCR3) && sum(goodF3) > 0 % If goodF and goodCR are empty, pause the update
                CRm3 = (1 - c) * CRm3 + c * mean(goodCR3);
                Fm3 = (1 - c) * Fm3 + c * sum(goodF3 .^ 2) / sum(goodF3); % Lehmer mean
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F3, CR3] = randFCR(popsize3, CRm3, 0.1, Fm3, 0.1);
            rot = (0:1:popsize3-1);
            ind = randperm(2);
            a1  = randperm(popsize3);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize3);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize3);
            a3  = a2(rt+1);
            pm1 = pop3(a1,:);             % shuffled population 1
            pm2 = pop3(a2,:);             % shuffled population 2
            pm3 = pop3(a3,:);             % shuffled population 3
            vi =pm1 + F3(:, ones(1, D)) .* (pm2 - pm3); %repmat(rand(popsize3,1),1,D) .* (pm1 - pop3) 
%             if sum(problem==[7,25]) ==0
%                 vi = boundConstraint(vi, pop3, lu);
%             end
            
            mask = rand(popsize3, n) > CR3(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize3)'; cols = floor(rand(popsize3, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize3 n], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop3(mask);
            
            %valOffspring3 = benchmark_func(ui, problem, o, A, M, a, alpha, b);
            valOffspring3 = cec14_func(ui',problem)'-100*problem;
            FES = FES + popsize3;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents3, I] = min([valParents3, valOffspring3], [], 2);
            popold3 = pop3;
%             archive = updateArchive(archive, popold(I == 2, :), valParents2(I == 2));
            popold3(I == 2, :) = ui(I == 2, :);
            goodCR3 = CR3(I == 2);
            goodF3 = F3(I == 2);            
            arrayGbestChange(3) = arrayGbestChange(3) + sum(prevalParents3- valParents3);
            
            if prevalParents3(I2(1)) == valParents3(I2(1))
                popold3(I2(1),:) =  mixPop(arrayThird(I2(1)),:);
                valParents3(I2(1)) = mixVal(arrayThird(I2(1)));
            end                     
            mixPop(arrayThird,:) = popold3;
            mixVal(arrayThird) = valParents3;        
           end
           if everbestfitall > min(mixVal)
               everbestfitall = min(mixVal);
           end
           if  rem(FES,n * 100 * times ) == 0
               curdiv = 0.0;
               for x =1 : n
                   midpoint(x) = median(mixPop(:,x));
               end
               distobest = 1 : xxpopsize;
               for x = 1: xxpopsize
                   distobest (x)= 0;
                   for y = 1 : n
                       distobest(x) = distobest(x) + abs((mixPop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                   end
                   distobest (x) = distobest (x) / n;
                   curdiv = curdiv + distobest (x);
               end
               curdiv = curdiv / xxpopsize;
               fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(mixVal), min(mixVal), everbestfitall);
           end
        end
        %%
        %problem
        Solve = [Solve, everbestfitall];
        time = time + 1;
        %min(mixVal)
    end
    
    disp(['F',num2str(problem),...
        ' | ',num2str(mean(Solve)),...
        ' | ',num2str(std(Solve))]);
%     disp(sort(Solve));
    fprintf(fp,'%e %e: ', mean(Solve), std(Solve));
      fprintf(fp2,'%e %e: ', mean(Solve), std(Solve));
    for x = 1 : totalTime
        fprintf(fp,'%s ', num2str(Solve(x)));
           fprintf(fp2,'%s ', num2str(Solve(x)));
    end
    fclose(fp);    
     fprintf(fp2,'\n');
end
toc;
% save DESequentialEnsemble-50D-25RUN