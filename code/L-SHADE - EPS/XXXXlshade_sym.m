%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of L-SHADE which is an improved version of SHADE 1.1.
%% Note that this source code is transferred from the C++ source code version.
%% About L-SHADE, please see following papers:
%%
%% * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%
%% For this package, we downloaded JADE's source code from Dr. Q. Zhang's website (http://dces.essex.ac.uk/staff/qzhang) and modified it.
%%
%% Update
%% 9/Oct/2014: incompatibilities-fix between Octave and Matlab and bug-fix of population size reduction procedure (thanks to Dr. Elsayed)
%%%%%%%%%%%%%%%%%%%

clc;
clear all;

format long;
format compact;

problem_size = 30;
beishu = 10;
max_nfes = 10000 * problem_size * beishu;

rand('seed', sum(100 * clock));

val_2_reach = 10^(-8);
max_region = 100.0;
min_region = -100.0;
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
fhd=@cec14_func;
for func = [21:30,1:20]
    filename = strcat(strcat('F',num2str(func)),'_1f1_09.txt');
    fp = fopen(filename,'a+');
    optimum = func * 100.0;
    thediv = 0.0;
    
    %% Record the best results
    outcome = [];
    
    fprintf('\n-------------------------------------------------------\n')
    %fprintf(fp,'Function = %d, Dimension size = %d\r\n', func, problem_size);
    
    for run_id = 1 : 30
        %%  parameter settings for L-SHADE
        xxx = 0;
        p_best_rate = 0.11;
        arc_rate = 1.4;
        memory_size = 5;
        pop_size = 18 * problem_size;
        
        max_pop_size = pop_size;
        min_pop_size = 4.0;
        
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        
        fitness = feval(fhd,pop',func);
        fitness = fitness';
        
        nfes = 0;
        renfes = 0;
        %fprintf(fp,'%d %e\r\n',nfes,min(fitness)-optimum);
        bsf_fit_var = 1e+30;
        bsf_solution = zeros(1, problem_size);
        
        curdiv = 1;
        everbestfitall = realmax();
        everbestfitinrun = realmax();
        %noprogressgen = 0;
        curbestfit = fitness(1)-optimum;
        prebestfit = realmax();
        %curbestchrom = popold(1, :);
        sign = 0;
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1 : pop_size
            nfes = nfes + 1;
            renfes =renfes + 1;
            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
            end
            
            %% if mod(nfes, 1000) == 0
            %% bsf_error_var = bsf_fit_var  - optimum;
            %% if bsf_error_var < val_2_reach; bsf_error_var = 0; end;
            %%	fprintf(sprintf('%1.16e \n', bsf_error_var));
            %%    fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
            %% end
            
            %%      if nfes > max_nfes; exit(1); end
            if nfes > max_nfes; break; end
        end
        %*************
        everbestfitall = min(fitness)-optimum;
        curdiv = 0.0;
        for x =1 : problem_size
            midpoint(x) = median(pop(:,x));
        end
        distobest = 1 : pop_size;
        for x = 1: pop_size
            distobest (x)= 0;
            for y = 1 : problem_size
                distobest(x) = distobest(x) + abs((pop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
            end
            distobest (x) = distobest (x) / problem_size;
            curdiv = curdiv + distobest (x);
        end
        curdiv = curdiv / pop_size;
        fprintf(fp,'%d %e %e %e %e\r\n', nfes, curdiv, mean(fitness)-optimum, min(fitness)-optimum, everbestfitall);
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_cr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;
        
        archive.NP = arc_rate * pop_size; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        %% main loop
        while nfes < max_nfes
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
            
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            
            %% for generating crossover rate
            cr = normrnd(mu_cr, 0.1);
            term_pos = find(mu_cr == -1);
            cr(term_pos) = 0;
            cr = min(cr, 1);
            cr = max(cr, 0);
            
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            pos = find(sf <= 0);
            
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            
            sf = min(sf, 1);
            
            r0 = [1 : pop_size];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
            [r3, r4] = gnR1R2(pop_size, size(popAll, 1), r0);
            
            pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
            
            prebestfit = curbestfit;
            %            prebestchrom = curbestchrom;
            %preind = ind;
            curbestfit = min(fitness)-optimum;
            %curbestchrom = popold(ind, :);
            if prebestfit < everbestfitinrun
                everbestfitinrun = prebestfit;
                %               everbestchrominrun = prebestchrom;
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
                    if curdiv > thediv
                        migrate = 1.0;
                    else
                        migrate = 0.9;
                    end
                    everbestfitinrun = realmax();
                    archive.pop = zeros(0, problem_size); % the solutions stored in te archive
                    archive.funvalues = zeros(0, 1); % the function value of the archived solutions                    
                end
            else
                curdiv = 0.0;
                for x =1 : problem_size
                    midpoint(x) = median(pop(:,x));
                end
                distobest = 1 : pop_size;
                for x = 1: pop_size
                    distobest (x)= 0;
                    for y = 1 : problem_size
                        distobest(x) = distobest(x) + abs((pop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                    end
                    distobest (x) = distobest (x) / problem_size;
                    curdiv = curdiv + distobest (x);
                end
                curdiv = curdiv / pop_size;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%sign = 0;%表示下一轮恢复正常
            end
            
            if sign == 0
                vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            else
                vi = pop + sf(:, ones(1, problem_size)) .* (pop(r3, :) - popAll(r4, :));               
            end
            vi = boundConstraint(vi, pop, lu);
            
            mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);
            
            if sign ==0 || (sign == 1 && curdiv > thediv) 
                children_fitness = feval(fhd, ui', func);
                children_fitness = children_fitness';
                if sign == 1 && curdiv > thediv%%%%%%%%%%%%%%%%%%%%% 
                    [d1,d2] = size(ui);
                    ui_sym = repmat(lu(2, :) + lu(1, :), d1, 1) - ui;
                    S = rand(pop_size,1) < migrate;
                    ui(S,:) = ui_sym(S,:);
                    children_fitness = feval(fhd, ui', func);
                    children_fitness = children_fitness';
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            for i = 1 : pop_size
                nfes = nfes + 1;
                renfes = renfes + 1;
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i, :);
                end
                
                %% if mod(nfes, 1000) == 0
                %% bsf_error_var = bsf_fit_var  - optimum;
                %% if bsf_error_var < val_2_reach; bsf_error_var = 0; end;
                %%       fprintf(sprintf('%1.16e \n', bsf_error_var));
                %%       fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
                %%end
                
                %%	if nfes > max_nfes; exit(1); end
                if nfes > max_nfes; break; end
            end
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            if sign == 0
                dif = abs(fitness - children_fitness);
                %% I == 1: the parent is better; I == 2: the offspring is better
                I = (fitness > children_fitness);
                goodCR = cr(I == 1);
                goodF = sf(I == 1);
                dif_val = dif(I == 1);
                
                %      isempty(popold(I == 1, :))
                archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
                
                [fitness, I] = min([fitness, children_fitness], [], 2);
                
                popold = pop;
                popold(I == 2, :) = ui(I == 2, :);
                
                num_success_params = numel(goodCR);
                
                if num_success_params > 0
                    sum_dif = sum(dif_val);
                    dif_val = dif_val / sum_dif;
                    
                    %% for updating the memory of scaling factor
                    memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                    
                    %% for updating the memory of crossover rate
                    if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                        memory_cr(memory_pos)  = -1;
                    else
                        memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                    end
                    
                    memory_pos = memory_pos + 1;
                    if memory_pos > memory_size;  memory_pos = 1; end
                end
                
                %% for resizing the population size
                plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * renfes) + max_pop_size);%%%%here
                
                if pop_size > plan_pop_size
                    reduction_ind_num = pop_size - plan_pop_size;
                    if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
                    
                    pop_size = pop_size - reduction_ind_num;
                    for r = 1 : reduction_ind_num
                        [valBest indBest] = sort(fitness, 'ascend');
                        worst_ind = indBest(end);
                        popold(worst_ind,:) = [];
                        pop(worst_ind,:) = [];
                        fitness(worst_ind,:) = [];
                    end
                    
                    archive.NP = round(arc_rate * pop_size);
                    
                    if size(archive.pop, 1) > archive.NP
                        rndpos = randperm(size(archive.pop, 1));
                        rndpos = rndpos(1 : archive.NP);
                        archive.pop = archive.pop(rndpos, :);
                    end
                end
                
            else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sign=1
                
                if pop_size < our_mem_pop_size
                    fitness = [fitness; children_fitness];%确实扩张了。
                    popold = [pop; ui];%%%%%%%
                    temp = size(fitness);%%%%%%%
                    pop_size = temp(1);%%%%%%%
                else
                    fitness = children_fitness;                
                    popold = pop;
                    popold = ui;
                end
            end
            if sign == 1 && curdiv > thediv%mycount > 3000 + upcount
                sign = 0;%表示这一代开始前进
                renfes = 0;%%%%%%%%%%%%%%%
            end
            
%            dif = abs(fitness - children_fitness);
%            
%             if sign == 0
%                 %% I == 1: the parent is better; I == 2: the offspring is better
%                 I = (fitness > children_fitness);
%                 goodCR = cr(I == 1);
%                 goodF = sf(I == 1);
%                 dif_val = dif(I == 1);
%                 
%                 %      isempty(popold(I == 1, :))
%                 archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
%                 
%                 [fitness, I] = min([fitness, children_fitness], [], 2);
%                 
%                 popold = pop;
%                 popold(I == 2, :) = ui(I == 2, :);
%                 
%                 num_success_params = numel(goodCR);
%             else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% selection for restart
%                 %%  parameter settings for L-SHADE
%                 p_best_rate = 0.11;
%                 arc_rate = 1.4;
%                 memory_size = 5;
%                 pop_size = 18 * problem_size;
%                 
%                 max_pop_size = pop_size;
%                 min_pop_size = 4.0;
%                 
%                 %% Initialize the main population
%                 popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
%                 pop = popold; % the old population becomes the current population
%                 
%                 fitness = feval(fhd,pop',func);
%                 fitness = fitness';
%                 
%                 bsf_fit_var = 1e+30;
%                 bsf_solution = zeros(1, problem_size);
%                 
%                 curdiv = 1;
%                 everbestfitinrun = realmax();
%                 %noprogressgen = 0;
%                 curbestfit = fitness(1);
%                 prebestfit = realmax();
%                 %curbestchrom = popold(1, :);
%                 sign = 0;
%                 %%%%%%%%%%%%%%%%%%%%%%%% for out
%                 for i = 1 : pop_size
%                     nfes = nfes + 1;
%                     renfes = renfes + 1;
%                     if fitness(i) < bsf_fit_var
%                         bsf_fit_var = fitness(i);
%                         bsf_solution = pop(i, :);
%                     end
%                     
%                     %% if mod(nfes, 1000) == 0
%                     %% bsf_error_var = bsf_fit_var  - optimum;
%                     %% if bsf_error_var < val_2_reach; bsf_error_var = 0; end;
%                     %%	fprintf(sprintf('%1.16e \n', bsf_error_var));
%                     %%    fprintf(sprintf('%d %1.16e \n', nfes, bsf_error_var));
%                     %% end
%                     
%                     %%      if nfes > max_nfes; exit(1); end
%                     if nfes > max_nfes; break; end
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%% for out
%                 
%                 memory_sf = 0.5 .* ones(memory_size, 1);
%                 memory_cr = 0.5 .* ones(memory_size, 1);
%                 memory_pos = 1;
%                 
%                 archive.NP = arc_rate * pop_size; % the maximum size of the archive
%                 archive.pop = zeros(0, problem_size); % the solutions stored in te archive
%                 archive.funvalues = zeros(0, 1); % the function value of the archived solutions
%                 
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%selection for restart
            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;
                
                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                
                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                    memory_cr(memory_pos)  = -1;
                else
                    memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end
                
                memory_pos = memory_pos + 1;
                if memory_pos > memory_size;  memory_pos = 1; end
            end
            
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * renfes) + max_pop_size);
            
            if pop_size > plan_pop_size
                reduction_ind_num = pop_size - plan_pop_size;
                if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
                
                pop_size = pop_size - reduction_ind_num;
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
                
                archive.NP = round(arc_rate * pop_size);
                
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (nfes / max_nfes) >= (xxx / 100)
                if  xxx > 0
                    curdiv = 0.0;
                    for x =1 : problem_size
                        midpoint(x) = median(pop(:,x));
                    end
                    distobest = 1 : pop_size;
                    for x = 1: pop_size
                        distobest (x)= 0;
                        for y = 1 : problem_size
                            distobest(x) = distobest(x) + abs((pop(x,y) - midpoint(y))/(lu(2, y) - lu(1, y)));
                        end
                        distobest (x) = distobest (x) / problem_size;
                        curdiv = curdiv + distobest (x);
                    end
                    curdiv = curdiv / pop_size;
                    fprintf(fp,'%d %e %e %e %e\r\n', nfes, curdiv, mean(fitness)-optimum, min(fitness)-optimum, everbestfitall);
                end
                xxx = xxx + 1;
            end
        end
        
        
        bsf_error_val = min(everbestfitall,min(fitness));%bsf_fit_var - optimum;这句被我改了
        %fprintf('%f  %f',min(fitness)-optimum,bsf_error_val);
        %if bsf_error_val < val_2_reach
        %    bsf_error_val = 0;
        %end
        
        %fprintf(fp,'%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
        outcome = [outcome bsf_error_val];
    end %% end 1 run
    
    fprintf(fp,'%e (%e): \r\n', mean(outcome), std(outcome));
    for x = 1 : 30
        fprintf(fp,'%s ', num2str(outcome(x)));
    end
    
end %% end 1 function run
