%%%%%%%%%%%%%%%%%%%
%%



clc;
clear all;

format long;
format compact;

problem_size = 30;

%%% change freq
freq_inti = 0.5;
    
beishu = 1;    
max_nfes = 10000 * problem_size * beishu;


rand('seed', sum(100 * clock));

val_2_reach = 10^(-8);
max_region = 100.0;
min_region = -100.0;
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
fhd=@cec14_func;
pb = 0.4;
ps = .5;

S.Ndim = problem_size;
S.Lband = ones(1, S.Ndim)*(-100);
S.Uband = ones(1, S.Ndim)*(100);

GenMaxSelected = 250; %%% For local search

%%%% Count the number of maximum generations before as NP is dynamically
%%%% decreased 
G_Max = 0;
if problem_size == 10
    G_Max = 2163;
end
if problem_size == 30
    G_Max = 2745;
end
if problem_size == 50
    G_Max = 3022;
end
if problem_size == 100
    G_Max = 3401;
end

num_prbs = 30;
runs = 30;
run_funcvals = [];
RecordFEsFactor = ...
	[0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, ...
	0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
progress = numel(RecordFEsFactor);

allerrorvals = zeros(progress, runs, num_prbs);
result=zeros(num_prbs,5);

for func = [1:30]
    filename = strcat(strcat('F',num2str(func)),'_Eps_30_2014.txt');
    fp = fopen(filename,'a+');
    filename2 = strcat(strcat('results','_Eps_30_2014.txt'));
    fp2 = fopen(filename2,'a+');
    optimum = func * 100.0;
    S.FuncNo = func;
  
  %% Record the best results
  outcome = []; 

%   fprintf('\n-------------------------------------------------------\n')
  %fprintf('Function = %d, Dimension size = %d\n', func, problem_size) 

  for run_id = 1 : runs
      xxx = 0;
      run_funcvals = [];
      col=1;              %% to print in the first column in all_results.mat
     
    %%  parameter settings for L-SHADE
    p_best_rate = 0.11;    %0.11
    arc_rate = 1.4;
    memory_size = 5;
    pop_size = 18 * problem_size;   %18*D
    SEL = round(ps*pop_size);

    max_pop_size = pop_size;
    min_pop_size = 4.0;

     nfes = 0;
    %% Initialize the main population
    popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
    pop = popold; % the old population becomes the current population

    fitness = feval(fhd,pop',func);
    fitness = fitness';
    
    curdiv = 1;
    everbestfitall = realmax();
    everbestfitinrun = realmax();
    curbestfit = fitness(1);
    prebestfit = realmax();
        
    bsf_fit_var = 1e+30;
    bsf_index = 0;
    bsf_solution = zeros(1, problem_size);
    
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
        nfes = nfes + 1;
        
        if fitness(i) < bsf_fit_var
            bsf_fit_var = fitness(i);
            bsf_solution = pop(i, :);
            bsf_index = i;
        end
        
        if nfes > max_nfes; break; end
    end
    %%%%%%%%%%%%%%%%%%%%%%%% for out
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
    fprintf(fp,'%d %e %e %e %e\r\n', nfes, curdiv, mean(fitness)-optimum, min(fitness)-optimum, everbestfitall);    %输出格式
    disp(['nfes:' ,num2str(nfes),'  curdiv:',num2str( curdiv), '  mean:',num2str(mean(fitness)-optimum), '  min:',num2str(min(fitness)-optimum),'  everbestfitall:' num2str(everbestfitall)] );
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    
    %%% Initialize LS population
    %locacl search based on gaussian walk
    flag_LS = false;
    counter = 0;
    popsize_LS = 10;
    
    %%% Initialize LS population for re-start 
    %10*30
    popLS = repmat(lu(1, :), popsize_LS, 1) + rand(popsize_LS, problem_size) .* (repmat(lu(2, :) - lu(1, :), popsize_LS, 1));
    fitness_LS = feval(fhd,popLS',func);
    fitness_LS = fitness_LS';
    nfes = nfes + popsize_LS;
    %%%%%%%%%%%%%
    [Sorted_FitVector, Indecis] = sort(fitness_LS);
    popLS = popLS(Indecis,:);%sorting the points based on obtaind result
    %==========================================================================
     %Finding the Best point in the group=======================================
    BestPoint = popLS(1, :);
    F = Sorted_FitVector(1);%saving the first best fitness
    %%%%%%%%%%%%%
    
  
 
    run_funcvals = [run_funcvals;fitness];
    run_funcvals = [run_funcvals;fitness_LS]; 
    
    memory_sf = 0.5 .* ones(memory_size, 1);
    memory_cr = 0.5 .* ones(memory_size, 1);

    memory_freq = freq_inti*ones(memory_size, 1);
    memory_pos = 1;

    archive.NP = arc_rate * pop_size; % the maximum size of the archive
    archive.pop = zeros(0, problem_size); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions

    %% main loop
    gg=0;  %%% generation counter used For Sin
    igen =1;  %%% generation counter used For LS
 
    flag1 = false;
    flag2 = false;
    while nfes < max_nfes
      gg=gg+1;
         
      pop = popold; % the old population becomes the current population
      [temp_fit, sorted_index] = sort(fitness, 'ascend');

      mem_rand_index = ceil(memory_size * rand(pop_size, 1));
      mu_sf = memory_sf(mem_rand_index);
      mu_cr = memory_cr(mem_rand_index);
      mu_freq = memory_freq(mem_rand_index);

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
      
      
      freq = mu_freq + 0.1 * tan(pi*(rand(pop_size, 1) - 0.5));
      pos_f = find(freq <=0);
      while ~ isempty(pos_f)
        freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(length(pos_f), 1) - 0.5));
        pos_f = find(freq <= 0);
      end

      sf = min(sf, 1);
      freq = min(freq, 1);
      
      if(nfes <= max_nfes/2)
          c=rand;
          if(c<0.5)
              sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(pop_size,problem_size);
          else
              sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
          end
      end
      
      r0 = [1 : pop_size];
      popAll = [pop; archive.pop];
      [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
      
      pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
      randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
      randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
      pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

      prebestfit = curbestfit;     
      [curbestfit,ind] = min(fitness);
      curbestchrom = popold(ind, :);
      if prebestfit < everbestfitinrun
          everbestfitinrun = prebestfit;
      end      
      if prebestfit < everbestfitall
          everbestfitall = prebestfit;
      end
     
      
      vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
      vi = boundConstraint(vi, pop, lu);
      
      mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
      rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
      jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
      ui = vi; ui(mask) = pop(mask);

      children_fitness = feval(fhd, ui', func);
      children_fitness = children_fitness';

     
      %%%% To check stagnation
      flag = false;
      bsf_fit_var_old = bsf_fit_var;
      %%%%%%%%%%%%%%%%%%%%%%%% for out
      for i = 1 : pop_size
          nfes = nfes + 1;
          
          if children_fitness(i) < bsf_fit_var
              bsf_fit_var = children_fitness(i);
              bsf_solution = ui(i, :);
              bsf_index = i; 
          end
          
          if nfes > max_nfes; break; end
      end      
      %%%%%%%%%%%%%%%%%%%%%%%% for out

      dif = abs(fitness - children_fitness);


      %% I == 1: the parent is better; I == 2: the offspring is better
      I = (fitness > children_fitness);
      goodCR = cr(I == 1);  
      goodF = sf(I == 1);
      goodFreq = freq(I == 1);
      dif_val = dif(I == 1);

%      isempty(popold(I == 1, :))   
      archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));

      [fitness, I] = min([fitness, children_fitness], [], 2);
      
      run_funcvals = [run_funcvals; fitness];
      
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
          
          %% for updating the memory of freq
          if max(goodFreq) == 0 || memory_freq(memory_pos)  == -1
              memory_freq(memory_pos)  = -1;
          else
              memory_freq(memory_pos) = (dif_val' * (goodFreq .^ 2)) / (dif_val' * goodFreq);
          end
          
          memory_pos = memory_pos + 1;
          if memory_pos > memory_size;  memory_pos = 1; end
      end

      %% for resizing the population size
      plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);

      if pop_size > plan_pop_size
          reduction_ind_num = pop_size - plan_pop_size;
          if pop_size - reduction_ind_num <  min_pop_size; 
              reduction_ind_num = pop_size - min_pop_size;
          end
          
          pop_size = pop_size - reduction_ind_num;
          SEL = round(ps*pop_size);
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
       %%%%%%%%%%%%%%% Call LS based on Gaussian works when NP is less than 20 for the first time  %%%%%
      if pop_size <= 20
          counter = counter + 1;
      end

      if counter == 1
          flag_LS = true;
      else
          flag_LS = false;
      end

      if flag_LS == true
          r_index = randi([1 pop_size],1,popsize_LS); %1*10的1-20随机数
          %%% Pick 10 random individuals from L-SHADE pop
          for gen_LS = 0 : GenMaxSelected
              New_Point = [];%creating new point
              FitVector = [];%creating vector of fitness functions

              for i = 1 : popsize_LS
                  [NP, fit] = LS_Process(popLS(i,:),S,gg,BestPoint);
                  New_Point = [New_Point;NP];
                  FitVector = [FitVector,fit];
              end
           
              %%%%
              fittemp = FitVector;
              for i = 1 : popsize_LS
                  %%% Update those 10 random individuals from pop L-SHADE
                  if FitVector(i) < fitness(r_index(i))    %比较新产生的个体和原来群体的适应值 若新产生的好，就替换
                      fitness (r_index(i)) = FitVector(i);
                      pop(r_index(i),:) = New_Point(i,:);
                      
                  else
                      fittemp(i) =  fitness (r_index(i));  % fittemp记录所有差解的适应值 fitness记录好的 FitVector记录新的
                  end

                  %%%% Update best individual L-SHADE
                  if FitVector(i) < bsf_fit_var
                      bsf_fit_var = FitVector(i);
                      bsf_solution = New_Point(i,:);
                  end

                  nfes = nfes + 1;
                  if nfes > max_nfes; break; end
              end
              
              %%%%%% To recored those changes
              fittemp = fittemp';
              run_funcvals = [run_funcvals; fittemp];
              
              %%%%%%%%%%%%%%

              [SortedFit,SortedIndex] = sort(FitVector);
              New_Point = New_Point(SortedIndex,:);
              BestPoint = New_Point(1,:);%first point is the best
              BestFitness = SortedFit(1,1);
              popLS = New_Point;
          end
      end
      if everbestfitall > min(fitness)-optimum
          everbestfitall = min(fitness)-optimum;
      end
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
               disp(['nfes:' ,num2str(nfes),'  curdiv:',num2str( curdiv), '  mean:',num2str(mean(fitness)-optimum), '  min:',num2str(min(fitness)-optimum),'  everbestfitall:' num2str(everbestfitall)] );
          end
          xxx = xxx + 1;
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %%%%%%%%nfes
    
    bsf_error_val = bsf_fit_var - optimum;
    outcome = [outcome bsf_error_val];
  end %% end 1 run 
  fprintf(fp,'%e (%e): \r\n', mean(outcome), std(outcome));   %平均数：（标准差）
  fprintf(fp2,'F%d   %e  %e:  ',func , mean(outcome), std(outcome));
  for x = 1 : 30
      fprintf(fp,'%s ', num2str(outcome(x)));   %依次输出每一次运行的结果
      fprintf(fp2,'%s ', num2str(outcome(x)));
  end
  fprintf(fp2,'\n');
end %% end 1 function run


