function [x, xold, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq,current_eval,res_det ] = ...
    EBO( x,xold, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,archive_T,archive_freq, xmin, xmax,  n,...
    PopSize,  current_eval, I_fno,res_det,Printing,Max_FES, G_Max, gg)


vi=zeros(PopSize,n);

%% calc CR and F
mem_rand_index = ceil(memory_size * rand(PopSize, 1));
mu_sf = archive_f(mem_rand_index);
mu_cr = archive_Cr(mem_rand_index);
mu_T  = archive_T(mem_rand_index);
mu_freq = archive_freq(mem_rand_index);
%% ========================= generate CR ==================================
cr = (mu_cr + 0.1*sqrt(pi)*(asin(-rand(1,PopSize))+asin(rand(1,PopSize))))';
cr(mu_cr == -1) = 0;
cr = min(cr, 1);
cr = max(cr, 0);


%% ========================= generate F ===================================
F = mu_sf + 0.1 * tan(pi * (rand( 1,PopSize) - 0.5));
pos = find(F <= 0);
while ~ isempty(pos)
    F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand( 1,length(pos)) - 0.5));
    pos = find(F <= 0);
end
F = min(F, 1)';

%% ========================= generate T ===================================
T = mu_T + 0.05*(sqrt(pi)*(asin(-rand(1, PopSize))+asin(rand(1, PopSize))));
T = max(T,0)'; T = min(T,0.5)';
l = floor(n*rand(1,PopSize))+1;
CR = [];
if n == 1
    CR = cr;
else
for i = 1:PopSize
    if rem(n,2) == 0
       mm = exp(-T(i)/n*(0:n/2-1));
       ll = cr(i).*[mm fliplr(mm)];
       CR(i,[l(i):n (1:l(i)-1)]) =ll;
    else
       mm = exp(-T(i)/n*(0:floor(n/2-1)));
       mm1 = exp(-T(i)/n*floor(n/2));
       ll = cr(i).*[mm mm1 fliplr(mm)];
       CR(i,[l(i):n (1:l(i)-1)]) =ll;
    end
end
end

%% ========================= genrate freq =================================
freq = mu_freq + 0.1 * tan(pi*(rand(1, PopSize) - 0.5));
pos_f = find(freq <=0);
while ~ isempty(pos_f)
      freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(1,length(pos_f)) - 0.5));
      pos_f = find(freq <= 0);
end
freq = min(freq, 1)';
if(current_eval <= Max_FES/2)
   c=rand;
   if(c<0.5)
      F = 0.5.*( tan(2.*pi.*0.5.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(PopSize,1);
   else
      F = 0.5 *( tan(2*pi .* freq(:, ones(1, 1)) .* gg) .* (gg/G_Max) + 1 ) .* ones(PopSize, 1);
   end
end
 
%% ======================== generate new x =================================
popAll = [x;archive.pop]; %% set archive
r0 = 1 : PopSize;
%% generate random integer numbers
[r1, r2,r3] = gnR1R2(PopSize, size(popAll, 1), r0);

%% mutation
bb= rand(PopSize, 1);
probiter = prob(1,:);
l2= sum(prob(1:2));
op_1 = bb <=  probiter(1)*ones(PopSize, 1);
op_2 = bb > probiter(1)*ones(PopSize, 1) &  bb <= (l2*ones(PopSize, 1)) ;



pNP = max(round(0.1 * PopSize), 2); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
randindex = bestt(PopSize,n);
phix = x(randindex, :); 
%% crisscross modification

vi(op_1==1,:) = x(op_1==1,:)+ F(op_1==1, ones(1, n)) .*(x(r1(op_1==1),:) - x(op_1==1,:) + x(r3(op_1==1), :) - popAll(r2(op_1==1), :));

%% towards-best modification
vi(op_2==1,:) =  x(op_2==1,:)+ F(op_2==1, ones(1, n)) .*(phix(op_2==1,:) - x(op_2==1,:) + x(r1(op_2==1), :) - x(r3(op_2==1), :));%+w.*( x(op_2==1,:)- xold(op_2==1,:));



%% handle boundaries
vi = han_boun(vi, xmax, xmin, x,PopSize,1);
%% crossover
mask = rand(PopSize, n) > CR; 
% mask = rand(PopSize, n) > cr(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
% mask = expo_mex(PopSize,n,cr,T);
rows = (1 : PopSize)'; cols = floor(rand(PopSize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
jrand = sub2ind([PopSize n], rows, cols); mask(jrand) = false;
ui = vi; ui(mask) = x(mask);
% ui = x; ui(mask) = vi(mask);
%% evaluate
fitx_new = cec14_func(ui',I_fno);
%% update FITNESS EVALUATIONS
current_eval =current_eval+PopSize;

%% calc. imprv. for Cr and F
diff = abs(fitx - fitx_new);
I =(fitx_new < fitx);
goodCR = cr(I == 1);
goodF = F(I == 1);
goodT = T(I == 1)';
goodFreq = freq(I == 1);

%% ========================= update archive ===============================
 archive = updateArchive(archive, x(I == 1, :), fitx(I == 1)');
%% ==================== update Prob. of each DE ===========================
diff2 = max(0,(fitx - fitx_new))./abs(fitx);
count_S(1)=max(0,mean(diff2(op_1==1)));
count_S(2)=max(0,mean(diff2(op_2==1)));
% count_S(3)=max(0,mean(diff2(op_3==1)));

%% update probs.
%%% Althouth count_S~=0 may slow down the convergence, it gives more
%%% diversity. In case you look for a better convergence you can set it to
%%% sum(count_S)~=0 
if count_S~=0 
prob= max(0.1,min(0.9,count_S./(sum(count_S))));
else
    prob=1/2 * ones(1,2);
end
%% ==================== update x and fitx =================================
fitx(I==1)= fitx_new(I==1); xold(I == 1, :) = x(I == 1, :);
x(I == 1, :) = ui(I == 1, :);
%% =================== update memory cr and F =============================
num_success_params = numel(goodCR);
if num_success_params > 0
    weightsDE = diff(I == 1)./ sum(diff(I == 1));
    %% for updating the memory of scaling factor
    archive_f(hist_pos) = (weightsDE * (goodF .^ 2))./ (weightsDE * goodF);
    
    %% for updating the memory of crossover rate
    if max(goodCR) == 0 || archive_Cr(hist_pos)  == -1
        archive_Cr(hist_pos)  = -1;
    else
        archive_Cr(hist_pos) = (weightsDE * (goodCR .^ 2)) / (weightsDE * goodCR);
    end
    
    hist_pos= hist_pos+1;
    if hist_pos > memory_size;  hist_pos = 1; end
    
    %% for updating the memory of T
    archive_T(hist_pos) = (weightsDE * (goodT .^ 2)) ./ (weightsDE * goodT);
    
    %% for updating the memory of freq
    if max(goodFreq) == 0 || archive_freq(hist_pos)  == -1
       archive_freq(hist_pos)  = -1;
    else
       archive_freq(hist_pos) = (weightsDE * (goodFreq .^ 2)) / (weightsDE * goodFreq);
    end
end

%% sort new x, fitness
[fitx, ind]=sort(fitx);
x=x(ind,:);xold = xold(ind,:);

%% record the best value after checking its feasiblity status
if fitx(1)<bestold  && min(x(ind(1),:))>=-100 && max(x(ind(1),:))<=100
    bestold=fitx(1);
     bestx= x(1,:);
end
%% check to print
if Printing==1
    res_det= [res_det repmat(bestold,1,PopSize)];
end

