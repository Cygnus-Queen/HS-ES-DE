function MLCC()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLCC:  Sheng Xin Zhang, Li Ming Zheng, Kit Sang Tang, Shao Yong Zheng, Wing Shing Chan,
%        Multi-Layer Competitive-Cooperative Framework for Performance Enhancement of Differential Evolution,
%        Information Sciences, https://doi.org/10.1016/j.ins.2018.12.065
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
format long;
format compact;
currentPath = cd;
cauchyPath = sprintf('%s/cauchy', currentPath);
addpath(genpath(cauchyPath));
rand('seed', sum(100 * clock));
val_2_reach = 10^(-8);
fbis = [100:100:3000];

D = 30;
times = 10;%老子翻了10倍
totalTime = 30;
NP = 150;
max_nfes = 10000 * D * times;
lb = ones(D,1).*(-100);
ub = ones(D,1).*(100);
genMax = times*D*10000/NP;
T = 1000*D/NP;
GT = 5*T;
X_init = ...
    repmat(lb, 1, NP) + ...
    repmat(ub - lb, 1, NP) .* lhsdesign(NP, D, 'iteration', 100)';

for func = [10:30,1:9]
    curbestfit = realmax();%%%%%%%%%%%%%%%%%%our parameters
    prebestfit = realmax();
    curdiv = 1.0;
    count = 0;
    sign = 0;
    everbestfitall = realmax();
    everbestfitinrun = realmax();%%%%%%%%%%%%%%%%%%our parameters
    MLCCDE_outcome = [];
    func = func;
    filename = strcat(strcat('f',num2str(func)),'_original.txt');
    fp = fopen(filename,'a+');
    % %   parfor irun = 1:totalTime
    for irun = 1:totalTime
        X = X_init;
        A = X;
        X = X';
        fx = cec14_func(X', func)';
        FES = NP;
        xxx = 0;
        gen = 1;
        Accumulation = 0;
        stage = 0;
        k = 1;
        A_size = 0;
        Chy = cauchyrnd(0, 0.1, NP + 10);
        H = NP;
        MF = 0.7 * ones(H, 1);
        MCR = 0.5 * ones(H, 1);
        sucIdx = zeros(2,NP);
        Prefer_M = ceil(rand(1,NP)*2);
        
        everbestfitall = min(fx)-fbis(func);
        curdiv = 0.0;
        for x =1 : D
            midpoint(x) = median(X(:,x));
        end
        distobest = 1 : NP;
        for x = 1: NP
            distobest (x)= 0;
            for y = 1 : D
                distobest(x) = distobest(x) + abs((X(x,y) - midpoint(y))/(ub(y) - lb(y)));
            end
            distobest (x) = distobest (x) / D;
            curdiv = curdiv + distobest (x);
        end
        curdiv = curdiv / NP;
        %fprintf(fp,'%d %e %e %e\r\n', FES, min(mixVal),curdiv,everbestfitall);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fx)-fbis(func), min(fx)-fbis(func), everbestfitall);
        xxx = 0;
        %everbestfitall = realmax();
        everbestfitinrun = realmax();
        noprogressgen = 0;
        %curbestfit = valParents(1);
        %prebestfit = realmax();	
        %curbestchrom = popold(1, :);        
        while FES < max_nfes
            
            S_CR = zeros(1, NP);
            S_F = zeros(1, NP);
            S_df = zeros(1, NP);
            [~, fidx] = sort(fx);
            rank = zeros(1,NP);
            rank(fidx) = 1:NP;
            gen = gen + 1;
            if gen > 2
                for i = 1 : NP
                    if sucIdx(2,i)==0
                        if sucIdx(1,i)==1
                            Prefer_M(i)=2;
                        elseif sucIdx(1,i)==2
                            Prefer_M(i)=1;
                        end
                    else
                        Prefer_M(i)=sucIdx(1,i);
                    end
                end
            end
            mbest_= ceil(rand*NP*0.05);
            MU = find(rank<=mbest_);
            shade_Idx = unique([find(Prefer_M==1),MU]);
            ide_Idx = unique([find(Prefer_M==2),MU]);
            %%%%%%%%%
            prebestfit = curbestfit;
            %prebestchrom = curbestchrom;
            %preind = ind;
            [curbestfit,ind] = min(fx);
            %curbestchrom = popold(ind, :);
            %           [curbestfit, index] = min(mixVal);
            %noprogressFES = noprogressFES + mixPopSizexx;
            if prebestfit < everbestfitinrun
                everbestfitinrun = prebestfit;
                %everbestchrominrun = prebestchrom;
            end
            
            if prebestfit < everbestfitall
                everbestfitall = prebestfit;
            end
            %%%%%%%
            [U1,F,CR] = shade(X,fx,MF,MCR,A_size,A,lb,ub,shade_Idx);
            U2 = ide(X,fx,stage,gen,lb,ub,genMax,ide_Idx);
            
            fu1 = fx;  fu2 = fx; fx_old = fx;
            mark = zeros(1,NP);
            nS = 0; e = zeros(1,NP);
            
            for i = 1 : NP
                if rank(i)<=mbest_
                    if FES >= max_nfes
                        break;
                    end
                    
                    fu1(i) = cec14_func(U1(i, :)',func);
                    fu2(i) = cec14_func(U2(i, :)',func);
                    FES = FES + 2;
                    [~,Prefer_M(i)] = min([fu1(i),fu2(i)]);
                    if fu1(i)<fx(i)
                        FLAG1 = 1;
                    else
                        FLAG1 = 0;
                    end
                    if fu2(i)<fx(i)
                        FLAG2 = 1;
                    else
                        FLAG2 = 0;
                    end
                else
                    if FES >= max_nfes
                        break;
                    end
                    FLAG1=0;FLAG2=0;
                    
                    if Prefer_M(i)==1
                        fu1(i) = cec14_func(U1(i, :)',func);
                    elseif Prefer_M(i)==2
                        fu2(i) = cec14_func(U2(i, :)',func);
                    end
                    FES = FES + 1;
                end
                if (FLAG1 == 1 && Prefer_M(i)~=1)
                    mark(i) = 1;
                    if fu1(i) < fx(i)
                        nS = nS + 1;
                        S_df(nS)	= abs(fu1(i) - fx(i));
                        S_CR(nS)	= CR(i);
                        S_F(nS)		= F(i);
                    end
                end
                %         if  (FLAG2 == 1 && Prefer_M(i)~=2)
                %             mark(i) = 2;
                %             if fu2(i) < fx(i)
                %             end
                %         end
                if Prefer_M(i) == 1
                    mark(i) = 1;
                    e(i) = fu1(i);
                    if fu1(i) < fx(i)
                        sucIdx(1,i) = 1;
                        sucIdx(2,i) = 1;
                        X(i,:) = U1(i,:);
                        nS = nS + 1;
                        S_df(nS)	= abs(fu1(i) - fx(i));
                        fx(i) = fu1(i);
                        S_CR(nS)	= CR(i);
                        S_F(nS)		= F(i);
                        if A_size < NP
                            A_size = A_size + 1;
                            A(:,A_size) = X(i, :)';
                        else
                            ri = floor(1 + NP * rand);
                            A(:, ri) = X(i, :)';
                        end
                    else
                        sucIdx(2,i) = 0;
                    end
                elseif Prefer_M(i) == 2
                    mark(i) = 2;
                    e(i) = fu2(i);
                    if fu2(i) < fx(i)
                        sucIdx(1,i) = 2;
                        sucIdx(2,i) = 1;
                        X(i,:) = U2(i,:);
                        fx(i) = fu2(i);
                        
                        if A_size < NP
                            A_size = A_size + 1;
                            A(:,A_size) = X(i, :)';
                        else
                            ri = floor(1 + NP * rand);
                            A(:, ri) = X(i, :)';
                        end
                    else
                        sucIdx(2,i) = 0;
                    end
                end
            end   %end of NP;
            
            SC = 0;
            for i = 1:NP
                if(fx_old(i)-e(i) > 1e-8)
                    SC = SC + 1;
                end
            end
            SR = SC/NP;
            
            if gen < GT
                if SR <= 0
                    Accumulation = Accumulation + 1;
                else
                    Accumulation = 0;
                end
                if(Accumulation >= T)
                    if(stage==0)
                        stage = 1;
                    end
                end
            else
                if SR <= 0.1
                    Accumulation = Accumulation + 1;
                else
                    Accumulation = 0;
                end
                if(Accumulation >= T)
                    if(stage==0)
                        stage = 1;
                    end
                end
            end
            
            
            if nS > 0
                w = S_df(1 : nS) ./ sum(S_df(1 : nS));
                MCR(k) = sum(w .* S_CR(1 : nS));
                MF(k) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
                k = k + 1;
                if k > H
                    k = 1;
                end
            end
            if min(fx)-fbis(func) < everbestfitall 
                everbestfitall = min(fx)-fbis(func);
            end
            %everbestfitall = min(fx)-fbis(func);
            if (FES/(D * 10000 * times)) >= (xxx / 100)
                if xxx > 0
                    curdiv = 0.0;
                    for x =1 : D
                        midpoint(x) = median(X(:,x));
                    end
                    distobest = 1 : NP;
                    for x = 1: NP
                        distobest (x)= 0;
                        for y = 1 : D
                            distobest(x) = distobest(x) + abs((X(x,y) - midpoint(y))/(ub(y) - lb(y)));
                        end
                        distobest (x) = distobest (x) / D;
                        curdiv = curdiv + distobest (x);
                    end
                    curdiv = curdiv / NP;
                    fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fx)-fbis(func), min(fx)-fbis(func), everbestfitall);
                end
                xxx = xxx + 1;
            end
            
        end
        minfx = min(fx)-fbis(func);
        %if minfx <= val_2_reach
        %    minfx = 0;
        %end
        MLCCDE_outcome = [MLCCDE_outcome,minfx];
        
        %fprintf('MLCC_%d th run, best-so-far error value = %1.8e\n', irun , minfx)
    end
    outcome = MLCCDE_outcome;
    %fprintf('\n')
    fprintf(fp,'%1.3e %1.3e:\n', mean(MLCCDE_outcome), std(MLCCDE_outcome));
    for x = 1 : totalTime
        fprintf(fp,'%s ', num2str(outcome(x)));
    end
    fclose(fp);
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [U,F,CR] = shade(X,fx,MF,MCR,A_size,A,lb,ub,shade_Idx)
X = X';
[D,NP] = size(X);
H = NP;
[~, fidx] = sort(fx);
CR = zeros(1, NP);
r = zeros(1, NP);
p = zeros(1, NP);

for i = shade_Idx
    r(i) = floor(1 + H * rand);
    CR(i) = MCR(r(i)) + 0.1 * randn;
end

CR(CR > 1) = 1;
CR(CR < 0) = 0;

F = zeros(1, NP);
for i = shade_Idx
    while F(i) <= 0
        F(i) = MF((i)) + 0.1 * tan(pi * (rand - 0.5));
    end
    
    if F(i) > 1
        F(i) = 1;
    end
end
pmin = 2 / NP;
for i = shade_Idx
    p(i) = pmin + rand * (0.2 - pmin);
end

XA = [X, A];
rt = zeros(1, NP);
r1 = zeros(1, NP);
r2 = zeros(1, NP);

for i = shade_Idx
    rt(i) = i;
    
    r1(i) = floor(1 + NP * rand);
    while rt(i) == r1(i)
        r1(i) = floor(1 + NP * rand);
    end
    
    r2(i) = floor(1 + (NP + A_size) * rand);
    while rt(i) == r2(i) || r1(i) == r2(i)
        r2(i) = floor(1 + (NP + A_size) * rand);
    end
end

V = X;
for i = shade_Idx
    pbest = fidx(floor(1 + round(p(rt(i)) * NP) * rand));
    
    V(:, i) = X(:, rt(i)) + F(rt(i)) .* (X(:, pbest) - X(:, rt(i))) ...%的确是 JADE的变异方式
        + F(rt(i)) .* (X(:, r1(i)) - XA(:, r2(i)));
end
U = X;
for i = shade_Idx
    jrand = floor(1 + D * rand);
    for j = 1 : D
        if rand < CR(i) || j == jrand
            U(j, i) = V(j, i);
        else
            U(j, i) = X(j, rt(i));
        end
    end
end

for i = shade_Idx
    for j = 1 : D
        if U(j, i) < lb(j)
            U(j, i) = 0.5 * (lb(j) + X(j, rt(i)));
        elseif U(j, i) > ub(j)
            U(j, i) = 0.5 * (ub(j) + X(j, rt(i)));
        end
    end
end
U = U';
end

function U = ide(X,fx,stage,gen,lb,ub,genMax,ide_Idx)
X = X';
[D,NP] = size(X);
[~, index] = sort(fx);
rank = zeros(1,NP);
rank(index) = 1:NP;
V = X;
U = X;
r1 = zeros(1,NP);
r2 = zeros(1,NP);
r3 = zeros(1,NP);
rt = zeros(1,NP);
for i = ide_Idx
    rt(i) = i;
    
    r1(i) = floor(1 + NP * rand);
    while rt(i) == r1(i)
        r1(i) = floor(1 + NP * rand);
    end
    
    r2(i) = floor(1 + NP * rand);
    while rt(i) == r2(i) || r1(i) == r2(i)
        r2(i) = floor(1 + NP * rand);
    end
    
    r3(i) = floor(1 + NP * rand);
    while rt(i) == r3(i) || r1(i) == r3(i) || r2(i) == r3(i)
        r3(i) = floor(1 + NP * rand);
    end
end

ps = 0.1 + 0.9 * 10.^(5*(gen/genMax-1));


if stage==0
    o = [1:NP];
else
    o = r1;
end

Fo = zeros(1,NP);
for i = ide_Idx
    Fo(i) = rank(o(i))/NP + 0.1 * randn;
    if Fo(i) < 0 || Fo(i) > 1
        Fo(i) = rank(o(i))/NP;
    end
end

CR = zeros(1,NP);
for i = ide_Idx
    CR(i) = rank(i)/NP + 0.1 * randn;
    if CR(i) < 0 || CR(i) > 1
        CR(i) = rank(i)/NP;
    end
end

pd = 0.1 * ps;
for i = ide_Idx
    d = X(:,r3(i));
    for j = 1 : D
        if rand < pd
            d(j,1) = lb(j) + rand*(ub(j) - lb(j));
        end
    end
    if rank(i)/NP < ps
        V(:,i) = X(:, o(i)) + Fo(i) * (X(:, r1(i)) - X(:, o(i))) + ...
            Fo(i)  * (X(:, r2(i)) - d);
    else
        V(:,i) = X(:, o(i)) + Fo(i) * (X(:, index(ceil(rand*ps*NP))) - X(:, o(i))) + ...
            Fo(i)  * (X(:, r2(i)) - d);
    end
end

for i = ide_Idx
    jrand = floor(1 + D * rand);
    for j = 1 : D
        if rand < CR(i) || j == jrand
            U(j, i) = V(j, i);
        else
            U(j, i) = X(j, i);
        end
    end
end

for i = ide_Idx
    for j = 1 : D
        if U(j, i) < lb(j)
            U(j, i) = lb(j) + rand*(ub(j) - lb(j));
        elseif U(j, i) > ub(j)
            U(j, i) = lb(j) + rand*(ub(j) - lb(j));
        end
    end
end
U = U';
end
