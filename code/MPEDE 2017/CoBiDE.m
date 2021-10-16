%**********************************************************************************
% Author: Yong Wang and Long Li
% Last Edited: 21/11/2013
% Email: ywang@csu.edu.cn

%************************************ Reference ***********************************
% * Y. Wang, H.-X. Li, T. Huang, and L. Li. Differential evolution based on
% * covariance matrix learning and bimodal distribution parameter setting. 
% * Applied Soft Computing, vol. 18, pp. 232-247, 2014.
%**********************************************************************************

clc;
clear all;

format long
format compact
addpath data

% dimensionality
n  = 30;
% population size
NP = 60;

%%

disp([' ------------------------------------------------------------------------------------------------- ']);
disp(['                     Differential Evolution Based on Covariance Matrix Learning                    ']);
disp(['                            and Bimodal Distribution Parameter Setting                             ']);
disp(['                                             base                                                   ']);
disp([' ------------------------------------------------------------------------------------------------- ']);
problem_set = [1 : 30];
for problem_index = [1 :30]
    
    tic;
    
    problem = problem_set(problem_index);
    filename = strcat(strcat('f',num2str(problem)),'_10_14.txt');
    fp = fopen(filename,'a+');
    pb = .4;
    ps = .5;
    SEL = round(ps*NP);
    
    % load data in CEC2005
    lu = [-100 * ones(1, n); 100 * ones(1, n)];
    LB = repmat(lu(1,:), NP, 1);
    UB = repmat(lu(2,:), NP, 1);
    
    TIMES = 30;
    Solve = [];
    
    for time = 1:TIMES
        curbestfit = realmax();
        curbestfitr = realmax();
        prebestfit = realmax();
        prebestfitr = realmax();
        recfit = realmax();
        curdiv = 1.0;
        prediv = 1.0;
        preavefit = realmax();
        curavefit = realmax();
        recdiv = 1.0;
        pro1 = 0.0;
        pro2 = 0.0;
        pro3 = 0.0;%diversity pro
        pro5 = 1.0;
        dp = 0.0;
        dp2 = 0.0;
        dp3 = 1.0;
        gap = 1;
        gap2 = 1;
		flag = false;
        dis = 1 : NP;
        distobest = 1 : NP;
        dis_ = 1 : NP; 
        t1 = 0;
        t2 = 0;
        count = 0;
        X = LB + rand(NP,n).*(UB-LB);
        U = X;
        %%%%%%
        times = 10;
        sign = 0;
        countloop = 0;
        discount = 0;
        insertsign = 0;
        %%%%%%
        X = LB + rand(NP,n).*(UB-LB);
        U = X;
        
        %%
        fit=cec14_func(X',problem)'-100*problem;
        FES = NP;        
        %
        igen =1;
        recordgen = 1;
        SavedX = X(1, :);
        Savedfit = fit(1);         
        %% Initialization for F
        F = zeros(NP, 1);
        sel = rand(NP, 1)<.5;
        num = sum(sel);
        F(sel)  = .65 + .1*tan(pi*(rand(num, 1) - .5));
        F(~sel) = 1 + .1*tan(pi*(rand(NP-num, 1) - .5));
        F  = min(1, F);
        BLF = F<=0;
        while BLF(BLF)
            F(BLF&sel)  = .65 + .1*tan(pi*(rand(sum(BLF&sel), 1) - .5));
            F(BLF&~sel) = 1 + .1*tan(pi*(rand(sum(BLF&~sel), 1) - .5));
            F(BLF)  = min(1, F(BLF));
            BLF_ = F(BLF)<=0;
            BLF(BLF) = BLF_;
        end
        
        %% Initialization for CR
        CR = zeros(NP, 1);
        sel = rand(NP, 1)<.5;
        num = sum(sel);
        CR(sel) = .1 + .1*tan(pi*(rand(num, 1) - .5));
        CR(~sel)= .95 + .1*tan(pi*(rand(NP-num, 1) - .5));
        BLC = CR<=0; BUC = CR>=1;
        CR(BLC) = 0.1; CR(BUC) = 1;
        
        S = true(NP, 1);
        KK=0;
        %%
        div=0.1;
        while FES < NP * 5000 * times
            
            %% random numbers
            [r1, r2, r3] = select3rands_(NP);
        
            %% Bimodal Distributions for F
            oldF = F;
            F = zeros(NP, 1);
            sel = rand(NP, 1)<.5;
            num = sum(sel);
            F(sel)  = .65 + .1*tan(pi*(rand(num, 1) - .5));
            F(~sel) = 1 + .1*tan(pi*(rand(NP-num,1) - .5));
            F  = min(1 ,F);
            BLF = F<=0;
            while BLF(BLF)
                F(BLF&sel)  = .65 + .1*tan(pi*(rand(sum(BLF&sel), 1) - .5));
                F(BLF&~sel) = 1 + .1*tan(pi*(rand(sum(BLF&~sel), 1) - .5));
                F(BLF)  = min(1 ,F(BLF));
                BLF_ = F(BLF)<=0;
                BLF(BLF) = BLF_;
            end
            chgF = rand(NP, 1) < 1;
            chgF = chgF&(~S);
            F = chgF.*F + (~chgF).*oldF;
            
            %% Bimodal Distributions for CR
            oldCR = CR;
            CR = zeros(NP, 1);
            sel = rand(NP, 1)<.5;
            num = sum(sel);
            CR(sel) = .1 + .1*tan(pi*(rand(num, 1) - .5));
            CR(~sel)= .95 + .1*tan(pi*(rand(NP-num,1) - .5));
            BLC = CR<=0; BUC = CR>=1;
            CR(BLC) = 0.1; CR(BUC) = 1;
            
            chgCR = rand(NP, 1) < 1;
            chgCR = chgCR&(~S);
            CR = chgCR.*CR + (~chgCR).*oldCR;
            
            %% rand/1 mutation strategy
            prebestfit = curbestfit;
            curbestfit = min(fit);

                V = X(r1,:) + F(:, ones(1,n)).*(X(r2,:)-X(r3,:));
            %%
%            if problem ~=7 && problem ~=25
                BL = V<LB; V(BL) = 2*LB(BL) - V(BL);
                BLU = V(BL)>UB(BL); BL(BL) = BLU; V(BL) = UB(BL);
                BU = V>UB; V(BU) = 2*UB(BU) - V(BU);
                BUL = V(BU)<LB(BU); BU(BU) = BUL; V(BU) = LB(BU);
%            end
            
            %% crossover
            J_= mod(floor(rand(NP, 1)*n), n) + 1;
            J = (J_-1)*NP + (1:NP)';
            crs = rand(NP, n) < CR(:, ones(1, n));
            %%
            if rand<pb
                %% coordinate ratation
                    [nouse, seq] = sort(fit);
                    Xsel = X(seq(1:SEL), :);
                    xmean = mean(Xsel);
                    % covariance matrix calculation
                    C =  1/(SEL-1)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
                    C = triu(C) + transpose(triu(C,1)); % enforce symmetry
                    [R,D] = eig(C);
                    % limit condition of C to 1e20 + 1
                    if max(diag(D)) > 1e20*min(diag(D))
                        tmp = max(diag(D))/1e20 - min(diag(D));
                        C = C + tmp*eye(n);
                        [R, D] = eig(C);
                    end
                TM = R;
                TM_=R';
                Xr = X*TM;
                Vr = V*TM;
                %% crossover according to the Eigen coordinate system
                Ur = Xr;
                Ur(J) = Vr(J);
                Ur(crs) = Vr(crs);
                %%
                U = Ur*TM_;
                %%
%                if problem ~=7 && problem ~=25
                    BL  = U<LB; U(BL) = 2*LB(BL) - U(BL);
                    BLU = U(BL)>UB(BL); BL(BL) = BLU; U(BL) = UB(BL);
                    BU  = U>UB; U(BU) = 2*UB(BU) - U(BU);
                    BUL = U(BU)<LB(BU); BU(BU) = BUL; U(BU) = LB(BU);
%                end
                
            else
                
                U = X;
                U(J) = V(J);
                U(crs) = V(crs);
                
            end
            
            fit_ = cec14_func(U',problem)'-100*problem;
            %% selection
                S = fit_ <= fit;
                X(S,:) = U(S,:);
                fit(S) = fit_(S);

                %%
                F(~S) = oldF(~S);
                CR(~S) = oldCR(~S);



               


            if rem(igen,50 * times) == 0
                curdiv = 0.0;   
                for x =1 : n
                   midpoint(x) = median(X(:,x));
                end
                distobest = 1 : NP;
                for x = 1: NP
                    distobest (x)= 0;
                    for y = 1 : n
                        distobest(x) = distobest(x) + abs((X(x,y) - midpoint(y))/(UB(y) - LB(y)));
                    end
                    distobest (x) = distobest (x) / n;
                    curdiv = curdiv + distobest (x);
                end                   
                curdiv = curdiv / NP;
            end     
            %% selection end
            igen = igen +1;
            FES = FES + NP;
            if rem(igen,50 * times) == 0
                disp(FES);
                disp(igen);
                fprintf(fp,'%d %g %g %g\r\n', igen, min(fit), curdiv, Savedfit);
            end
        end
        disp('end');

            Solve = [Solve min(fit)];
   
    end
    disp(['F',num2str(problem),...
        ' | ',num2str(mean(Solve)),...
        ' | ',num2str(std(Solve))]);
%     disp(sort(Solve));
    fprintf(fp,'%s %s: ', num2str(mean(Solve)), num2str(std(Solve)));
    for x = 1 : TIMES
        fprintf(fp,'%s ', num2str(Solve(x)));
    end
    fclose(fp);    
    toc;
    
end


