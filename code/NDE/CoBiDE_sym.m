%for all algorithms: 10倍评估次数，500&2检测提高，1f5提高率
clc;
clear all;

format long
format compact
addpath data
times = 10;
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
for problem_index = [17:30,1:16]
    
    tic;
    
    problem = problem_set(problem_index);
    filename = strcat(strcat('f',num2str(problem)),'_8f3_09.txt');
    fp = fopen(filename,'a+');
    pb = .4;
    ps = .5;
    SEL = round(ps*NP);
    thediv = 8.0e-03;%只有这个参数变化最大
    % load data in CEC2014
    lu = [-100 * ones(1, n); 100 * ones(1, n)];
    LB = repmat(lu(1,:), NP, 1);
    UB = repmat(lu(2,:), NP, 1);
    
    TIMES = 30;
    Solve = [];
    
    for time = 1:TIMES
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
        X = LB + rand(NP,n).*(UB-LB);          
        U = X;
        %U_init_record = U;
        U_sym = U;
        curbestchrom = X(1, :);%%%%%%%%%%%%%%%%%%our parameters
        %%
        fit=cec14_func(X',problem)'-100*problem;
        everbestfitall = min(fit);

        FES = NP;
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
        disp(curdiv);
        fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fit), min(fit), everbestfitall);     

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
            %%%%%%%%%%%%%%%%%%%
            prebestfit = curbestfit;
            [curbestfit,ind] = min(fit);
            if prebestfit < everbestfitinrun
                everbestfitinrun = prebestfit;
            end
            
            if prebestfit < everbestfitall
                everbestfitall = prebestfit;
            end
            if sign == 0
                if (curbestfit >= everbestfitinrun) ||  (curbestfit < everbestfitinrun && (everbestfitinrun - curbestfit) / everbestfitinrun < 1e-5)
                    mycount = mycount +1;
                else
                    mycount = 0;
                end
                upcount =500;
                if abs(everbestfitinrun - everbestfitall) < 1e-10
                    upcount = upcount * 2;
                end
                if mycount >= upcount
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
                    %if curdiv > thediv 
                        sign = 1;%全体移动
                    %else
                        %sign = 1;%%%%%%%%%%%sign = 1;%表示这一代开始回退
                    %end
                    mycount = 0;
                    everbestfitinrun = realmax();
                end
            else
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
            %%%%%%%%%%%%%%%%%%%
            %% rand/1 mutation strategy
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
            if sign ==0 
                fit_ = cec14_func(U',problem)'-100*problem;
            end
            if sign == 1 && curdiv > thediv %&& Issym == 1%%%%%%%
                U_sym = LB + UB - U;
                S = rand(NP,1) < 0.9;
                U(S,:) = U_sym(S,:);
                fit_ = cec14_func(U',problem)'-100*problem;
            end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% selection
            if sign == 0%正常选择
                S = fit_ <= fit;
                X(S,:) = U(S,:);
                fit(S) = fit_(S);
                
                %%
                F(~S) = oldF(~S);
                CR(~S) = oldCR(~S);
            else
                for i = 1 : NP
                    S(i) = 1;
                end
                X(S,:) = U(S,:);
                fit(S) = fit_(S);
                %%
                F = oldF;
                CR = oldCR;         
            end                                
            %% selection end
            if min(fit) < everbestfitall 
                everbestfitall = min(fit);
            end

            if sign ==0 
                FES = FES + NP; 
            end
            if sign == 1 && curdiv > thediv%%%%%%%%%%%%%%%%%%%%%%                    
                FES = FES + NP;                    
            end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

            if rem(FES,NP * 5000 * times / 100) == 0
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
                fprintf(fp,'%d %e %e %e %e\r\n', FES, curdiv, mean(fit), min(fit), everbestfitall);                
            end
            if sign == 1 && curdiv > thediv%mycount > 3000 + upcount
                sign = 0;%表示这一代开始前进
            end
        end
        disp('end');
        
        Solve = [Solve everbestfitall];
        
    end
    disp(['F',num2str(problem),...
        ' | ',num2str(mean(Solve)),...
        ' | ',num2str(std(Solve))]);
    %     disp(sort(Solve));
    fprintf(fp,'%e %e: ', mean(Solve), std(Solve));
    for x = 1 : TIMES
        fprintf(fp,'%s ', num2str(Solve(x)));
    end
    fclose(fp);
    toc;
    
end


