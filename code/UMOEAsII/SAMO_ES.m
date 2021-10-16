%% ============United Multi-Operator Evolutionary AlgorithmsII ============
% Should you have any queries, please contact
% Dr. Saber Elsayed. University of New South Wales at Canberra
% s.elsayed@adfa.edu.au
% www.saberelsayd.net or
% https://sites.google.com/site/saberelsayed3/home 
% =========================================================================

%% Parts of this code were taken from  (https://www.lri.fr/~hansen/cmaes_inmatlab.html#matlab).

function[ x, fitx, setting,bestold,bestx,bnd,fitness,current_eval,res_det] = ...
    SAMO_ES( x, ~, setting, iter,bestold,bestx,fitness,bnd,xmin,xmax,n,PopSize,current_eval,I_fno,res_det,Printing,Max_FES)

stopOnWarnings=0;
noiseReevals = 0;
fitness.raw = NaN(1, PopSize + noiseReevals);
fitness.raw(PopSize + find(isnan(fitness.raw(1:noiseReevals)))) = NaN;

arz = randn(n,PopSize);
arx = repmat(setting.xmean, 1, PopSize) + setting.sigma * (setting.BD * arz);

%% ignore handling the boundaries constraints during the first 50% evolutionary process
%%-this is based on our earlier analysis carried out on UMOEAs in 2014.
handle_limit=0.5;
if current_eval >=handle_limit*Max_FES
    arxvalid =han_boun(arx', xmax, xmin, x,PopSize,2);
    arxvalid=arxvalid';
else
    arxvalid=arx;
end
%% evaluate and update cfe
fitness.raw = cec14_func(arxvalid,I_fno);
current_eval=current_eval+PopSize; %% increase the fitness evaluations

fitness.sel= fitness.raw ;
[fitness.sel, fitness.idxsel] = sort(fitness.sel);

fitness.raw= fitness.raw(fitness.idxsel);
arxvalid= arxvalid(:, fitness.idxsel);
arx= arx(:, fitness.idxsel);
arz=arz(:, fitness.idxsel);
[~,pos_ro]=min(fitness.raw);

%% record the best value after checking its feasiblity status
if fitness.raw(pos_ro) < bestold && (min(arxvalid(:,pos_ro)))>=-100 && (max(arxvalid(:,pos_ro)))<=100
    bestold=fitness.raw(pos_ro);
    bestx= arxvalid(:,pos_ro)';
end
if Printing==1
    res_det= [res_det repmat(bestold,1,PopSize)];
end


% Calculate new setting.xmean, this is selection and recombination
setting.xold = setting.xmean; % for speed up of Eq. (2) and (3)
cmean =1;% 1/min(max((PopSize-1*n)/2, 1), n);  % == 1/kappa
setting.xmean = (1-cmean) * setting.xold + cmean * arx(:,(1:setting.mu))*setting.weights;
if  current_eval >=handle_limit*Max_FES
    % setting.xmean = xintobounds(setting.xmean, xmin', xmax');
    setting.xmean =han_boun(setting.xmean', xmax, xmin, x(1,:),1,2);
    setting.xmean=setting.xmean';
end
zmean = arz(:,(1:setting.mu))*setting.weights;%==D^-1*setting.B'*(setting.xmean-setting.xold)/setting.sigma
% Cumulation: update evolution paths
setting.ps = (1-setting.cs)*setting.ps + sqrt(setting.cs*(2-setting.cs)*setting.mueff) * (setting.B*zmean);          % Eq. (4)
hsig = norm(setting.ps)/sqrt(1-(1-setting.cs)^(2*iter))/setting.chiN < 1.4 + 2/(n+1);

setting.pc = (1-setting.cc)*setting.pc ...
    + hsig*(sqrt(setting.cc*(2-setting.cc)*setting.mueff)/setting.sigma/cmean) * (setting.xmean-setting.xold);     % Eq. (2)
if hsig == 0
    % disp([num2str(iter) ' ' num2str(counteval) ' setting.pc update stalled']);
end
% Adapt covariance matrix
neg.ccov = 0;  % TODO: move parameter setting upwards at some point
if setting.ccov1 + setting.ccovmu > 0                                                    % Eq. (3)
    if setting.flgDiagonalOnly % internal linear(?) complexity
        setting.diagC = (1-setting.ccov1_sep-setting.ccovmu_sep+(1-hsig)*setting.ccov1_sep*setting.cc*(2-setting.cc)) * setting.diagC ... % regard old matrix
            + setting.ccov1_sep * setting.pc.^2 ...               % plus rank one update
            + setting.ccovmu_sep ...                      % plus rank setting.mu update
            * (setting.diagC .* (arz(:,(1:setting.mu)).^2 * setting.weights));
        %             * (repmat(setting.diagC,1,setting.mu) .* arz(:,(1:setting.mu)).^2 * setting.weights);
        setting.diagD = sqrt(setting.diagC); % replaces eig(setting.C)
    else
        arpos = (arx(:,(1:setting.mu))-repmat(setting.xold,1,setting.mu)) / setting.sigma;
        setting.C = (1-setting.ccov1-setting.ccovmu+(1-hsig)*setting.ccov1*setting.cc*(2-setting.cc)) * setting.C ... % regard old matrix
            + setting.ccov1 * setting.pc*setting.pc' ...     % plus rank one update
            + setting.ccovmu ...             % plus rank setting.mu update
            * arpos * (repmat(setting.weights,1,n) .* arpos');
        % is now O(setting.mu*n^2 + setting.mu*n), was O(setting.mu*n^2 + setting.mu^2*n) when using diag(setting.weights)
        %   for setting.mu=30*n it is now 10 times faster, overall 3 times faster
        
        setting.diagC = diag(setting.C);
    end
end


% Adapt setting.sigma
setting.sigma = setting.sigma * exp(min(1, (sqrt(sum(setting.ps.^2))/setting.chiN - 1) * setting.cs/setting.damps));             % Eq. (5)
% disp([iter norm(setting.ps)/setting.chiN]);

if 11 < 3   % testing with optimal step-size
    setting.sigma = 0.04 * setting.mueff * sqrt(sum(setting.xmean.^2)) / n; % 20D,lam=1000:25e3
    setting.sigma = 0.3 * setting.mueff * sqrt(sum(setting.xmean.^2)) / n; % 20D,lam=(40,1000):17e3
    %      75e3 with def (1.5)
    %      35e3 with setting.damps=0.25
end
if 11 < 3
    
    setting.xmean = ones(n,1);
end

% Update setting.B and D from setting.C

if ~setting.flgDiagonalOnly && (setting.ccov1+setting.ccovmu+neg.ccov) > 0 && mod(iter, 1/(setting.ccov1+setting.ccovmu+neg.ccov)/n/10) < 1
    setting.C=triu(setting.C)+triu(setting.C,1)'; % enforce symmetry to prevent complex numbers
    [setting.B,tmp] = eig(setting.C);     % eigen decomposition, setting.B==normalized eigenvectors
    % effort: approx. 15*n matrix-vector multiplications
    setting.diagD = diag(tmp);
    
    % limit condition of setting.C to 1e14 + 1
    if min(setting.diagD) <= 0
        
        setting.diagD(setting.diagD<0) = 0;
        tmp = max(setting.diagD)/1e14;
        setting.C = setting.C + tmp*eye(n,n); setting.diagD = setting.diagD + tmp*ones(n,1);
        
    end
    if max(setting.diagD) > 1e14*min(setting.diagD)
        
        tmp = max(setting.diagD)/1e14 - min(setting.diagD);
        setting.C = setting.C + tmp*eye(n,n); setting.diagD = setting.diagD + tmp*ones(n,1);
        
    end
    
    setting.diagC = diag(setting.C);
    setting.diagD = sqrt(setting.diagD); % D contains standard deviations now
    % setting.diagD = setting.diagD / prod(setting.diagD)^(1/n);  setting.C = setting.C / prod(setting.diagD)^(2/n);
    setting.BD = setting.B.*repmat(setting.diagD',n,1); % O(n^2)
end % if mod

% Align/rescale order of magnitude of scales of setting.sigma and setting.C for nicer output
% TODO: interference with sigmafacup: replace 1e10 with 2*sigmafacup
% not a very usual case
if 1 < 2 && setting.sigma > 1e10*max(setting.diagD) && setting.sigma > 8e14 * max(setting.insigma)
    fac = setting.sigma; % / max(setting.diagD);
    setting.sigma = setting.sigma/fac;
    setting.pc = fac * setting.pc;
    setting.diagD = fac * setting.diagD;
    if ~setting.flgDiagonalOnly
        setting.C = fac^2 * setting.C; % disp(fac);
        setting.BD = setting.B .* repmat(setting.diagD',n,1); % O(n^2), but repmat might be inefficient todo?
    end
    setting.diagC = fac^2 * setting.diagC;
end

if setting.flgDiagonalOnly > 1 && iter > setting.flgDiagonalOnly
    % full covariance matrix from now on
    setting.flgDiagonalOnly = 0;
    setting.B = eye(n,n);
    setting.BD = diag(setting.diagD);
    setting.C = diag(setting.diagC); % is better, because correlations are spurious anyway
end

% ----- numerical error management -----
% Adjust maximal coordinate axis deviations
if any(setting.sigma*sqrt(setting.diagC) > setting.maxdx')
    setting.sigma = min(setting.maxdx ./ sqrt(setting.diagC'));
    %warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
    %         'deviation at upper limit of ' num2str(setting.maxdx)]);
    % stopflag(end+1) = {'maxcoorddev'};
end
% Adjust minimal coordinate axis deviations
if any(setting.sigma*sqrt(setting.diagC) < setting.mindx)
    setting.sigma = max(setting.mindx ./ sqrt(setting.diagC)) * exp(0.05+setting.cs/setting.damps);
    %warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
    %         'deviation at lower limit of ' num2str(setting.mindx)]);
    % stopflag(end+1) = {'mincoorddev'};;
end
% Adjust too low coordinate axis deviations
if any(setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC))
    if stopOnWarnings
        %         stopflag(end+1) = {'warnnoeffectcoord'};
    else
        %         warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
        %             'deviation too low' ]);
        if setting.flgDiagonalOnly
            setting.diagC = setting.diagC + (setting.ccov1_sep+setting.ccovmu_sep) * (setting.diagC .* ...
                (setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC)));
        else
            setting.C = setting.C + (setting.ccov1+setting.ccovmu) * diag(setting.diagC .* ...
                (setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC)));
        end
        setting.sigma = setting.sigma * exp(0.05+setting.cs/setting.damps);
    end
end
% Adjust step size in case of (numerical) precision problem
if setting.flgDiagonalOnly
    tmp = 0.1*setting.sigma*setting.diagD;
else
    tmp = 0.1*setting.sigma*setting.BD(:,1+floor(mod(iter,n)));
end
if all(setting.xmean == setting.xmean + tmp)
    %     ii = 1+floor(mod(iter,n));
    if stopOnWarnings
    else
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end
% Adjust step size in case of equal function values (flat fitness)
% isequalfuncvalues = 0;
if fitness.sel(1) == fitness.sel(1+ceil(0.1+PopSize/4))
    % isequalfuncvalues = 1;
    if setting.stopOnEqualFunctionValues
        setting.arrEqualFunvals = [iter setting.arrEqualFunvals(1:end-1)];
        % stop if this happens in more than 33%
        if setting.arrEqualFunvals(end) > iter - 3 * length(setting.arrEqualFunvals)
            %             stopflag(end+1) = {'equalfunvals'};
        end
    else
        if flgWarnOnEqualFunctionValues
            %             warning(['Iteration ' num2str(iter) ...
            %                 ': equal function values f=' num2str(fitness.sel(1)) ...
            %                 ' at maximal main axis setting.sigma ' ...
            %                 num2str(setting.sigma*max(setting.diagD))]);
        end
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end
% Adjust step size in case of equal function values
if iter > 2 && myrange([fitness.hist fitness.sel(1)]) == 0
    if stopOnWarnings
        % 	stopflag(end+1) = {'warnequalfunvalhist'};
    else
        %         warning(['Iteration ' num2str(iter) ...
        %             ': equal function values in history at maximal main ' ...
        %             'axis setting.sigma ' num2str(setting.sigma*max(setting.diagD))]);
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end

%% print out final results
x= arxvalid';
fitx= fitness.raw;

function res=myrange(x)
res = max(x) - min(x);



