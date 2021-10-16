%% ============United Multi-Operator Evolutionary AlgorithmsII ============
% Should you have any queries, please contact
% Dr. Saber Elsayed. University of New South Wales at Canberra
% s.elsayed@adfa.edu.au
% www.saberelsayd.net or
% https://sites.google.com/site/saberelsayed3/home 
% =========================================================================

function [setting,fitness,bnd]= init_cmaes_par(setting,fitness,bnd,EA_obj2, EA_2, n, n2, xmin, xmax)
% Initialize dynamic internal state parameters

setting.flgActiveCMA=0;
setting.flgresume=0;
setting.flgDiagonalOnly=0;
%% remember that EA_2 was uniformly generated between xmin and xmax. 
%% So, mean(EA_2) does not violate the initialization condition in the competition
setting.xmean = mean(EA_2);
setting.xmean=setting.xmean';
setting.insigma=0.3;
setting.sigma = setting.insigma;
setting.maxdx = Inf;
setting.mindx = 0;

setting.sigma = max(setting.insigma);              % overall standard deviation
setting.pc = zeros(n,1); setting.ps = zeros(n,1);  % evolution paths for setting.C and setting.sigma


setting.idx = (xmin > -Inf) | (xmax < Inf);
if length(setting.idx) == 1
    setting.idx = setting.idx * ones(n,1);
end
if length(setting.insigma) == 1
    setting.insigma = setting.insigma * ones(n,1) ;
end
setting.diagD = setting.insigma/max(setting.insigma);      % diagonal matrix D defines the scaling
setting.diagC = setting.diagD.^2;
if setting.flgDiagonalOnly ~= 1            % use at some point full covariance matrix
    setting.B = eye(n,n);                      % setting.B defines the coordinate system
    setting.BD = setting.B.*repmat(setting.diagD',n,1);        % setting.B*D for speed up only
    setting.C = diag(setting.diagC);                   % covariance matrix == setting.BD*(setting.BD)'
end

if setting.flgDiagonalOnly
    setting.B = 1;
end
  setting.D = ones(n,1); 
  setting.invsqrtC =   setting.B * diag(  setting.D.^-1) *   setting.B';

fitness.hist=NaN*ones(1,10+ceil(3*10*n/n2)); % history of fitness values
fitness.histsel=NaN*ones(1,10+ceil(3*10*n/n2)); % history of fitness values
fitness.histbest=[]; % history of fitness values
fitness.histmedian=[]; % history of fitness values

% Initialize boundary handling
bnd.isactive = any(xmin > -Inf) || any(xmax < Inf);
if bnd.isactive
    bnd.weights = zeros(n,1);         % setting.weights for bound penalty
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero
    if bnd.flgscale ~= 0
        bnd.scale = setting.diagC/mean(setting.diagC);
    else
        bnd.scale = ones(n,1);
    end
    setting.idx = (xmin > -Inf) | (xmax < Inf);
    if length(setting.idx) == 1
        setting.idx = setting.idx * ones(n,1);
    end
    setting.idx = (xmin > -Inf) | (xmax < Inf);
    if length(setting.idx) == 1
        setting.idx = setting.idx * ones(n,1);
    end
    bnd.isbounded = zeros(n,1);
    bnd.isbounded(setting.idx) = 1;
    setting.maxdx = min(setting.maxdx, (xmax - xmin)/2);
    if any(setting.sigma*sqrt(setting.diagC) > setting.maxdx')
         setting.fac = min(setting.maxdx ./ sqrt(setting.diagC))/setting.sigma;
        setting.sigma = min(setting.maxdx ./ sqrt(setting.diagC));
    end
    
    setting.dd = setting.diagC;
    bnd.dfithist = 1;              % delta fit for setting setting.weights
    bnd.aridxpoints = [];          % remember complete outside points
    bnd.arfitness = [];            % and their fitness
    bnd.validfitval = 0;
    bnd.iniphase = 1;
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero
    if bnd.flgscale ~= 0
        bnd.scale = setting.diagC/mean(setting.diagC);
    else
        bnd.scale = ones(n,1);
    end
end

% Initialize dynamic (internal) strategy parameters and constants
% usually close to 1
% Initialize dynamic (internal) strategy parameters and constants


setting.endinvsqrtC = 0;%setting.B * diag(D.^-1) * setting.B';    % setting.C^-1/2
setting.eigeneval = 0;                      % track update of setting.B and D
setting.chiN=n^0.5*(1-1/(4*n)+1/(21*n^2));  % expectation of
%   ||n(0,I)|| == norm(randn(n,1))

setting.mu = ceil(n2/2);               % number of parents/points for recombination
setting.weights = log(max(setting.mu, n/2) + 1/2)-log(1:setting.mu)'; % muXone array for weighted recombination setting.mu = floor(setting.mu);
setting.mueff=sum(setting.weights)^2/sum(setting.weights.^2); % variance-effective size of setting.mu
setting.weights = setting.weights/sum(setting.weights);     % normalize recombination setting.weights array

% Strategy parameter setting: Adaptation
setting.cc = (4 + setting.mueff/n) / (n+4 + 2*setting.mueff/n); % time constant for cumulation for setting.C
setting.cs = (setting.mueff+2) / (n+setting.mueff+3);  % t-const for cumulation for setting.sigma control
setting.ccov1 = 2 / ((n+1.3)^2+setting.mueff);    % learning rate for rank-one update of setting.C
setting.ccovmu = 2 * (setting.mueff-2+1/setting.mueff) / ((n+2)^2+setting.mueff);  % and for rank-setting.mu update
setting.damps = 0.5 + 0.5*min(1, (0.27*n2/setting.mueff-1)^2) + 2*max(0,sqrt((setting.mueff-1)/(n+1))-1) + setting.cs; % damping for setting.sigma

if setting.flgDiagonalOnly
    setting.ccov1_sep = min(1, setting.ccov1 * (n+1.5) / 3);
    setting.ccovmu_sep = min(1-setting.ccov1_sep, setting.ccovmu * (n+1.5) / 3);
else
    setting.ccov1_sep=0;
    setting.ccovmu_sep=0;
end

setting.stopOnEqualFunctionValues= 2 + n/3;
setting.arrEqualFunvals = zeros(1, 10+n);
if isempty(setting.insigma)
    if all(size(( setting.xmean)) > 1)
        setting.insigma = std( setting.xmean, 0, 2);
    else
    end
end


fitness.hist(1)= EA_obj2(1);%cec14_func( setting.xmean,I_fno);
fitness.histsel(1)=fitness.hist(1);

setting.xold =  setting.xmean;
