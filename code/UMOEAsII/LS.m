%% ============United Multi-Operator Evolutionary AlgorithmsII ============
% Should you have any queries, please contact
% Dr. Saber Elsayed. University of New South Wales at Canberra
% s.elsayed@adfa.edu.au
% www.saberelsayd.net or 
% https://sites.google.com/site/saberelsayed3/home
% =========================================================================
function [x,f,current_eval,succ] = LS (bestx,f,Par,current_eval,I_fno,Max_FES,xmin,xmax)


Par.LS_FE=ceil(20.0000e-003*Max_FES ); %% Max FFEs_LS

options=optimset('Display','off','algorithm','interior-point',...
    'UseParallel','never','MaxFunEvals',Par.LS_FE) ;

[Xsqp, FUN , ~ , details]=fmincon(@cec14_func, bestx(1,:)',[],[],[],[],xmin,xmax, [],options,I_fno);

%% check if there is an improvement in the fitness value and update P_{ls}
if (f-FUN)>0
    succ=1;
    f = FUN;
    x(1,:)=Xsqp;
else
    succ=0;
    x=bestx;
   
end
%% update FFEs
current_eval=current_eval+details.funcCount;
end


