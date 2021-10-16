%% ============United Multi-Operator Evolutionary AlgorithmsII ============
% Should you have any queries, please contact
% Dr. Saber Elsayed. University of New South Wales at Canberra
% s.elsayed@adfa.edu.au
% www.saberelsayd.net or
% https://sites.google.com/site/saberelsayed3/home
% =========================================================================
function [Par] = Introd_Par(I_fno)

%% loading

Par.n_opr=3;  %% number of operators in MODE
Par.n=50;     %% number of decision vriables

if Par.n==10
    Par.CS=50; %% cycle
elseif Par.n==30
    Par.CS=100; %% cycle
elseif Par.n==50
    Par.CS=150; %% cycle
else
    Par.CS=150; %% cycle
end
opt= 100:100:3000;       %% define the optimal solution as shown in the TR
Par.xmin= -100*ones(1,Par.n);
Par.xmax= 100*ones(1,Par.n);
Par.Max_FES=Par.n*10000;
Par.f_optimal=opt(I_fno);
Par.PopSize=18*Par.n; %% population size
Par.MinPopSize=4;
Par.prob_ls=0.1;
%% printing the detailed results- this will increase the computational time
Par.Printing=1; %% 1 to print; 0 otherwise

end