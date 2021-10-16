%This function is used as a local search, and creates some
%new points based on Gaussian Walks.

%**************************************************************************
%The input function is:                                                   %
%Point: the input point which is going to be diffused                     %
%S: structure of problem information                                      %
%g: generation number                                                     %
%BestPoint: the best point in group                                       %                               %
%==========================================================================
%The output function is:                                                  %
%createPoint: the new points created by Diffusion process                 %
%扩散过程产生的新点
%fitness: the value of fitness function                                   %
%**************************************************************************

function [createPoint, fitness] = LS_Process(Point,S,g,BestPoint)
   
    fhd=@cec14_func;

    GeneratePoint = normrnd(BestPoint, (log(g)/g)*(abs((Point - BestPoint))), [1 size(Point,2)]) + ...
        (randn*BestPoint - randn*Point);
    
    %check bounds of generated point
    GeneratePoint = Bound_Checking(GeneratePoint,S.Lband,S.Uband);
    
%     size(GeneratePoint)  
%     for i=1:size(Point,2)
%         if GeneratePoint(1,i) > S.Uband
%             fprintf('violate upper');
%         end
%         if GeneratePoint(1,i) < S.Lband
%              fprintf('violate lower');
%         end  
%     end
    
    fitness = feval(fhd,GeneratePoint',S.FuncNo);

    createPoint = GeneratePoint;
    %======================================================================
end