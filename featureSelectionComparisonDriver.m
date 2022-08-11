% Comparison_driver.m
% 
% Tests a variety of algorrithms on one benchmark function benchmark is 
% defined seperatley and the handle is passed in, paramaters are set inside
% files for some algorithms but others are set in this file

close all
clear
clc

%% Initalization

Algorithm_Names = ["I-GWO" "PSO" "DE" "BRO" "GA" "MVO" "SCA" "DFO" "ACO" "AHA" "SMA" "MSA" "HGS" "RUN" "HHO" "HGSO"];

Pop = 10; % Number of search agents
Max_iteration = 100; % Maximum number of iterations
NumRun = 25; % run for this many times

% Ratio of validation data
ho = 0.2;
% Common parameter settings 
opts.N = 10;     % number of solutions
opts.T = 100;    % maximum number of iterations
opts.method = 'knn'; % knn, rf, dt
% Number of k in K-nearest neighbor
opts.k = 5; 
% Number of trees in random forest
opts.numtrees = 200; 

% Load dataset
% works with any dataset in .mat format with categories in label and 
% parameters in feat, can take input in other dataset formats (.csv, .xlsx)
% but will need to format them in this manner and then the code can be run
% as is from the partician section. The equivalencies of label and feat
% would be y and x respectivley (this is how they are fed into the models)
load ionosphere.mat;

% Divide data into training and validation sets
HO = cvpartition(label,'HoldOut',ho); 
opts.Model = HO; 
% Perform feature selection 
% fobj = @jFitnessFunction;
FS    =  jfs('igwo',feat,label,opts);
% Define index of selected features
sf_idx = FS.sf;
% Accuracy  
Acc    = jknn(feat(:,sf_idx),label,opts); 
% Plot convergence
plot(FS.c); grid on;
xlabel('Number of Iterations'); 
ylabel('Fitness Value'); 
title('IGWO');

Function_name = 'knn'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
    
    % Load details of the selected benchmark function
    %[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

    %% move datasets here in updated version to test many at once
    %for i=1:23

    %% IGWO
    %for j = 1:NumRun
    [IGWO_bestfit,IGWO_bestpos,IGWO_curve] = IGWO(dim,Pop,Max_iteration,lb,ub,fobj);
    display(['The best solution obtained by I-GWO for function ', Function_name ' is : ', num2str(IGWO_bestpos)]);
    display(['The best optimal value of the objective funciton found by I-GWO is : ', num2str(IGWO_bestfit), newline]);
    
        
    %end
    %% PSO

    [Bestfit, Bestpos, PSO_CVNG] = PSO(problem, params);
    display(['The best solution obtained by PSO for function ', Function_name ' is : ', num2str(Bestpos)]);
    display(['The best optimal value of the objective funciton found by PSO is : ', num2str(Bestfit), newline]);
    
    %% DE

    [DE_curve,de_bestsol,de_bestfitness] = DE(dim,Pop,Max_iteration,lb,ub,fobj);
    display(['The best solution obtained by DE for function ', Function_name ' is : ', num2str(de_bestsol)]);
    display(['The best optimal value of the objective funciton found by DE is : ', num2str(de_bestfitness), newline]);
    
    %% BRO 

    [BRO_bestfit, BRO_bestpos,BRO_cg_curve] = BRO_Fun(Pop,Max_iteration,fobj,lb,ub,dim);
    display(['The best solution obtained by BRO for function ', Function_name ' is : ', num2str((BRO_bestpos)')]);
    display(['The best optimal value of the objective funciton found by BRO is : ', num2str(BRO_bestfit), newline]);
     
    %% GA
                 
    gaDat2=ga(gaDat2);
    display(['The best solution obtained by GA for function ', Function_name ' is : ', num2str(gaDat2.xmin)]);
    display(['The best optimal value of the objective funciton found by GA is : ', num2str(gaDat2.fxmin), newline]);

    %% MVO
    
    [MVO_Best_score,MVO_Best_pos,MVO_cg_curve] = MVO(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by MVO for function ', Function_name ' is : ', num2str(MVO_Best_pos)]);
    display(['The best optimal value of the objective funciton found by MVO is : ', num2str(MVO_Best_score), newline]);


    %% SCA

    [SCA_Best_score,SCA_Best_pos,SCA_cg_curve]=SCA(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by SCA for function ', Function_name ' is : ', num2str(SCA_Best_pos)]);
    display(['The best optimal value of the objective funciton found by SCA is : ', num2str(SCA_Best_score), newline]);

    
    %% DFO
    
    [DFO_curve,DFO_best_pos,DFO_best_fit]=DFO(dim,Pop,Max_iteration,lb,ub,fobj);
    display(['The best solution obtained by DFO for function ', Function_name ' is : ', num2str(DFO_best_pos)]);
    display(['The best optimal value of the objective funciton found by DFO is : ', num2str(DFO_best_fit), newline]);
    
    %% ACO

    [ACO_curve, ACO_best_sol, ACO_best_fit] = ACO(Pop,Max_iteration,fobj,lb,ub,dim);
    display(['The best solution obtained by ACO for function ', Function_name ' is : ', num2str(ACO_best_sol)]);
    display(['The best optimal value of the objective funciton found by ACO is : ', num2str(ACO_best_fit), newline]);

    %% AHA
    
    [AHA_BestX,AHA_BestF,AHA_HisBestFit,AHA_VisitTable] = AHA(Max_iteration,Pop,dim,lb,ub,fobj);
    display(['The best solution obtained by AHA for function ', Function_name ' is : ', num2str(AHA_BestX)]);
    display(['The best optimal value of the objective funciton found by AHA is : ', num2str(AHA_BestF), newline]);    

    %% SMA

    [SMA_Destination_fitness, SMA_bestPosition, SMA_Convergence_curve] = SMA(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by SMA for function ', Function_name ' is : ', num2str(SMA_bestPosition)]);
    display(['The best optimal value of the objective funciton found by SMA is : ', num2str(SMA_Destination_fitness), newline]); 

    %% MSA

    [MSA_Best_pos,MSA_Best_score,MSA_Convergence_curve] = MSA(Pop,Max_iteration,ub,lb,dim,fobj);
    display(['The best solution obtained by MSA for function ', Function_name ' is : ', num2str(MSA_Best_pos)]);
    display(['The best optimal value of the objective funciton found by MSA is : ', num2str(MSA_Best_score), newline]); 

    %% HGS

    [HGS_Destination_fitness,HGS_bestPositions,HGS_Convergence_curve] = HGS(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by HGS for function ', Function_name ' is : ', num2str(HGS_bestPositions)]);
    display(['The best optimal value of the objective funciton found by HGS is : ', num2str(HGS_Destination_fitness), newline]); 

    %% RUN

    [RUN_Best_Cost,RUN_Best_X,RUN_Convergence_curve] = RUN(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by RUN for function ', Function_name ' is : ', num2str(RUN_Best_X)]);
    display(['The best optimal value of the objective funciton found by RUN is : ', num2str(RUN_Best_Cost), newline]); 

    %% HHO

    [HHO_Rabbit_Energy,HHO_Rabbit_Location,HHO_CNVG] = HHO(Pop,Max_iteration,lb,ub,dim,fobj);
    display(['The best solution obtained by HHO for function ', Function_name ' is : ', num2str(HHO_Rabbit_Location)]);
    display(['The best optimal value of the objective funciton found by HHO is : ', num2str(HHO_Rabbit_Energy), newline]); 

    %% HGSO

    [HGSO_pos, HGSO_fit, HGSO_CNVG] = HGSO(fobj, dim, lb, ub, Max_iteration, Pop, 5);
    display(['The best solution obtained by HGSO for function ', Function_name ' is : ', num2str(HGSO_pos)]);
    display(['The best optimal value of the objective funciton found by HGSO is : ', num2str(HGSO_fit), newline]); 

    % end

    %% convert everything to averages and find stddev

    
    % divide all curves by number of runs to get averages

    % find stddev



    %% graph creation

    % add more colours to plot
    colours = hsv(16);

    fh = figure();

    semilogy(IGWO_curve, 'LineWidth', 2,'color',colours(1,:));
    hold on;
    semilogy(PSO_CVNG, 'LineWidth', 2,'color',colours(2,:));
    semilogy(DE_curve, 'LineWidth', 2,'color',colours(3,:));
    semilogy(BRO_cg_curve, 'LineWidth', 2,'color',colours(4,:));
    semilogy(gaDat2.fxmingen, 'Linewidth', 2,'color',colours(5,:));
    semilogy(MVO_cg_curve, 'LineWidth', 2,'color',colours(6,:));
    semilogy(SCA_cg_curve, 'LineWidth', 2,'color',colours(7,:));
    semilogy(DFO_curve, 'LineWidth', 2,'color',colours(8,:));
    semilogy(ACO_curve, 'LineWidth', 2,'color',colours(9,:));
    semilogy(AHA_HisBestFit, 'LineWidth', 2,'color',colours(10,:));
    semilogy(SMA_Convergence_curve, 'LineWidth', 2,'color',colours(11,:));
    semilogy(MSA_Convergence_curve, 'LineWidth', 2,'color',colours(12,:));
    semilogy(HGS_Convergence_curve, 'LineWidth', 2,'color',colours(13,:));
    semilogy(RUN_Convergence_curve, 'LineWidth', 2,'color',colours(14,:));
    semilogy(HHO_CNVG, 'LineWidth', 2,'color',colours(15,:));
    semilogy(HGSO_CNVG, 'LineWidth', 2,'color',colours(16,:));
    
    xlabel('Iterations');
    ylabel('Best Cost');
    title(Function_name);
    grid on;
    axis tight
    box on
    legend(Algorithm_Names);
    fh.WindowState = 'maximized';

    % uncomment to save figure
    %saveas(gcf,strcat(Algorithm_Name,Function_name)); %needs changes

%end

