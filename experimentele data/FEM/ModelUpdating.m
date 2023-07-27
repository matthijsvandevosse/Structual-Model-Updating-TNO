clc;clear
close all
warning('off')

model



%% Load Sensitivity matrix

n_modes = 3; % Number of measured modes
modeIndex = 1:n_modes ; % Indexes of these measured modes

%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

structModel.K_j{1} = sparse(Kdiff1);
structModel.K_j{2} = sparse(Kdiff2);
structModel.K_j{3} = sparse(Kdiff4);   
structModel.K_j{4} = sparse(Kdiff5); 
%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = ones(1,length(structModel.K_j));
n_alpha = length(alpha_act);

%%

N = size(structModel.K0,1);

% measurement points -- sensor_loc = [55 155 255 355 455]/1000;
measDOFs = N/2+[6; 8; 10; 14; 5; 7; 11; 16;];

% measDOFs = 1:325


nonmeasDOFs =  N/2+[12];

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-6;
optimzOpts.tolX = 1e-6;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Simulate "experimental data"


load("eigenFrequencies.mat")
freqExp = eigenFrequencies(1:n_modes);
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")

for i = 1: 5
    [~,index] = max(abs(L(:,i)));
    L(:,i) = L(:,i) / (L(index,i));
end


psiExpAll = zeros(length(measDOFs),1);
for i = 1:n_modes
psiExpAll([1], i) = L(1,i);
psiExpAll([2], i) = L(2,i);
psiExpAll([3], i) = L(3,i);
psiExpAll([4], i) = L(5,i);
psiExpAll([5], i) = 0.99*L(1,i);
psiExpAll([6], i) = 0.99*L(2,i);
psiExpAll([7], i) = 0.99*L(3,i);
psiExpAll([8], i) = 0.99*L(5,i);
end

psi_m = psiExpAll;

%%
expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);
% expModes.lambdaWeights = [1 1 1];
expModes.lambdaWeights = 100*[1 1 1];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [10000 100 1000];
expModes.resWeights = ones(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = [1:(n_modes)];
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -1*ones(n_alpha,1);
    updatingOpts.x_ub =  5*ones(n_alpha,1);
    
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% MultiStart optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating

x(1:n_alpha)
%%
Ksolved = full(structModel.K0);
for n = 1:n_alpha
Ksolved = Ksolved + x(n)*full(structModel.K_j{n});
end

[Psi_solvedxy, lambdasolved] = eig(full(structModel.M0)\Ksolved);
 


[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
lambdasolved = lambdasolved(1:(n_modes));
Psi_solved = Psi_solvedxy(N/2+1:end, dummyInd(1:(n_modes)));
%%



naturalfrequency = sqrt(lambdasolved)/(2*pi)
freqExp
%%
for i = 1:size(Psi_solved,2)
    [index] = find(abs(psi_m(:,i)) == 1);
    Psi_solved(:,i) = Psi_solved(:,i) / (Psi_solved(measDOFs(index)-N/2,i));
end



Psi_solved(measDOFs(1:4)-N/2,:)
L([1 2 3 5],1:n_modes)

Psi_solved(measDOFs-N/2,:) = psi_m;

Psi_solved(nonmeasDOFs-N/2,:)
Psi_measured =  L(4,1:n_modes)


%%
figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,1))
axis equal

if n_modes >1
figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,2))
axis equal
end
if n_modes >2
    figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,3))
axis equal
end
if n_modes >3
    figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,4))
axis equal
end