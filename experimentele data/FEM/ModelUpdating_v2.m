clc;clear
close all
warning('off')

model_v2

pause(0.5)

%% Load Sensitivity matrix

modeIndex = 1:4; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

structModel.K_j = Kdiff;

%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act);

%%

N = size(structModel.K0,1);

% measurement points -- sensor_loc = [55 155 255 355 455]/1000;
measDOFs = N/2+[6; 8; 10; 14; 5; 7; 9; 16;];

% measDOFs = 1:325


nonmeasDOFs =  N/2+[12];

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-10;
optimzOpts.tolX = 1e-10;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
% optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Simulate "experimental data"


load("eigenFrequencies.mat")
freqExp = eigenFrequencies(modeIndex);
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")
freqExp(4) = eigenFrequencies(5);
L(:,4) = L(:,5);
for i = modeIndex
    [m,index] = max(abs(L([1 2 3 5],i)));
    L(:,i) = L(:,i) / m;
end


psiExpAll = zeros(length(measDOFs),1);
n = 1;
for i = modeIndex
psiExpAll([1], n) = L(1,i);
psiExpAll([2], n) = L(2,i);
psiExpAll([3], n) = L(3,i);
psiExpAll([4], n) = L(5,i);
psiExpAll([5], n) = 0.9999*L(1,i);
psiExpAll([6], n) = 0.9999*L(2,i);
psiExpAll([7], n) = 0.9999*L(3,i);
psiExpAll([8], n) = 0.9999*L(5,i);
n = n+1;
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
expModes.lambdaWeights = [10 10 10 1];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [2000 200 400 10];
expModes.resWeights = ones(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -2*ones(n_alpha,1);
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
% lambdasolved = lambdasolved(modeIndex);
Psi_solved = Psi_solvedxy(N/2+1:end, dummyInd(modeIndex));
%%



naturalfrequency = sqrt(lambdasolved(modeIndex))/(2*pi)
freqExp
%%
sign = [-1 1 1 -1];
for i = modeIndex
    [~, index] = max(abs(Psi_solved(measDOFs(1:4)-N/2,i)));
    Psi_solved(:,i) = sign(i)* Psi_solved(:,i) / (Psi_solved(measDOFs(index)-N/2,i));
end

Psi_solved(measDOFs(1:4)-N/2,:)
L([1 2 3 5],modeIndex)

Psi_solved(nonmeasDOFs-N/2,:)
Psi_measured =  L(4,modeIndex)


%%
figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,1))
axis equal

figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,2))
axis equal

figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,3))
axis equal
%%
figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,4))
axis equal
%%
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,1), LineWidth=2)
hold on
plot([1 2 3 5], psi_m(1:4,1), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,1), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'Measured', 'Validation')
ylim([-1 1])
title('First Mode')

if n_modes >1
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,2), LineWidth=2)
hold on
plot([1 2 3 5], psi_m(1:4,2), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,2), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'Measured', 'Validation')
ylim([-1 1])
title('Second Mode')

end
if n_modes >2
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,3), LineWidth=2)
hold on
plot([1 2 3 5], psi_m(1:4,3), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,3), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'Measured', 'Validation')
ylim([-1 1])
title('Third Mode')
end
if n_modes >3
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,4), LineWidth=2)
hold on
plot([1 2 3 5], psi_m(1:4,4), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,4), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'Measured', 'Validation')
ylim([-1.2 1.2])
title('Fourth Mode')
end


