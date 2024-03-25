%% Load Sensitivity matrix

modeIndex = [1 2 3 4 5 6]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);
structModel.M_j = [];
structModel.K_j = Kdiff;


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act);


%% Optimization structure parameter;
optimzOpts.tolFun = 1e-5;
optimzOpts.tolX = 1e-5;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
% optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Load experimental mode shapes and resonace frequencies

% Modal model fitted using Least Squared Optimalization
load("G_ss_modal_fitted.mat")

% load resonance frequencies 
lambdaExp = diag(-G_ss_Modal.A(9:end,1:8));
freqExp =sqrt(lambdaExp)/2/pi;

% Convert rad/s to hz
% lambdaExp = freqExp.^2;

% Load sensor mode shapes
L = G_ss_Modal.C(:,1:8);
R = G_ss_Modal.B(9:end,:)';

% Normalize the sensor mode shapes, note that 4 sensor is used for
% validation and not used to find maximum entry. 
for i = 1:length(L)
    [m,index] = max(abs(L([1 2 3 4 5],i)));
    L(:,i) = L(:,i) / m;
    [m,index] = max(abs(R(:,i)));
    R(:,i) = R(:,i) / m;
end

% Save mode shapes correctly to psi_m
psiExpAll = zeros(length(measDOFs),1);
n = 1;
for i = modeIndex
psiExpAll([1], n) = L(1,i);
psiExpAll([2], n) = L(2,i);
psiExpAll([3], n) = L(3,i);
psiExpAll([4], n) = L(5,i);
n = n+1;
end
psi_m = [ psiExpAll(1,:); psiExpAll(2,:); psiExpAll(3,:); psiExpAll(4,:)];

%%
expModes.lambdaExp = lambdaExp(1:n_modes);
expModes.psiExp = psi_m(:,1:n_modes);
expModes.measDOFs = measDOFs;

unmeasDOFs = setdiff(1 : N, measDOFs(:,1));
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);


expModes.lambdaWeights = [50 50 30 30 0 10];
expModes.lambdaWeights = expModes.lambdaWeights(1:n_modes);

expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [2000 400 400 200 0 100 ];
expModes.psiWeights = expModes.psiWeights(1:n_modes);
expModes.resWeights = ones(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1.4;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
                                 % 1.3: Zero Matching
                                 % 1.4: Delta Psi Matching

updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -2*ones(n_alpha,1);
    updatingOpts.x_ub =  2*ones(n_alpha,1);
    
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% Start optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];


n_x = length(updatingOpts.x_lb);
optimzOpts.x0 = zeros(n_x, 1);

optimzOpts.x0 = 0.1*rand(n_x, 1);



updtResults = StructModelUpdating(structModel, expModes, updatingOpts, optimzOpts);

x = updtResults.xOpt;
fval = updtResults.fvalOpt;
optmzSolvOutput = updtResults.output;
    
