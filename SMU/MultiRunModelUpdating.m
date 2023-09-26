% Script: MultiRunModelUpdating 
%
%   Yang Wang, Xinjun Dong, Dan Li, Yu Otsuki
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.1
%
% The script can perform multiple runs of model updating by starting the
% search from randomized starting points. 
%
% Must have following variables before running the script. See
% documentation of StructModelUpdating.m for descriptions.
%   structModel, expModes, updatingOpts, optimzOpts, filename
%
% If following variables do not exist, default actions will be
% assigned.
%    numRuns: number of starting points.  Default = 1
%    randSeed: random seed value for rng.  By default the random seed will
%       NOT be fixed.
%
% If multiple runs are requested through "numRuns" variable and the runs
% are disrupted half way, the script can pick up and continue the runs by
% loading previous run results from the *.mat named by the "filename"
% variable. IF A USER WANTS TO START ALL THE RUNS AFRESH, THE *.MAT FILE
% MUST BE DELETED FIRST.
%

if ~exist('numRuns', 'var')
    numRuns = 1;
end

if exist('randSeed', 'var')
    rng (randSeed);
end

n_x = length(updatingOpts.x_lb);

if (exist(filename, 'file') ~= 2)
    % If the result file does not exist, start the runs afresh by
    % initializing result storage.
    fval = zeros(numRuns, 1); % objective function values.
    exit_flag = zeros(numRuns, 1); % exit flags
    gradient = zeros(n_x, numRuns); % gradients
    x = zeros(n_x, numRuns); % optimal solutions
    
    % Generate starting points based on uniform distrubution
    x0 = zeros(n_x, numRuns); % starting points
    % t: computing time spend on feasible in-bounds results
    % t_waste: computing time wasted on infeasible out-of-bounds results
    t = zeros(1, numRuns);
    t_waste = 0;
    
    % runNum: number of valid runs with feasible in-bounds results
    % numWasteRuns: number of discarded result sets outside feasible bounds.
    runNum = 1;
    numWasteRuns = 1;
    % eval(['save ' filename ' x x0 fval exit_flag gradient t ' ...
    %     't_waste structModel expModes updatingOpts optimzOpts']);
    
else
    % If the result file exists, pick up from the previous runs that may
    % have been interrupted in the middle.
    prompt = sprintf(['\nWarning: A mat file ("%s") \n'...
        'containing results from a previous run exists in the current folder. \n'...
        'The file may have results from previously finished starting points. \n'...
        'If you intend for a fresh new run, the file needs to be deleted first.\n'...
        'Would you like to keep results in the mat file, and continue at the next\n'...
        'random starting point?\nY/N [Y]: '], filename);
    userInput = input(prompt, 's');
    
    if isempty (userInput)
        userInput = 'Y';
    end
    if( upper(userInput) ~= 'N')
        fprintf('\nOptimization continues at the next random starting point.\n')
    else
        error('\nPlease delete the file "%s" with previous results and run again.\n', filename);
    end
    
    % load(filename);
    runNum = length(find(t ~= 0)) + 1;
    numWasteRuns = size(t_waste, 2);
    
    % If the random seed is fixed, re-generated previously used starting
    % points and discard, so that future rand provides new starting points.
    if exist('randSeed', 'var')
        for i = 1 : (runNum + numWasteRuns - 1)
            rand(length(updatingOpts.x_lb), 1);
        end
    end
end

while(runNum <= numRuns)
    tic
    fprintf(1, '\n*************************************************************');
    fprintf(1, '\nSearching from #%d among the %d starting points.', runNum, numRuns);
    fprintf(1, '\n*************************************************************\n');
    % load( filename );
    x0(:,runNum) = 0.05*updatingOpts.x_lb + ...
            .05*(updatingOpts.x_ub - updatingOpts.x_lb) .* rand(n_x, 1);

    % For modal dynamic residual approach
    if(updatingOpts.formID == 3)
        K_ini = structModel.K0;
        for i = 1 : n_alpha
            K_ini = K_ini + x0(i, runNum) * structModel.K_j{i};
        end
        [psiSim,lambdaSim] = eig(full(structModel.M0)\full(K_ini));
        [lambdaSim,dummyInd] = sort((diag(lambdaSim)),'ascend') ;
        lambdaSim = lambdaSim(modeIndex);
        psiSim = psiSim(:,dummyInd(modeIndex));
        psiSim_m = psiSim(measDOFs,:);
        psiSim_u = psiSim(unmeasDOFs,:);
        % Normalize the mode shape vectors by maximum entry
        for i = 1:n_modes
            [~,index] = max(abs(psiSim_m(:,i)));
            psiSim_u(:,i) = psiSim_u(:,i) / psiSim_m(index,i);
        end
        x0(n_alpha + 1 : end, runNum) = reshape(psiSim_u, ...
            length(unmeasDOFs) * n_modes, 1);
    end
    
    optimzOpts.x0 = x0(:, runNum);
    optimzOpts.x0 = 0.*x0;
    % if  exist('x_2','var')
    %     optimzOpts.x0 = [x_2; zeros(n_beta, 1 )];
    % end

    updtResults = StructModelUpdating(structModel, expModes, updatingOpts, optimzOpts);
    x(:,runNum) = updtResults.xOpt;
    fval(runNum) = updtResults.fvalOpt;
    exit_flag(runNum) = updtResults.exitFlag;
    gradient(:,runNum) = updtResults.gradient;
    optmzSolvOutput{runNum} = updtResults.output;
    
    % Decide whether the results are in-bounds feasible; if not, discard
    % the results and re-run at another randomized starting point.
    if(isempty( find(x(:,runNum) > updatingOpts.x_ub, 1)) && ...
            isempty(find(x(:,runNum) < updatingOpts.x_lb, 1)))
        % Results are in-bounds
        fprintf(1, 'Successfully found a solution point within bounds.\n');
        t(runNum) = toc;
        % eval(['save ' filename ' x x0 fval exit_flag gradient t ' ...
        %     't_waste structModel expModes updatingOpts optimzOpts optmzSolvOutput']);
        runNum = runNum + 1;
    else
        % Some result entries are out of bounds.
        fprintf(1, ['The solution point is out of bounds and discarded. ' ...
            'Restart now from next random starting point.\n']);
        t_waste(numWasteRuns) = toc;
        % eval(['save ' filename ' x x0 fval exit_flag gradient t ' ...
        %     't_waste structModel expModes updatingOpts optimzOpts optmzSolvOutput']);
        numWasteRuns = numWasteRuns + 1;
    end
end
