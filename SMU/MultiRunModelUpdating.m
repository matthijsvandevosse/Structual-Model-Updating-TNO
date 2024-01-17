
    n_x = length(updatingOpts.x_lb);
    optimzOpts.x0 = zeros(n_x, 1);


    updtResults = StructModelUpdating(structModel, expModes, updatingOpts, optimzOpts);
    x(:,runNum) = updtResults.xOpt;
    fval(runNum) = updtResults.fvalOpt;
    exit_flag(runNum) = updtResults.exitFlag;
    gradient(:,runNum) = updtResults.gradient;
    optmzSolvOutput{runNum} = updtResults.output;
    
