function [jac] = JacMPDLsqnonlinMassDeltaPsi(structModel, expModes, simModes, ...
    eigFreqOpt, normOpt, objOpt)
% function [jac] = JacMPDLsqnonlin(structModel, expModes, simModes, ...
%     eigFreqOpt, normOpt, objOpt)
%
%   Yang Wang, Xinjun Dong, Dan Li, Yu Otsuki
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.1
%
% For implementation with MATLAB lsqnonlin, this function calculates
% the Jacobian matrix for various forms of the modal property
% difference approaches.
%
% Input:
%   structModel - a structure array with following fields of structural
%   model information:
%       M0 (N x N)- mass matrix (assumed accurate enough and no need to
%          update in current revision). Here N refers to the number of
%          degrees of freedom of the finite element model
%       K0 (N x N) - nominal stiffness matrix constructed with nominal
%          parameter values
%       K_j {N x N x n_alpha} - influence matrix corresponding to updating
%          variables (Note: the third dimension of K_j should be
%          equal to the number of updating variables). Here n_alpha refers
%          the number of stiffness updating variables
%       K (N x N) - stiffness matrix constructed with the current alpha
%         values, using K0 and K_j
%
%   expModes - a structure array with experimental modal properties for
%     model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs.Here n_meas refers to the number of measured DOFs
%       measDOFs (n_meas x 1) - measured DOFs
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
%
%   simModes - a structure array with simulated modal properties for
%     model updating:
%       lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at
%          measured DOFs
%       psi    (N x n_modes) - simulated mode shape vector at all DOFs

%   eigFreqOpt
%       0 - eigenvalue difference
%       1 - angular frequency difference (rad/s)
%       2 - ordinary frequency differnce (Hz)
%
%   normOpt
%       1 - normalize eigenvecotr qi-th entry equal to 1
%       2 - normalize eigenvecotr norm equal to 1
%
%   objOpt
%       1 - MAC value formulation
%       2 - eigenvector difference formulation
%
% Output:
%   jac: the Jacobian of the objective function (matrix)

n_alpha = length( structModel.K_j ) ;
n_beta = length(structModel.M_j);
N = size( structModel.K0, 1 );

omegaSim = sqrt( simModes.lambda );
omegaExp = sqrt( expModes.lambdaExp );

n_meas = expModes.n_meas;
n_modes = expModes.n_modes;

for i = 1 : n_modes
    if (normOpt == 1)
        simModes.psi_m(:,i) = simModes.psi_m(:,i) / simModes.psi_m(expModes.qm(i), i);
        simModes.psi(:,i) = simModes.psi(:,i) / simModes.psi(expModes.q(i), i);

    elseif (normOpt == 2)
        simModes.psi_m(:,i) = simModes.psi_m(:,i) / norm( simModes.psi(:, i) );
        simModes.psi(:,i) = simModes.psi(:,i) / norm( simModes.psi(:, i) );

    else
        error('\nWrong nomalization option (normOpt) input for JacMPDLsqnonlin.');
    end
end

modalMass = zeros( n_modes, 1 );
for i = 1 : n_modes
    modalMass(i) = simModes.psi(:, i)' * structModel.M * simModes.psi(:, i);
end

% dLambda will store values for d_lambda_i / d_alpha_j
dLambda = zeros( n_modes, n_alpha+n_beta );

% d_ri_eigFreqTerm will store the first part of d_r_i / d_alpha, i.e.
% the part involving the derivative of eigenvalue | angular
% frequency " ordinary frequency over alpha. For example, when eiegenvalue
% is used, the formulation is:
%       -weight_Lambda_i * D_alpha(Lambda_i) / Lambda_i^EXP
d_ri_eigFreqTerm = zeros( n_modes, n_alpha );

for i = 1 : n_modes
    for j = 1 : n_alpha
        dLambda(i,j) = simModes.psi(:,i)' * structModel.K_j{j} *...
            simModes.psi(:,i) / modalMass(i);
        if eigFreqOpt == 0
            d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                / expModes.lambdaExp(i) ;
        elseif eigFreqOpt == 1 || eigFreqOpt == 2
            d_ri_eigFreqTerm(i,j) = - dLambda(i,j) * expModes.lambdaWeights(i) ...
                / (omegaExp(i) * 2 * omegaSim(i));
        end
    end
    for j = 1 : n_beta
        dLambda(i,j+n_alpha) = simModes.psi(:,i)' * -simModes.lambda(i) * structModel.M_j{j} *...
            simModes.psi(:,i) / modalMass(i);
        if eigFreqOpt == 0
            d_ri_eigFreqTerm(i,j+n_alpha) = - dLambda(i,j+n_alpha) * expModes.lambdaWeights(i) ...
                / expModes.lambdaExp(i) ;
        elseif eigFreqOpt == 1 || eigFreqOpt == 2
            d_ri_eigFreqTerm(i,j+n_alpha) = - dLambda(i,j+n_alpha) * expModes.lambdaWeights(i) ...
                / (omegaExp(i) * 2 * omegaSim(i));
        end
    end
end

% The 3rd dimension corresponds to i in Psi_i^m -- the mode index
dPsi_m = zeros(n_meas, n_alpha + n_beta, n_modes);
dPsi_m2 = zeros(n_meas, n_alpha + n_beta, n_modes);
dPsi_m1 = zeros(n_meas, n_alpha + n_beta, n_modes);

for i = 1 : n_modes
    dPsi_dAlpha_j = zeros(N, 1);
    B = sparse( structModel.K - simModes.lambda(i) * structModel.M );

    if (normOpt == 1)
        % % % % % The maximum entry of Psi_m is normalized to 1
        % % % % P_i = setdiff(1 : N, expModes.q(i));
        % % % P_i = 1 : N;
        % Cross out q_i-th row in B and b, and q_i-th column in B
        % % B = B(P_i, P_i);

        % factorization
        dcompB = decomposition(B);

        for j = 1 : n_alpha
            b = dLambda(i,j) * structModel.M * simModes.psi(:, i) -...
                structModel.K_j{j} * simModes.psi(:, i);
            b = sparse(b);
            
            dPsi_dAlpha_j = dcompB \ b;
            
            dPsi_m1(:, j, i) = full(dPsi_dAlpha_j(expModes.measDOFs(:,1)));

            dPsi_m2(any(expModes.measDOFs(:,2),2), j, i) = full(dPsi_dAlpha_j(nonzeros(expModes.measDOFs(:,2))));
            
            dPsi_m(:, j, i)  = dPsi_m1(:, j, i)  - dPsi_m2(:, j, i);

            % Largest mode shape DOF is always normalized to 1, so no
            % change. 
%             dPsi_m(expModes.qm(i), j, i) = 0;

        end
        for j = 1 : n_beta
            b = dLambda(i,j+n_alpha) * structModel.M * simModes.psi(:, i) +...
                simModes.lambda(i) * structModel.M_j{j} * simModes.psi(:, i);
            b = sparse(b);

            dPsi_dAlpha_j = dcompB \ b;
  
            dPsi_m1(:, j+n_alpha, i) = full(dPsi_dAlpha_j(expModes.measDOFs(:,1)));

            dPsi_m2(any(expModes.measDOFs(:,2),2), j+n_alpha, i) = full(dPsi_dAlpha_j(nonzeros(expModes.measDOFs(:,2))));
            
            dPsi_m(:, j+n_alpha, i)  = dPsi_m1(:, j+n_alpha, i)  - dPsi_m2(:, j+n_alpha, i);


            % The q_i-th entry of dPsi_m remains as 0 
%             dPsi_m(expModes.qm(i), j, i) = 0;

        end

    else
        % Length of Psi is normalized to 1.
        % factorization
        dcompB = decomposition(B);

        for j = 1 : n_alpha
            b = dLambda(i,j) * structModel.M * simModes.psi(:, i) -...
                structModel.K_j{j} * simModes.psi(:, i);
            b = sparse(b);
            
            v = dcompB \ b;

            % R. B. Nelson, "Simplified calculation of eigenvector
            % derivatives," AIAA journal, vol. 14, pp. 1201-1205, 1976.
            c = -simModes.psi(:,i)' *  v ;
            dPsi_dAlpha_j = v + c * simModes.psi(:,i);
            dPsi_m(:, j, i) = dPsi_dAlpha_j(1 : n_meas);
            
            dPsi_m1(:, j, i) = full(dPsi_dAlpha_j(expModes.measDOFs(:,1)));

            dPsi_m2(any(expModes.measDOFs(:,2),2), j, i) = full(dPsi_dAlpha_j(nonzeros(expModes.measDOFs(:,2))));
            
            dPsi_m(:, j, i)  = dPsi_m1(:, j, i)  - dPsi_m2(:, j, i);


        end
        for j = 1 : n_beta
            b = dLambda(i,j+n_alpha) * structModel.M * simModes.psi(:, i) +...
                simModes.lambda(i) * structModel.M_j{j} * simModes.psi(:, i);
            b = sparse(b);

            v = dcompB \ b;

            % R. B. Nelson, "Simplified calculation of eigenvector
            % derivatives," AIAA journal, vol. 14, pp. 1201-1205, 1976.
            c = -simModes.psi(:,i)' *  v ;
            dPsi_dAlpha_j = v + c * simModes.psi(:,i);
  
            dPsi_m1(:, j+n_alpha, i) = full(dPsi_dAlpha_j(expModes.measDOFs(:,1)));

            dPsi_m2(any(expModes.measDOFs(:,2),2), j+n_alpha, i) = full(dPsi_dAlpha_j(nonzeros(expModes.measDOFs(:,2))));
            
            dPsi_m(:, j+n_alpha, i)  = dPsi_m1(:, j+n_alpha, i)  - dPsi_m2(:, j+n_alpha, i);
        end
    end
end

if (objOpt == 1)
    jac_p = zeros(2 * n_modes, n_alpha + n_beta);
else
    if(normOpt == 1)
        jac_p = zeros(n_meas * n_modes, n_alpha + n_beta);
    else
        jac_p = zeros((n_meas + 1) * n_modes, n_alpha + n_beta);
    end
end

for i = 1 : n_modes
    psiExp = expModes.psiExp(:,i);
    psiSim = simModes.psi_m(:,i);

    if (objOpt == 1)
        % For the MAC value formulation, d_ri_MAC is the second part of
        % d_r_i/d_alpha that involves MAC value:
        MACValue = MAC(psiExp, psiSim);
        d_ri_MAC = (-expModes.psiWeights(i) / sqrt(MACValue)) * ...
            (psiExp' / (psiExp' * psiSim) - psiSim' / norm(psiSim)^2) ...
            * dPsi_m(:, :, i);
        jac_p((i - 1) * 2 + 1 : i * 2, :) = [d_ri_eigFreqTerm(i,:); d_ri_MAC];
    else
        % For the eigenvector difference formulation, d_ri_psi is the
        % second part of d_r_i/d_alpha that involves eigenvectors:
        %      -Qi * D_alpha(Psi_i^m) * weight_Psi_i
        if(normOpt == 1)
            Q_i = setdiff(1 : n_meas, expModes.q(i));
            d_ri_psi = -dPsi_m(Q_i, :, i) * expModes.psiWeights(i);
            jac_p((i - 1) * n_meas + 1 : i * n_meas, :) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
        else
            d_ri_psi = -expModes.psiWeights(i) * ...
                ((psiExp' * psiSim* eye(n_meas) + psiSim *...
                psiExp') / norm(psiSim)^2 - 2 * (psiSim * psiExp' * (psiSim * psiSim'))...
                / norm(psiSim)^4) * dPsi_m(:, :, i);
            jac_p((i - 1) * (n_meas + 1) + 1 : i * (n_meas + 1),:) = ...
                [d_ri_eigFreqTerm(i,:); d_ri_psi];
        end
    end
end

jac = jac_p;

end
