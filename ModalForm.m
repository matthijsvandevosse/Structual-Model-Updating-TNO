
structModel.D0 = sparse(10e-6*eye(N));


gain = 800;

% LeastSquaredOpt

% Derived using LeastSquaredOpt (takes a long time!!)

structModel.D = sparse(9e-6*eye(N));


nModes =100;

[Psi_solved, lambdasolved] = eigs(structModel.K, structModel.M0,nModes,1e3);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);

for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*full(structModel.M0)*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(i,i) = Psi_solved(:,i)'*structModel.D*Psi_solved(:,i);
end

Dm_fitted = -diag(G_ss_Modal.A(9:end,9:end));
Dm(1:length(Dm_fitted),1:length(Dm_fitted)) = diag(Dm_fitted);
Dm(16:18,16:18) = 10;
%%
clear R G

freq = ffrf;

for ii = 1:length(freq)
    w = freq(ii);
    for m = 1:nModes%length(lambda)
        R(:,:,m) = Psi_solved([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],m)*Psi_solved(actDOFs,m)';
        if m == 1
            G(:,:,ii) = gain* R(:,:,m)./( -w^2 + 1i*Dm(m,m)*w + lambdasolved(m) );
        else
            G(:,:,ii) = G(:,:,ii) + gain* R(:,:,m)./( -w^2 + 1i*Dm(m,m)*w + lambdasolved(m) );  
        end
    end
end


G_mag = (abs(G));
    
G_pha = angle(G);


%% Initial model

[Psi_solved, lambdainit] = eigs(structModel.K0, structModel.M0,nModes,1e3);

[lambdainit,dummyInd] = sort(diag(lambdainit), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);


for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*full(structModel.M0)*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(i,i) = Psi_solved(:,i)'*structModel.D*Psi_solved(:,i);

end
Dm(1:length(Dm_fitted),1:length(Dm_fitted)) = diag(Dm_fitted);
clear R G0

for ii = 1:length(freq)
    w = freq(ii);
    for m = 1:nModes%length(lambda)
        R(:,:,m) = Psi_solved([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],m)*Psi_solved(actDOFs,m)';
        if m == 1
            G0(:,:,ii) = gain* R(:,:,m)./( -w^2 + 1i*Dm(m,m)*w + lambdainit(m) );
        else
            G0(:,:,ii) = G0(:,:,ii) + gain* R(:,:,m)./( -w^2 + 1i*Dm(m,m)*w + lambdainit(m) );  
        end
    end
end


G0_mag = (abs(G0));
    
G0_pha = angle(G0);

