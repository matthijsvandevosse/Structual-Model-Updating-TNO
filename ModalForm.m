
nModes =30;


structModel.K = Ksolved;

[Psi_solved, lambdasolved] = eigs(Ksolved, structModel.M0,nModes,1e3);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);

strut_nodes1 = nodes(coord(:,1) > .25);
strut_nodes2 = nodes(coord(:,2)<0);

structModel.D0(strut_nodes1,strut_nodes1) = strut_damping1*structModel.D0(strut_nodes1,strut_nodes1);
structModel.D0(strut_nodes2,strut_nodes2) = strut_damping2*structModel.D0(strut_nodes2,strut_nodes2);

for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*full(structModel.M0)*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(i,i) = Psi_solved(:,i)'*structModel.D0*Psi_solved(:,i);
end

% Dm(1,1) = Dm(1,1)*3.25;
% Dm(2,2) = Dm(2,2)*1.5;
%%
clear G

for i = 1:nModes

    Psi_squared = gain*Psi_solved(:,i)*Psi_solved(:,i)';

    if i > 1
        G = G+ Psi_squared([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);
    else
          G = Psi_squared([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);
    end
end

w = logspace(0,4,800);
%%%
[G_mag,G_pha, w] = bode(G, w);

%%
for loop = 1:5
for i = 1:5
    for j = 1:3
        if G_pha(i,j,1) > 200
            G_pha(i,j,:) = G_pha(i,j,:) - 360;
        end
    end
end
end

%% Initial model
nModes =30;




[Psi_solved, lambdasolved] = eigs(structModel.K0, structModel.M0,nModes,1e3);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);


for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*full(structModel.M0)*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(i,i) = Psi_solved(:,i)'*structModel.D0*Psi_solved(:,i);
end

Dm(1,1) = Dm(1,1)*3.5;
Dm(2,2) = Dm(2,2)*1.5;
%%
clear G

for i = 1:nModes

    Psi_squared = gain*Psi_solved(:,i)*Psi_solved(:,i)';

    if i > 1
        G = G+ Psi_squared([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);
    else
          G = Psi_squared([measDOFs(1) measDOFs(2) measDOFs(3) nonmeasDOFs(1) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);
    end
end

w = logspace(0,4,800);
%%%
[G0_mag,G0_pha, ffrf] = bode(G, w);

%%
for loop = 1:5
for i = 1:5
    for j = 1:3
        if G0_pha(i,j,1) > 200
            G0_pha(i,j,:) = G0_pha(i,j,:) - 360;
        end
    end
end
end
