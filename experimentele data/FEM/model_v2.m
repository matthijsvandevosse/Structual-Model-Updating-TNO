clear all
close all


E = 200e9;
nu = 0.27;
rho = 8e3;


structuralmodelinit = createpde(structural = "modal-planestrain");
structuralmodelspring1 = createpde(structural = "modal-planestrain");
structuralmodelspring2 = createpde(structural = "modal-planestrain");
structuralmodelspringmotor = createpde(structural = "modal-planestrain");
structuralmodelE = createpde(structural = "modal-planestrain");
structuralmodelnu = createpde(structural = "modal-planestrain");
structuralmodelE1 = createpde(structural = "modal-planestrain");
structuralmodelE2 = createpde(structural = "modal-planestrain");
structuralmodelE3 = createpde(structural = "modal-planestrain");
structuralmodelE4 = createpde(structural = "modal-planestrain");
structuralmodelE5 = createpde(structural = "modal-planestrain");
structuralmodelE6 = createpde(structural = "modal-planestrain");

l = 500/1000;
h = 22/1000;
t = 2/1000;

% sensor_loc = [50 150 250 350 450]/1000;
sensor_loc = [55 155 255 355 455]/1000;
sensor_loc = [45 145 245 345 445]/1000;

% For a rectangle, the first row contains 3, and the second row contains 4. 
% The next four rows contain the x-coordinates of the starting points of 
% the edges, and the four rows after that contain the y-coordinates of 
% the starting points of the edges.


%% Geometry
make_geometry

%%
pdegplot(geometry, "VertexLabels",'on')
%%
pdegplot(geometry, "EdgeLabels","on", "FaceLabels","on")

%%
structuralmodelinit.Geometry = geometry;
structuralmodelspring1.Geometry = geometry;
structuralmodelspring2.Geometry = geometry;
structuralmodelspringmotor.Geometry = geometry;
structuralmodelE.Geometry= geometry;
structuralmodelE1.Geometry= geometry;
structuralmodelE2.Geometry= geometry;
structuralmodelE3.Geometry= geometry;
structuralmodelE4.Geometry= geometry;
structuralmodelE5.Geometry= geometry;
structuralmodelE6.Geometry= geometry;
structuralmodelnu.Geometry =geometry;

%% Meshing
hmax = t*4-0.000001;
msh = generateMesh(structuralmodelinit,Hmax=hmax);
generateMesh(structuralmodelspring1,Hmax=hmax);
generateMesh(structuralmodelspring2,Hmax=hmax);
generateMesh(structuralmodelspringmotor,Hmax=hmax);
generateMesh(structuralmodelE,Hmax=hmax);
generateMesh(structuralmodelE1,Hmax=hmax);
generateMesh(structuralmodelE2,Hmax=hmax);
generateMesh(structuralmodelE3,Hmax=hmax);
generateMesh(structuralmodelE4,Hmax=hmax);
generateMesh(structuralmodelE5,Hmax=hmax);
generateMesh(structuralmodelE6,Hmax=hmax);

generateMesh(structuralmodelnu,Hmax=hmax);
figure
pdeplot(msh,"NodeLabels","on")

%% Create alternative models

structuralProperties(structuralmodelinit,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu)
structuralProperties(structuralmodelspring1,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);
structuralProperties(structuralmodelspring2,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);

structuralProperties(structuralmodelspringmotor,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);


structuralProperties(structuralmodelE,YoungsModulus=.25.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);
structuralProperties(structuralmodelE1,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 1);
structuralProperties(structuralmodelE1,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 2:8);
structuralProperties(structuralmodelE2,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 2);
structuralProperties(structuralmodelE2,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= [1 3:8]);
structuralProperties(structuralmodelE3,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 3);
structuralProperties(structuralmodelE3,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= [1:2 4:8]);
structuralProperties(structuralmodelE4,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 4);
structuralProperties(structuralmodelE4,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= [1:3 5:8]);

structuralProperties(structuralmodelE5,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 5);
structuralProperties(structuralmodelE5,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= [1:4 6:8]);

structuralProperties(structuralmodelE6,YoungsModulus=.6.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= 6);
structuralProperties(structuralmodelE6,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu, ...
                                     Face= [1:5 7:8]);

structuralProperties(structuralmodelnu,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=.9*nu);


%% Boundaries

stiffness = [0; 1e4];
structuralBoundaryLoad(structuralmodelinit, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelinit, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)

structuralBoundaryLoad(structuralmodelE, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE1, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE1, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE2, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE2, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE3, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE3, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE4, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE4, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE5, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE5, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelE6, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE6, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelspringmotor, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
% structuralBoundaryLoad(structuralmodelspringmotor, ...
%                               "Edge",[3; 5; 13], ...
%                               "TranslationalStiffness",0.01*stiffness)
structuralBoundaryLoad(structuralmodelnu, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelnu, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)

structuralBoundaryLoad(structuralmodelspring1, ...
                              "Edge",[2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring1, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)
structuralBoundaryLoad(structuralmodelspring2, ...
                              "Edge",1, ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring2, ...
                              "Edge",[3; 5; 13], ...
                              "TranslationalStiffness",10*stiffness)



structuralBoundaryLoad(structuralmodelinit, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE1, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE2, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE3, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE4, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE5, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE6, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelnu, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelspring1, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelspring2, ...
                              "Edge",15, ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelspringmotor, ...
                              "Edge",15, ...
                              "TranslationalStiffness",[1e20; 0])




%% Results

modalresultsinit = solve(structuralmodelinit,FrequencyRange=[1,1e6]);
% modalresultsinit = solve(structuralmodelinit,FrequencyRange=[1,1e6]);
% figure
% pdeplot(modalresultsinit.Mesh,XYData=modalresultsinit.ModeShapes.y(:,1));
% axis equal
% pdeplot3D(modalresultsinit.Mesh, ColorMapData = modalresultsinit.ModeShapes.y(:,1))
% pdeplot3D(msh)
% title(['First Mode with Frequency ', ...
        % num2str(modalresultsinit.NaturalFrequencies(1)/(2*pi)),' Hz'])

%    2.6902     5.1000   32.6888   51.8581   92.9761



%% Init Matrices
FEMs = assembleFEMatrices(structuralmodelinit, 'stiff-spring');
Kinit = full(FEMs.Ks);
Minit = full(FEMs.M);

Minitx = Minit(1:size(Minit,1)/2,1:size(Minit,1)/2);
Myinit = Minit(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:end);
Kyinit = Kinit(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
%%
% 
% 
% 
% [Psi_solved, lambdasolved] = eig(Minit\Kinit);
% 
% 
% 
% [lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
% lambdasolved = lambdasolved(1:5);
% Psi_solved = Psi_solved(:, dummyInd(1:5));
% 
% for i = 1:5
%     [~,index] = max(abs(Psi_solved(:,i)));
%     Psi_solved(:,i) = Psi_solved(:,i) / Psi_solved(index,i);
% end


%% Spring1 Matrices
FEMs = assembleFEMatrices(structuralmodelspring1,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
Myspring1 = M(size(Minit,1)/2+1:end,size(Minit,1)/2+1:end);
Kyspring1 = K(size(Kinit,1)/2+1:end,size(Kinit,1)/2+1:end);
Mspring1 = M;
Kspring1 = K;
%% Spring2 Matrices
FEMs = assembleFEMatrices(structuralmodelspring2,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
Myspring2 = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
Kyspring2 = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
Mspring2 = M;
Kspring2 = K;

%% E Matrices
FEMs = assembleFEMatrices(structuralmodelE,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
MyE = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
KyE = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
ME = M;
KE = K;

%% E1 Matrices
FEMs = assembleFEMatrices(structuralmodelE1,'stiff-spring');
KE1= full(FEMs.Ks);
ME1 = full(FEMs.M);

%% E2 Matrices
FEMs = assembleFEMatrices(structuralmodelE2,'stiff-spring');
KE2= full(FEMs.Ks);
ME2 = full(FEMs.M);
%% E3 Matrices
FEMs = assembleFEMatrices(structuralmodelE3,'stiff-spring');
KE3= full(FEMs.Ks);
ME3 = full(FEMs.M);
%% E4 Matrices
FEMs = assembleFEMatrices(structuralmodelE4,'stiff-spring');
KE4= full(FEMs.Ks);
ME4 = full(FEMs.M);
%% E5 Matrices
FEMs = assembleFEMatrices(structuralmodelE5,'stiff-spring');
KE5= full(FEMs.Ks);
ME5 = full(FEMs.M);
%% E6 Matrices
FEMs = assembleFEMatrices(structuralmodelE6,'stiff-spring');
KE6= full(FEMs.Ks);
ME6 = full(FEMs.M);

%% nu Matrices
FEMs = assembleFEMatrices(structuralmodelnu,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
Mynu = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
Kynu = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
Mnu = M;
Knu = K;

%%
FEMs = assembleFEMatrices(structuralmodelspringmotor,'stiff-spring');
Kmotor= full(FEMs.Ks);
Mmotor = full(FEMs.M);



Kdiff{1}= sparse(Kinit - Kspring1);
Kdiff{2} = sparse(Kinit - Kspring2);
Kdiff{3} = sparse(Kinit - KE1);
Kdiff{4} = sparse(Kinit - KE2);
Kdiff{5} = sparse(Kinit - KE3);
Kdiff{6} = sparse(Kinit - KE4);
Kdiff{7} = sparse(Kinit - KE5);
Kdiff{8} = sparse(Kinit - KE6);
Kdiff{9} = sparse(Kinit - Knu); 
Kdiff{10} = sparse(Kinit - KE);
% Kdiff{11} = sparse( Kinit - Kmotor);

% [Psi_solved, lambdasolved] = eig(Minit\Kinit);
% 
% [lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
% lambdasolved = lambdasolved(1:6);
% Psi_solved = Psi_solved(:, dummyInd(1:6));

%% Nodes

TopNodes = structuralmodelinit.Mesh.Nodes(1,:) == 0.022;
