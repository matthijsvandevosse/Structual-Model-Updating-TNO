clear all
close all


E = 100e9;
nu = 0.29;
rho = 7e3;


structuralmodelinit = createpde(structural = "modal-planestrain");
structuralmodelspring1 = createpde(structural = "modal-planestrain");
structuralmodelspring2 = createpde(structural = "modal-planestrain");
structuralmodelspring3 = createpde(structural = "modal-planestrain");
structuralmodelspring4 = createpde(structural = "modal-planestrain");
structuralmodelE = createpde(structural = "modal-planestrain");
structuralmodelnu = createpde(structural = "modal-planestrain");

l = 500/1000;
h = 22/1000;
t = 2/1000;


sensor_loc = [55 155 255 355 455]/1000;
% sensor_loc = sensor_loc-10/1000;
% sensor_loc = [125 250 375]/1000;
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
structuralmodelspring3.Geometry = geometry;
structuralmodelspring4.Geometry = geometry;
structuralmodelE.Geometry= geometry;
structuralmodelnu.Geometry =geometry;

%% Meshing
hmax = t*4-0.000001;
msh = generateMesh(structuralmodelinit,Hmax=hmax);
generateMesh(structuralmodelspring1,Hmax=hmax);
generateMesh(structuralmodelspring2,Hmax=hmax);
generateMesh(structuralmodelspring3,Hmax=hmax);
generateMesh(structuralmodelspring4,Hmax=hmax);
generateMesh(structuralmodelE,Hmax=hmax);
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
structuralProperties(structuralmodelspring3,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);
structuralProperties(structuralmodelspring4,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);


structuralProperties(structuralmodelE,YoungsModulus=.1.*E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=nu);
structuralProperties(structuralmodelnu,YoungsModulus=E, ...
                                     MassDensity=rho, ... 
                                     PoissonsRatio=.9*nu);


%% Boundaries

stiffness = [0; 5e4];
structuralBoundaryLoad(structuralmodelinit, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelE, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelnu, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring1, ...
                              "Edge",[2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring2, ...
                              "Edge",1, ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring3, ...
                              "Edge",[1; 2], ...
                              "TranslationalStiffness",stiffness)
structuralBoundaryLoad(structuralmodelspring3, ...
                              "Edge",[1;], ...
                              "TranslationalStiffness",0.5*stiffness)


structuralBoundaryLoad(structuralmodelinit, ...
                              "Edge",[15], ...
                              "TranslationalStiffness",[1e20; 0])
structuralBoundaryLoad(structuralmodelE, ...
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
structuralBoundaryLoad(structuralmodelspring3, ...
                              "Edge",[15], ...
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



[Psi_solved, lambdasolved] = eig(Myinit\Kyinit);



[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
lambdasolved = lambdasolved(1:5);
Psi_solved = Psi_solved(:, dummyInd(1:5));

for i = 1:5
    [~,index] = max(abs(Psi_solved(:,i)));
    Psi_solved(:,i) = Psi_solved(:,i) / Psi_solved(index,i);
end


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
%% Spring3 Matrices
FEMs = assembleFEMatrices(structuralmodelspring3,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
Myspring3 = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
Kyspring3 = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
Mspring3 = M;
Kspring3 = K;

%% E Matrices
FEMs = assembleFEMatrices(structuralmodelE,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
MyE = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
KyE = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
ME = M;
KE = K;
%% nu Matrices
FEMs = assembleFEMatrices(structuralmodelnu,'stiff-spring');
K = full(FEMs.Ks);
M = full(FEMs.M);
Mynu = M(size(Minit,1)/2+1:2*size(Minit,1)/2,size(Minit,1)/2+1:2*size(Minit,1)/2);
Kynu = K(size(Kinit,1)/2+1:2*size(Kinit,1)/2,size(Kinit,1)/2+1:2*size(Kinit,1)/2);
Mnu = M;
Knu = K;

Kdiff1 = Kinit - Kspring1;
Kdiff2 = Kinit - Kspring2;
Kdiff3 = Kspring3 - Kinit;
Kdiff4 = Kinit - KE;
Kdiff5 = Kinit - Knu;


[Psi_solved, lambdasolved] = eig(Minit\Kinit);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
lambdasolved = lambdasolved(1:6);
Psi_solved = Psi_solved(:, dummyInd(1:6));

% for i = 1:3
%     [~,index] = max(abs(Psi_solved(:,i)));
%     Psi_solved(:,i) = Psi_solved(:,i) / abs(Psi_solved(index,i));
% end
%%
naturalfrequency = sqrt(lambdasolved)/(2*pi)

%%
figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit,1)/2+1:end,2))
axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit)/2+1:end,2))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit)/2+1:end,3))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit)/2+1:end,4))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit)/2+1:end,5))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(size(Kinit)/2+1:end,6))
% axis equal

% 
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,1))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,2))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,3))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,4))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,5))
% axis equal
% figure
% pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:size(Kinit)/2,6))
% axis equal