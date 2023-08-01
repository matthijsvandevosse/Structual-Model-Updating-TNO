clc;clear
close all
warning('off')

model_v2

pause(0.5)

%% Load Sensitivity matrix

modeIndex = [1 2 3 4 5]; % Indexes of these measured modes
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
measDOFs = N/2+[6; 8; 10; 14];
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
freqExp(4) = eigenFrequencies(5);
freqExp(5) = 170;
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")
L(:,4) = L(:,5);
for i = 1:modeIndex(end)
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
% psiExpAll([5], n) = 0.9999*L(1,i);
% psiExpAll([6], n) = 0.9999*L(2,i);
% psiExpAll([7], n) = 0.9999*L(3,i);
% psiExpAll([8], n) = 0.9999*L(5,i);
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
expModes.lambdaWeights = 2*[10 10 10 10 1];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [1000 200 400 400 0.1];
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
 
[Psi_o, lambda_o] = eig(full(structModel.M0)\structModel.K0);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');
[lambda_o,dummyInd_o] = sort(diag(lambda_o), 'ascend');
% lambdasolved = lambdasolved(modeIndex);
Psi_solved = Psi_solvedxy(N/2+1:end, dummyInd(1:5));
Psi_o = Psi_o(N/2+1:end, dummyInd_o(1:5));
%%
naturalfrequency = sqrt(lambdasolved(1:5))/(2*pi)
freqExp
%%
sign = [-1 1 1 -1 -1 1 1];
for i = 1:modeIndex(end)
    [~, index] = max(abs(Psi_solved(measDOFs(1:4)-N/2,i)));
    Psi_solved(:,i) = sign(i)* Psi_solved(:,i) / (Psi_solved(measDOFs(index)-N/2,i));
        
    [~, index] = max(abs(Psi_o(measDOFs(1:4)-N/2,i)));
    Psi_o(:,i) = sign(i)* Psi_o(:,i) / (Psi_o(measDOFs(index)-N/2,i));
end

Psi_solved(measDOFs(1:4)-N/2,modeIndex)
L([1 2 3 5],modeIndex)

Psi_solved(nonmeasDOFs-N/2,modeIndex)
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

figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,4))
axis equal

figure
pdeplot(modalresultsinit.Mesh,XYData=Psi_solved(1:end,5))
axis equal
%%
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,1), LineWidth=2)
hold on
plot(Psi_o([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,1), ':k', LineWidth=2)
plot([1 2 3 5], psi_m(1:4,1), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,1), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1 1])
title('First Mode')

if n_modes >1
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,2), LineWidth=2)
hold on
plot(Psi_o([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,2), ':k', LineWidth=2)
plot([1 2 3 5], psi_m(1:4,2), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,2), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1 1])
title('Second Mode')

end
if n_modes >2
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,3), LineWidth=2)
hold on
plot(Psi_o([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,3), ':k', LineWidth=2)
plot([1 2 3 5], psi_m(1:4,3), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,3), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1 1])
title('Third Mode')
end
if n_modes >3
figure
plot(Psi_solved([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,4), LineWidth=2)
hold on
plot(Psi_o([measDOFs(1:3); nonmeasDOFs(1);measDOFs(4)]-N/2,4), ':k', LineWidth=2)
plot([1 2 3 5], psi_m(1:4,4), '*', 'MarkerSize',15, LineWidth=2)
plot(4, L(4,4), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1.2 1.2])
title('Fourth Mode')
end

pause(0.5)

%% x = [position; velocity]
Minv = inv(structModel.M0);
C = [zeros(5,2*N)];
actDOFs = [237+N/2 427+N/2 188+N/2];
C(1,measDOFs(1)) = 30000;
C(2,measDOFs(2)) = 30000;
C(3,measDOFs(3)) = 30000;
C(4,nonmeasDOFs(1)) = 30000;
C(5,measDOFs(4)) = 30000;
A = [zeros(N,N) eye(N); -inv(structModel.M0)*Ksolved zeros(N,N)];
A_o = [zeros(N,N) eye(N); -inv(structModel.M0)*structModel.K0 zeros(N,N)];
B = [zeros(N,3); Minv(:,actDOFs)];
load("Modal_Systeem_V1.mat")
sys1 = ss(A,B,C,0);
syso = ss(A_o,B,C,0);
w = logspace(0,4,400);
h = freqresp(sys1, w);
h_o = freqresp(syso, w);
hm = freqresp(Gmod7, w);
%%
h_frf = frd(h, w);
ho_frf = frd(h_o, w);
hm_frf = frd(hm, w);

opts.PhaseWrapping = 'on';
[h_mag,h_pha   ,wfrf]   = bode(h_frf,w);
[ho_mag,ho_pha   ,wfrf]   = bode(ho_frf);
[hm_mag,hm_pha   ,wfrf]   = bode(hm_frf);

h_pha(h_pha>1) = h_pha(h_pha>1) - 360;

%%
load('Modal_Systeem_V1.mat')
load('BeamModal_O12345I123.mat')
Inp = [1,2,3];
Out = [1,2,3,4,5];
data{1} = ModalFitData;
G_ref = data{1}.G_ref;

[Gfrf_mag,Gfrf_pha   ,wfrf]   = bode(G_ref);
ffrf = wfrf/2/pi;
%%
f = figure(23);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:5
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-50,0,50])


            plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
            plot(w/(2*pi),mag2db(squeeze(h_mag(i,j,:))), 'LineWidth',1.2)
            % plot(w/(2*pi),mag2db(squeeze(hm_mag(i,j,:))), ':k', 'LineWidth',1.2)
        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
            xticklabels({})
        end
        if j==1 && i==2
            ylabel('Mag [dB]')
        elseif j~=1

            yticklabels({})
        end
        xlim([1,200])
        ylim([-75,80])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        fig = get(groot,'CurrentFigure');
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'Units', 'centimeters');
        x_width=8.85;          %x_width of the figure, set at template column width
        y_width=0.75*x_width;   %y_width of the figure, free to choose
        set(fig, 'PaperPosition', [0 0 x_width y_width]); 
        set(fig, 'PaperSize', [x_width y_width]); 
        set(fig, 'InnerPosition', [0 0 x_width y_width]);
        set(gca,'FontSize',9)

    end
end

f = figure(24);
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:4
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        set(gca,'ytick',[-360,-180,0])


            plot(ffrf,squeeze((Gfrf_pha( Out(i),Inp(j),:)))+40/100*ffrf, 'b',   'LineWidth', 1.2)
            plot(w/(2*pi),(squeeze(h_pha(i,j,:))), 'LineWidth',1.2)
            % plot(w/(2*pi),(squeeze(-rad2deg(angle(h(i,j,:))))), 'LineWidth',1.2)
            % plot(w/(2*pi),(squeeze(hm_pha(i,j,:))), ':k', 'LineWidth',1.2)

        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
%             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==2
            ylabel('Phase [$^\circ$]')
        elseif j~=1
%             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
        xlim([1,500])
        ylim([-400,0])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        fig = get(groot,'CurrentFigure');
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'Units', 'centimeters');
        x_width=8.85;          %x_width of the figure, set at template column width
        y_width=0.75*x_width;   %y_width of the figure, free to choose
        set(fig, 'PaperPosition', [0 0 x_width y_width]); 
        set(fig, 'PaperSize', [x_width y_width]); 
        set(fig, 'InnerPosition', [0 0 x_width y_width]);
        set(gca,'FontSize',9)

    end
end