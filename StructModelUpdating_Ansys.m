clc;clear all
% close all
warning('off')

Retrieve_K_M_Ansys

pause(0.5)
%%

N = size(Kinit,1);

% measurement points -- sensor_loc = [45 145 245 345 445]/1000;
measDOFs = N/2+[6; 8; 10; 14]; 

actDOFs = [21+N/2 141+N/2 233+N/2];

nonmeasDOFs =  N/2+[12];

%%
PoleUpdating
%%
% PoleZeroUpdating

%%

[lambdasolved_sorted,dummyInd] = sort(diag(lambdasolved), 'ascend');



[lambda_o,dummyInd_o] = sort(diag(lambda_o), 'ascend');

Psi_solvedxy = Psi_solvedxy(:,dummyInd);
Psi_solved = Psi_solvedxy(N/2+1:end,:);
Psi_o = Psi_o(:, dummyInd_o);
Psi_o = Psi_o(N/2+1:end, :);


naturalfrequency = sqrt(lambdasolved_sorted(modeIndex))/(2*pi)
freqExp
%%
sign = [-1 1 1 1 1 1 -1 1 1 1 1];
for i = modeIndex 

    [~, index] = max(abs(Psi_solved(measDOFs(1:4)-N/2,i)));
    Psi_solved(:,i) = sign(i)* Psi_solved(:,i) / (Psi_solved(measDOFs(index)-N/2,i));

    [~, index] = max(abs(Psi_o(measDOFs(1:4)-N/2,i)));
    Psi_o(:,i) = sign(i)* Psi_o(:,i) / (Psi_o(measDOFs(index)-N/2,i));
end


Psi_solved(nonmeasDOFs-N/2,modeIndex)
Psi_measured =  L(4,modeIndex)



%%
figure
plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,1), '*' , LineWidth=2)
hold on
plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,1), '.k', LineWidth=2)
plot(sensor_loc([1 2 3 5]), psi_m(1:4,1), '*', 'MarkerSize',15, LineWidth=2)
plot(sensor_loc(4), L(4,1), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1.2 1.2])
title('First Mode')

if n_modes >1
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,2), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,2), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,2), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,2), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.2 1.2])
    title('Second Mode')

end
if n_modes >2
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,3), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,3), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,3), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,3), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.2 1.2])
    title('Third Mode')
end
if n_modes >3
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,4), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,4), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,4), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,4), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Fourth Mode')
end
if n_modes >4
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,5), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,5), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,5), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,5), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Fifth Mode')
end
if n_modes >5
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,6), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,6), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,6), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,6), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Sixth Mode')
end
if n_modes >6
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,7), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,7), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,7), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,7), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Seventh Mode')
end
if n_modes >7
    figure
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_solved(:,8), '*' , LineWidth=2)
    hold on
    plot([structuralmodelinit.Mesh.Nodes(1,:)], Psi_o(:,8), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), L([1 2 3 5],8), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,8), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Eigth Mode')
end


pause(0.5)

%%

x_1_K = [x_2];
x_1_M = 0;


structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);
Ksolved = full(structModel.K0);
Msolved = full(structModel.M0);

% for n = 1:n_alpha
%     Ksolved = Ksolved + (x_1_K(n))*full(structModel.K_j{n});
% end
% for n = 1:n_beta
%     Msolved = Msolved + (x_1_M(n))*full(structModel.M_j{n});
% end

structModel.K = Ksolved;
structModel.M = Msolved;

Kp = Kinit;
for n = 1:n_alpha
    Kp = Kp + (x_2(n))*full(structModel.K_j{n});
end

%% x = [position; velocity]

Minv_o = inv(full(structModel.M0));
Minv = inv(structModel.M);
C_0 = [zeros(5,2*N)];
C_0(1,measDOFs(1)) =      35000;
C_0(2,measDOFs(2)) =      31000;
C_0(3,measDOFs(3)) =      31000;
C_0(4,nonmeasDOFs(1)) =   31000;
C_0(5,measDOFs(4)) =      31000;

C = [zeros(5,2*N)];
C(1,measDOFs(1)) =      30000+6000*x_1_M(1);
C(2,measDOFs(2)) =      30000+6000*x_1_M(1);
C(3,measDOFs(3)) =      30000+6000*x_1_M(1);
C(4,nonmeasDOFs(1)) =   30000+6000*x_1_M(1);
C(5,measDOFs(4)) =      30000+6000*x_1_M(1);


%%
structModel.D0 = 8e-2.*eye(N,N);
structModel.D = 8e-2.*eye(N,N);
% structModel.D0 = 0e-2.*eye(N,N);

%%
A = [zeros(N,N) eye(N); -inv(structModel.M)*structModel.K -inv(structModel.M)*structModel.D];

A_o = [zeros(N,N) eye(N); -inv(structModel.M0)*structModel.K0 zeros(N,N)];

A_p = [zeros(N,N) eye(N); -inv(structModel.M0)*Kp -inv(structModel.M0)*structModel.D0];


B = [zeros(N,3); Minv(:,actDOFs)];
B_o = [zeros(N,3); Minv_o(:,actDOFs)];
load("Modal_Systeem_V1.mat")
sys1 = ss(A,B,C,0);
syso = ss(A_o,B_o,C_0,0);
sysp = ss(A_p,B_o,C_0,0);

%%
w = logspace(0,4,800);
h = freqresp(sys1, w);
% h_o = freqresp(syso, w);
h_p = freqresp(sysp, w);
hm = freqresp(Gmod7, w);
%%
load('Modal_Systeem_V1.mat')
load('BeamModal_O12345I123.mat')

data{1} = ModalFitData;
G_ref = data{1}.G_ref;

[Gfrf_mag,Gfrf_pha   ,wfrf]   = bode(G_ref);
ffrf = wfrf;


load('ResponseData0_04')
% Remove phase delay
for i = 1:length(G_bla)
    new_angle = angle(G_bla.ResponseData(:,:,i))+(1.95/1000)*ones(5,3).*G_bla.Frequency(i)*2*pi;
    G_bla.ResponseData(:,:,i) = abs(G_bla.ResponseData(:,:,i)).*exp(1i*new_angle);
end

[Gfrf_mag,Gfrf_pha   ,ffrf]   = bode(G_bla);

%%
data_opt = data{1,1};
%%
% LeastSquaredOpt
% structModel.D = 5.8e-2*eye(N);
% structModel.D = 5e-3.*eye(N,N);
%%
% A_new = [zeros(N) eye(N); -inv(structModel.M)*structModel.K -inv(structModel.M)*(structModel.D)];
% 
% sys3 = ss(A_new,B,C,0);
% 
% h2 = freqresp(sys3, w);

%%

% h2_frf = frd(h2, w);
% h4_frf = frd(h4, w);

h_frf = frd(h, w);
hp_frf = frd(h_p, w);
% ho_frf = frd(h_o, w);
hm_frf = frd(hm, w);
%%

% [h2_mag,h2_pha   ,wfrf]   = bode(h2_frf,w);
% [h4_mag,h4_pha   ,wfrf]   = bode(h4_frf,w);
[hp_mag,hp_pha   ,wfrf]   = bode(hp_frf,w);
[h_mag,h_pha   ,wfrf]   = bode(h_frf,w);
% [ho_mag,ho_pha   ,wfrf]   = bode(ho_frf);
[hm_mag,hm_pha   ,wfrf]   = bode(hm_frf);


%%
f = figure(24);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        set(gca,'ytick',[-360,-180,0])


        plot(ffrf,squeeze((Gfrf_pha(i,j,:)))+40/100*ffrf/2/pi, 'b',   'LineWidth', 1.5)
        % plot(w,(squeeze(ho_pha(i,j,:))), 'k--', 'LineWidth',1.5)
        % plot(w,(squeeze(hm_pha(i,j,:))), 'LineWidth',1.2)
        % plot(w,(squeeze(h_pha(i,j,:))), 'g', 'LineWidth',1.5)
        % plot(w,(squeeze(h2_pha(i,j,:))), 'y', 'LineWidth',2)
        plot(w,(squeeze(hp_pha(i,j,:))), 'r', 'LineWidth',1.5)


        if i ==1 & j == 1
        lgd = legend(["FRF", "Mass updated", "Updated",],Location="southwest");
        end
        lgd.FontSize = 18;  % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i==3 && j==2;
            xlabel('f [Hz]')
        elseif i~=3
            %             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==3
            ylabel('Phase [^\circ]')
        elseif j~=1
            %             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
              xlim([1,2000])
        ylim([-380,0])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        % fig = get(groot,'CurrentFigure');
        % set(fig, 'PaperUnits', 'centimeters');
        % set(fig, 'Units', 'centimeters');
        % x_width=8.85;          %x_width of the figure, set at template column width
        % y_width=0.75*x_width;   %y_width of the figure, free to choose
        % set(fig, 'PaperPosition', [0 0 x_width y_width]);
        % set(fig, 'PaperSize', [x_width y_width]);
        % set(fig, 'InnerPosition', [0 0 x_width y_width]);
        % set(gca,'FontSize',9)

    end
end
%%
Out = [1 1 2 2 3 3];
In = [1 1 2 2 3 3];
f = figure(23);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-50,0,50])


        plot(ffrf,squeeze(mag2db(Gfrf_mag(i,j,:))), 'b',   'LineWidth', 1.5)
        % plot(w,mag2db(squeeze(ho_mag(Out(i),Inp(j),:))), 'k--', 'LineWidth',1.5)

        % plot(w,mag2db(squeeze(h2_mag(Out(i),Inp(j),:))), 'LineWidth',1.5)
        plot(w,mag2db(squeeze(hp_mag(i,j,:))), 'r','LineWidth',1.5)
        % plot(w,mag2db(squeeze(h_mag(i,j,:))),'g', 'LineWidth',2 )



        % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i ==1 & j == 1
        lgd = legend(["FRF", "Init", "Updated",],Location="southwest");
        end
        lgd.FontSize = 18;
        % legend(["FRF", "Updated Zeros and Poles",],Location="southwest")
        if i==3 && j==2;
            xlabel('f [Hz]')
        elseif i~=3
            xticklabels({})
        end
        if j==1 && i==3
            ylabel('Mag [dB]')
        elseif j~=1
            yticklabels({})
        end
        xlim([1,2000])
        ylim([-75,80])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on

    end
end

