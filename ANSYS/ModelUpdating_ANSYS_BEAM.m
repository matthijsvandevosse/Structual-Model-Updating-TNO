clc;clear all; close all
% close all
warning('off')
addpath("/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/BEAM FEM")
addpath(genpath("/Users/matthijsvandevosse/Structual-Model-Updating_matthijs"))

%%

Retrieve_K_M_Ansys

%%

N = size(Kinit,1);

% measurement points --
sensor_loc = [45 145 245 345 445]/1000;
update_sensor_loc = [45 145 245 445]/1000;
val_sensor_loc_x = 345/1000;




% act_loc_x = [50  250  450]/1000;
act_loc_x = [45  245  445]/1000;

[~, idx] = min(abs( coord(output_nodes,1)+0.25 - update_sensor_loc(:)'));
measDOFs = node_mapping.z(output_nodes(idx))
[~, idx] = min(abs( coord(output_nodes,1)+0.25 - act_loc_x(:)'));
actDOFs = node_mapping.z(output_nodes(idx))
[~, idx] = min(abs( coord(output_nodes,1)+0.25 - val_sensor_loc_x(:)'));
nonmeasDOFs = node_mapping.z(output_nodes(idx));


%%
PoleUpdating
%%
% PoleZeroUpdating

%%

[lambdasolved_sorted,dummyInd] = sort(diag(lambdasolved), 'ascend');
[lambda_o,dummyInd_o] = sort(diag(lambda_o), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);
Psi_o = Psi_o(:, dummyInd_o);

naturalfrequency = sqrt(lambdasolved_sorted(1:8))/(2*pi)
naturalfrequency_0 =  sqrt(lambda_o(1:8))/(2*pi)
freqExp
%%
for i = 1:length(modeIndex)

    [~, index] = max(abs(Psi_solved(measDOFs(1:4),i)));
    Psi_solved(:,i) = sign(dot(Psi_solved(measDOFs,i), psi_m(:,i)))* Psi_solved(:,i) / abs(Psi_solved(measDOFs(index),i));

    [~, index] = max(abs(Psi_o(measDOFs(1:4),i)));
    Psi_o(:,i) = sign(dot(Psi_o(measDOFs,i), psi_m(:,i)))* Psi_o(:,i) / abs(Psi_o(measDOFs(index),i));
end


Psi_solved(nonmeasDOFs,modeIndex)
Psi_measured =  L(4,modeIndex)



%%
plot_modeshapes

pause(0.5)

%%

x_1_K = [x_2];
x_1_M = 0;


structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

structModel.K = Ksolved;
structModel.M = Minit;

Kp = Kinit;
for n = 1:n_alpha
    Kp = Kp + (x_2(n))*full(structModel.K_j{n});
end

%%
structModel.D0 = 3.5e-6*eye(N);
gain = 800;
% StateSpaceRep
ModalForm

hp_mag = G_mag;
hp_pha = G_pha;
h0_mag = G0_mag;
h0_pha = G0_pha;
%
%%
load('Modal_Systeem_V1.mat')
load('BeamModal_O12345I123.mat')

data{1} = ModalFitData;
G_ref = data{1}.G_ref;

[Gfrf_mag,Gfrf_pha   ,wfrf]   = bode(G_ref);
ffrf = wfrf;

% Use new measurements
load('ResponseData0_04_corrected_2.mat')
% Remove phase delay
for i = 1:length(G_bla.Frequency)
    new_angle = angle(G_bla.ResponseData(:,:,i))+(1.95/1000)*ones(5,3).*G_bla.Frequency(i)*2*pi;
    G_bla.ResponseData(:,:,i) = abs(G_bla.ResponseData(:,:,i)).*exp(1i*new_angle);
end
%%

[Gfrf_mag,Gfrf_pha   ,ffrf]   = bode(G_bla);

[hm_mag , hm_pha ] = bode(G_ss_Modal, ffrf);

%%
% data_opt = data{1,1};
%%
% LeastSquaredOpt
% structModel.D = 5.8e-2*eye(N);
% structModel.D = 5e-3.*eye(N,N);

%%
f = figure(24);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        set(gca,'ytick',[-360,-180,0])


        plot(ffrf,squeeze((wrapTo180(Gfrf_pha(i,j,:)))), 'b',   'LineWidth', 1)
        % plot(w,(squeeze(ho_pha(i,j,:))), 'k--', 'LineWidth',1.5)
        plot(w,(squeeze(h0_pha(i,j,:))),':k', 'LineWidth',1)
        % plot(w,(squeeze(h_pha(i,j,:))), 'g', 'LineWidth',1.5)
        % plot(w,(squeeze(h2_pha(i,j,:))), 'y', 'LineWidth',2)
        plot(w,(squeeze(wrapTo180(hp_pha(i,j,:)))), 'r', 'LineWidth',1)


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
        ylim([-200,200])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
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


        plot(ffrf,squeeze(mag2db(Gfrf_mag(i,j,:))), 'b',   'LineWidth', 1)
        % plot(w,mag2db(squeeze(ho_mag(Out(i),Inp(j),:))), 'k--', 'LineWidth',1.5)

        plot(w,mag2db(squeeze(h0_mag(i,j,:))), ':k', 'LineWidth',1)
        plot(w,mag2db(squeeze(hp_mag(i,j,:))), 'r','LineWidth',1)
        % plot(ffrf,mag2db(squeeze(hm_mag(i,j,:))),'g', 'LineWidth',1 )

        % xline(sqrt(lambdasolved_sorted))

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
        xlim([5,5000])
        ylim([-75,80])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on

    end
end
