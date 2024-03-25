%% Normalize Mode shapes
n_alpha = length(structModel.K_j);
n_beta = length(structModel.M_j);

structModel.K = structModel.K0;
for i = 1 : n_alpha
    structModel.K = structModel.K + x(i) * structModel.K_j{i};
end

structModel.M = structModel.M0;
for i = 1 : n_beta
    structModel.M = structModel.M + x(i + n_alpha) * structModel.M_j{i};
end

% Updated mode shape and resonance frequencies
[Psi_solved, lambdasolved]  = eigs(structModel.K,  structModel.M, max(modeIndex),  expModes.lambdaExp(1)-100, 'IsSymmetricDefinite', 1, 'Tolerance', 1e-3);

% Initial mode shape and resonance frequencies
[Psi_o, lambda_o]  = eigs(structModel.K0,  structModel.M0, max(modeIndex),  expModes.lambdaExp(1)-100, 'IsSymmetricDefinite', 1, 'Tolerance', 1e-3);

% Sign matching of the simulated and measured mode shapes
for i = 1:length(modeIndex)
    [~, index] = max(abs(Psi_solved(measDOFs(:,1),i)));
    Psi_solved(:,i) = sign(dot(Psi_solved(measDOFs(:,1),i), psi_m(:,i)))* Psi_solved(:,i) / abs(Psi_solved(measDOFs(index,1),i));

    [~, index] = max(abs(Psi_o(measDOFs(:,1),i)));
    Psi_o(:,i) = sign(dot(Psi_o(measDOFs(:,1),i), psi_m(:,i)))* Psi_o(:,i) / abs(Psi_o(measDOFs(index,1),i));
end

%% Plotting for thesis

figure('Renderer', 'painters', 'Position', [10 10 400 225])
tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');

for i = [1 2 3 4 6]
  nexttile

hold on
plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),i), '.' , LineWidth=1, color = "#0072BD", MarkerSize = 10)
plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),i), '.k', LineWidth=1, MarkerSize =5)


plot(sensor_loc([1 2 3 5]), psi_m(1:4,i), '*', 'MarkerSize',10, LineWidth=2)
plot(sensor_loc(4), L(4,i), '*', 'MarkerSize',10, LineWidth=2)
% legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1.2 1.2])
title(num2str(round(freqExp(i),1)))
xlim([0 0.5])
        set(gca, 'Fontsize', 10)
        box on
        set(gca,'linewidth',.25)
end



%% Large Plotting
return

if n_modes >1
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),2), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),2), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,2), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,2), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.2 1.2])
    title('Second Mode')

end
if n_modes >2
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),3), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),3), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,3), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,3), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.2 1.2])
    title('Third Mode')
end
if n_modes >3
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),4), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),4), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,4), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,4), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Fourth Mode')
end

% 5 fifth mode not plotted as it is the out of plane mode
% if n_modes >4
%     figure
%     plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),5), '*' , LineWidth=2)
%     hold on
%     plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),5), '.k', LineWidth=2)
%     plot(sensor_loc([1 2 3 5]), psi_m(1:4,5), '*', 'MarkerSize',15, LineWidth=2)
%     plot(sensor_loc(4), L(4,5), '*', 'MarkerSize',15, LineWidth=2)
%     legend('Solved', 'init', 'Measured', 'Validation')
%     ylim([-1.5 1.5])
%     title('Fifth Mode')
% end

if n_modes >5
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),6), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),6), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,6), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,6), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Sixth Mode')
end
if n_modes >6
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),7), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),7), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,7), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,7), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Seventh Mode')
end
if n_modes >7
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),8), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),8), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), L([1 2 3 5],8), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,8), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Eigth Mode')
end


