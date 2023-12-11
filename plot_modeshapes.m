figure
plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),1), '*' , LineWidth=2)
hold on
plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),1), '.k', LineWidth=2)
plot(sensor_loc([1 2 3 5]), psi_m(1:4,1), '*', 'MarkerSize',15, LineWidth=2)
plot(sensor_loc(4), L(4,1), '*', 'MarkerSize',15, LineWidth=2)
legend('Solved', 'init', 'Measured', 'Validation')
ylim([-1.2 1.2])
title('First Mode')


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
if n_modes >4
    figure
    plot([coord(output_nodes,1)+.25], Psi_solved(node_mapping.z(output_nodes),5), '*' , LineWidth=2)
    hold on
    plot([coord(output_nodes,1)+.25], Psi_o(node_mapping.z(output_nodes),5), '.k', LineWidth=2)
    plot(sensor_loc([1 2 3 5]), psi_m(1:4,5), '*', 'MarkerSize',15, LineWidth=2)
    plot(sensor_loc(4), L(4,5), '*', 'MarkerSize',15, LineWidth=2)
    legend('Solved', 'init', 'Measured', 'Validation')
    ylim([-1.5 1.5])
    title('Fifth Mode')
end
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
