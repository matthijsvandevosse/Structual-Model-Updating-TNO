%% Find locations of file.full
% Location of nodes file
dir_FEM = '/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/IRTF_DM_V2';
addpath(dir_FEM)
% Location of Ansys work files
dir_FEM = '/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/IRTF_DM_V2/v12_files/dp0';
addpath(dir_FEM)

% Find all subfolders
dir_folders = dir(fullfile(dir_FEM));

% Find all file.full in subfolders
FILES = [];
for ii = 1:length(dir_folders)
    ii;
    f =dir(fullfile([dir_folders(ii).folder '/' dir_folders(ii).name ],'/MECH/file.full'));
    if ~isempty(f)
        FILES{length(FILES)+1} = f;
    end
end


%% Find input and output nodes

filename = 'IRTF_DM_nodes.txt';
[nodes,coord] = ansysnodes(filename);
[nodes, sort_i] =sort(nodes, 'ascend');
coord = coord(sort_i,:);

% Rotate coordinate system to measurement setup
coord(:,1:2) = -flip(coord(:,1:2),2);
nodes = nodes-min(nodes)+1;

% Nodes number at the front of the mirror
filename = 'Mirror_nodes.txt';
[output_nodes_nr, coord_output_nodes] = ansysnodes(filename);

% Find in which row of the matrix the nummber corresponds to
for ii = 1:length(output_nodes_nr)
    output_nodes(ii) = find(output_nodes_nr(ii)==nodes);
end

% Nodes number at the top of the actuators
filename = 'Actuator_nodes.txt';
[input_nodes_top_nr, ~] = ansysnodes(filename);

% Find in which row of the matrix the nummber corresponds to
for ii = 1:length(input_nodes_top_nr)
    input_nodes_top(ii) = find(input_nodes_top_nr(ii)==nodes);
end

% -4.39480e-2

% Nodes number at the top of the actuators
filename = 'Backplate_nodes.txt';
[input_nodes_back_nr,~] = ansysnodes(filename);

% Find in which row of the matrix the nummber corresponds to
for ii = 1:length(input_nodes_back_nr)
    input_nodes_back(ii) = find(input_nodes_back_nr(ii)==nodes);
end

load('INFO_sensors_actuators.mat')

Actuator_pos = Actuator_pos(:,1:2)/1000;

for ii = Acts
    [min_check_b(ii), act_loc_b(ii)] = min( sum(abs(coord(input_nodes_back,1:2) - Actuator_pos(ii,1:2)),2));
    [min_check(ii), act_loc_t(ii)] = min( sum(abs(coord(input_nodes_top,1:2) - Actuator_pos(ii,1:2)),2));
end

input_nodes_back = input_nodes_back(act_loc_b);
input_nodes_top = input_nodes_top(act_loc_t);

Kdiff = [];
Mdiff = [];
%% Retrieve Stiffness and Mass matrix
if load_matrices
    for ii = 1:length(FILES)
        % Only retreive Mass matrix first time
        if ii ==1 || ii == 8
            retreive_mass = 1;
        else
            retreive_mass = 0;
        end
        
        tic
        [Stiff, Mass, node_mapping] = ansys2Structural([FILES{ii}.folder '/' 'file.full'], retreive_mass);
        toc
        
        Stiff_cleaned = Stiff;
        Mass_cleaned = Mass;
        
        % remove "lose" nodes
        lose_DOFS = find(abs(diag(Stiff_cleaned))==0);
        
        Mass_cleaned(abs(diag(Stiff))==0,:) = [];
        Mass_cleaned(:,abs(diag(Stiff))==0) = [];
        
        
        Stiff_cleaned(abs(diag(Stiff))==0,:) = [];
        Stiff_cleaned(:,abs(diag(Stiff))==0) = [];
        
        % Number of DOF removed
        DOFS_removed(ii) = length(lose_DOFS);
        
        for jj = 1:DOFS_removed(ii)
            
            % Row and Column of lose node is removed and will not be mapped
            % anymore
            
            node_mapping.x(node_mapping.x == lose_DOFS(jj)) = "deleted node";
            node_mapping.y(node_mapping.y == lose_DOFS(jj)) = "deleted node";
            node_mapping.z(node_mapping.z == lose_DOFS(jj)) = "deleted node";
            
            % Because the Row and Columns are removed, the nodes after the
            % removed column and row will need to point to an index lower.
            
            node_mapping.x(node_mapping.x>lose_DOFS(jj)) = node_mapping.x(node_mapping.x>lose_DOFS(jj))-1;
            node_mapping.y(node_mapping.y>lose_DOFS(jj)) = node_mapping.y(node_mapping.y>lose_DOFS(jj))-1;
            node_mapping.z(node_mapping.z>lose_DOFS(jj)) = node_mapping.z(node_mapping.z>lose_DOFS(jj))-1;
            
            % The remaining DOFs also map now to an index lower
            if jj < DOFS_removed(ii)-1
                lose_DOFS(jj+1:end) = lose_DOFS(jj+1:end)-1;
            end
        end
        
        
        % Save initial stiffness and mass matrix and the delta stiffness matrix
        if ii == 1
            Kinit = Stiff_cleaned;
            Minit = Mass_cleaned;
        else
            if any(any(Stiff_cleaned - Kinit))
                Kdiff{length(Kdiff)+1} = Stiff_cleaned - Kinit;
            end
            
            if retreive_mass
                Mdiff{length(Mdiff)+1} = Mass_cleaned - Minit;
            end
        end
    end
    save((["IRTF_DM_K_M_Matrices_" + date + ".mat"]), "Kdiff", "Kinit", "Mdiff", "Minit", "node_mapping" )

end



return
ModalForm_IRTF

return
%% Find resonant frequencies and modeshapes
nModes =28;
Mode_hz = [526; 1398.7; 1420.8; 1476.8; 1592.6]; %hz
Mode_nr = [8; 16; 18; 21; 28];
first_mode = 150;

tic
[Psi_solved, lambdasolved] = eigs(Kinit, Minit, nModes, (first_mode*2*pi)^2, 'IsSymmetricDefinite', "Display", 0);
% [Psi_solved, lambdasolved] = eigs(Kinit, Minit, nModes, lambdasolved(1), 'IsSymmetricDefinite', 1, 'Tolerance',1e-4, "Display", 0, "StartVector",Psi_solved(:,1));

lambdasolved = diag(lambdasolved);
toc
omega_hz = real(sqrt(lambdasolved)/2/pi)

% Kdiff{2} = 1e2*Kdiff{2};
% %%
% figure
% plot(Psi_solved(node_mapping.z(output_nodes),1))

% clearvars -except Kdiff Kinit Minit output_nodes node_mapping coord

%%
close all
displaymode = 24;

clear x y

% omega_hz(displaymode)
figure
for i = 1:length(output_nodes)
    x(i) = coord(( output_nodes(i)),1);
    y(i) =coord(( output_nodes(i)),2);
end
plot(x,y,'*')

T = delaunay(x,y);
T = unique(T,'rows')
trisurf(T, x, y, -Psi_solved(node_mapping.z(output_nodes), displaymode))
hold on
xlabel("X")
ylabel("Y")
% title(['mode at ' num2str(omega_hz(displaymode)) 'hz' ])

%%

clear x y

omega_hz(displaymode)
figure
for i = 1:length(input_nodes_back)
    x(i) = coord(( input_nodes_back(i)),1);
    y(i) =coord(( input_nodes_back(i)),2);
end
plot(x,y,'*')

%%
T = delaunay(x,y);
T = unique(T,'rows')
trisurf(T, x, y, -Psi_solved(node_mapping.z(input_nodes_top), displaymode) + -Psi_solved(node_mapping.z(input_nodes_back), displaymode))
hold on
xlabel("X")
ylabel("Y")
% title(['mode at ' num2str(omega_hz(displaymode)) 'hz' ])

%%
displaymode = 4;
x = Actuator_pos(:,1);
y = Actuator_pos(:,2);
T = delaunay(x,y);
T = unique(T,'rows');
trisurf(T, x, y, -simModes.psi_m(:, displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['mode at ' num2str(sqrt(simModes.lambda(displaymode))/2/pi) 'hz' ])
%%

clear x y

omega_hz(displaymode)
figure
for i = 1:length(input_nodes_back)
    x(i) = coord(( input_nodes_back(i)),1);
    y(i) =coord(( input_nodes_back(i)),2);
end
plot(x,y,'*')

T = delaunay(x,y);
T = unique(T,'rows')
trisurf(T, x, y,  -Psi_solved(node_mapping.z(input_nodes_back), displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['mode at ' num2str(omega_hz(displaymode)) 'hz' ])
