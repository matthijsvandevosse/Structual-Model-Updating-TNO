clear all
%% Find locations of file.full
% Location of Ansys work files
dir_FEM = '/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/BEAM FEM/Beam_FEM_files/dp0';
addpath(dir_FEM)

% Find all subfolders
dir_folders = dir(fullfile(dir_FEM))

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
filename = 'beam_nodes.txt';
[nodes,coord] = ansysnodes(filename);
[nodes, sort_i] =sort(nodes, 'ascend');
coord = coord(sort_i,:);

% find nodes at the front of the beam
output_nodes = nodes(and(coord(:,2)== [8.5e-3], coord(:,3)== [0]));

% Remove nodes with same coordinates
[~,sort_i] = unique(coord(output_nodes,:), 'rows');
output_nodes = output_nodes(sort_i)

% Sort from left to right
[~,sort_i] = sort(coord(output_nodes(:,1),1), 'ascend');
output_nodes = output_nodes(sort_i,:);


%% Retrieve Stiffness and Mass matrix

for ii = 1:length(FILES)-6
    % Only retreive Mass matrix first time
    if ii ==1
        retreive_mass = 1;
    else
        retreive_mass = 0;
    end
    tic
    [Stiff, Mass, node_mapping] = ansys2Structural([FILES{ii}.folder '/' 'file.full'],nodes, retreive_mass);
    toc

    % remove "lose" nodes
    lose_DOFS = find(abs(diag(Stiff))==0);

    Mass_cleaned = Mass;
    Mass_cleaned(abs(diag(Stiff))==0,:) = [];
    Mass_cleaned(:,abs(diag(Stiff))==0) = [];

    Stiff_cleaned = Stiff;
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
        Kdiff{ii-1} = Stiff_cleaned - Kinit;
    end
end

%% Find resonant frequencies and modeshapes
tic
[Psi_solved, lambdasolved] = eigs(Kinit, Minit, 10, 1e3);
toc
omega_hz = real(sqrt(diag(lambdasolved))/2/pi)

Kdiff{2} = 1e2*Kdiff{2};
% %%
% figure 
% plot(Psi_solved(node_mapping.z(output_nodes),1))

% clearvars -except Kdiff Kinit Minit output_nodes node_mapping coord

