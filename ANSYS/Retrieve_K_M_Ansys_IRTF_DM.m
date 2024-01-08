clear all
%% Find locations of file.full
% Location of nodes file
dir_FEM = '/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/IRTF_DM';
addpath(dir_FEM)
% Location of Ansys work files
dir_FEM = '/Users/matthijsvandevosse/Library/CloudStorage/OneDrive-TUEindhoven/IRTF_DM/v12_files/dp0';
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
nodes = nodes-min(nodes)+1;

% find nodes at the front of the mirror
output_nodes = nodes(coord(:,3) == 0.021585);

input_nodes_mirror = nodes(coord(:,3)== 0.018335);




for ii = 1:length(input_nodes_mirror)
[min_d,idx] = min(sum(abs(coord - [  coord(find(nodes==input_nodes_mirror(ii)),1:2)  0.015085 ]),2));
input_nodes_backplate(ii,:) =    nodes(idx);
end

for ii = 1:length(input_nodes_mirror)
[min_d,idx] = min(sum(abs(coord - [  coord(find(nodes==input_nodes_mirror(ii)),1:2)  021585 ]),2));
frf_output(ii) = nodes(idx);
end


% % Remove nodes with same coordinates
% [~,sort_i] = unique(coord(output_nodes,:), 'rows');
% output_nodes = output_nodes(sort_i)
% 
% % Sort from left to right
% [~,sort_i] = sort(output_nodes, 'ascend');
% output_nodes = output_nodes(sort_i,:);


%% Retrieve Stiffness and Mass matrix

for ii = 1:length(FILES)
    % Only retreive Mass matrix first time
    if ii ==1
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
        Kdiff{ii-1} = Stiff_cleaned - Kinit;
    end
end
%%
deletednodes = find(isnan(node_mapping.z(:)));
output_nodes_cleaned = output_nodes;
for ii = 1:length(deletednodes)

   output_nodes_cleaned(find(output_nodes_cleaned == deletednodes(ii))) = [];
   % input_nodes_mirror(find(input_nodes_mirror == deletednodes(ii))) = [];

end



%% Find resonant frequencies and modeshapes
nModes =30;
first_mode = 150; %hz

tic
[Psi_solved, lambdasolved] = eigs(Kinit, Minit, nModes, (first_mode*2*pi)^2, 'IsSymmetricDefinite', 1, 'Tolerance',1e-4, "Display", 0);
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
displaymode = 1;
 

nodes_mirror = sqrt((coord(output_nodes_cleaned,[1])  +0.14688).^2 +(coord(output_nodes_cleaned,[2])  -0.16367).^2) < .0801;

output_nodes_mirror = output_nodes_cleaned(nodes_mirror);
output_nodes_cleaned = output_nodes_mirror;
clear x y

omega_hz(displaymode)
figure
for i = 1:length(output_nodes_cleaned)
x(i) = coord(find(nodes == output_nodes_cleaned(i)),1);
y(i) =coord(find(nodes == output_nodes_cleaned(i)),2);
end
plot(x,y,'*')

T = delaunay(x,y);

trisurf(T, x, y, -Psi_solved(node_mapping.z(output_nodes_cleaned), displaymode))
hold on
xlabel("X")
ylabel("Y")



%%
 clear Dm m Y Q FRF
for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*Minit*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(1:size(Psi_solved,2),1:size(Psi_solved,2)) = 5e-1*sqrt(lambdasolved(i));
end
%%
freq = logspace(0,5,10000);

Y = Psi_solved(node_mapping.z(frf_output),:);
Q = Psi_solved(node_mapping.z(input_nodes_mirror),:) - Psi_solved(node_mapping.z(input_nodes_backplate),:);

for ii = 1:length(freq)
    w = freq(ii);
    Hinv = ( -w^2*eye(nModes) + 1i*Dm*w + diag(lambdasolved) )\Q';
    FRF(:,:,ii) = Y*Hinv;
end
%%
figure

hold on
for i  = 1:40
semilogx(freq/2/pi,mag2db(abs(squeeze(FRF(1,i,:)))), 'LineWidth', 1)
end
% semilogx(freq/2/pi,mag2db(abs(squeeze(FRF(1,1,:)))), 'r', 'LineWidth', 3)
%%





