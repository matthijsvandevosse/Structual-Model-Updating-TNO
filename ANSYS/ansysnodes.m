function [nodes,coord] = ansysnodes(filename)

% ANSYSNODES Extract relevant node numbers from an ANSYS .txt file.
%
%   NODES = ANSYSNODES(filename)
%     filename = name of the .txt file with ANSYS nodes
%
%   [NODES,COORD] = ANSYSNODES(filename)
%     The optional second argument COORD returns the correspoding
%     coordinates of the nodes, provided that these are included in
%     filename.txt.
%
% Implemented by Gert Witvoet

fid = fopen(filename);
% read header
header = textscan(fid,'%s',1,'delimiter',';');
hpiece = textscan(header{1}{1},'%s','delimiter','\t');
hpiece = hpiece{1};

% determine header length and data format
n = length(hpiece);
form = '%u';
for i=2:n
    form = [form,'%s'];
end
node_data = textscan(fid,form);
nodes = node_data{1};

if isempty(nodes)
    warning('ansys2bode:readnodes:nonodes',...
        'No node numbers stored in %s',filename);
    return
end

% if available, determine nodal coordinates X Y Z
if nargout > 1
    try
        coord = zeros(length(nodes),3);
        idx(1) = find(~cellfun('isempty',regexp(hpiece,'X Location')));
        idx(2) = find(~cellfun('isempty',regexp(hpiece,'Y Location')));
        idx(3) = find(~cellfun('isempty',regexp(hpiece,'Z Location')));
        for i=1:3
            coord(:,i) = str2double(strrep(node_data{idx(i)},',','.'));
        end
    catch
        warning('ansys2bode:readnodes:nocoord',...
            'File does not contain nodal coordinates XYZ.');
        coord = [];
    end
end

fclose(fid);

