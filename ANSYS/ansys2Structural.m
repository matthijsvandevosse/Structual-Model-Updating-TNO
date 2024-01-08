function [Stiff, Mass, Node_mapping] = ansys2Structural(filename_full, retreive_mass)
%ANSYS2KMATRIX rewrite the K matrix from ansys based on the solver ordning
%to the specified user nodes
%
%   Does not accept distributed files (so no file0.full)
%
%   Inputs:
%     filename_full  = the "file.full" file as used by ANSYS in the
%       mechanical environment. It can be found in the same folder as the
%       "file.mode" file
%     User_nodes_in  = The user nodes input reference
%     User_nodes_out = The user nodes output reference
%     d              = The degree of freedom.
%                           1 = x;
%                           2 = y;
%                           3 = z;
%   Outputs:
%     Stiff     = stiffness matrix
%     idx_in    = DOF reference at input stiffness matrix
%     idx_out   = DOF reference at output stiffness matrix
%
%   Notes:
%   * filename should be a file with extension *.full
%     don't forget to add the extension
%
%   Check the ANSYS Mechanical APDL Programmer's Reference
%
% Implemented by Stijn Langedijk,
% Edited by Matthijs van de Vosse (04-12-2023)


%%%%% Initialization %%%%%
disp(' ');     % to new line in output window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Open ANSYS Binary file %%%%%
tic;
try
    [fid,message]=fopen(filename_full,'r',type);    % options: real only, binary file ('b' may not be omitted)
    % N Rijnveld, Oct 2011: for ANSYS 13, file works with omitting of 'b'
catch
    [fid,message]=fopen(filename_full,'r');
end

if fid==-1
    disp(message);
    return
else
    str=['ANSYS Binary file ',filename_full,' opened succesfully'];
    disp(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Read Sets from file %%%%%

% Set appropriate format
format1='int32';         % integer of 32 bits, or 4 bytes
fac1=4;                % 4 bytes per integer; fac is used in fseek command
format2='double';      % float of 64 bits, or 8 bytes
fac2=8;                % 8 bytes per integer; fac is used in fseek command


% Set 1: Standard ANSYS file header
reclength=fread(fid,1,format1);                 % length of ANSYS header (100)
dum=fread(fid,1,format1);                       % signed 32-bit integer: -2^31
ANSYS_file_header=fread(fid,reclength,format1); % read ANSYS file header
dum=fread(fid,1,format1);                       % length of ANSYS header (100)


% Set 2: K file header
reclength=fread(fid,1,format1) ;                % length of Mode header (60 or 100)
dum=fread(fid,1,format1) ;                      % signed 32-bit integer: -2^31
Stiff_file_header=fread(fid,reclength,format1);  % read Mode file header
dum=fread(fid,1,format1);                       % length of Mode header (60 or 100)

fun04=Stiff_file_header(1);
neqn=Stiff_file_header(2);
nmrow=Stiff_file_header(3);          % number of (master) degrees of freedom (=lenbac*numdof)
nmatrx=Stiff_file_header(4);          % number of modes extracted: #eigenvalues
numdof=Stiff_file_header(8);         % number of degrees of freedom per node
lenbac=Stiff_file_header(7);         % number of nodes

ptfSTF = Stiff_file_header(19);      %pointer to stiffness matrix
keyuns = Stiff_file_header(14);      %pointer to stiffness matrix
ptrMASl = Stiff_file_header(27);      %pointer to mass matrix
ptrDOFl = Stiff_file_header(36);
ptrRHSl = Stiff_file_header(38);
if keyuns == 1
    disp('there exists an unsymetric matrix in the file, probz better stop')
end

%%
% c fun04, neqn, nmrow, nmatrx, kan,
% c wfmax, lenbac, numdof,jcgtrmL,jcgtrmH,
% c lumpm, jcgeqn, jcgtrm, keyuns, extopt,
% c keyse, sclstf, nxrows,ptrIDXl,ptrIDXh,
% c ncefull,ncetrm,ptrENDl,ptrENDh, 0,

% Display relevant ANSYS data
%ansysinfo(ANSYS_file_header,Mode_file_header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set 3: Analysis information

% Group 1: Degrees of freedom per node
dum=fread(fid,1,format1);           % number of DOF per node = numdof
dum=fread(fid,1,format1);           % signed 32-bit integer: -2^31
dof=fread(fid,numdof,format1);      % read DOF per node
dum=fread(fid,1,format1);           % number of DOF per node = numdof

% Group 2: Nodal equivalence table:
% possibly not agreeing with the stiffness matrix
dum=fread(fid,1,format1);           % length of nodal eq.table = lenbac (number of nodes)
dum=fread(fid,1,format1);           % signed 32-bit integer: -2^31
nodeq=fread(fid,lenbac,format1);    % read Nodal equivalence table
dum=fread(fid,1,format1);           % length of nodal eq.table = lenbac (number of nodes)

% Group: DOF information

fseek(fid,ptrDOFl*fac1,'bof')       %
l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);           % tells you about the data type used
Nod_ext=fread(fid,l1,format1);      % DOF at node
dum =fread(fid,1,format1);          % length of thing

l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);           % tells you about the data type used
DOF_vec=fread(fid,l1,format1);      % content
dum =fread(fid,1,format1);

l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);           % tells you about the data type used
DOF_im=fread(fid,l1,format1);       % DOF that have imposed values
dum =fread(fid,1,format1);

l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);           % tells you about the data type used
DOF_im_val=fread(fid,l1/2,format2); % The imposed values
dum =fread(fid,1,format1);

% - The nodal equivalence table nodeq relates the user nodal numeration to
% the solver numeration
% Using the Nodal extent vector, the degrees of freedom per node can be
% found, relating to the total degrees of freedom as used in the stiffness
% matrix. This is not necesarily the same as used in Modal analysis.
%
% In order to truncate the K matrix to the values as requested by the in-
% and output, the DOF is taken relating these.
% [nodeq_s, idx] = sort(nodeq);
% Nod_ext_s = Nod_ext(idx);
clear idx
% step one: find position in nodeq:

X_DOF = find(abs(DOF_vec)==1);
Y_DOF = find(abs(DOF_vec)==2);
Z_DOF = find(abs(DOF_vec)==3);
[~, index_sort] = sort(nodeq, 'ascend');

Node_mapping.x = X_DOF(index_sort);
Node_mapping.y = Y_DOF(index_sort);
Node_mapping.z = Z_DOF(index_sort);
Node_mapping.ReadMea = "NODE: -> ROW/COLUMN. Maps usernode to the row and column of M and K";



% Group 3: Stifness matrix
str=['  ... Reading Stiffness matrix from binary file ',filename_full];
disp(str);

fseek(fid,ptfSTF*fac1,'bof');    % goto stiffness record
perc = -1;
t1 = clock;
idx = cell(1,neqn);
cont = cell(1,neqn);
for ii = 1:neqn

    perc1 = round(ii/neqn*100);
    %give update of proces
    if perc ~= perc1
        t2 = clock;
        elapsed_time = (t2 - t1)*[0,0,0,60,1,1/60]';
        time_left = neqn/ii*elapsed_time - elapsed_time;
        perc = perc1;
        disp(['Retrieving of stiffness matrix: ',num2str(perc1),'%, Time left: ',num2str(time_left), ' Minutes'])
    end
    %     fseek(fid,ptfSTF*fac1,'bof')    % goto stiffness record
    %     fseek(fid,ii*fac1,'cof')        % goto stiffness record
    % indexing
    l1=fread(fid,1,format1);        % length of content
    dum=fread(fid,1,format1);       % info about the data type used
    idx{ii}=[fread(fid,l1,format1)'; repelem(ii,1,l1)];     % content
    dum =fread(fid,1,format1);      % length of content


    % load content
    l1=fread(fid,1,format1);        % length of content
    dum=fread(fid,1,format1);       % tells you about the data type used
    cont{ii}=fread(fid,l1/2,format2);   % content
    dum =fread(fid,1,format1);      % length of content

    % if (ismember(ii, DOF_rot))
    %     cont{ii}= 0*cont{ii};
    %     ii
    % end
    % populate the symmetric matrix
    %Stiff(ii,idx1) = idx2;
end
disp('Writing stiffness sparse matrix')
N = [0 cumsum(cellfun('size', idx,2))]; % same for column
idx1 = zeros(2,N(end));
content = zeros(1,N(end));
for ii = 1:length(N)-1
    % try
        idx1(:,N(ii)+1:N(ii+1)) = idx{ii};
        content(N(ii)+1:N(ii+1)) = cont{ii};
    % catch
    %     disp('Error in stiffness matrix')
    % end
end
Stiff = sparse([idx1(1,:)],[idx1(2,:)],[content],neqn,neqn);
Stiff = Stiff+Stiff'-diag(diag(Stiff));

fseek(fid,ptrRHSl*fac1,'bof');      %
%Load vector
l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);           % tells you about the data type used
dum=fread(fid,l1/2,format2);        % content
dum =fread(fid,1,format1);          % length of thing

%Stiffness matrix diagonal vector
l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);            % tells you about the data type used
stiff_diag=fread(fid,l1/2,format2);             % content
dum =fread(fid,1,format1);            % length of thing

% norm(full(diag(Stiff))-stiff_diag) % check if done correctly

%Stiffness matrix diagonal scaling vector
l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);            % tells you about the data type used
stiff_scale=fread(fid,l1/2,format2);             % content
dum =fread(fid,1,format1);            % length of thing

% DOF marker array
l1=fread(fid,1,format1);            % length of thing
dum=fread(fid,1,format1);            % tells you about the data type used
DOF_m_ar=fread(fid,l1/2,format2);             % content
dum =fread(fid,1,format1);            % length of thing

% Constrain nodes fixed to ground (by default no stiffness assigned in ansys)
idx = find(DOF_vec < 0);
for i = 1:length(idx)
    Stiff(idx(i),idx(i))  =1e20;
    % Mass(idx,idx) = 1e20;
end

%% Extracting the mass matrix
% w = waitbar(0,'Extracting the mass matrix');     
Mass = sparse(neqn,neqn);
if retreive_mass
    fseek(fid,ptrMASl*fac1,'bof');       % goto mass record

    display("Extracting the mass matrix")

    for ii = 1:neqn
        % try
            % first indexing
            l1=fread(fid,1,format1);        % length of content
            dum=fread(fid,1,format1);       % info about the data type used
            idx1=fread(fid,l1,format1);     % content
            dum =fread(fid,1,format1);      % length of content

            % load content
            l1=fread(fid,1,format1);        % length of content
            dum=fread(fid,1,format1);       % tells you about the data type used
            idx2=fread(fid,l1/2,format2);   % content
            dum =fread(fid,1,format1);      % length of content

            % populate the symmetric matrix
            Mass(ii,idx1) = rmmissing(idx2);
            Mass(idx1,ii) = rmmissing(idx2);

        % w = waitbar(ii/neqn,w,"Extracting the mass matrix");
        % catch
        %     display("error in Mass matrix")
        % end
    end
end

%Stiff = full(Stiff(idx_in,idx_out));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Close ANSYS Binary file %%%%%
status=fclose(fid);
if status==0
    str=['ANSYS Binary file ',filename_full,' closed succesfully'];
    disp(str);
else
    str=['Could not close ANSYS Binary file ',filename_full];
    disp(str);
end

disp(['   (Data loaded in ',num2str(toc),' seconds)'])



