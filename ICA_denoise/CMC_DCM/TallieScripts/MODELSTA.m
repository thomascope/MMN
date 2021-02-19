function [M,C,L,B,labs] = MODELSTA
% For use with ERP_DCM_ALEX.m
% Model architectures from Holly/Garrido work - See Phillips JNeurosc
%
% M = matrix of directional connectios (upper = forward).  % TA: changed to
%     'A' coz I prefer it.. for now
% B = B-matrix of condition specific changes. 
% C = vector of inputs (length of sources). 
% L = Lateral connections
%
% For Alex MMN, order of sensors / sources = LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
%
% AS2016

labs = {'LIFG','LSTG','LAUD','RIFG','LSTG','LAUD' };%,'RSTG'}; warning('!!!! ---- CHANGE MODEL ---- !!!!') %,'RAUD'};

% M{1} = [...
%   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
%     1    1    3    0    0    0;
%     2    1    0    3    0    0;
%     3    0    1    1    0    0;
%     0    3    2    1    0    0
%     0    0    0    0    0    0
%     0    0    0    0    0    0];

% % or perhaps.. no lat, no LIFG, input Rifg
% M{14} = ...
% [...
% 0     0     0     0     0     0;
% 0     0     1     0     0     0;
% 0     1     1     0     0     0;
% 0     0     0     0     1     0;
% 0     0     0     1     0     1;
% 0     0     0     0     1     1];
% C{14} = [0 0 1 1 0 1]';
% L{14} = sparse(6,6);
% B{14} = M{14} + L{14};

% % % % A{1}{1} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     0    0    0    0    0    0;
% % % %     0    0    1    0    0    0;
% % % %     0    0    1    0    0    0;
% % % %     0    0    0    0    1    0
% % % %     0    0    0    0    0    1
% % % %     0    0    0    0    0    1];
% % % % 
% % % % A{1}{2} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     0    0    0    0    0    0;
% % % %     0    0    0    0    0    0;
% % % %     0    1    0    0    0    0;
% % % %     0    0    0    0    0    0
% % % %     0    0    0    1    0    0
% % % %     0    0    0    0    1    0];
% % % % 
% % % % L{1} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     0    0    0    0    0    0;
% % % %     0    0    0    0    0    0;
% % % %     0    0    0    0    0    0;
% % % %     0    0    0    0    0    0
% % % %     0    0    0    0    0    0
% % % %     0    0    0    0    0    0];
% % % % 
% % % % ord  = [1 2 3 4 5 6];      % order of LFP chans
% % % % A{1}{1} = A{1}{1}(ord,ord);  % adjacency matrix
% % % % A{1}{2} = A{1}{2}(ord,ord);  % adjacency matrix
% % % % %L{1} = A{1}==3;        % laterals
% % % % C{1} = [1 1 0 0 0 0]';     % inputs
% % % % %A{1} = A{1}.*(A{1}~=3);% remove laterals from M
% % % % B{1}{1} = logical(A{1}{1}) + logical(L{1});
% % % % 
% % % % % same as 1 but with inputs on LpSTS and RFFA
% % % % A{2}{1} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     1    1    0    0    0    0;
% % % %     0    1    0    0    0    0;
% % % %     0    0    1    1    0    0;
% % % %     0    0    0    1    0    0
% % % %     0    0    0    0    0    0
% % % %     0    0    0    0    0    0];
% % % % 
% % % % A{2}{2} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     0    0    0    0    0    0;
% % % %     1    0    0    0    0    0;
% % % %     0    0    0    0    0    0;
% % % %     0    0    1    0    0    0
% % % %     0    0    0    0    0    0
% % % %     0    0    0    0    0    0];
% % % % 
% % % % L{2} = [...
% % % %   %LIFG,LSTG,LAUD,RIFG,RSTG,RAUD
% % % %     0    0    1    0    0    0;
% % % %     0    0    0    1    0    0;
% % % %     1    0    0    0    0    0;
% % % %     0    1    0    0    0    0
% % % %     0    0    0    0    0    0
% % % %     0    0    0    0    0    0];
% % % % 
% % % % ord  = [1 2 3 4 5 6];      % order of LFP chans
% % % % A{2}{1} = A{2}{1}(ord,ord);  % adjacency matrix
% % % % A{2}{2} = A{2}{2}(ord,ord);  % adjacency matrix
% % % % %L{1} = A{1}==3;        % laterals
% % % % C{2} = [1 1 1 1 0 0]';     % inputs
% % % % %A{1} = A{1}.*(A{1}~=3);% remove laterals from M
% % % % B{2}{1} = logical(A{2}{1}) + logical(L{2});



% % 
% % labs = {'lV1','rV1','lpSTS','rFFA'};
% % A{1}{1} = ...
% % [...
% % % lV1   rV1   lpSTS rFFA
% %   1     0     1     1     ;
% %   0     1     1     1     ;
% %   0     0     1     1     ;
% %   0     0     1     1     ];
% % 
% % A{1}{2} = ...
% % [...
% % % lV1   rV1   lpSTS rFFA
% %   0     0     0     0     ;
% %   0     0     0     0     ;
% %   1     1     0     0     ;
% %   1     1     0     0     ];
% % 
% % %M{1} = ones(4,4);
% % 
% % C{1} = [1 1 0 0]';
% % L{1} = sparse(4,4);
% % 
% % %L{1}(1,2) = 1;
% % %L{1}(2,1) = 1;
% % 
% % B{1} = A{1}{1} + L{1};
% % 




M{1} = ...
[...
1     1     1     1     1     1;
1     1     1     1     1     1;
1     1     1     1     1     1;
1     1     1     1     1     1;
1     1     1     1     1     1;
0     0     0     0     0     0];
C{1} = [0 0 1 0 1 0 ]'; %1
L{1} = sparse(6,6);
B{1} = M{1} + L{1};



M{2} = ...
[...
0     0     0     0     0     0;
0     0     0     0     0     0;
0     0     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1];
C{2} = [0 0 1 0 0 1]';
L{2} = sparse(6,6);
B{2} = M{2} + L{2};


M{3} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     0     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     0];
C{3} = [0 0 1 0 0 1]';
L{3} = sparse(6,6);
B{3} = M{3} + L{3};


M{4} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     1];
C{4} = [0 0 1 0 0 1]';
L{4} = sparse(6,6);
B{4} = M{4} + L{4};

% R IFG models
M{5} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     0     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     0];
C{5} = [0 0 1 0 0 1]';
L{5} = sparse(6,6);
B{5} = M{5} + L{5};


M{6} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{6} = [0 0 1 0 0 1]';
L{6} = sparse(6,6);
B{6} = M{6} + L{6};


% L IFG model
M{7} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     1];
C{7} = [0 0 1 0 0 1]';
L{7} = sparse(6,6);
B{7} = M{7} + L{7};

% bilat IFG models
M{8} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{8} = [0 0 1 0 0 1]';
L{8} = sparse(6,6);
B{8} = M{8} + L{8};


% with lateral conn STGs + rIFG
M{9} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{9} = [0 0 1 0 0 1]';
L{9} = sparse([5 2],[2 5],1,6,6);
B{9} = M{9} + L{9};


% with lateral conn STGs + LIFG
M{10} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     1];
C{10} = [0 0 1 0 0 1]';
L{10} = sparse([5 2],[2 5],1,6,6);
B{10} = M{10} + L{10};

% with lateral conn STGs + both IFGs
M{11} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{11} = [0 0 1 0 0 1]';
L{11} = sparse([5 2],[2 5],1,6,6);
B{11} = M{11} + L{11};


% with lateral conn IFGs (not lat STG)
M{12} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{12} = [0 0 1 0 0 1]';
L{12} = sparse([4 1],[1 4],1,6,6);
B{12} = M{12} + L{12};


% with lateral conn IFGs + lat STGs
M{13} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{13} = [0 0 1 0 0 1]';
L{13} = sparse([4 1 2 5],[1 4 5 2],1,6,6);
B{13} = M{13} + L{13};


% no lat, no LIFG, input Rifg
M{14} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{14} = [0 0 1 1 0 1]';
L{14} = sparse(6,6);
B{14} = M{14} + L{14};


% no lat, no rIFG, input Lifg
M{15} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     1];
C{15} = [1 0 1 0 0 1]';
L{15} = sparse(6,6);
B{15} = M{15} + L{15};


% no lat, , input both IFGs
M{16} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{16} = [1 0 1 1 0 1]';
L{16} = sparse(6,6);
B{16} = M{16} + L{16};


% Lat STG, RIFG + input RIFG
M{17} = ...
[...
0     0     0     0     0     0;
0     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{17} = [0 0 1 1 0 1]';
L{17} = sparse([5 2],[2 5],1,6,6);
B{17} = M{17} + L{17};


% Lat STG, LIFG + input LIFG
M{18} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     0     0;
0     0     0     0     0     1;
0     0     0     0     1     1];
C{18} = [1 0 1 0 0 1]';
L{18} = sparse([5 2],[2 5],1,6,6);
B{18} = M{18} + L{18};


% Lat STG, bilat inputs to IFG
M{19} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{19} = [1 0 1 1 0 1]';
L{19} = sparse([5 2],[2 5],1,6,6);
B{19} = M{19} + L{19};


% Lat IFG, bilat inputs to IFG
M{20} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{20} = [1 0 1 1 0 1]';
L{20} = sparse([4 1],[1 4],1,6,6);
B{20} = M{20} + L{20};


% Lat IFG, lat STG + inputs to both IFG
M{21} = ...
[...
0     1     0     0     0     0;
1     0     1     0     0     0;
0     1     1     0     0     0;
0     0     0     0     1     0;
0     0     0     1     0     1;
0     0     0     0     1     1];
C{21} = [1 0 1 1 0 1]';
L{21} = sparse([4 1 5 2],[1 4 2 5],1,6,6);
B{21} = M{21} + L{21};


