%% data reformat 
clear
indir = 'C:\Users\fuq01\Documents\GitHub\Sinai-UG-EM-FIT\example_data';
infile = 'em_leap_online_baseline.csv';
T = readtable(fullfile(indir, infile));

%% Healthy group
T_hc = T(strcmp(T.group, 'Healthy'), :);
sub = unique(T_hc.participant);
% Preallocate HC_beh as a cell array
beh = cell(1, length(sub));
for i = 1:length(sub)
    keepStruct = struct();
    
    % Use strcmp for cell array comparison
    keepStruct.offer = T_hc.offer(strcmp(T_hc.participant, sub{i}), :);
    keepStruct.choice = T_hc.choice_acc(strcmp(T_hc.participant, sub{i}), :);
    
    % Assign keepStruct to the cell array
    beh{i} = keepStruct;
end
HC = struct();
HC.beh = beh;
HC.expname = 'hc';
HC.ID = sub;
HC.em = {};

%% Depression group
T_dep = T(strcmp(T.group, 'Depression'), :);
sub = unique(T_dep.participant);
% Preallocate HC_beh as a cell array
beh = cell(1, length(sub));
for i = 1:length(sub)
    keepStruct = struct();
    
    % Use strcmp for cell array comparison
    keepStruct.offer = T_dep.offer(strcmp(T_dep.participant, sub{i}), :);
    keepStruct.choice = T_dep.choice_acc(strcmp(T_dep.participant, sub{i}), :);
    
    % Assign keepStruct to the cell array
    beh{i} = keepStruct;
end
DEP = struct();
DEP.beh = beh;
DEP.expname = 'dep';
DEP.ID = sub;
DEP.em = {};

%% Anhedonia group
T_anh = T(strcmp(T.group, 'Anhedonia'), :);
sub = unique(T_anh.participant);
% Preallocate HC_beh as a cell array
beh = cell(1, length(sub));
for i = 1:length(sub)
    keepStruct = struct();
    
    % Use strcmp for cell array comparison
    keepStruct.offer = T_anh.offer(strcmp(T_anh.participant, sub{i}), :);
    keepStruct.choice = T_anh.choice_acc(strcmp(T_anh.participant, sub{i}), :);
    
    % Assign keepStruct to the cell array
    beh{i} = keepStruct;
end
ANH = struct();
ANH.beh = beh;
ANH.expname = 'anh';
ANH.ID = sub;
ANH.em = {};

%% Both group
T_both = T(strcmp(T.group, 'Both'), :);
sub = unique(T_both.participant);
% Preallocate HC_beh as a cell array
beh = cell(1, length(sub));
for i = 1:length(sub)
    keepStruct = struct();
    
    % Use strcmp for cell array comparison
    keepStruct.offer = T_both.offer(strcmp(T_both.participant, sub{i}), :);
    keepStruct.choice = T_both.choice_acc(strcmp(T_both.participant, sub{i}), :);
    
    % Assign keepStruct to the cell array
    beh{i} = keepStruct;
end
BOTH = struct();
BOTH.beh = beh;
BOTH.expname = 'both';
BOTH.ID = sub;
BOTH.em = {};

%% Final structure for 4 groups seperatly 
% Put into into the s structure
s = struct();
s.hc = HC;
s.dep = DEP;
s.anh = ANH;
s.both = BOTH;

%% All group combines
T_all = T(strcmp(T.group, 'Healthy'), :);
sub = unique(T.participant);
% Preallocate All_beh as a cell array
beh = cell(1, length(sub));
for i = 1:length(sub)
    keepStruct = struct();
    
    % Use strcmp for cell array comparison
    keepStruct.offer = T.offer(strcmp(T.participant, sub{i}), :);
    keepStruct.choice = T.choice_acc(strcmp(T.participant, sub{i}), :);
    
    % Assign keepStruct to the cell array
    beh{i} = keepStruct;
end
All = struct();
All.beh = beh;
All.expname = 'all';
All.ID = sub;
All.em = {};

%% Final structure for all groups together
% Put into into the s structure
s = struct();
s.all = All;


%% Save structure
save(fullfile(indir, 'DATA_LEAP_online_baseline_2024_4G.mat'),'s');
clear
load("example_data\DATA_LEAP_online_baseline_2024_4G.mat")
