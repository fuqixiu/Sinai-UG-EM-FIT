% Convert UG2 to a format for EM procedure
% Created by BRK Shevlin, May 2023
clear;
%indir = '~example_data';
indir = 'C:\Users\fuq01\Documents\GitHub\Sinai-UG-EM-FIT\example_data';
infile_bed = 'em_pre_data.csv';
T_bed = readtable(fullfile(indir, infile_bed));
infile_bn = 'em_post_data.csv';
T_bn = readtable(fullfile(indir, infile_bn));
load(fullfile(indir,"beh001.mat"));

%%
% Number of subejcts and extract participant IDS for each group
possIDs = string(cell2mat(ID)); % All IDS

% Only BED matched HC
possIDs_hc_bed = cell2mat(T_bed.prolific_pid(T_bed.group == "HC",:));

% Only BED
possIDs_bed = cell2mat(T_bed.prolific_pid(T_bed.group == "BED",:));

% Only BN matched HC
possIDs_hc_bn = cell2mat(T_bn.prolific_pid(T_bn.group == "HC",:));

% Only BN
possIDs_bn = cell2mat(T_bn.prolific_pid(T_bn.group == "BN",:));

%%
% BED
tot_bed_IC=0;
for i = 1:size(possIDs_bed,1)
    this_ID = possIDs_bed(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(IC_Choice(:,id_idx))==30 || sum(IC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_bed(i,:),i);
        continue
    end
    tot_bed_IC = tot_bed_IC + 1;

    keepStruct_IC = struct();
    keepStruct_IC.offer = IC_offer(:,id_idx);
    keepStruct_IC.choice = IC_Choice(:,id_idx);
    IC_beh_bed(:,tot_bed_IC) = {keepStruct_IC};
end

tot_bed_NC=0;
for i = 1:size(possIDs_bed,1)
    this_ID = possIDs_bed(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(NC_Choice(:,id_idx))==30 || sum(NC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_bed(i,:),i);
        continue
    end
    tot_bed_NC = tot_bed_NC + 1;

    keepStruct_NC = struct();
    keepStruct_NC.offer = NC_offer(:,id_idx);
    keepStruct_NC.choice = NC_Choice(:,id_idx);
    NC_beh_bed(:,tot_bed_NC) = {keepStruct_NC};
end

% BN
tot_bn_IC=0;
for i = 1:size(possIDs_bn,1)
    this_ID = possIDs_bn(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(IC_Choice(:,id_idx))==30 || sum(IC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_bn(i,:),i);
        continue
    end
    tot_bn_IC = tot_bn_IC + 1;

    keepStruct_IC = struct();
    keepStruct_IC.offer = IC_offer(:,id_idx);
    keepStruct_IC.choice = IC_Choice(:,id_idx);
    IC_beh_bn(:,tot_bn_IC) = {keepStruct_IC};

end

tot_bn_NC=0;
for i = 1:size(possIDs_bn,1)
    this_ID = possIDs_bn(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(NC_Choice(:,id_idx))==30 || sum(NC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_bn(i,:),i);
        continue
    end
    tot_bn_NC = tot_bn_NC + 1;

    keepStruct_NC = struct();
    keepStruct_NC.offer = NC_offer(:,id_idx);
    keepStruct_NC.choice = NC_Choice(:,id_idx);
    NC_beh_bn(:,tot_bn_NC) = {keepStruct_NC};

end

% BED-matched HCS
tot_bedhc_IC=0;
for i = 1:size(possIDs_hc_bed,1)
    this_ID = possIDs_hc_bed(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(IC_Choice(:,id_idx))==30 || sum(IC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_hc_bed(i,:),i);
        continue
    end
    tot_bedhc_IC = tot_bedhc_IC + 1;

    keepStruct_IC = struct();
    keepStruct_IC.offer = IC_offer(:,id_idx);
    keepStruct_IC.choice = IC_Choice(:,id_idx);
    
    IC_beh_bedhc(:,tot_bedhc_IC) = {keepStruct_IC};

end

tot_bedhc_NC=0;
for i = 1:size(possIDs_hc_bed,1)
    this_ID = possIDs_hc_bed(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(NC_Choice(:,id_idx))==30 || sum(NC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_hc_bed(i,:),i);
        continue
    end
    tot_bedhc_NC = tot_bedhc_NC + 1;

    keepStruct_NC = struct();
    keepStruct_NC.offer = NC_offer(:,id_idx);
    keepStruct_NC.choice = NC_Choice(:,id_idx);
    NC_beh_bedhc(:,tot_bedhc_NC) = {keepStruct_NC};

end

% BN-matched HCs
tot_bnhc_IC=0;
for i = 1:size(possIDs_hc_bn,1)
    this_ID = possIDs_hc_bn(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(IC_Choice(:,id_idx))==30 || sum(IC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_hc_bn(i,:),i);
        continue
    end
    tot_bnhc_IC = tot_bnhc_IC + 1;

    keepStruct_IC = struct();
    keepStruct_IC.offer = IC_offer(:,id_idx);
    keepStruct_IC.choice = IC_Choice(:,id_idx);  
    IC_beh_bnhc(:,tot_bnhc_IC) = {keepStruct_IC};
end

tot_bnhc_NC=0;
for i = 1:size(possIDs_hc_bn,1)
    this_ID = possIDs_hc_bn(i,:);

    [~,id_idx] = ismember(this_ID,possIDs);

    if sum(NC_Choice(:,id_idx))==30 || sum(NC_Choice(:,id_idx))==0
        fprintf('\nSubject %s (%i) was a flat responder\n',possIDs_hc_bn(i,:),i);
        continue
    end
    tot_bnhc_NC = tot_bnhc_NC + 1;
    keepStruct_NC = struct();

    keepStruct_NC.offer = NC_offer(:,id_idx);
    keepStruct_NC.choice = NC_Choice(:,id_idx);
    NC_beh_bnhc(:,tot_bnhc_NC) = {keepStruct_NC};

end


%%
% Create structures

% IC condition
BED_IC = struct();
BED_IC.beh = IC_beh_bed;
BED_IC.expname = 'BED_IC';
BED_IC.ID = possIDs_bed;
BED_IC.em = {};

BEDHC_IC = struct();
BEDHC_IC.beh = IC_beh_bedhc;
BEDHC_IC.expname = 'BEDHC_IC';
BEDHC_IC.ID = possIDs_hc_bed;
BEDHC_IC.em = {};

BN_IC = struct();
BN_IC.beh = IC_beh_bn;
BN_IC.expname = 'BN_IC';
BN_IC.ID = possIDs_bn;
BN_IC.em = {};

BNHC_IC = struct();
BNHC_IC.beh = IC_beh_bnhc;
BNHC_IC.expname = 'BNHC_IC';
BNHC_IC.ID = possIDs_hc_bn;
BNHC_IC.em = {};

% NC condition
BED_NC = struct();
BED_NC.beh = NC_beh_bed;
BED_NC.expname = 'BED_NC';
BED_NC.ID = possIDs_bed;
BED_NC.em = {};

BEDHC_NC = struct();
BEDHC_NC.beh = NC_beh_bedhc;
BEDHC_NC.expname = 'BEDHC_NC';
BEDHC_NC.ID = possIDs_hc_bed;
BEDHC_NC.em = {};

BN_NC = struct();
BN_NC.beh = NC_beh_bn;
BN_NC.expname = 'BN_NC';
BN_NC.ID = possIDs_bn;
BN_NC.em = {};

BNHC_NC = struct();
BNHC_NC.beh = NC_beh_bnhc;
BNHC_NC.expname = 'BNHC_NC';
BNHC_NC.ID = possIDs_hc_bn;
BNHC_NC.em = {};

% Put into into the s structure
s = struct();
s.BN_IC = BN_IC;
s.BN_NC = BN_NC;
s.BED_IC = BED_IC;
s.BED_NC = BED_NC;
s.BEDHC_IC = BEDHC_IC;
s.BEDHC_NC = BEDHC_NC;
s.BNHC_IC = BNHC_IC;
s.BNHC_NC = BNHC_NC;

%%

save(fullfile(indir, 'DATA_ED_Round1_May23_2023.mat'),'s');


