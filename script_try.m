IC_Choice= IC_Choice(:, 1:5);
IC_offer = IC_offer(:, 1:5);
IC_reward = IC_reward(:, 1:5);
IC_pc = IC_pc(1:5,:);
IC_pc_rt = IC_pc_rt(1:5,:);
IC_RT = IC_RT(:,1:5);
ID = ID(1:5,:);

%% data reformat 
clear
%indir = 'C:\Users\fuq01\Documents\GitHub\Sinai-UG-EM-FIT\example_data';

indir = 'C:\Users\fuq01\Documents\Scripts\leap_ddm\Data\';
infile_all = 'ddm_or_data.csv';
T_bed = readtable(fullfile(indir, infile_all));

%possIDs_hc_bed = cell2mat(T_bed.participant(T_bed.stim == "pre",:));

%infile_pre = 'em_pre_data.csv';
T%_pre = readtable(fullfile(indir, infile_pre));
%infile_post = 'em_post_data.csv';
%T_post = readtable(fullfile(indir, infile_post));

preIDs = T_pre.participant;

sub = unique(T_bed.participant);
for s = unique(T_bed.participant)
    disp(s);
    offer = T_bed.offer(T_bed.participant == s,:);
    %disp(offer)
end

T_pre.offer(T_pre.participant == 701,:)


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

% Put into into the s structure
s = struct();
s.BED_NC = BED_NC;
s.BEDHC_NC = BEDHC_NC;


%%
save(fullfile(indir, 'DATA_ED_Round1_May23_2023.mat'),'s');S
