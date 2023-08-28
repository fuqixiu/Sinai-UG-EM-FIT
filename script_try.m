%% data reformat 
clear
indir = 'C:\Users\fuq01\Documents\GitHub\Sinai-UG-EM-FIT\example_data';
infile_pre = 'em_pre_data.csv';
T_pre = readtable(fullfile(indir, infile_pre));
infile_post = 'em_post_data.csv';
T_post = readtable(fullfile(indir, infile_post));

sub = unique(T_pre.participant);
for i = 1:length(sub)
    keepStruct_pre = struct();
    keepStruct_pre.offer = T_pre.offer(T_pre.participant == sub(i),:);
    keepStruct_pre.choice = T_pre.choice_b(T_pre.participant == sub(i), :);
    pre_beh(:,i) = {keepStruct_pre};
end

for i = 1:length(sub)
    keepStruct_post = struct();
    keepStruct_post.offer = T_post.offer(T_post.participant == sub(i),:);
    keepStruct_post.choice = T_post.choice_b(T_post.participant == sub(i), :);
    post_beh(:,i) = {keepStruct_post};
end

% pre_stim session 
pre_stim = struct();
pre_stim.beh = pre_beh;
pre_stim.expname = 'pre_stim';
pre_stim.ID = num2str(sub);
pre_stim.em = {};


% post_stim session 
post_stim = struct();
post_stim.beh = post_beh;
post_stim.expname = 'post_stim';
post_stim.ID = num2str(sub);
post_stim.em = {};

% Put into into the s structure
s = struct();
s.pre = pre_stim;
s.post = post_stim;


%%
save(fullfile(indir, 'DATA_LEAP_Aug28_2023.mat'),'s');
