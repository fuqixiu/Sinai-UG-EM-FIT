% Convert UG to a format for EM procedure
% Created by BRK Shevlin, April 2023
clear
addpath('simulation');

global offer0
global eta
global n

offer0 = 5;
n = 30;             % number of trials
eta = 0.8;          % fixed parameter

% Number of subejcts
nP = 1500;
IDs = 1:nP;

gen_params = rand(nP,5);
gen_params(:,2) = 15*gen_params(:,2);
gen_params(:,3) = 10;%20*gen_params(:,3);  
gen_params(:,5) = 2*gen_params(:,5);  

% Set floor for beta
for i = 1:nP
    if gen_params(i,2) < 1; gen_params(i,2) = 1; end
end

%%

% Simulate UG0
offer_cell = {};
choice_cell = {};
for i = 1:nP

    [offer, choice] = simulate_0step_ic(gen_params(i,:));
    
    offer_cell{i} = offer;
    choice_cell{i} = choice;
end
% Behavior cells
sim_beh = struct();
sim_beh.offer = offer_cell;
sim_beh.choice = choice_cell; 


simUG0 = struct();
simUG0.beh = sim_beh;
simUG0.expname = 'simUG0';
simUG0.ID = IDs;
simUG0.params = gen_params;
simUG0.em = {};


%%
% Simulate UG1
offer_cell = {};
choice_cell = {};
for i = 1:nP

    [offer, choice] = simulate_1step_ic(gen_params(i,:));

    offer_cell{i} = offer;
    choice_cell{i} = choice;
end
% Behavior cells
sim_beh = struct();
sim_beh.offer = offer_cell;
sim_beh.choice = choice_cell; 


simUG1 = struct();
simUG1.beh = sim_beh;
simUG1.expname = 'simUG1';
simUG1.ID = IDs;
simUG1.params = gen_params;
simUG1.em = {};
%%
% Simulate UG2
offer_cell = {};
choice_cell = {};
for i = 1:nP

    [offer, choice] = simulate_2step_ic(gen_params(i,:));

    offer_cell{i} = offer;
    choice_cell{i} = choice;
end
% Behavior cells
sim_beh = struct();
sim_beh.offer = offer_cell;
sim_beh.choice = choice_cell; 

simUG2 = struct();
simUG2.beh = sim_beh;
simUG2.expname = 'simUG2';
simUG2.ID = IDs;
simUG2.params = gen_params;
simUG2.em = {};
%%
% Simulate UG3
offer_cell = {};
choice_cell = {};
for i = 1:nP

    [offer, choice] = simulate_3step_ic(gen_params(i,:));

    offer_cell{i} = offer;
    choice_cell{i} = choice;
end
% Behavior cells
sim_beh = struct();
sim_beh.offer = offer_cell;
sim_beh.choice = choice_cell; 


simUG3 = struct();
simUG3.beh = sim_beh;
simUG3.expname = 'simUG3';
simUG3.ID = IDs;
simUG3.params = gen_params;
simUG3.em = {};

%%
% Put into into the s structure
s = struct();
s.simUG0 = simUG0;
s.simUG1 = simUG1;
s.simUG2 = simUG2;
s.simUG3 = simUG3;

save('DATA_Simulated_1500Subj_30Trials.mat', 's');


