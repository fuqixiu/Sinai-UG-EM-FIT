function [ out ] = EMmc_ms(rootfile,modnames,omitBMS)
% Model comparison for EM fitted models
% MK Wittmann, Nov 2017
%
% INPUT:        - rootfile: root of a specific experiment
%               - modnames: which model to looks at
%               - omitBMS: you can omit plotting the bayesian model comparison, because it requres you to install spm
%
%%


expname     = rootfile.expname;
ID          = rootfile.ID; nm = numel(ID);
moncolset   = [    1.0000    1.0000    0.6980
    0.9961    0.8000    0.3608
    0.9922    0.5529    0.2353
    0.8902    0.1020    0.1098 ];

%%%% get info
lme= zeros(1,)
for imod=1:numel(modnames)
   lme(:,imod)    = rootfile.em.(modnames{imod}).fit.lme;
   bicint(imod)   = rootfile.em.(modnames{imod}).fit.bicint;
end


% I. Plot LME and BIC: --------------------------------------------------
% 1) LME sum
h=figure('name', expname);
suptitle('Bayesian Model Comparison');

subplot(2,2,1);   bar(sum(lme)); xtickrot = 25;
set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
ylabel('Summed log evidence (more is better)','FontWeight','bold');

% 2) BICint
subplot(2,2,2);
bar(bicint);
set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames);
set(gca,'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
ylabel('BICint (less is better)','FontWeight','bold');
%------------------------------------------------------------------------


% II Calculate exceedance probability and compare log model evidence ----
% Compare LME of two models:
compM   = {'ms_RL_cl_cs','ms_RL_cl_cs_rt'};                                    % models to compare specifically; only executed if these models have been fitted
if isfield(rootfile.em,compM{1})&&isfield(rootfile.em,compM{2})
   subplot(2,2,3);
   for imod=1:numel(compM)
      lmepair(:,imod)    = rootfile.em.(compM{imod}).fit.lme;
   end
   lmediff = lmepair(:,2) - lmepair(:,1);
   [lmediffsort,scode] = sort(lmediff);
   IDsort = ID(scode);
   MKIDS = sort(unique(ID));
   for imon = 1:numel(MKIDS)
      curID = MKIDS(imon);
      plotlmes = zeros(nm,1); plotlmes(IDsort==curID) = lmediffsort(IDsort==curID);
      hbar = barh(1:nm, plotlmes); hold all;
      set(hbar,'Facecolor',moncolset(curID,:));
   end
   set(gca,'ytick',1:numel(IDsort),'yticklabel',IDsort);
   ylabel('Sessions (sorted)');
   xlabel([compM{1} ' < ' compM{2} ]);
   title('Log Evidence Differences / Log Bayes factor');
end

% this is coming from SPM you need to have it somewhere in your paths
if omitBMS == 0
   [~,~,BMS.xp] = spm_BMS(lme);                                                   % log model evidence
   subplot(2,2,4)
   bar(BMS.xp);
   set(gca,'XTick',1:numel(modnames),'XTickLabel',modnames,'XTickLabelRotation',xtickrot);
   rl1 = refline(0,.95);     set(rl1,'linestyle','--','Color','r');
   ylabel('Exceedance Probability','FontWeight','bold');
end
setfp(gcf)
figname=['EM_BMC_' expname '_'  date];
saveas(gcf,figname);



