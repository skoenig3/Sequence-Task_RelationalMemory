% written by Seth Konig December 2014
% code runs all the other code preprocess then process behavioral only data 
% when task was run on cortex and not during a recording

datadir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Cortex Data\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Figures\';

% PW pre-lesion
% %may or may not want to drop sessions for SeqL13.itm, PW150410.2 got
% %through 93.6% of all trials with 81.1% of all trials being valid
cch25_files = {'PW150324.1','PW150325.1','PW150326.1','PW150327.1','PW150330.1',...
               'PW150331.1','PW150401.1','PW150402.1','PW150403.1','PW150406.1',...
               'PW150407.1','PW150409.1','PW150413.1','PW150414.1',...
               'PW150415.1'};%'PW150410.1'
SeqL_files = {'PW150324.2','PW150325.2','PW150326.2','PW150327.2','PW150330.2',...
              'PW150331.2','PW150401.2','PW150402.2','PW150403.2','PW150406.2',...
              'PW150407.2','PW150409.2','PW150413.2','PW150414.2',...
              'PW150415.2'}; %'PW150410.2'
item_sets =  {'SeqL01.itm','SeqL02.itm','SeqL03.itm','SeqL04.itm','SeqL05.itm',...
              'SeqL06.itm','SeqL07.itm','SeqL08.itm','SeqL09.itm','SeqL10.itm',...
              'SeqL11.itm','SeqL12.itm','SeqL14.itm','SeqL15.itm',...
              'SeqL16.itm'}; %'SeqL13.itm'

% PW post-lesion
% cch25_files = {'PW160125.1','PW160126.1','PW160127.1','PW160128.1','PW160129.1',...
%                'PW160201.1','PW160202.1','PW160208.1','PW160209.1','PW160210.1',...
%                'PW160212.1','PW160216.1','PW160217.1','PW160218.1','PW160219.1'};
% SeqL_files =  {'PW160125.2','PW160126.2','PW160127.2','PW160128.2','PW160129.2',...
%                'PW160201.2','PW160202.2','PW160208.2','PW160209.2','PW160210.2',...
%                'PW160212.2','PW160216.2','PW160217.2','PW160218.2','PW160219.2'};
% item_sets =   {'SeqL21.itm','SeqL22.itm','SeqL23.itm','SeqL24.itm','SeqL25.itm',...
%                'SeqL26.itm','SeqL27.itm','SeqL30.itm','SeqL31.itm','SeqL32.itm',...
%              'SeqL34.itm','Seql35.itm','SeqL36.itm','Seql37.itm','SeqL38.itm'};
min_rt = 156;%mininum reaction time varies by monkey. Taken as 5 percentile for random sequences
% 
% pre_sets = {'PW150324.2','PW150325.2','PW150326.2','PW150327.2','PW150330.2',...
%               'PW150331.2','PW150401.2','PW150402.2','PW150403.2','PW150406.2',...
%               'PW150407.2','PW150409.2','PW150413.2','PW150414.2','PW150415.2'};
% post_sets = {'PW160125.2','PW160126.2','PW160127.2','PW160128.2','PW160129.2',...
%                'PW160201.2','PW160202.2','PW160208.2','PW160209.2','PW160210.2',...
%                'PW160212.2','PW160216.2','PW160217.2','PW160218.2','PW160219.2'};

%TT
%  cch25_files = {'TT150402.1','TT150403.1','TT150406.1','TT150407.1','TT150408.1',...
%                 'TT150409.1','TT150413.1','TT150414.1','TT150415.1','TT150416.1',...
%                 'TT150417.1','TT150420.1','TT150421.1','TT150422.1','TT150423.1',};%TT150424.1
% SeqL_files = {'TT150402.2','TT150403.2','TT150406.2','TT150407.2','TT150408.2',...
%               'TT150409.2','TT150413.2','TT150414.2','TT150415.2','TT150416.2',...
%               'TT150417.2','TT150420.2','TT150421.2','TT150422.2','TT150423.2',};%TT150424.2
% item_sets =  {'SeqL01.itm','SeqL02.itm','SeqL03.itm','SeqL04.itm','SeqL05.itm',...
%               'SeqL06.itm','SeqL07.itm','SeqL08.itm','SeqL09.itm','SeqL10.itm',...
%               'SeqL11.itm','SeqL12.itm','SeqL13.itm','SeqL14.itm','SeqL15.itm'};
%min_rt = 172;
%SeqL11.itm run with SeqL01.cnd NOT SeqL11.cnd

%RR pre-lesion
% cch25_files = {'RR150324.1','RR150325.1','RR150326.1','RR150327.1','RR150330.1',...
%                'RR150331.1','RR150401.1','RR150408.1','RR150409.1','RR150410.1',...
%                'RR150413.1','RR150414.1','RR150415.1','RR150417.1','RR150422.1'};
% SeqL_files = {'RR150324.2','RR150325.2','RR150326.2','RR150327.2','RR150330.2',...
%               'RR150331.2','RR150401.2','RR150408.2','RR150409.2','RR150410.2',...
%               'RR150413.2','RR150414.2','RR150415.2','RR150417.2','RR150422.2'};
% item_sets =  {'SeqL21.itm','SeqL22.itm','SeqL23.itm','SeqL24.itm','SeqL25.itm',...
%               'SeqL26.itm','SeqL27.itm','SeqL31.itm','SeqL32.itm','SeqL33.itm',...
%               'SeqL34.itm','SeqL35.itm','SeqL36.itm','SeqL38.itm','SeqL39.itm'};

%RR post-lesion
% cch25_files = {'RR160210.2','RR160211.2','RR160212.2','RR160217.2','RR160218.2',...
%               'RR160219.2','RR160222.2','RR160223.2','RR160224.2','RR160225.2',...
%               'RR160226.2','RR160229.2','RR160301.2','RR160302.2','RR160303.2'};
% SeqL_files = {'RR160210.1','RR160211.1','RR160212.1','RR160217.1','RR160218.1',...
%               'RR160219.1','RR160222.1','RR160223.1','RR160224.1','RR160225.1',...
%               'RR160226.1','RR160229.1','RR160301.1','RR160302.1','RR160303.1'};
% item_sets =  { 'SeqL06.itm','SeqL08.itm','SeqL07.itm','SeqL09.itm','SeqL10.itm',...
%               'SeqL11.itm','SeqL12.itm','SeqL13.itm','SeqL14.itm','SeqL15.itm',...
%               'SeqL16.itm','SeqL17.itm','SeqL18.itm','SeqL19.itm','SeqL20.itm'};
%            also have SeqL40.itm RR160304.1. Collected Just In case
% min_rt = 109; 

% pre_sets = {'RR150324.2','RR150325.2','RR150326.2','RR150327.2','RR150330.2',...
%               'RR150331.2','RR150401.2','RR150408.2','RR150409.2','RR150410.2',...
%               'RR150413.2','RR150414.2','RR150415.2','RR150417.2','RR150422.2'};
% post_sets = {'RR160210.1','RR160211.1','RR160212.1','RR160217.1','RR160218.1',...
%               'RR160219.1','RR160222.1','RR160223.1','RR160224.1','RR160225.1',...
%               'RR160226.1','RR160229.1','RR160301.1','RR160302.1','RR160303.1'};

%TO pre-lesion
% cch25_files = {'TO150421.1','TO150422.2','TO150423.1','TO150424.1','TO150427.2',...
%                'TO150428.1','TO150429.1','TO150430.1','TO150501.1','TO150504.1',...
%                'TO150505.1','TO150506.1','TO150507.1','TO150508.1','TO150512.1'};
% SeqL_files = {'TO150421.3','TO150422.3','TO150423.2','TO150424.2','TO150427.3',...
%               'TO150428.2','TO150429.2','TO150430.2','TO150501.2','TO150504.2',...
%               'TO150505.2','TO150506.2','TO150507.2','TO150508.2','TO150512.2'};
% item_sets = {'SeqL22.itm','SeqL23.itm','SeqL24.itm','SeqL25.itm','SeqL26.itm',...
%              'SeqL27.itm','SeqL28.itm','SeqL29.itm','SeqL30.itm','SeqL31.itm',...
%              'SeqL32.itm','SeqL33.itm','SeqL34.itm','SeqL35.itm','SeqL37.itm'}; 

%TO post-lesion
%  cch25_files = {'TO170228.1','TO170301.1','TO170306.1','TO170307.1','TO170308.1',...
%                 'TO170309.1','TO170310.1','TO170313.1','TO170314.1','TO170315.1',...
%                 'TO170316.1','TO170317.1','TO170320.1','TO170321.1','TO170322.1'};
%  
% SeqL_files = {'TO170228.2','TO170301.2','TO170306.2','TO170307.2','TO170308.2',...
%               'TO170309.2','TO170310.2','TO170313.2','TO170314.2','TO170315.2',...
%               'TO170316.2','TO170317.2','TO170320.2','TO170321.2','TO170322.2'};
% 
% 
% item_sets = {'SeqL01.itm','SeqL02.itm','SeqL03.itm','SeqL04.itm','SeqL05.itm',...
%              'SeqL06.itm','SeqL07.itm','SeqL08.itm','SeqL09.itm','SeqL10.itm',...
%              'SeqL11.itm','SeqL12.itm','SeqL13.itm','SeqL14.itm','SeqL15.itm'};
% min_rt = 135;
% 

%MF pre-lesion
% cch25_files ={'MF161209.1','MF161212.1','MF161213.1','MF161214.1',...
%               'MF161215.1','MF161219.1','MF161220.1','MF161221.1',...
%               'MF161229.1','MF161230.1','MF170103.1','MF170104.2',...
%               'MF170105.1','MF170106.1','MF170109.1','MF170110.1'};
% SeqL_files = {'MF161209.2','MF161212.2','MF161213.2','MF161214.2',...
%               'MF161215.2','MF161219.2','MF161220.2','MF161221.2',...
%               'MF161229.2','MF161230.2','MF170103.2','MF170104.3',...
%               'MF170105.2','MF170106.2','MF170109.2','MF170110.2'};
% item_sets = {'SeqL02.itm','SeqL03.itm','SeqL04.itm','SeqL05.itm',...
%              'SeqL06.itm','SeqL08.itm','SeqL09.itm','SeqL10.itm',...
%              'SeqL11.itm','SeqL12.itm','SeqL13.itm','SeqL14.itm',...
%              'SeqL15.itm','SeqL16.itm','SeqL17.itm','SeqL18.itm'}; 
%should skip set SeqL01.itm since had calibration issue half way through
%and looks like creates weird RTs, and skip SeqL07.itm since restarted
%after a few trials
% min_rt = 112;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Preprocess all the SeqL data---%%%
% for ls =length(SeqL_files) %7,9,12-15
%     ImportSeqLData(cch25_files{ls},SeqL_files{ls},item_sets{ls})
%     close all
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---Autmoatically analyze Reaction times in Sequence Task---%%%
for ls = 1:length(SeqL_files)
    DetermineSeqLTime2Fixation([datadir SeqL_files{ls}(1:8) '_'  SeqL_files{ls}(end) '-fixation.mat'])
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---Autmoatically Plot Eye Position Density---%%%
% for ls = length(SeqL_files)
%     PlotEyeDensity([datadir SeqL_files{ls}(1:8) '_'  SeqL_files{ls}(end) '-fixation.mat'],figure_dir)
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Autmoatically analyze Fixation Durations in List task---%%%
% for ls = length(SeqL_files)
%     DetermineSeqLFixationDurations([datadir SeqL_files{ls}(1:8) '_'  SeqL_files{ls}(end) '-fixation.mat'])
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%---Automatically analyze KL Divergencitem_sets{ls}e Data for List task---%%%
% for ls = 1:length(SeqL_files)
%    DetermineSeqLKLdiverence([datadir SeqL_files{ls}(1:8) '_'  SeqL_files{ls}(end) '-fixation.mat']);
% end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%---Plot mean reaction time across sessions---%%%
block_prop = NaN(length(SeqL_files),21);
sess_means = NaN(1,length(SeqL_files)); 
all_rts = [];
n_trials = [];
bad_trials = []; 
mean_rt_diff = [];
mean_rt = [];
rt_diff = [];
cortex_break = [];
means = NaN(1,4);
num_trials = NaN(1,length(SeqL_files));
for ls =1:length(SeqL_files)
    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-SQRTs.mat'],'reaction_time','cortex_rt','cortexbreak');
    
    num_trials(ls) = sum(~isnan(reaction_time(:,1)));
    [~,~,sequence_locations] = read_SeqL_itm_and_cnd_files(item_sets{ls},SeqL_files{ls}(1:2)); %load item and condition files
    overlap = find(sum(sequence_locations{1} == sequence_locations{2}) == 2);

    mean_rt = [mean_rt; nanmean(reaction_time)];
    rt_diff = [rt_diff; reaction_time-cortex_rt];
    mean_rt_diff(ls) = nanmean(nanmean(reaction_time-cortex_rt));
    [~,~,sequence_locations] = read_SeqL_itm_and_cnd_files(item_sets{ls},SeqL_files{ls}(1:2)); %load item and condition files
    n_trials(ls) = sum(~isnan(reaction_time(:,1)));
    means(ls,:) =100*sum(reaction_time(21:end,:) < min_rt)./sum(~isnan(reaction_time(21:end,:)));
    cortex_break = [cortex_break; cortexbreak];
    if overlap ~= 1
        all_rts = [all_rts; reaction_time];
    else
        temp = reaction_time;
        temp(:,1) = NaN;
        all_rts = [all_rts; temp];
    end
    sess_means(ls) = mean(100*sum(reaction_time(21:end,:) < min_rt)./sum(~isnan(reaction_time(21:end,:))));
    nb = floor(size(reaction_time,1)/20);
    for n = 1:nb
        temp = reaction_time(20*(n-1)+1:20*n,2:4);
        temp = temp(1:end);
        block_prop(ls,n) = 100*sum(temp < min_rt)./sum(~isnan(temp));
    end
end

n_sess = sampsizepwr('z',[mean(sess_means) std(sess_means)],5);
n_item2 = sampsizepwr('z',[mean(means(:,2)) std(means(:,2))],5);
n_item3 = sampsizepwr('z',[mean(means(:,3)) std(means(:,3))],5);
n_item4 = sampsizepwr('z',[mean(means(:,4)) std(means(:,4))],5);

[~,p_item(1)] = ztest(means(:,1),5,std(means(:,1)),'tail','right');
[~,p_item(2)] = ztest(means(:,2),5,std(means(:,2)),'tail','right');
[~,p_item(3)] = ztest(means(:,3),5,std(means(:,3)),'tail','right');
[~,p_item(4)] = ztest(means(:,4),5,std(means(:,4)),'tail','right');

[~,p_item2(1)] = ttest2(means(:,1),means(:,2));
[~,p_item2(2)] = ttest2(means(:,1),means(:,3));
[~,p_item2(3)] = ttest2(means(:,1),means(:,4));

figure
hold on
bar(nanmean(means));
errorb(nanmean(means),nanstd(means)./sqrt(sum(~isnan(means))))
set(gca,'Xtick',[1:4])
xlabel('Item Number')
% ylabel('Mean Reacion times (time to fixation) (ms)')
ylabel(['Percentage of Reaction times < ' num2str(min_rt) ' ms'])
title(['Relation Memory Sequence Task Reaction times for ' ...
    SeqL_files{1}(1:2) ' n = ' num2str(size(means,1))])

for ii = 1:4
    if p_item(ii) < 0.05
        plot(ii,mean(means(:,ii))+1.5*std(means(:,ii))./sqrt(15),'*k')
    end
end
for ii = 1:3
    if p_item2(ii) < 0.05
        plot(ii+1,mean(means(:,ii+1))+2.5*std(means(:,ii+1))./sqrt(15),'*r')
    end
end
hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot mean reaction time across sessions < 50---%%%
block_prop = NaN(length(SeqL_files),21);
sess_means = NaN(1,length(SeqL_files)); 
n_trials = [];
bad_trials = []; 
means = NaN(1,4);
for ls =1:length(SeqL_files)
    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-SQRTs.mat'],'reaction_time');
    n_trials(ls) = sum(~isnan(reaction_time(:,1)));
    means(ls,:) =100*sum(reaction_time(21:end,:) < 50)./sum(~isnan(reaction_time(21:end,:)));
    sess_means(ls) = mean(100*sum(reaction_time(21:end,:) < 50)./sum(~isnan(reaction_time(21:end,:))));
    nb = floor(size(reaction_time,1)/20);
    for n = 1:nb
        temp = reaction_time(20*(n-1)+1:20*n,2:4);
        temp = temp(1:end);
        block_prop(ls,n) = 100*sum(temp < 50)./sum(~isnan(temp));
    end
end

n_sess = sampsizepwr('z',[mean(sess_means) std(sess_means)],5);
n_item2 = sampsizepwr('z',[mean(means(:,2)) std(means(:,2))],5);
n_item3 = sampsizepwr('z',[mean(means(:,3)) std(means(:,3))],5);
n_item4 = sampsizepwr('z',[mean(means(:,4)) std(means(:,4))],5);

% [~,p_item(1)] = ztest(means(:,1),0,std(means(:,1)),'tail','right');
% [~,p_item(2)] = ztest(means(:,2),0,std(means(:,2)),'tail','right');
% [~,p_item(3)] = ztest(means(:,3),0,std(means(:,3)),'tail','right');
% [~,p_item(4)] = ztest(means(:,4),0,std(means(:,4)),'tail','right');

[~,p_item2(1)] = ttest2(means(:,1),means(:,2));
[~,p_item2(2)] = ttest2(means(:,1),means(:,3));
[~,p_item2(3)] = ttest2(means(:,1),means(:,4));

figure
hold on
bar(nanmean(means));
errorb(nanmean(means),nanstd(means)./sqrt(sum(~isnan(means))))
set(gca,'Xtick',[1:4])
xlabel('Item Number')
% ylabel('Mean Reacion times (time to fixation) (ms)')
ylabel(['Percentage of Reaction times < 50 ms'])
title(['Relation Memory Sequence Task Reaction times for ' ...
    SeqL_files{1}(1:2) ' n = ' num2str(size(means,1))])

% for ii = 1:4
%     if p_item(ii) < 0.05
%         plot(ii,mean(means(:,ii))+1.5*std(means(:,ii))./sqrt(15),'*k')
%     end
% end
for ii = 1:3
    if p_item2(ii) < 0.05
        plot(ii+1,mean(means(:,ii+1))+2.5*std(means(:,ii+1))./sqrt(15),'*r')
    end
end
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot mean Reaction times ignoring 1st item overlapping for stats---%
sess_means = NaN(length(SeqL_files),4); 
sess_means50 = NaN(length(SeqL_files),4); 
n_trials = zeros(1,4);
overlap = NaN(1,length(SeqL_files));
for ls =1:length(SeqL_files)
    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-SQRTs.mat'],'reaction_time');
    
    [~,~,sequence_locations] = read_SeqL_itm_and_cnd_files(item_sets{ls},SeqL_files{ls}(1:2)); %load item and condition files
    overlap(ls) = find(sum(sequence_locations{1} == sequence_locations{2}) == 2);

    sess_means(ls,:) = 100*sum(reaction_time(21:end,:) < min_rt)./sum(~isnan(reaction_time(21:end,:)));
    sess_means50(ls,:) = 100*sum(reaction_time(21:end,:) < 50)./sum(~isnan(reaction_time(21:end,:)));
    
    if overlap(ls) == 1;
        n_trials(2:4) = n_trials(2:4)+1;
        sess_means(ls,1) = NaN;
        sess_means50(ls,:)= NaN;
    else
        n_trials = n_trials+1;
    end
end

[~,~,predict_ci]= ztest(sess_means(:,1),5,nanstd(sess_means(:,1)));
[~,~,predict_ci50]= ztest(sess_means50(:,1),0,nanstd(sess_means(:,1)));

figure
subplot(1,2,1)
hold on
bar(nanmean(sess_means));
errorb(nanmean(sess_means),nanstd(sess_means)./sqrt(n_trials))
set(gca,'Xtick',[1:4])
xlabel('Item Number')
plot([0 5],[predict_ci(2) predict_ci(2)],'k--')
ylabel(['Percentage of Reaction times < ' num2str(min_rt) ' ms'])
title(['Relation Memory Sequence Task Reaction times for ' ...
    SeqL_files{1}(1:2) ' n = ' num2str(size(means,1))])

subplot(1,2,2)
hold on
bar(nanmean(sess_means50));
errorb(nanmean(sess_means50),nanstd(sess_means50)./sqrt(n_trials))
set(gca,'Xtick',[1:4])
xlabel('Item Number')
plot([0 5],[predict_ci50(2) predict_ci50(2)],'k--')
ylabel(['Percentage of Reaction times < 50 ms'])
title(['Relation Memory Sequence Task Reaction times for ' ...
    SeqL_files{1}(1:2) ' n = ' num2str(size(means,1))])


%%
prop_predict = NaN(1,length(SeqL_files)); %proportion of items predicted
prop_predict50 = NaN(1,length(SeqL_files)); %proportion of items predicted
prop_predict_10 = NaN(1,length(SeqL_files)); %proportion of items predicted
prop_predict50_10 = NaN(1,length(SeqL_files)); %proportion of items predicted
sess_overlap =  NaN(1,length(SeqL_files)); %percentage of predicted for overlapping item
sess_overlap50 =  NaN(1,length(SeqL_files)); %percentage of predicted for overlapping item
sess_max =  NaN(1,length(SeqL_files)); %max percentage per session predicted
sess_max50 =  NaN(1,length(SeqL_files)); %max percentage per session predicted
overlap = NaN(1,length(SeqL_files));
which_sequence = []; %Matlab confused about variable name or function name so place holder
for ls =1:length(SeqL_files)
    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-SQRTs.mat'],'reaction_time','which_sequence');
    
    [~,~,sequence_locations] = read_SeqL_itm_and_cnd_files(item_sets{ls},SeqL_files{ls}(1:2)); %load item and condition files
    overlap(ls) = find(sum(sequence_locations{1} == sequence_locations{2}) == 2)-1;
    
    seq10 = reaction_time(which_sequence(21:end) == 1,:);
    seq10 = seq10(11:end,2:4);
    seq1 = 100*sum(seq10 < min_rt)./sum(~isnan(seq10));
    seq150 = 100*sum(seq10 < 50)./sum(~isnan(seq10));
    
    seq20 = reaction_time(which_sequence(21:end) == 2,:);
    seq20 = seq20(11:end,2:4);
    seq2 = 100*sum(seq20 < min_rt)./sum(~isnan(seq20));
    seq250 = 100*sum(seq20 < 50)./sum(~isnan(seq20));
    
    prop_predict(ls) = sum([seq1 > predict_ci(2)  seq2 >= predict_ci(2)])/6;
    prop_predict50(ls) = sum([seq150 > predict_ci50(2)  seq250 >= predict_ci50(2)])/6;
    
    prop_predict_10(ls) = sum([seq1 > 10  seq2 >10])/6;
    prop_predict50_10(ls) = sum([seq150 > 10  seq250 > 10])/6;
    
    if overlap(ls) ~= 0 %aka 1st item is overlapping
        sess_overlap(ls) = mean([seq1(overlap(ls)) seq2(overlap(ls))]);
        sess_overlap50(ls) = mean([seq150(overlap(ls)) seq250(overlap(ls))]);
    end
    
    sess_max(ls) = max([seq1 seq2]);
    sess_max50(ls) = max([seq150 seq250]);
end


figure
subplot(2,2,1)
hold on
bar([mean(prop_predict) mean(prop_predict_10)])
errorb([mean(prop_predict) mean(prop_predict_10)],...
    [std(prop_predict) std(prop_predict_10)]./sqrt(length(prop_predict)))
hold off
ylabel(['Proportion of Items Predicted (w/ RTs < ' num2str(min_rt) ' ms)'])
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'> CI','>10%'})

subplot(2,2,2)
hold on
bar([mean(prop_predict50) mean(prop_predict50_10)])
errorb([mean(prop_predict50) mean(prop_predict50_10)],...
    [std(prop_predict50) std(prop_predict50_10)]./sqrt(length(prop_predict50)))
hold off
ylabel('Proportion of Items Predicted (w/  RTs < 50 ms)')
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'> CI','>10%'})

subplot(2,2,3)
hold on
bar([nanmean(sess_overlap) nanmean(sess_overlap50)])
errorb([nanmean(sess_overlap) nanmean(sess_overlap50)],...
    [nanstd(sess_overlap) nanstd(sess_overlap50)]./sqrt(sum(~isnan((sess_overlap)))))
hold off
ylabel('% of Predictive Eye Movements on Overlapping Item')
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{['< '  num2str(min_rt) ' ms'],'< 50 ms'})

subplot(2,2,4)
hold on
bar([mean(sess_max) mean(sess_max50)])
errorb([mean(sess_max) mean(sess_max50)],...
    [std(sess_max) std(sess_max50)]./sqrt(length(sess_max)))
hold off
ylabel('Max % of Predictive Eye Movements')
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{['< '  num2str(min_rt) ' ms'],'< 50 ms'})

subtitle(SeqL_files{1}(1:2))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot mean Fixation Durations and Saccade Amplitudes Across Sets---%%%
sac_means =cell(1,2); 
means = cell(1,2);
n = cell(1,2); 
for ls =1:length(SeqL_files)
    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-fixdurs.mat'],...
        'fixation_durations','saccade_amplitudes');
    means{1} = [means{1}; nanmean(fixation_durations{1}(:,2:21))];
    means{2} = [means{2}; nanmean(fixation_durations{2}(:,2:21))];
    n{1} = [n{1} sum(~isnan(fixation_durations{1}'))];
    n{2} = [n{2} sum(~isnan(fixation_durations{2}'))];
    sac_means{1} = [sac_means{1}; nanmean(saccade_amplitudes{1}(:,1:20))];
    sac_means{2} = [sac_means{2}; nanmean(saccade_amplitudes{2}(:,1:20))];
end
figure
hold on
errorbar(nanmean(means{1}),nanstd(means{1})./sqrt(sum(~isnan(means{1}))),'b')
errorbar(nanmean(means{2}),nanstd(means{2})./sqrt(sum(~isnan(means{2}))),'r')
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat')

figure
hold on
errorbar(nanmean(sac_means{1}),nanstd(sac_means{1})./sqrt(sum(~isnan(sac_means{1}))),'b')
errorbar(nanmean(sac_means{2}),nanstd(sac_means{2})./sqrt(sum(~isnan(sac_means{2}))),'r')
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (ms)')
legend('Novel','Repeat')


nn_90 = [];
nn_80 = [];
nr_90 = [];
nr_80 = [];
p = [];
for f = 1:20;
    nn_80(f) = sampsizepwr('t',[mean(means{1}(:,f)) std(means{1}(:,f))], mean(means{2}(:,f)),0.8);
    nn_90(f) = sampsizepwr('t',[mean(means{1}(:,f)) std(means{1}(:,f))], mean(means{2}(:,f)),0.9);
    nr_80(f) = sampsizepwr('t',[mean(means{2}(:,f)) std(means{2}(:,f))], mean(means{1}(:,f)),0.8);
    nr_90(f) = sampsizepwr('t',[mean(means{2}(:,f)) std(means{2}(:,f))], mean(means{1}(:,f)),0.9);
   [~,p(f)] = ttest2( means{1}(:,f),means{2}(:,f));
end

figure
hold on
errorbar(nanmean(sac_means{1}),nanstd(sac_means{1})./sqrt(sum(~isnan(sac_means{1}))),'b')
errorbar(nanmean(sac_means{2}),nanstd(sac_means{2})./sqrt(sum(~isnan(sac_means{2}))),'r')
xlabel('Ordinal Sacccade #')
ylabel('Saccade Amplitude (dva)')
legend('Novel','Repeat')
%% Plot KL divergence 
% 
% KLnorm_mean = NaN(length(SeqL_files),5);
% KLshuff_mean = NaN(length(SeqL_files),5);
% means = NaN(1,4);
% for ls =1:length(SeqL_files)
%    load([datadir SeqL_files{ls}(1:8) '_' SeqL_files{ls}(end) '-KLdivergence.mat'],'KLnorm','KLshuff');
%    KLnorm_mean(ls,:) = nanmean(KLnorm);
%    KLshuff_mean(ls,:) = nanmean(KLshuff);
% end
%%
%% Pre-vs-post lesion Analysis
% all rts
pre_rts = [];
post_rts = [];

% average rts
sess_means = cell(1,2); 
sess_means{1} = NaN(length(pre_sets),4); 
sess_means{2} = NaN(length(post_sets),4);

% percent less than minimum reaction time
sess_minrt = cell(1,2); 
sess_minrt{1} = NaN(length(pre_sets),4); 
sess_minrt{2} = NaN(length(post_sets),4);

% percent less than 50 ms
sess_50 = cell(1,2); 
sess_50{1} = NaN(length(pre_sets),4); 
sess_50{2} = NaN(length(post_sets),4);


for ls =1:length(pre_sets)
    load([datadir pre_sets{ls}(1:8) '_' pre_sets{ls}(end) '-SQRTs.mat'],'reaction_time');
    
    these_rts = reaction_time(21:end,2:4);
    pre_rts = [pre_rts these_rts(1:end)];
    
    sess_means{1}(ls,:) = nanmean(reaction_time(21:end,:));
    sess_minrt{1}(ls,:) = 100*sum(reaction_time(21:end,:) < min_rt)./sum(~isnan(reaction_time(21:end,:)));
    sess_50{1}(ls,:) = 100*sum(reaction_time(21:end,:) < 50)./sum(~isnan(reaction_time(21:end,:))); 
end

for ls =1:length(post_sets)
    load([datadir post_sets{ls}(1:8) '_' post_sets{ls}(end) '-SQRTs.mat'],'reaction_time');
    
    these_rts = reaction_time(21:end,2:4);
    post_rts = [post_rts these_rts(1:end)];
    
    sess_means{2}(ls,:) = nanmean(reaction_time(21:end,:));
    sess_minrt{2}(ls,:) = 100*sum(reaction_time(21:end,:) < min_rt)./sum(~isnan(reaction_time(21:end,:)));
    sess_50{2}(ls,:) = 100*sum(reaction_time(21:end,:) < 50)./sum(~isnan(reaction_time(21:end,:))); 
end

sess_minrt{1} = [sess_minrt{1} mean(sess_minrt{1}(:,2:4),2)];
sess_minrt{2} = [sess_minrt{2} mean(sess_minrt{2}(:,2:4),2)];
sess_means{1} = [sess_means{1} mean(sess_means{1}(:,2:4),2)];
sess_means{2} = [sess_means{2} mean(sess_means{2}(:,2:4),2)];
sess_50{1} = [sess_50{1} mean(sess_50{1}(:,2:4),2)];
sess_50{2} = [sess_50{2} mean(sess_50{2}(:,2:4),2)];

numps = [size(sess_means{1},1)*ones(1,5); size(sess_means{2},1)*ones(1,5)]';

figure
subplot(2,2,1)
hold on
bar([mean(sess_means{1});mean(sess_means{2})]')
errorb([mean(sess_means{1});mean(sess_means{2})]',[std(sess_means{1});std(sess_means{2})]'...
     ./numps);
 for g = 1:5
     [~,p] = ttest2(sess_means{1}(:,g),sess_means{2}(:,g));
     if p < 0.05
         plot(g,mean(sess_means{1}(:,g))+25,'k*');
     end
 end
 
hold off
yl = ylim;
ylim([125 yl(2)]);
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'1','2','3','4','2-4'})
xlim([0.5 5.5])
xlabel('Item #')
ylabel('Reaction time (ms)')

subplot(2,2,2)
hold on
bar([mean(sess_minrt{1});mean(sess_minrt{2})]')
errorb([mean(sess_minrt{1});mean(sess_minrt{2})]',[std(sess_minrt{1});std(sess_minrt{2})]'...
    ./numps);
for g = 1:5
    [~,p] = ttest2(sess_minrt{1}(:,g),sess_minrt{2}(:,g));
    if p < 0.05
        plot(g,mean(sess_minrt{1}(:,g))+5,'k*');
    end
end
hold off
yl = ylim;
xlim([0.5 5.5])
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'1','2','3','4','2-4'})
xlabel('Item #')
ylabel('% of RTs < Min RT')

subplot(2,2,3)
hold on
bar([mean(sess_50{1});mean(sess_50{2})]')
errorb([mean(sess_50{1});mean(sess_50{2})]',[std(sess_50{1});std(sess_50{2})]'...
     ./numps);
 for g = 1:5
    [~,p] = ttest2(sess_50{1}(:,g),sess_50{2}(:,g));
    if p < 0.05
        plot(g,mean(sess_50{1}(:,g))+2,'k*');
    end
end
hold off
yl = ylim;
xlim([0.5 5.5])
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'1','2','3','4','2-4'})
xlabel('Item #')
ylabel('% of RTs < 50')
legend('Pre-Lesion','Post-Lesion','Location','NorthEastOutside')


subtitle(['Monkey: ' pre_sets{1}(1:2)])


figure
subplot(2,2,1)
[n1,x1] = hist(pre_rts,250);
hist(pre_rts,250);
title('Pre-Lesion')
xlim([-250 600])
xlabel('Reaction Time (ms)')
ylabel('Count')

n1 = n1/sum(n1);

subplot(2,2,2)
[n2,x2] = hist(post_rts,250);
hist(post_rts,250)
title('Post-Lesion')
xlim([-250 600])
xlabel('Reaction Time (ms)')
ylabel('Count')

n2 = n2/sum(n2);

subplot(2,2,[3 4])
hold on
plot(x1,n1)
plot(x2,n2,'r')
hold off
xlim([-250 600])
xlabel('Reaction Time (ms)')
ylabel('Fraction')
legend('Pre-Lesion','Post-Lesion')

subtitle(['Monkey: ' pre_sets{1}(1:2)])
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot mean Fixation Durations and Saccade Amplitudes Across Sets---%%%
sac_means =cell(2,2); 
means = cell(2,2);
n = cell(2,2); 
for ls =1:length(pre_sets)
    load([datadir pre_sets{ls}(1:8) '_' pre_sets{ls}(end) '-fixdurs.mat'],...
        'fixation_durations','saccade_amplitudes');
    
    means{1,1} = [means{1,1}; nanmean(fixation_durations{1}(:,1:20))];
    means{1,2} = [means{1,2}; nanmean(fixation_durations{2}(:,1:20))];
    n{1,1} = [n{1,1} sum(~isnan(fixation_durations{1}'))];
    n{1,2} = [n{1,2} sum(~isnan(fixation_durations{2}'))];
    sac_means{1,1} = [sac_means{1,1}; nanmean(saccade_amplitudes{1}(:,1:20))];
    sac_means{1,2} = [sac_means{1,2}; nanmean(saccade_amplitudes{2}(:,1:20))];
end

for ls =1:length(post_sets)
    load([datadir post_sets{ls}(1:8) '_' post_sets{ls}(end) '-fixdurs.mat'],...
        'fixation_durations','saccade_amplitudes');
    
    means{2,1} = [means{2,1}; nanmean(fixation_durations{1}(:,1:20))];
    means{2,2} = [means{2,2}; nanmean(fixation_durations{2}(:,1:20))];
    n{2,1} = [n{2,1} sum(~isnan(fixation_durations{1}'))];
    n{2,2} = [n{2,2} sum(~isnan(fixation_durations{2}'))];
    sac_means{2,1} = [sac_means{2,1}; nanmean(saccade_amplitudes{1}(:,1:20))];
    sac_means{2,2} = [sac_means{2,2}; nanmean(saccade_amplitudes{2}(:,1:20))];
end

figure
hold on
errorbar(nanmean(means{1,1}),nanstd(means{1,1})./sqrt(sum(~isnan(means{1,1}))),'b')
errorbar(nanmean(means{1,2}),nanstd(means{1,2})./sqrt(sum(~isnan(means{1,2}))),'r')
errorbar(nanmean(means{2,1}),nanstd(means{2,1})./sqrt(sum(~isnan(means{2,1}))),'k')
errorbar(nanmean(means{2,2}),nanstd(means{2,2})./sqrt(sum(~isnan(means{2,2}))),'g')
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel Pre','Repeat Pre','Novel Post','Repeat Post')

figure
hold on
errorbar(nanmean(sac_means{1,1}),nanstd(sac_means{1,1})./sqrt(sum(~isnan(sac_means{1,1}))),'b')
errorbar(nanmean(sac_means{1,2}),nanstd(sac_means{1,2})./sqrt(sum(~isnan(sac_means{1,2}))),'r')
errorbar(nanmean(sac_means{2,1}),nanstd(sac_means{2,1})./sqrt(sum(~isnan(sac_means{2,1}))),'k')
errorbar(nanmean(sac_means{2,2}),nanstd(sac_means{2,2})./sqrt(sum(~isnan(sac_means{2,2}))),'g')
xlabel('Ordinal Fixation #')
ylabel('Saccade Amplitude (dva)')
legend('Novel Pre','Repeat Pre','Novel Post','Repeat Post')