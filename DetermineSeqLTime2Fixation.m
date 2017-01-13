function DetermineSeqLTime2Fixation(ListSQdatafile)
%written by Seth Konig December, 2014
%code determines how long it takes to fixation an item in the sequence
%portion of the ListSQ task. For cortex/behavioral only days. Updated
%2/15/2015 by Seth Konig. Old code may have had some bugs and less
%efficient. Main analysis replaced by code that works across all tasks that
% use sequence trials including recording days.

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Figures\';%where to save figures

load(ListSQdatafile)

[itmlist,sequence_items,sequence_locations] = read_SeqL_itm_and_cnd_files(item_set,ListSQdatafile(end-22:end-21)); %load item and condition files

%convert item locations back into dva.
sequence_locations{1} = (sequence_locations{1}-[400;300]*ones(1,4))/24;
sequence_locations{2} = (sequence_locations{2}-[400;300]*ones(1,4))/24;

num_trials = length(fixationstats);
fixwin = 5;
reaction_time = NaN(length(fixationstats),4); %how long after a stimulus disappears do they fixate
cortex_rt = NaN(length(fixationstats),4); %how long after a stimulus disappears do they fixate according to cortex
fixation_accuracy = NaN(length(fixationstats),4); %how far off
extrafixations = NaN(length(fixationstats),4); %how many extra fixations to get to next item
fixation_numbers = NaN(length(fixationstats),4);%which fixation number is associated with which item
which_sequence = NaN(length(fixationstats),1); %sequence 1 or sequence 2
cortex_predicted = NaN(length(fixationstats),4);%was eye movement truly predictive according to cortex
cortexbreak =NaN(length(fixationstats),4);%monkey broke fixation then refixated
% figure
% pos = zeros(600,800); 
for trial = 1:num_trials
    if ~any(trial == not_first_showing) %if they've seen any part of the 1st trial don't analyze it
        if any(itmlist(per(trial).cnd) == sequence_items)
            if itmlist(per(trial).cnd) == sequence_items(1)
                sequence = 1;
            else
                sequence = 2;
            end
            
            if ~isempty(find(per(trial).allval == 30))%had a small error with compiling save file for red in
                % which 2-3 trial occured before if fixed it but she got
                % rewarded anyway.
                 
%                 xy = fixationstats{trial}.XY;
%                 xy(1,:) = round(24*xy(1,:)+400);
%                 xy(2,:) = round(24*xy(2,:)+300);
%                 xy(xy < 1) = NaN;
%                 xy(1,xy(1,:) > 800) = NaN;
%                 xy(2,xy(2,:) > 600) = NaN;
%                 ind = sub2ind([600 800],xy(2,:),xy(1,:));
%                 ind(isnan(ind)) = [];
%                 for in = 1:length(ind);
%                     pos(ind) = pos(ind)+1;
%                 end
%                 if sequence == 1
%                     subplot(1,2,1)
%                     hold on
%                     plot(xy(1,:),xy(2,:))
%                     hold off
%                 else
%                     subplot(1,2,2)
%                     hold on
%                     plot(xy(1,:),xy(2,:))
%                     hold off
%                 end
                
                trialdata = analyze_sequence_trial(fixationstats{trial},sequence_locations{sequence},...
                    fixwin,per(trial).allval,per(trial).alltim,0);
                
                fixation_numbers(trial,:) = trialdata.fixationnums;
                reaction_time(trial,:) = trialdata.t2f;
                cortex_rt(trial,:) = trialdata.cortext2f;
                fixation_accuracy(trial,:) = trialdata.accuracy;
                extrafixations(trial,:) = trialdata.extrafixations;
                cortex_predicted(trial,:) = trialdata.cortexpredict;
                which_sequence(trial) = sequence;
                cortexbreak(trial,:) = trialdata.cortexbreak;
            end
        end
    end
end

% subplot(1,2,1)
% axis equal
% xlim([-13 13])
% ylim([-10 10])
% subplot(1,2,2)
% axis equal
% xlim([-13 13])
% ylim([-10 10])

%remove image trials
image_trials = 21:5:size(reaction_time,1);
fixation_numbers(image_trials,:) = [];
reaction_time(image_trials,:) = [];
fixation_accuracy(image_trials,:) = [];
extrafixations(image_trials,:) = [];
cortex_predicted(image_trials,:) = [];
which_sequence(image_trials,:) = [];
cortex_rt(image_trials,:) = [];
cortexbreak(image_trials,:) = [];


%for figures below ignoring the 1st 20 trials since these are the
%familizarization block

figure
%reaction times
subplot(2,2,1)
hold on
bar(nanmean(reaction_time(21:end,:)),'r')
errorb(nanmean(reaction_time(21:end,:)),nanstd(reaction_time(21:end,:))...
    ./sqrt(sum(~isnan(reaction_time(21:end,:)))));
hold off
title('Average Reaction time by item')
xlabel('Item #')
ylabel('Reaction time/time to fixation (ms)')
[y_lims] = ylim;
if y_lims(1) > 0;
    ylim([0 y_lims(2)]);
end
set(gca,'Xtick',1:4)

%cortex predicted
subplot(2,2,2)
predict = cortex_predicted(21:end,:);
bar(100*nansum(predict)./sum(~isnan(predict)),'r')
xlabel('Item #')
ylabel('Percentage')
title('Percent of trials with "Truly" Predictive Saccades')

%reacion times < 150 ms
subplot(2,2,3)
rts = reaction_time(21:end,:);
bar(100*nansum(rts < 150)./sum(~isnan(rts)),'r')
xlabel('Item #')
ylabel('Percentage')
title('Percent of trials with reaction times less than 150 ms')

%fixation accuracy
subplot(2,2,4)
hold on
plot(nanmean(fixation_accuracy(which_sequence == 1,:),2),'.')
plot(nanmean(fixation_accuracy(which_sequence == 2,:),2),'r.')
hold off
xlabel('"Trial #"')
ylabel('Fixation Accuracy (dva)')
title('Fixation Accuracy by Trial and Sequence (by color)')

save_and_close_fig(figure_dir,[ListSQdatafile(end-22:end-13) '-SessionDataSummary'])

save([ListSQdatafile(1:end-13) '-SQRTs.mat'],'reaction_time','fixwin',...
    'fixation_accuracy','extrafixations','which_sequence','fixation_numbers','cortex_predicted','cortex_rt','cortexbreak')