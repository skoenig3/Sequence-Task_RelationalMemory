function PlotEyeDensity(ListSQdatafile,figure_dir)

load(ListSQdatafile) 
[itmlist,sequence_items,sequence_locations] = read_SeqL_itm_and_cnd_files(item_set,ListSQdatafile(end-22:end-21)); %load item and condition files

num_trials = length(fixationstats);

seq1_pos = zeros(600,800);
seq2_pos = zeros(600,800);

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
                
                xy = fixationstats{trial}.XY;
                xy(1,:) = round(24*xy(1,:)+400);
                xy(2,:) = round(24*xy(2,:)+300);
                xy(xy < 1) = NaN;
                xy(1,xy(1,:) > 800) = NaN;
                xy(2,xy(2,:) > 600) = NaN;
                ind = sub2ind([600 800],xy(2,:),xy(1,:));
                ind(isnan(ind)) = [];
                
                if sequence == 1
                    seq1_pos(ind) = seq1_pos(ind)+1;
                else
                    seq2_pos(ind) = seq2_pos(ind)+1;
                end
            end
        end
    end
end

shape = {'kx','ko'};

figure
subplot(1,2,1)
hold on  
imagesc(log10(seq1_pos))
for s = 1:2
    for item = 1:4
        plot(sequence_locations{s}(1,item),sequence_locations{s}(2,item),shape{s},'markersize',15)
    end
end
box off
axis off
title('Sequence 1')

subplot(1,2,2)
hold on  
imagesc(log10(seq2_pos))
for s = 1:2
    for item = 1:4
        plot(sequence_locations{s}(1,item),sequence_locations{s}(2,item),shape{s},'markersize',15)
    end
end
box off
axis off
title('Sequence 2')

subtitle('ListSQdatafile')
save_and_close_fig(figure_dir,[ListSQdatafile(end-22:end-13) '-EyeDensity'])