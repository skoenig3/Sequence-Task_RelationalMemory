% fixation duration by image #

all_nov_fix_durs = [];
all_rep_fix_durs = [];
for set = 1:size(first_fix,2)
    all_nov_fix_durs = [all_nov_fix_durs; first_fix{set}(1,:)];
    all_rep_fix_durs = [all_rep_fix_durs; first_fix{set}(2,:)]; 
end

figure
hold on
plot(nanmean(all_nov_fix_durs))
plot(nanmean(all_rep_fix_durs),'r')
plot([1 1],[130 220],'k--')   %first image block 1
plot([13 13],[130 220],'k--') %first image block 2
hold off
xlabel('Image #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat','Block Start')

%% fixation duration averaged over 12 images by block #
an = [];
ar = [];

for block = 1:2;
    an = [an; nanmean(nanmean(all_nov_fix_durs(:,12*(block-1)+1:block*12)))];
    ar = [ar; nanmean(nanmean(all_rep_fix_durs(:,12*(block-1)+1:block*12)))];
end

figure
hold on
plot(an)
plot(ar,'r')
hold off
xlabel('Block #')
ylabel('Fixation Duration (ms)')

%% fixation duration by image within block
anb = cell(1,12);
arb = cell(1,12);

for img = 1:12;
    anb{img} = [anb{img}; all_nov_fix_durs(:,img:12:end)];
    arb{img} = [arb{img}; all_rep_fix_durs(:,img:12:end)];
end

anbm = NaN(1,12);
arbm = NaN(1,12);
for img = 1:12
    anbm(img) = nanmean(anb{img}(1:end));
    arbm(img) = nanmean(arb{img}(1:end));
end

figure
hold on
plot(anbm)
plot(arbm,'r')
hold off
xlabel('Image #')
ylabel('Fixation Duration (ms)')


%% by Set number
set_means = [];
for set = 1:length(first_fix)
    set_means(1,set) = nanmean(first_fix{set}(1,:));
    set_means(2,set) = nanmean(first_fix{set}(2,:));
end
figure
subplot(3,1,1)
plot(set_means(:,1:15)')
xlabel('Data Session')
ylabel('Fixation Duration (ms)')
title('Vivian')

subplot(3,1,2)
plot(set_means(:,16:30)')
xlabel('Data Session')
ylabel('Fixation Duration (ms)')
title('Red')

subplot(3,1,3)
plot(set_means(:,31:end)')
xlabel('Data Session')
ylabel('Fixation Duration (ms)')
title('Tobii')

%% get fixation pdfs for first fixation
nov_map = zeros(200,800);
rep_map = zeros(200,800);
for set = 1:length(first_fix_location);
    for imgnum = 1:92
        if ~isnan(first_fix_location{1,set}(1,imgnum))
            xy = floor(first_fix_location{1,set}(:,imgnum));
            xy(xy < 1) = 1;
            xy(2,xy(2) > 200) = 200;
            xy(1,xy(1) > 800) = 800;
            nov_map(xy(2),xy(1)) = nov_map(xy(2),xy(1))+1;
            
        end
        
        if ~isnan(first_fix_location{2,set}(1,imgnum))
            xy = floor(first_fix_location{2,set}(:,imgnum));
            xy(xy < 1) = 1;
            xy(2,xy(2) > 200) = 200;
            xy(1,xy(1) > 800) = 800;
            rep_map(xy(2),xy(1)) = rep_map(xy(2),xy(1))+1;
        end
    end
end
