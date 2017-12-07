%% Pre-vs-post lesion Analysis
clar


datadir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Cortex Data\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Figures\';

%---PW---%
% pre_sets = {'PW150324.2','PW150325.2','PW150326.2','PW150327.2','PW150330.2',...
%               'PW150331.2','PW150401.2','PW150402.2','PW150403.2','PW150406.2',...
%               'PW150407.2','PW150409.2','PW150413.2','PW150414.2','PW150415.2'};
% post_sets = {'PW160125.2','PW160126.2','PW160127.2','PW160128.2','PW160129.2',...
%                'PW160201.2','PW160202.2','PW160208.2','PW160209.2','PW160210.2',...
%                'PW160212.2','PW160216.2','PW160217.2','PW160218.2','PW160219.2'};
% min_rt = 156;%mininum reaction time varies by monkey. Taken as 5 percentile for random sequences


%---RR---%
% pre_sets = {'RR150324.2','RR150325.2','RR150326.2','RR150327.2','RR150330.2',...
%               'RR150331.2','RR150401.2','RR150408.2','RR150409.2','RR150410.2',...
%               'RR150413.2','RR150414.2','RR150415.2','RR150417.2','RR150422.2'};
% post_sets = {'RR160210.1','RR160211.1','RR160212.1','RR160217.1','RR160218.1',...
%               'RR160219.1','RR160222.1','RR160223.1','RR160224.1','RR160225.1',...
%               'RR160226.1','RR160229.1','RR160301.1','RR160302.1','RR160303.1'};
% min_rt = 109; 

%---TO---%
% pre_sets = {'TO150421.3','TO150422.3','TO150423.2','TO150424.2','TO150427.3',...
%               'TO150428.2','TO150429.2','TO150430.2','TO150501.2','TO150504.2',...
%               'TO150505.2','TO150506.2','TO150507.2','TO150508.2','TO150512.2'};
% post_sets = {'TO170228.2','TO170301.2','TO170306.2','TO170307.2','TO170308.2',...
%               'TO170309.2','TO170310.2','TO170313.2','TO170314.2','TO170315.2',...
%               'TO170316.2','TO170317.2','TO170320.2','TO170321.2','TO170322.2'};
% min_rt = 135;


%MF pre-lesion

pre_sets = {'MF161209.2','MF161212.2','MF161213.2','MF161214.2',...
    'MF161215.2','MF161219.2','MF161220.2','MF161221.2',...
    'MF161229.2','MF161230.2','MF170103.2','MF170104.3',...
    'MF170105.2','MF170106.2','MF170109.2','MF170110.2'};
post_sets = {'MF161209.2','MF161212.2','MF161213.2','MF161214.2',...
    'MF161215.2','MF161219.2','MF161220.2','MF161221.2',...
    'MF161229.2','MF161230.2','MF170103.2','MF170104.3',...
    'MF170105.2','MF170106.2','MF170109.2','MF170110.2'};
min_rt = 112;




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
%%

npre = length(pre_sets);
npost = length(post_sets);
%%


ids = [[ones(npre,1); 2*ones(npre,1); 3*ones(npre,1); ones(npost,1); 2*ones(npost,1); 3*ones(npost,1)] ...
      [ones(npre,1); ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1); 2*ones(npost,1)]];
  
  
vals = [sess_means{1}(:,2); sess_means{1}(:,3); sess_means{1}(:,4); ...
    sess_means{2}(:,2); sess_means{2}(:,3); sess_means{2}(:,4)];
[P_ANOVA_rt] = anovan(vals,ids,'model','interaction','varnames',{'Item #','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [sess_minrt{1}(:,2); sess_minrt{1}(:,3); sess_minrt{1}(:,4); ...
    sess_minrt{2}(:,2); sess_minrt{2}(:,3); sess_minrt{2}(:,4)];
[P_ANOVA_minrt] = anovan(vals,ids,'model','interaction','varnames',{'Item #','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [sess_50{1}(:,2); sess_50{1}(:,3); sess_50{1}(:,4); ...
    sess_50{2}(:,2); sess_50{2}(:,3); sess_50{2}(:,4)];
[P_ANOVA_rt50] = anovan(vals,ids,'model','interaction','varnames',{'Item #','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures
%%

figure
subplot(2,2,1)
errorb([mean(sess_means{1});mean(sess_means{2})]',[std(sess_means{1});std(sess_means{2})]'...
     ./numps);
 hold on
 for g = 1:5
     [~,p] = ttest2(sess_means{1}(:,g),sess_means{2}(:,g));
     if p < 0.05
         plot(g,mean(sess_means{1}(:,g))+25,'k*');
     end
 end
hold off
yl = ylim;
ylim([100 yl(2)]);
set(gca,'Xtick',1:5)
set(gca,'XtickLabel',{'1','2','3','4','2-4'})
xlim([0.5 5.5])
xlabel('Item #')
ylabel('Reaction time (ms)')
box off
title(sprintf(['ANOVA Items 2-4: \n p_{lesion} = ' num2str(P_ANOVA_rt(2),2) ...
    ' , p_{Item} = ' num2str(P_ANOVA_rt(1),2) ', p_{inter} = ' num2str(P_ANOVA_rt(3),2)]))

subplot(2,2,2)
errorb([mean(sess_minrt{1});mean(sess_minrt{2})]',[std(sess_minrt{1});std(sess_minrt{2})]'...
    ./numps);
hold on
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
box off
title(sprintf(['ANOVA Items 2-4: \n p_{lesion} = ' num2str(P_ANOVA_minrt(2),2) ...
    ' , p_{Item} = ' num2str(P_ANOVA_minrt(1),2) ', p_{inter} = ' num2str(P_ANOVA_minrt(3),2)]))

subplot(2,2,3)
errorb([mean(sess_50{1});mean(sess_50{2})]',[std(sess_50{1});std(sess_50{2})]'...
     ./numps);
 hold on
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
box off
title(sprintf(['ANOVA Items 2-4: \n p_{lesion} = ' num2str(P_ANOVA_rt50(2),2) ...
    ' , p_{Item} = ' num2str(P_ANOVA_rt50(1),2) ', p_{inter} = ' num2str(P_ANOVA_rt50(3),2)]))

[n1,x1] = hist(pre_rts,-250:5:400);
n1 = n1/sum(n1);
[n2,x2] = hist(post_rts,-250:5:400);
n2 = n2/sum(n2);

subplot(2,2,4)
hold on
plot(x1,n1)
plot(x2,n2,'r')
hold off
xlim([-250 400])
xlabel('Reaction Time (ms)')
ylabel('Fraction')
legend('Pre-Lesion','Post-Lesion')


subtitle(['Pre vs Post: ' pre_sets{1}(1:2)]) 