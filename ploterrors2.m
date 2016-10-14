% Calculate and plot absolute evacuation errors (# evacuated in model - #
% evacuated in experiment) for Ind50 and Ind25 trials

load data/mle_Ind50 P_Ind50 evac
evac50 = evac;
load data/mle_Ind25_test2 P_Ind25 evac
evac25 = evac;
clear evac;

res50 = P_Ind50 - evac50;
res25 = P_Ind25 - evac25;

% Plot Ind50 errors
figure('position',[0 0 1200 800])
s1 = subplot(2,2,1);
plot(1:60,zeros(60,1),'color','k') % x-axis
hold on
for i = 1:7
    plot(1:60,res50(:,i),'LineStyle',':','color',[1 0 0 0.4])
    hold on
end
sh50 = shadedErrorBar(1:60,nanmean(res50,2),nanstd(res50,0,2),{'color','r'},1,1);
set(sh50.edge(1),'visible','off'); set(sh50.edge(2),'visible','off');
r50 = plot(1:60,mean(res50,2),'color','r','LineWidth',4);
xlabel('Time','FontSize',26)
ylabel('Error','FontSize',26)
set(gca,'ytick',-50:20:50);
set(gca,'fontsize',22)  
xlim([1 60])
ylim([-50 50])
title('Number evacuated error, Ind50 trials','FontSize',26)

% Plot Ind25 errors
s2 = subplot(2,2,2);
plot(1:60,zeros(60,1),'color','k') % x-axis
hold on
for i = 1:13
    plot(1:60,res25(:,i),'LineStyle',':','color',[0 0.8 0 0.6])
    hold on
end
sh25 = shadedErrorBar(1:60,nanmean(res25,2),nanstd(res25,0,2),{'color','g'},1,1);
set(sh25.edge(1),'visible','off'); set(sh25.edge(2),'visible','off');
r25 = plot(1:60,nanmean(res25,2),'color','g','LineWidth',4);
xlabel('Time','FontSize',26)
ylabel('Error','FontSize',26)
set(gca,'ytick',-25:10:25);
set(gca,'fontsize',22)  
xlim([1 60])
ylim([-25 25])
title('Number evacuated error, Ind25 trials','FontSize',26)

% Plot Ind50 final errors
s3 = subplot(2,2,3);
load data/mle_Ind50 rP_hits
bar(res50(end,:));
xlabel('Trial','FontSize',26)
ylabel('Error','FontSize',26)
ybuff=2.5;
XDATA=get(get(gca,'Children'),'XData');
YDATA=get(get(gca,'Children'),'YData');
for j=1:size(XDATA,2)
    x=XDATA(j);
    if YDATA(j)>0
        y=YDATA(j)+ybuff;
    else
        y = YDATA(j)-ybuff;
    end
    if rP_hits(end,j) == 1
        t = 'H';
    else
        t = 'M';
    end
    text(x,y,t,'Color','k','HorizontalAlignment','center','fontsize',22)
end
xlim([0.25 16.75])
ylim([-15 35])
set(gca,'xtick',2:2:16)
set(gca,'fontsize',22)  
title('Final evacuation error, Ind50 trials','FontSize',26)

% Plot Ind25 final errors
s4 = subplot(2,2,4);
load data/mle_Ind25 rP_hits
bar(res25(end,:),'m');
xlabel('Trial','FontSize',26)
ylabel('Error','FontSize',26)
ybuff=2.5;
XDATA=get(get(gca,'Children'),'XData');
YDATA=get(get(gca,'Children'),'YData');
for j=1:size(XDATA,2)
    x=XDATA(j);
    if YDATA(j)>0
        y=YDATA(j)+ybuff;
    else
        y = YDATA(j)-ybuff;
    end
    if rP_hits(end,j) == 1
        t = 'H';
    else
        t = 'M';
    end
    text(x,y,t,'Color','k','HorizontalAlignment','center','fontsize',22)
end
xlim([0.25 13.75])
ylim([-15 35])
set(gca,'fontsize',22)  
title('Final evacuation error, Ind25 trials', 'FontSize',26)
tightfig;

ps1 = s1.Position;
ps2 = s2.Position;
ps3 = s3.Position;
ps4 = s4.Position;

annotation('textbox',[ps1(1) ps1(2)+ps1(4)-0.1 0.1 0.1],'string','A','LineStyle','none','FontSize',30,'FontWeight','bold')
annotation('textbox',[ps2(1) ps2(2)+ps2(4)-0.1 0.1 0.1],'string','B','LineStyle','none','FontSize',30,'FontWeight','bold')
annotation('textbox',[ps3(1) ps3(2)+ps3(4)-0.1 0.1 0.1],'string','C','LineStyle','none','FontSize',30,'FontWeight','bold')
annotation('textbox',[ps4(1) ps4(2)+ps4(4)-0.1 0.1 0.1],'string','D','LineStyle','none','FontSize',30,'FontWeight','bold')

% savefig('newplots/Ind_errors')
% print('newplots/Ind_errors','-dpdf','-r300')
% print('newplots/Ind_errors','-dsvg','-r300')


%%
load data/mle_Ind50 rP_hits
figure()
final_error_hit = mean(res50(end,rP_hits(end,:)==1));
final_error_miss = mean(res50(end,rP_hits(end,:)==0));
h2 = bar([1 2],[final_error_hit final_error_miss],'r');
hold on
errorbar([final_error_hit final_error_miss],[std(res50(end,rP_hits(end,:)==1)),std(res50(end,rP_hits(end,:)==0))],'linestyle','none','color','k')
xlim([0.5 2.5])
ylim([-10 15])
set(gca,'XTickLabel',{'H','M'},'fontsize',50)
ylabel('Mean final error','fontsize',50)
tightfig;
set(gcf,'position',[0 0 535 450])
tightfig;
% savefig('figures/meanfinalerror_50')
print('newplots/meanfinalerror_50','-dsvg','-r300')
%%
load data/mle_Ind25_test2 rP_hits
figure()
final_error_hit = mean(res25(end,rP_hits(end,:)==1));
final_error_miss = mean(res25(end,rP_hits(end,:)==0));
bar([1 2],[final_error_hit final_error_miss],'r');
hold on
errorbar([final_error_hit final_error_miss],[std(res25(end,rP_hits(end,:)==1)),std(res25(end,rP_hits(end,:)==0))],'linestyle','none','color','k')
xlim([0.5 2.5])
ylim([-10 15])
set(gca,'XTickLabel',{'H','M'},'fontsize',50) 
ylabel('Mean final error','fontsize',50) 
tightfig;
set(gcf,'position',[0 0 535 450])
tightfig;
% savefig('figures/meanfinalerror_25')
print('newplots/meanfinalerror_25','-dsvg','-r300')


