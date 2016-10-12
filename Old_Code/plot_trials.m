function [] = plot_trials(plotfun)
% this function plots an array of single-trial trajectories
% it takes the parameters for figure array size and 

figure()

% plot large figure
i = 16;
    subplot(5,4,[15 16 19 20])
    hold on
    dataa = area(1:60,evac(:,i),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')

    confinterval=shadedErrorBar_2(T_power{i},mProbs_power{i},3*stdevs{i},{'color','k','LineWidth',5,'LineStyle',':'},1);
    set(confinterval.edge(1),'visible','off')
    set(confinterval.edge(2),'visible','off')
    [ax,loop,php]=plotyy(1:60,P_LOO_power(:,i),1:60,rP_hits(:,i));
    set(loop,'color',[0 1 1 0.4]);
    set(loop,'LineWidth',5)
    set(php,'color',[0 1 0])
    set(php,'LineWidth',1)
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50;
    ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50];
    ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(i)];
    ax(2).XLim = [1 endTimes(i)];
    ax(1).FontSize = 13;
    ax(2).FontSize = 13;
    ax(2).YLabel.String = 'likelihood';
    ax(2).YLabel.FontSize=18;
    xlabel('time','FontSize',18)
    ylabel('evacuations','FontSize',18)
    title(['Trial Ind50-' num2str(i)],'FontSize',14)
    legend([dataa, confinterval.mainLine,confinterval.edge(1),loop,php],{'   data','   model average','   99.7% confidence','   LOOCV','   P_{hit}'},'location','northwest','fontsize',17)
    
for i = 1:length(z)-2
    subplot(5,4,i) 
    hold on
    dataa = area(1:60,evac(:,i),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')    
    confinterval=shadedErrorBar_2(T_power{i},mProbs_power{i},3*stdevs{i},{'color','k','LineWidth',5,'LineStyle',':'},1);
    set(confinterval.edge(1),'visible','off')
    set(confinterval.edge(2),'visible','off')
    [ax,loop,php]=plotyy(1:60,P_LOO_power(:,i),1:60,rP_hits(:,i));
    set(loop,'color',[0 1 1 0.4]);
    set(loop,'LineWidth',5)
    set(php,'color',[0 1 0])
    set(php,'LineWidth',1)
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50;
    ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50];
    ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(i)];
    ax(2).XLim = [1 endTimes(i)];
    ax(1).FontSize = 11;
    ax(2).FontSize = 11;
    title(['Trial Ind50-' num2str(i)],'FontSize',14)
end

i = 15;
    subplot(5,4,18) 
    hold on
    dataa = area(1:60,evac(:,i),'FaceColor',[192/255 192/255 192/255]);
    set(gca,'Layer','top')    
    confinterval=shadedErrorBar_2(T_power{i},mProbs_power{i},3*stdevs{i},{'color','k','LineWidth',5,'LineStyle',':'},1);
    set(confinterval.edge(1),'visible','off')
    set(confinterval.edge(2),'visible','off')
    [ax,loop,php]=plotyy(1:60,P_LOO_power(:,i),1:60,rP_hits(:,i));
    set(loop,'color',[0 1 1 0.4]);
    set(loop,'LineWidth',5)
    set(php,'color',[0 1 0])
    set(php,'LineWidth',1)
    set(ax(1),'box','on')
    set(ax(2),'YColor','g')
    ax(1).YTick = 0:10:50;
    ax(2).YTick = 0:0.2:1;
    ax(1).YColor = 'k';
    ax(1).YLim = [0 50];
    ax(2).YLim = [0 1];
    ax(1).XLim = [1 endTimes(i)];
    ax(2).XLim = [1 endTimes(i)];
    ax(1).FontSize = 11;
    ax(2).FontSize = 11;
    title(['Trial Ind50-' num2str(i)],'FontSize',14)

set(gcf,'units','pixels')
set(gcf,'position',[0 0 1200 800])
tightfig;
pos=get(gcf,'position');
set(gcf,'position',[pos(1:3) 800])
tightfig;
set(gcf,'units','centimeters')
set(gcf,'position',[pos(1:2) 40 28])
set(gcf,'papersize',[40 28])
savefig('individualtrials')

end