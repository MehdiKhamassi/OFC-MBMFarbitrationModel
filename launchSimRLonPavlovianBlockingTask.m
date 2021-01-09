%% launch simRLonPavlovianBlockingTask (with OFC inactivation)
nbExperiments = 1000;

%% LAUNCH X experiments
totalprob = [];
totalY = [];
for ttt=1:nbExperiments
    % experiment starts
    simRLonPavlovianBlockingTask
    % experiment finished
    % record logs:
    totalprob = [totalprob P4];
    totalY = [totalY Y4];
end

%% plot figure
% figure parameters
if (~ofclesion)
    if (condAB)
        figure
    end
    markerColor = 'w';
    condColor = 'k';
    offset = 0;
    styletrait = '--';
else
    markerColor = 'k';
    condColor = 'r';
    offset = 4;
    styletrait = '';
end
if (condAB)
    stylemarque = '-o';
else
    stylemarque = '--^';
end

if (condAB) % no need to replot the same beginning for stim C and D
    % Stage 1
    subplot(2,5,1:2)
    if (ofclesion||~condAB)
        hold on
    end
    plot(1:nbDayS1,mean(totalprob(1:nbDayS1,1:3:end-2)')',[styletrait 'b'],'LineWidth',2) % P4(1:nbDayS1,1)
    hold on
    plot(1:nbDayS1,mean(totalprob(1:nbDayS1,2:3:end-1)')',[styletrait 'r'],'LineWidth',2) % P4(1:nbDayS1,2)
    plot(1:nbDayS1,mean(totalprob(1:nbDayS1,3:3:end)')',[styletrait 'k'],'LineWidth',2) % P4(1:nbDayS1,3)
    plot(nbDayS1+2:nbDayS1+nbDayIn+1,mean(totalprob(nbDayS1+1:nbDayS1+nbDayIn,1:3:end-2)')',[styletrait 'b'],'LineWidth',2) % P4(nbDayS1+1:nbDayS1+nbDayIn,1)
    plot(nbDayS1+2:nbDayS1+nbDayIn+1,mean(totalprob(nbDayS1+1:nbDayS1+nbDayIn,2:3:end-1)')',[styletrait 'r'],'LineWidth',2) % P4(nbDayS1+1:nbDayS1+nbDayIn,2)
    plot(nbDayS1+2:nbDayS1+nbDayIn+1,mean(totalprob(nbDayS1+1:nbDayS1+nbDayIn,3:3:end)')',[styletrait 'k'],'LineWidth',2) % P4(nbDayS1+1:nbDayS1+nbDayIn,3)

    legend('MB','MF','EXP','Location','northwest')
    axis([0 nbDayS1+nbDayIn+2 0 1])
    xticks([1:nbDayS1 nbDayS1+2:nbDayS1+nbDayIn+1])
    xticklabels([1:nbDayS1 nbDayS1+1:nbDayS1+nbDayIn])
    yticks([0:0.2:1])
    xlabel('Day','FontSize',14)
    ylabel('Probability of Selection','FontSize',14)
    if (~ofclesion&&condAB)
        title('Stage 1','FontSize',18)
    end
    % Stage 2
    subplot(2,5,3:4)
    if (ofclesion||~condAB)
        hold on
    end
    plot(offset+1:offset+nbDayS2,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,1:3:end-2)')',[styletrait 'b'],'LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,1)
    hold on
    plot(offset+1:offset+nbDayS2,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,2:3:end-1)')',[styletrait 'r'],'LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,2)
    plot(offset+1:offset+nbDayS2,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,3:3:end)')',[styletrait 'k'],'LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,3)
    axis([0 8 0 1])
    xticks([1:3 5:7])
    xticklabels([12:14 12:14])
    yticks([0:0.2:1])
    xlabel('Day','FontSize',14)
    if (~ofclesion&&condAB)
        title('Stage 2','FontSize',18)
    end
    % Test
    subplot(2,5,5)
    if (ofclesion||~condAB)
        hold on
    end
    plot(1,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,1:3:end-2)')','ob','LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,1)
    hold on
    plot(1,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,2:3:end-1)')','or','LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,2)
    plot(1,mean(totalprob(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,3:3:end)')','ok','LineWidth',2) % P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,3)
    %legend('A','B','C','D','Location','east')
    axis([0 nbDayTe+1 0 1])
    %xticks([1 nbDayTe])
    yticks([0:0.2:1])
    if (~ofclesion&&condAB)
        xlabel('Day','FontSize',14)
        title('Test','FontSize',18)
    end
end


%% number of magazine entries for Saline condition
% Stage 1
if (condAB) % no need to replot the same beginning for stim C and D
    subplot(2,5,6:7)
    if (ofclesion||~condAB)
        hold on
    end
    plot(1:nbDayS1,mean(totalY(1:nbDayS1,:)')','-sk','LineWidth',2,'MarkerFaceColor',markerColor) % Y4(1:nbDayS1)
    hold on
    plot(121:122,[-2 -2],'-sk','LineWidth',2,'MarkerFaceColor',[0 0 0])
    plot(nbDayS1+2:nbDayS1+nbDayIn+1,mean(totalY(nbDayS1+1:nbDayS1+nbDayIn,:)')','-sk','LineWidth',2,'MarkerFaceColor',markerColor) % Y4(nbDayS1+1:nbDayS1+nbDayIn)
    legend('Saline: A','Muscimol: A','Location','southeast')
    axis([0 nbDayS1+nbDayIn+2 0 1])
    xticks([1:nbDayS1 nbDayS1+2:nbDayS1+nbDayIn+1])
    xticklabels([1:nbDayS1 nbDayS1+1:nbDayS1+nbDayIn])
    yticks([0:0.2:1])
    xlabel('Day','FontSize',14)
    ylabel('Probability of Magazine Entry','FontSize',14)
end
% Stage 2
subplot(2,5,8:9)
if (ofclesion||~condAB)
    hold on
end
plot(121:122,[-2 -2],'-ok','LineWidth',2,'MarkerFaceColor',markerColor)
hold on
plot(121:122,[-2 -2],'--^k','LineWidth',2,'MarkerFaceColor',markerColor)
plot(121:122,[-2 -2],'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0])
plot(121:122,[-2 -2],'--^k','LineWidth',2,'MarkerFaceColor',[0 0 0])
plot(offset+1:offset+nbDayS2,mean(totalY(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,:)')',[stylemarque 'k'],'LineWidth',2,'MarkerFaceColor',markerColor) % Y4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2)
legend('Saline: AB','Saline: CD','Muscimol: AB','Muscimol: CD','Location','southeast')
axis([0 8 0 1])
xticks([1:3 5:7])
xticklabels([12:14 12:14])
yticks([0:0.2:1])
xlabel('Day','FontSize',14)
% Test
subplot(2,5,10)
if (ofclesion||~condAB)
    hold on
end
if (condAB)
    bar(ofclesion*3+1,mean(totalY(nbDays,:)')','LineWidth',2,'EdgeColor',condColor,'FaceColor',condColor) % Y4(nbDays)
else
    bar(ofclesion*3+2,mean(totalY(nbDays,:)')','LineWidth',2,'EdgeColor',condColor,'FaceColor','w') % Y4(nbDays)
end
hold on
bar(1000,-2,'LineWidth',2,'EdgeColor',condColor,'FaceColor','w')
legend('B','D','Location','northwest')
axis([0 6 0 1])
xticks([1.5 4.5])
xticklabels({'Saline','Muscimol'})
yticks([0:0.2:1])

FigHandle = gcf;
set(FigHandle, 'Position', [100, 10, 700, 500]);
    