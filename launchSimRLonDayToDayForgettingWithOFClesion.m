%% launch simRLonDayToDayForgettingWithOFClesion
nbExperiments = 1000;
% LAUNCH X experiments
totalprob = [];
totalY = [];
for ttt=1:nbExperiments,
    % experiment starts
    simRLonDayToDayForgettingWithOFClesion
    % experiment finished
    % record logs:
    totalprob = [totalprob probaMBMF(1:end-1,:)];
    totalY = [totalY Y];
end

% PLOT FIGURE
figure
subplot(2,1,1)
totalpmb = totalprob(:,1:3:end-2);
totalpmf = totalprob(:,2:3:end-1);
totalpex = totalprob(:,3:3:end);
smoothpmb = [];
smoothpmf = [];
smoothpex = [];
for iii=1:nbTrial/4
    smoothpmb = [smoothpmb mean(mean(totalpmb((iii-1)*4+1:iii*4,:)))];
    smoothpmf = [smoothpmf mean(mean(totalpmf((iii-1)*4+1:iii*4,:)))];
    smoothpex = [smoothpex mean(mean(totalpex((iii-1)*4+1:iii*4,:)))];
end
plot(1,smoothpmb(36),'ob','LineWidth',2) % acquisition
hold on
plot(1,smoothpmf(36),'or','LineWidth',2) % acquisition
plot(1,smoothpex(36),'ok','LineWidth',2) % acquisition
plot(3:6,smoothpmb(37:40),'b','LineWidth',2) % day 1
plot(3:6,smoothpmf(37:40),'r','LineWidth',2) % day 1
plot(3:6,smoothpex(37:40),'k','LineWidth',2) % day 1
plot(8:11,smoothpmb(41:44),'b','LineWidth',2) % day 2
plot(8:11,smoothpmf(41:44),'r','LineWidth',2) % day 2
plot(8:11,smoothpex(41:44),'k','LineWidth',2) % day 2
plot(13:16,smoothpmb(45:48),'b','LineWidth',2) % day 3
plot(13:16,smoothpmf(45:48),'r','LineWidth',2) % day 3
plot(13:16,smoothpex(45:48),'k','LineWidth',2) % day 3
plot(18:21,smoothpmb(49:52),'b','LineWidth',2) % day 4
plot(18:21,smoothpmf(49:52),'r','LineWidth',2) % day 4
plot(18:21,smoothpex(49:52),'k','LineWidth',2) % day 4
plot(23:26,smoothpmb(53:56),'b','LineWidth',2) % day 5
plot(23:26,smoothpmf(53:56),'r','LineWidth',2) % day 5
plot(23:26,smoothpex(53:56),'k','LineWidth',2) % day 5
plot(28:31,smoothpmb(57:60),'b','LineWidth',2) % day 6
plot(28:31,smoothpmf(57:60),'r','LineWidth',2) % day 6
plot(28:31,smoothpex(57:60),'k','LineWidth',2) % day 6
ylabel('Probability of Selection')
xlabel('Blocks of 4 trials')
legend('MB','MF','EXP','Location','east')
xticks([1 3:6 8:11 13:16 18:21 23:26 28:31])
xticklabels({'Acq','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4'})

subplot(2,1,2)
smoothY = [];
for iii=1:nbTrial/4
    smoothY = [smoothY mean(mean(totalY((iii-1)*4+1:iii*4,:)))];
end
plot(100,smoothY(36)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % acquisition
hold on
plot(103:106,smoothY(37:40)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 1
plot(1,smoothY(36)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % acquisition
plot(3:6,smoothY(37:40)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 1
plot(8:11,smoothY(41:44)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 2
plot(13:16,smoothY(45:48)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 3
plot(18:21,smoothY(49:52)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 4
plot(23:26,smoothY(53:56)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 5
plot(28:31,smoothY(57:60)-1,'-ok','LineWidth',2,'MarkerFaceColor','w','MarkerSize',8) % day 6
ylabel('Probability of Magazine Entry')
xlabel('Blocks of 4 trials')
xticks([1 3:6 8:11 13:16 18:21 23:26 28:31])
xticklabels({'Acq','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4','1','2','3','4'})
legend('Saline','Muscimol','Location','northeast')
axis([0 35 0 1])


plot(1,smoothY(36)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % acquisition
hold on
plot(3:6,smoothY(37:40)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 1
plot(8:11,smoothY(41:44)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 2
plot(13:16,smoothY(45:48)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 3
plot(18:21,smoothY(49:52)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 4
plot(23:26,smoothY(53:56)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 5
plot(28:31,smoothY(57:60)-1,'-sk','LineWidth',2,'MarkerFaceColor','r','MarkerSize',8) % day 6
