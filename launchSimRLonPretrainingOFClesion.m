%% launch simRLonPavlovianBlockingTask (with OFC inactivation)
nbExperiments = 1000;

% LAUNCH X experiments
totalprob = [];
totalY = [];
for ttt=1:nbExperiments
    % experiment starts
    simRLonPretrainingOFClesion
    % experiment finished
    % record logs:
    totalprob = [totalprob probaMBMF(1:end-1,:)];
    totalY = [totalY Y];
end

% PLOT FIGURE
figure
subplot(2,1,1)
plot(mean(totalprob(:,1:3:end-2)')','LineWidth',2)
hold on
plot(mean(totalprob(:,2:3:end-1)')','r','LineWidth',2)
plot(mean(totalprob(:,3:3:end)')','k','LineWidth',2)
ylabel('Probability of Selection')
xlabel('Trial Number')
legend('MB','MF','EXP','Location','east')
subplot(2,1,2)
plot(mean(totalY')'-1,'LineWidth',2)
ylabel('Probability of Magazine Entry')
xlabel('Trial Number')
legend('MB-alone','MF-alone','Sham','Lesion','Location','southeast')
