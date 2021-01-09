% Experiment 4 - Marios Panayi, Mehdi Khamassi, Simon Killcross (2021) 
% lOFC muscimol inactivation during Pavlovian Blocking Task

%% MODEL PARAMETERS
ofclesion = false;
condAB = true; % if false, stim C and D are presented in Stage 2 and Test
nbState = 4; % A or C context x reward delivered/undelivered
nbAction = 2; % magazine entry vs. do nothing
nbStim = 4; % total number of stim used in the task (fixed to 4)
init = 0; % initial stimulus values (V_0)
lambda = 2; % inverse temperature for MB-MF arbitration
alpha = 0.01; % MF learning rate
kappa = 0.9; % MF forgetting rate
Qinit = ones(nbStim,nbAction) * init; % initial Q-values
gamma = 0; % discount factor (time horizon for reward prediction)
beta = 10; % exploration rate (inverse temperature)
deltaMax = 1; % max reliability
decisionRule = 'softmax'; % decision-rule
infusionLevel = 0; % increases step-by-step during infusion

%% TASK PARAMETERS
nbDayS1 = 4; % 1-4 Stage 1
nbDayIn = 6; % 5-10 Stage 1+infusion
nbDayPE = 1; % 11 pre-exposure
nbDayS2 = 3; % 12-14 Stage 2
nbDayTe = 1; % 15 Test
nbDays = nbDayS1 + nbDayIn + nbDayPE + nbDayS2 + nbDayTe;
nbTrialSta1 = 16; % nb trials (all with stim A) per day for Stage 1
nbTrialPreE = 4; % nb trials per stim during pre exposure (B or D) on day 11
nbTrialSta2 = 8; % nb trials for each compound (AB or CD) per day for Stage 2
nbTrialTest = 8; % nb trials per blocked cue (B or D) per day for Test
nbTrial = nbTrialSta1 * (nbDayS1 + nbDayIn) + nbTrialPreE * nbDayPE + nbTrialSta2 * nbDayS2 + nbTrialTest * nbDayTe;

%% INIT MF-RL MODEL
QMF = Qinit; % MF-values for each action
probaMF = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not
memQMF = QMF; % memory of QMF so as to monitore changes in Q-values

%% INIT MB-RL MODEL (transition and reward functions)
hatP = ones(nbState, nbAction, nbState) / nbState;
hatR = ones(nbState, nbAction) * 0; % / nbAction;
N    = ones(nbState, nbAction) * 0;
QMB = Qinit; % MB-values for each action
memQMB = QMB; % memory of QMB so as to monitore changes in Q-values
probaMB = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not

%% INIT RANDOM EXPLORATION MODEL
probaEXP = ones(1,nbAction) / nbAction;

%% INIT META-CONTROLLER (FOR MB-MF ARBITRATION)
probaMBMF = ones(nbTrial+1,3) / 3; % proba to rely on MB or MF or EXP system
choiceMBMF = zeros(nbTrial,1); % chosen system
Y = zeros(nbTrial,1); % chosen action
proba = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not
reliability = zeros(nbTrial+1,3) * deltaMax; % delta^2 of each system
    
%% IF WE WANT TO PLOT SOME FIGURES
plotFigure = true;

%% DATA INPUT STRUCTURE
% 1 - day
% 2 - trial
% 3 - nb magazine entries
% 4:7 - stim (format : 1010)
% 8 - reward
% 9 - model predicted value
% 10 - model reward prediction error
% 11-14 - model stim values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT DATA FIRST VERSION
DATA = [];
% data for Stage 1 + Pre-Exposure with B and D + Stage 2 with AB and CD +
% Test with B and D
for nd=1:nbDays % loop over days
    switch (nd)
        case {1,2,3,4,5,6,7,8,9,10}
            nbt = nbTrialSta1; % Stage 1
            Vmatrix = [ones(nbt,1) zeros(nbt,nbStim-1)]; % A
            reward = 1;
        case 11
            nbt = nbTrialPreE; % Pre-Exposure
            if (condAB)
                Vmatrix = [zeros(nbTrialPreE,1) ones(nbTrialPreE,1) zeros(nbTrialPreE,nbStim-2)]; % B
            else
                Vmatrix = [zeros(nbTrialPreE,nbStim-1) ones(nbTrialPreE,1)]; % D
            end
            reward = 0;
        case {12,13,14}
            nbt = nbTrialSta2;  % Stage 2
            if (condAB)
                Vmatrix = [ones(nbTrialSta2,2) zeros(nbTrialSta2,nbStim-2)]; % AB
            else
                Vmatrix = [zeros(nbTrialSta2,nbStim-2) ones(nbTrialSta2,2)]; % CD
            end
            reward = 1;
        otherwise
            nbt = nbTrialTest; % Test
            if (condAB)
                Vmatrix = [zeros(nbTrialTest,1) ones(nbTrialTest,1) zeros(nbTrialTest,nbStim-2)]; % B
            else
                Vmatrix = [zeros(nbTrialTest,nbStim-1) ones(nbTrialTest,1)]; % D
            end
            reward = 0;
    end
    DATA = [DATA ; zeros(nbt, 8)];
    DATA(end-nbt+1:end,1) = nd;
    DATA(end-nbt+1:end,2) = (1:nbt)';
    DATA(end-nbt+1:end,3) = round(DATA(end-nbt+1:end,8) * 10);
    DATA(end-nbt+1:end,4:7) = Vmatrix;
    DATA(end-nbt+1:end,8) = reward;
end

%% INIT EXPERIMENT
day = 1;
for nt=1:nbTrial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% trial initialization
    % Observe presented stimuli
    stim = DATA(nt,4:4+nbStim-1); % presented stimuli
    listStim = (1:nbStim) .* stim; % list of stim that are present
    listStim(listStim==0) = []; % removing zeros from the list
    % detect if main stim A present or not
    stimA = stim(1);
    % determine current state for MB
    if (stimA)
        s = 1; % reward not yet delivered
    else % stim A not present
        s = 3; % reward not yet delivered
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MB-MF-EXP arbitration (meta-control)
    % reliability-based coordination
    [choiceMBMF(nt), probaMBMF(nt+1,:)] = valueBasedDecision([reliability(nt,:)], decisionRule, lambda, 0);
    if (ofclesion&&(day>=5)&&(day<=10)) % muscimol injection
        infusionLevel = infusionLevel + 1;
        theBest = argmax(probaMBMF(nt+1,:));
        % 1ST SOLUTION: slow decrease of best and increase of others
        probaMBMF(nt+1,theBest) = max(1/3,probaMBMF(nt+1,theBest) - infusionLevel/100); % proba of choosing MB system
        probaMBMF(nt+1,mod(theBest+0,3)+1) = min(1/3,probaMBMF(nt+1,theBest) + infusionLevel/200); % proba of choosing MF system
        probaMBMF(nt+1,mod(theBest+1,3)+1) = min(1/3,probaMBMF(nt+1,theBest) + infusionLevel/200); % proba of choosing EXP system
        probaMBMF(nt+1,:) = probaMBMF(nt+1,:) ./ sum(probaMBMF(nt+1,:)); % normalization
        choiceMBMF(nt) = drand01(probaMBMF(nt+1,:)); % rolls a dice and chooses a system depending on its proba
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% decision-making
    % Ask the chosen system (MB or MF) which action to perform
    switch (choiceMBMF(nt))
        case 1 % MB
            % the MB system is state-based
            [Y(nt), probaMB(nt+1,:)] = valueBasedDecision(QMB(s,:), decisionRule, beta, 0);
            probaMF(nt+1,:) = probaMF(nt,:);
            proba(nt+1,:) = probaMB(nt+1,:);
            predictedV = QMB(s,:);
        case 2 % MF
            % The proba of the MF system to visit the magazine for a certain
            % amount of time depends on the sum over presented stim values
            [Y(nt), probaMF(nt+1,:)] = valueBasedDecision([sum(QMF(listStim,1)) sum(QMF(listStim,2))], decisionRule, beta, 0);
            probaMB(nt+1,:) = probaMB(nt,:);
            proba(nt+1,:) = probaMF(nt+1,:);
            predictedV = [sum(QMF(listStim,1)) sum(QMF(listStim,2))];
        case 3 % EXP
            proba(nt+1,:) = [sum(probaEXP(1,1)) sum(probaEXP(1,2))];
            proba(nt+1,:) = proba(nt+1,:) / sum(proba(nt+1,:));
            Y(nt) = drand01(proba(nt+1,:)); % rolls a dice and chooses an action depending on its proba
            probaMF(nt+1,:) = probaMF(nt,:);
            probaMB(nt+1,:) = probaMB(nt,:);
            predictedV = [0 0];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% observing consequence in the environment
    if (Y(nt) == 2)&&(DATA(nt,4+nbStim)) % if entered magazine and (Stage 1 or Stage 2)
        r = 1; % reward delivered
        % determine new state
        if (stimA)
            y = 2;
        else % not the context where stim A is present
            y = 4;
        end
    else % if (pre-exposure of B/D or test phase)
        % reward undelivered
        if (Y(nt)==2) % if action==enter magazine
            r = -1; % energy cost of movement
        else % no move
            r = 0;
        end
        % determine new state
        if (stimA)
            y = 1;
        else % not the context where stim A is present
            y = 3;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MB learning
    % computer MB reliability based on MB prediction error
    delta = temporalDifferenceError( r, 0, 1, Y(nt), QMB(s,:), QMB(y,:), 1, alpha, gamma, 1, 0, 0, 0 );
    reliability(nt+1,1) = reliability(nt,1) + alpha * (1 - delta(Y(nt))^2 / deltaMax^2 - reliability(nt,1)); % monitoring reliability
    % Update the number of visits for the current state-action pair
    N(s, Y(nt)) = N(s, Y(nt)) + 1;
    % Update transition matrix (stochastic)
    hatP(s, Y(nt), :) = (1 - 1/N(s, Y(nt))) * hatP(s, Y(nt), :) + reshape((1:nbState == y) / N(s, Y(nt)), 1, 1, nbState);
    % Update reward function (deterministic)
    hatR(s, Y(nt)) = r; %(1 - 1/N(s, Y(nt))) * hatR(s, Y(nt)) + r / N(s, Y(nt));
    if (choiceMBMF(nt) == 1) % MB system has been chosen
        % good old simple value iteration method to update QMB
        Qmax = max(QMB, [], 2);
        for sss = 1:nbState
            for aaa = 1:nbAction
                QMB(sss, aaa) = hatR(sss, aaa) + gamma * sum(reshape(hatP(sss, aaa, :), nbState, 1) .* Qmax);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MF learning
    % Update Q in a model-free manner
    delta = temporalDifferenceError( r, 0, 1, Y(nt), [sum(QMF(listStim,1)) sum(QMF(listStim,2))], 0, 1, alpha, gamma, 1, 0, 0, 0 );
    reliability(nt+1,2) = reliability(nt,2) + alpha * (1 - delta(Y(nt))^2 / deltaMax^2 - reliability(nt,2)); % monitoring reliability
    % updating QMF
    for aaa = 1:nbAction
        QMF(:,aaa) = QMF(:,aaa) + alpha * delta(aaa) * stim';
    end
    % forgetting over night
    if (DATA(nt,1) ~= day)
        day = DATA(nt,1);
        for sss = 1:nbState
            for aaa = 1:nbAction
                % if (ofclesion&&(day>=5)&&(day<=10)) % muscimol injection
                if (probaMBMF(nt+1,1)<0.5) % no forgetting if MB inference can compensate
                    QMF(sss, aaa) = QMF(sss, aaa) + kappa * (Qinit(sss, aaa) - QMF(sss, aaa));
                    % forgetting makes Q-values converge back to memorized Qinit
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXP
    % Update EXP reliability
    reliability(nt+1,3) = reliability(nt,3) + alpha * (1 - deltaMax^2 / deltaMax^2 - reliability(nt,3)); % monitoring reliability
    
    % The next state becomes the current state
    s = y;

    %% LOGS
    DATA(nt,4+nbStim+1) = predictedV(Y(nt));
    DATA(nt,4+nbStim+2) = delta(Y(nt));
    DATA(nt,4+nbStim+2+1:4+nbStim+2+nbStim) = sum(QMF'); % value of each stim
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA POST-PROCESSING
% organizing data to plot (V values) and compute distance to
% experimental data curves
Y4 = zeros(nbDays,1);
P4 = zeros(nbDays,3);
for nd=1:nbDays
    switch nd
        case {1,2,3,4,5,6,7,8,9,10} % Stage 1
            Y4(nd) = mean(Y((nd-1)*nbTrialSta1+1:nd*nbTrialSta1)-1);
            P4(nd,1) = mean(probaMBMF((nd-1)*nbTrialSta1+1:nd*nbTrialSta1,1));
            P4(nd,2) = mean(probaMBMF((nd-1)*nbTrialSta1+1:nd*nbTrialSta1,2));
            P4(nd,3) = mean(probaMBMF((nd-1)*nbTrialSta1+1:nd*nbTrialSta1,3));
        case 11 % Pre-Exposure
            Y4(nd) = mean(Y((nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn-1)*nbTrialPreE+1:(nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn)*nbTrialPreE)-1);
            P4(nd,1) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn-1)*nbTrialPreE+1:(nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn)*nbTrialPreE,1));
            P4(nd,2) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn-1)*nbTrialPreE+1:(nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn)*nbTrialPreE,2));
            P4(nd,3) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn-1)*nbTrialPreE+1:(nbDayS1+nbDayIn)*nbTrialSta1+(nd-nbDayS1-nbDayIn)*nbTrialPreE,3));
        case {12,13,14} % Stage 2
            Y4(nd) = mean(Y((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE-1)*nbTrialSta2+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE)*nbTrialSta2)-1);
            P4(nd,1) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE-1)*nbTrialSta2+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE)*nbTrialSta2,1));
            P4(nd,2) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE-1)*nbTrialSta2+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE)*nbTrialSta2,2));
            P4(nd,3) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE-1)*nbTrialSta2+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+(nd-nbDayS1-nbDayIn-nbDayPE)*nbTrialSta2,3));
        otherwise % Test
            Y4(nd) = mean(Y((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2-1)*nbTrialTest+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2)*nbTrialTest)-1);
            P4(nd,1) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2-1)*nbTrialTest+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2)*nbTrialTest,1));
            P4(nd,2) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2-1)*nbTrialTest+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2)*nbTrialTest,2));
            P4(nd,3) = mean(probaMBMF((nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2-1)*nbTrialTest+1:(nbDayS1+nbDayIn)*nbTrialSta1+nbDayPE*nbTrialPreE+nbDayS2*nbTrialSta2+(nd-nbDayS1-nbDayIn-nbDayPE-nbDayS2)*nbTrialTest,3));
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURES
% if (plotFigure)
%     
%     % figure parameters
%     if (~ofclesion)
%         if (condAB)
%             figure
%         end
%         markerColor = 'w';
%         condColor = 'k';
%         offset = 0;
%         styletrait = '--';
%     else
%         markerColor = 'k';
%         condColor = 'r';
%         offset = 4;
%         styletrait = '';
%     end
%     if (condAB)
%         stylemarque = '-o';
%     else
%         stylemarque = '--^';
%     end
% 
%     if (condAB) % no need to replot the same beginning for stim C and D
%         % Stage 1
%         subplot(2,5,1:2)
%         if (ofclesion||~condAB)
%             hold on
%         end
%         plot(1:nbDayS1,P4(1:nbDayS1,1),[styletrait 'b'],'LineWidth',2)
%         hold on
%         plot(1:nbDayS1,P4(1:nbDayS1,2),[styletrait 'r'],'LineWidth',2)
%         plot(1:nbDayS1,P4(1:nbDayS1,3),[styletrait 'k'],'LineWidth',2)
%         plot(nbDayS1+2:nbDayS1+nbDayIn+1,P4(nbDayS1+1:nbDayS1+nbDayIn,1),[styletrait 'b'],'LineWidth',2)
%         plot(nbDayS1+2:nbDayS1+nbDayIn+1,P4(nbDayS1+1:nbDayS1+nbDayIn,2),[styletrait 'r'],'LineWidth',2)
%         plot(nbDayS1+2:nbDayS1+nbDayIn+1,P4(nbDayS1+1:nbDayS1+nbDayIn,3),[styletrait 'k'],'LineWidth',2)
% 
%         legend('MB','MF','EXP','Location','northwest')
%         axis([0 nbDayS1+nbDayIn+2 0 1])
%         xticks([1:nbDayS1 nbDayS1+2:nbDayS1+nbDayIn+1])
%         xticklabels([1:nbDayS1 nbDayS1+1:nbDayS1+nbDayIn])
%         yticks([0:0.2:1])
%         xlabel('Day','FontSize',14)
%         ylabel('Probability of Selection','FontSize',14)
%         if (~ofclesion&&condAB)
%             title('Stage 1','FontSize',18)
%         end
%         % Stage 2
%         subplot(2,5,3:4)
%         if (ofclesion||~condAB)
%             hold on
%         end
%         plot(offset+1:offset+nbDayS2,P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,1),[styletrait 'b'],'LineWidth',2)
%         hold on
%         plot(offset+1:offset+nbDayS2,P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,2),[styletrait 'r'],'LineWidth',2)
%         plot(offset+1:offset+nbDayS2,P4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2,3),[styletrait 'k'],'LineWidth',2)
%         axis([0 8 0 1])
%         xticks([1:3 5:7])
%         xticklabels([12:14 12:14])
%         yticks([0:0.2:1])
%         xlabel('Day','FontSize',14)
%         if (~ofclesion&&condAB)
%             title('Stage 2','FontSize',18)
%         end
%         % Test
%         subplot(2,5,5)
%         if (ofclesion||~condAB)
%             hold on
%         end
%         plot(1,P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,1),'ob','LineWidth',2)
%         hold on
%         plot(1,P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,2),'or','LineWidth',2)
%         plot(1,P4(nbDayS1+nbDayIn+nbDayPE+nbDayS2+nbDayTe,3),'ok','LineWidth',2)
%         %legend('A','B','C','D','Location','east')
%         axis([0 nbDayTe+1 0 1])
%         %xticks([1 nbDayTe])
%         yticks([0:0.2:1])
%         if (~ofclesion&&condAB)
%             xlabel('Day','FontSize',14)
%             title('Test','FontSize',18)
%         end
%     end
% 
% 
%     %% number of magazine entries for Saline condition
%     % Stage 1
%     if (condAB) % no need to replot the same beginning for stim C and D
%         subplot(2,5,6:7)
%         if (ofclesion||~condAB)
%             hold on
%         end
%         plot(1:nbDayS1,Y4(1:nbDayS1),'-sk','LineWidth',2,'MarkerFaceColor',markerColor)
%         hold on
%         plot(121:122,[-2 -2],'-sk','LineWidth',2,'MarkerFaceColor',[0 0 0])
%         plot(nbDayS1+2:nbDayS1+nbDayIn+1,Y4(nbDayS1+1:nbDayS1+nbDayIn),'-sk','LineWidth',2,'MarkerFaceColor',markerColor)
%         legend('Saline: A','Muscimol: A','Location','southeast')
%         axis([0 nbDayS1+nbDayIn+2 0 1])
%         xticks([1:nbDayS1 nbDayS1+2:nbDayS1+nbDayIn+1])
%         xticklabels([1:nbDayS1 nbDayS1+1:nbDayS1+nbDayIn])
%         yticks([0:0.2:1])
%         xlabel('Day','FontSize',14)
%         ylabel('Probability of Magazine Entry','FontSize',14)
%     end
%     % Stage 2
%     subplot(2,5,8:9)
%     if (ofclesion||~condAB)
%         hold on
%     end
%     plot(121:122,[-2 -2],'-ok','LineWidth',2,'MarkerFaceColor',markerColor)
%     hold on
%     plot(121:122,[-2 -2],'--^k','LineWidth',2,'MarkerFaceColor',markerColor)
%     plot(121:122,[-2 -2],'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0])
%     plot(121:122,[-2 -2],'--^k','LineWidth',2,'MarkerFaceColor',[0 0 0])
%     plot(offset+1:offset+nbDayS2,Y4(nbDayS1+nbDayIn+nbDayPE+1:nbDayS1+nbDayIn+nbDayPE+nbDayS2),[stylemarque 'k'],'LineWidth',2,'MarkerFaceColor',markerColor)
%     legend('Saline: AB','Saline: CD','Muscimol: AB','Muscimol: CD','Location','southeast')
%     axis([0 8 0 1])
%     xticks([1:3 5:7])
%     xticklabels([12:14 12:14])
%     yticks([0:0.2:1])
%     xlabel('Day','FontSize',14)
%     % Test
%     subplot(2,5,10)
%     if (ofclesion||~condAB)
%         hold on
%     end
%     if (condAB)
%         bar(ofclesion*3+1,Y4(nbDays),'LineWidth',2,'EdgeColor',condColor,'FaceColor',condColor)
%     else
%         bar(ofclesion*3+2,Y4(nbDays),'LineWidth',2,'EdgeColor',condColor,'FaceColor','w')
%     end
%     hold on
%     bar(1000,-2,'LineWidth',2,'EdgeColor',condColor,'FaceColor','w')
%     legend('B','D','Location','northwest')
%     axis([0 6 0 1])
%     xticks([1.5 4.5])
%     xticklabels({'Saline','Muscimol'})
%     yticks([0:0.2:1])
% 
%     FigHandle = gcf;
%     set(FigHandle, 'Position', [100, 10, 700, 500]);
% end % end of if (plotFigure)
