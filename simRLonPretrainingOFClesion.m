% Experiment 1 - Marios Panayi, Mehdi Khamassi, Simon Killcross (2021)
% pre-training lOFC lesion

%% PARAMETERS
ofclesion = true;
nbSession = 21; % 21 days
nbTrial = nbSession * 16; % 16 trials per days (16 trials a day, 15s CS, 90s variable Inter-Trial-Interval)
nbState = 2; % reward undelivered vs. reward delivered
nbAction = 2; % magazine entry vs. do nothing
lambda = 2; % inverse temperature for MB-MF arbitration
alpha = 0.01; % MF learning rate
kappa = 0.1; % MF forgetting rate
Qinit = zeros(nbState,nbAction); % initial Q-values
gamma = 0; % discount factor (time horizon for reward prediction)
beta = 10; % exploration rate (inverse temperature)
deltaMax = 1; % max reliability
decisionRule = 'softmax'; % decision-rule

%% INIT MF-RL MODEL
QMF = Qinit; % MF-values for each action
probaMF = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not
memQMF = zeros(nbTrial+1,nbAction); % proba to visit magazine or not
memQMF(1,:) = QMF(1,:); % storing the initial Q-values in decision state
HMF = ones(nbTrial+1,1); % uncertainty in MF system

%% INIT MB-RL MODEL (transition and reward functions)
hatP = ones(nbState, nbAction, nbState) / nbState;
hatR = ones(nbState, nbAction) * 0; % / nbAction;
N    = ones(nbState, nbAction) * 10;
QMB = ones(nbState,nbAction) * 1; %Qinit; % MB-values for each action
probaMB = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not
HMB = ones(nbTrial+1,1); % uncertainty in MB system

%% INIT RANDOM EXPLORATION MODEL
probaEXP = ones(nbState,nbAction) / nbAction;
HEXP = zeros(nbState,1);
for sss=1:nbState
    HEXP(sss) = entropyProba(probaEXP(sss,:));
end

%% INIT META-CONTROLLER (FOR MB-MF ARBITRATION)
probaMBMF = ones(nbTrial+1,3) / 3; % proba to rely on MB or MF or EXP system
choiceMBMF = zeros(nbTrial,1); % chosen system
Y = zeros(nbTrial,1); % chosen action
proba = ones(nbTrial+1,nbAction) / nbAction; % proba to visit magazine or not
reliability = ones(nbTrial+1,3) * deltaMax; % delta^2 of each system

%% Run the experiment
for nt=1:nbTrial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% trial initialization
    s = 1; % reward not yet delivered
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MB-MF-EXP arbitration (meta-control)
    % measuring uncertainty in the MF system
    [~, probaMF(nt+1,:)] = valueBasedDecision(QMF(s,:), decisionRule, 1, 0);
    HMF(nt+1) = entropyProba(probaMF(nt+1,:));
    % measuring uncertainty in the MB system
    [~, probaMB(nt+1,:)] = valueBasedDecision(QMB(s,:), decisionRule, 1, 0);
    HMB(nt+1) = entropyProba(probaMB(nt+1,:));
    % choose which system (MB or MF) decides
    if (ofclesion) % MF wins because it is never inhibited
        probaMBMF(nt+1,2) = 1; % proba of choosing MF system
        probaMBMF(nt+1,1) = (1 - probaMBMF(nt+1,2)) / 2; % proba of choosing MB system
        probaMBMF(nt+1,3) = (1 - probaMBMF(nt+1,2)) / 2; % proba of choosing EXP system
        choiceMBMF(nt) = drand01(probaMBMF(nt+1,:)); % rolls a dice and chooses a system depending on its proba
    else % intact OFC (thus intact arbitration)
        % entropy-based coordination
        %[choiceMBMF(nt), probaMBMF(nt+1,:)] = valueBasedDecision(1-[HMB(nt+1) HMF(nt+1) HEXP(s)], decisionRule, lambda, 0);
        % reliability-based coordination
        [choiceMBMF(nt), probaMBMF(nt+1,:)] = valueBasedDecision([reliability(nt,:)], decisionRule, lambda, 0);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% decision-making
    % Ask the chosen system (MB or MF) which action to perform
    switch (choiceMBMF(nt))
        case 1 % MB
            [Y(nt), probaMB(nt+1,:)] = valueBasedDecision(QMB(s,:), decisionRule, beta, 0);
            probaMF(nt+1,:) = probaMF(nt,:);
            proba(nt+1,:) = probaMB(nt+1,:);
        case 2 % MF
            [Y(nt), probaMF(nt+1,:)] = valueBasedDecision(QMF(s,:), decisionRule, beta, 0);
            probaMB(nt+1,:) = probaMB(nt,:);
            proba(nt+1,:) = probaMF(nt+1,:);
        case 3 % EXP
            proba(nt+1,:) = probaEXP(s,:);
            Y(nt) = drand01(proba(nt+1,:)); % rolls a dice and chooses an action depending on its proba
            probaMF(nt+1,:) = probaMF(nt,:);
            probaMB(nt+1,:) = probaMB(nt,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% observing consequence in the environment
    if (Y(nt) == 2) % reward
        y = 2;
        r = 1;
    else % no reward
        y = 1;
        r = 0;
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
    hatR(s, Y(nt)) = (1 - 1/N(s, Y(nt))) * hatR(s, Y(nt)) + r / N(s, Y(nt));
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
    % Update Q in a model-free manner and update MF reliability
    [ delta, QMF(s,:) ] = temporalDifferenceError( r, 0, 1, Y(nt), QMF(s,:), QMF(y,:), 1, alpha, gamma, 1, 0, 0, 0 );
    %reliability(nt+1,2) = reliability(nt,2) + alpha * (1 - delta^2 / deltaMax^2 - reliability(nt,2)); % monitoring reliability
    reliability(nt+1,2) = reliability(nt,2) + alpha * (1 - delta(Y(nt))^2 / deltaMax^2 - reliability(nt,2)); % monitoring reliability
    % forgetting
    if (mod(nt,16)==0)
        for sss = 1:nbState
            for aaa = 1:nbAction
                if (probaMBMF(nt+1,1)<0.5) % no forgetting if MB inference can compensate
                    QMF(sss, aaa) = QMF(sss, aaa) + kappa * (Qinit(sss, aaa) - QMF(sss, aaa));
                    % forgetting makes Q-values converge back to memorized Qinit
                end
            end
        end
    end
    memQMF(nt,:) = QMF(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXP
    % Update EXP reliability
    reliability(nt+1,3) = reliability(nt,3) + alpha * (1 - deltaMax^2 / deltaMax^2 - reliability(nt,3)); % monitoring reliability
    
    % The next state becomes the current state
    s = y;
end

%% PLOT FIGURES

% % trial figure
% figure
% subplot(3,1,1)
% plot(probaMBMF,'LineWidth',2)
% ylabel('proba(MB)')
% xlabel('trials')
% legend('MB','MF')
% subplot(3,1,2)
% plot(Y-1,'LineWidth',2)
% ylabel('magazine entries')
% xlabel('trials')
% subplot(3,1,3)
% plot(proba(:,2),'LineWidth',2)
% ylabel('proba magazine entry')
% xlabel('blocks')
% 
% 
% 
% % block figure
% figure
% subplot(3,1,1)
% smoothprobaMBMF = [];
% for iii=1:nbSession
%     smoothprobaMBMF = [smoothprobaMBMF ; mean(probaMBMF((iii-1)*(nbTrial/nbSession)+1:iii*(nbTrial/nbSession),:))];
% end
% plot(smoothprobaMBMF,'LineWidth',2)
% ylabel('proba(MB)')
% xlabel('blocks')
% legend('MB','MF')
% subplot(3,1,2)
% smoothY = [];
% for iii=1:nbSession
%     smoothY = [smoothY mean(Y((iii-1)*(nbTrial/nbSession)+1:iii*(nbTrial/nbSession),:))];
% end
% plot(smoothY-1,'LineWidth',2)
% ylabel('magazine entries')
% xlabel('blocks')
% subplot(3,1,3)
% smoothproba = [];
% for iii=1:nbSession
%     smoothproba = [smoothproba mean(proba((iii-1)*(nbTrial/nbSession)+1:iii*(nbTrial/nbSession),2))];
% end
% plot(smoothproba,'LineWidth',2)
% ylabel('proba magazine entry')
% xlabel('blocks')
% 
%