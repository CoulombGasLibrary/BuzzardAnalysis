% Bootstrap error calculation for Coulomb gasses. A number of realisations
% for an interval of \beta must have been generated with CoulombGasSim.m
% and the distances calculated with CoulombGasLoad.m. These are compared to
% each other to get a sense of the uncertainty in the fitting algorithm.
% This is done with moving averages, but can be adapted to independent
% samples.
% To save time, the fitting here assumes a convex landscape for \beta. That
% means that we start the fit at the true \beta and investigate the
% neighbour around it. This method should only be used with caution on real
% data, as several local minima may exist there.

clear;
% Coulomb gas to be fitted
N = 200; % Number of particles
Nconf = 1e2; % Number of configurations
Nsteps = 1e2; % Number of iterations
ep = N^(-1/3); % Step size

% \beta interval to be fitted. Note that the highest values will be biased,
% so some margin should be included.
BetMin = 0.0;
BetMax = 0.3;
BetStep = 0.1;
BetVec = BetMin:BetStep:BetMax;
BetN = length(BetVec);

% Prelocate for the fits and the corresponding Kolmogorov distance and save
% them in a file with the names TEXT.
betaFits = NaN * ones(BetN,Nconf);
ksFits = NaN * ones(BetN,Nconf);
betaFitsNext = NaN * ones(BetN,Nconf);
ksFitsNext = NaN * ones(BetN,Nconf);

% Grouping size for the moving average
Group = 5;

TEXT = ['CoulombBetaFit_Group',num2str(Group),'_N',num2str(N),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_BetMin',num2str(BetMin),...
                '_BetMax',num2str(BetMax),'_BetStep',num2str(BetStep)];


% Plot a check of the fits. That is, the average value of the fit for each
% \beta_0 +/- the uncertainty. The average should of course be a straight
% line with slope 1. The bias of the highest values is apparent here. This
% is done for both symmetric and asymmetric errors, and both the NN and NNN
% distances.
PlotCheck = true;


% Load all the Coulomb gas distances first and saving them in the cells
% DistBet and DistBet2.
DistBet = cell(BetN,1); % Nearest Neighbour
DistBet2 = cell(BetN,1); % Next-to-nearest Neighbour
DistBetMat = cell(BetN,1); % Nearest Neighbour
DistBet2Mat = cell(BetN,1); % Next-to-nearest Neighbour

for iBet = 1:BetN
    
    betaNum = BetVec(iBet);
    
    % TEXT for the Coulomb gasses
    TEXTc = ['Coulomb_N',num2str(N),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(betaNum),'_ep',num2str(ep)];
    
    DistBetTemp = importdata(['data/',TEXTc,'_dists.txt']);
    DistBet{iBet} = sort(DistBetTemp(:));
    DistBetMat{iBet} = DistBetTemp;
    
    DistBet2Temp = importdata(['data/',TEXTc,'_NextDists.txt']);
    DistBet2{iBet} = sort(DistBet2Temp(:));
    DistBet2Mat{iBet} = DistBet2Temp;
    
end


%% The main fitting part
for iBet = 1:BetN
    
    % Pick the corresponding \beta
    betaNum = BetVec(iBet);
    % Display the time to keep track of the simulation
    disp(datetime('now'))
    disp(['Looking at \beta=',num2str(betaNum)])
    
    % TEXT for the Coulomb gasses
    TEXTc = ['Coulomb_N',num2str(N),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(betaNum),'_ep',num2str(ep)];
    
    % Load distances
    DistMat = DistBetMat{iBet};
    DistNextMat = DistBet2Mat{iBet};
    
    parfor iTrials = 1:Nconf
        
        % Display the time
        if mod(iTrials,Nconf/100)==0
            disp(datetime('now'))
            disp(['Looking at trial ',num2str(iTrials)]);
        end
        
        % Pick out the group. If independent entries are preferred, this is
        % where it should be changed.
        MinDistGroup = DistMat(:,iTrials:min(iTrials+Group-1,Nconf));
        MinDistNextGroup = DistNextMat(:,iTrials:min(iTrials+Group-1,Nconf));

        % To make the drawing of groups periodic.
        if iTrials+Group-1 > Nconf
            MinDistGroup = [MinDistGroup, DistMat(:,1:iTrials+Group-1-Nconf) ];
            MinDistNextGroup = [MinDistNextGroup, DistNextMat(:,1:iTrials+Group-1-Nconf) ];
        end
        
        
        % Start out with the true beta and make a Kolmogorov comparison.
        [~,~,ksDist] = kstest2( MinDistGroup(:) , DistBet{iBet} ); % NN
        [~,~,ksDist2] = kstest2( MinDistNextGroup(:) , DistBet2{iBet} ); % NNN
        
        % ----------------------- Nearest Neighbour -----------------------
        % Test the \beta below
        BetMoveMinus = iBet-1;
        % If we are at the smallest index, set the Kolmogorov distance to
        % 1, so this direction is not picked.
        if BetMoveMinus >= 1
            [~,~,ksDistMinus] = kstest2( MinDistGroup(:) , DistBet{BetMoveMinus} );
        else
            ksDistMinus = 1;
        end
        
        % Test the \beta above
        BetMovePlus = iBet+1;
        % If we are at the highest index, set the Kolmogorov distance to 1,
        % so this direction is not picked.
        if BetMovePlus <= BetN
            [~,~,ksDistPlus] = kstest2( MinDistGroup(:) , DistBet{BetMovePlus} );
        else
            ksDistPlus = 1;
        end
        
        % Find the lowest Kolmogorow distance and move in that direction
        % until the distance becomes larger again. If the lowest point is
        % at the original point, do nothing.
        if ksDistMinus < ksDist
            DirMove = -1;
            ksDist = ksDistMinus;
            iBetMove = BetMoveMinus;
        elseif ksDistPlus < ksDist
            DirMove = 1;
            ksDist = ksDistPlus;
            iBetMove = BetMovePlus;
        else
            DirMove = 0;
            iBetMove = iBet;
        end
        
        % If we are not at edge of the available \beta and if a direction
        % has been chosen, move in that direction until a larger distance
        % is found.
        if DirMove ~= 0 && iBetMove > 1 && iBetMove < BetN
            [~,~,ksDistTemp] = kstest2( MinDistGroup(:) , DistBet{iBetMove + DirMove} );
            while iBetMove >= 1 && iBetMove <= BetN && ksDistTemp < ksDist
                ksDist = ksDistTemp;
                iBetMove = iBetMove + DirMove;
                if iBetMove > 1 && iBetMove < BetN 
                    [~,~,ksDistTemp] = kstest2( MinDistGroup(:) , DistBet{iBetMove + DirMove} );
                end
            end
        end
        BetaFit = BetVec(iBetMove);
        
        
        % ------------------- Next-to-Nearest Neighbour -------------------
        % Test the \beta below
        BetMoveNextMinus = iBet-1;
        % If we are at the smallest index, set the Kolmogorov distance to
        % 1, so this direction is not picked.
        if BetMoveNextMinus >= 1
            [~,~,ksDistNextMinus] = kstest2( MinDistNextGroup(:) , DistBet2{BetMoveNextMinus} );
        else
            ksDistNextMinus = 1;
        end
        
        % Test the \beta above
        BetMovePlus = iBet+1;
        % If we are at the highest index, set the Kolmogorov distance to 1,
        % so this direction is not picked.
        if BetMovePlus <= BetN
            [~,~,ksDistNextPlus] = kstest2( MinDistNextGroup(:) , DistBet2{BetMovePlus} );
        else
            ksDistNextPlus = 1;
        end
        
        % Find the lowest Kolmogorow distance and move in that direction
        % until the distance becomes larger again. If the lowest point is
        % at the original point, do nothing.
        if ksDistNextMinus < ksDist2
            DirMoveNext = -1;
            ksDist2 = ksDistNextMinus;
            iBetMoveNext = BetMoveNextMinus;
        elseif ksDistNextPlus < ksDist2
            DirMoveNext = 1;
            ksDist2 = ksDistNextPlus;
            iBetMoveNext = BetMovePlus;
        else
            DirMoveNext = 0;
            iBetMoveNext = iBet;
        end
        
        % If we are not at edge of the available \beta and if a direction
        % has been chosen, move in that direction until a larger distance
        % is found.
        if DirMoveNext ~= 0 && iBetMoveNext > 1 && iBetMoveNext < BetN
            [~,~,ksDistTemp2] = kstest2( MinDistNextGroup(:) , DistBet2{iBetMoveNext + DirMoveNext} );
            while iBetMoveNext >= 1 && iBetMoveNext <= BetN && ksDistTemp2 < ksDist2
                ksDist2 = ksDistTemp2;
                iBetMoveNext = iBetMoveNext + DirMoveNext;
                if iBetMoveNext > 1 && iBetMoveNext < BetN 
                    [~,~,ksDistTemp2] = kstest2( MinDistNextGroup(:) , DistBet2{iBetMoveNext + DirMoveNext} );
                end
            end
        end
        BetaFit2 = BetVec(iBetMoveNext);
        
        % Add the fitted \beta to the file matrices.
        betaFits(iBet,iTrials) = BetaFit;
        ksFits(iBet,iTrials) = ksDist;
        
        betaFitsNext(iBet,iTrials) = BetaFit2;
        ksFitsNext(iBet,iTrials) = ksDist2;
        
    end
    
    % Save the found \beta.
    csvwrite(strcat('data/',TEXT,'_betaFits.txt'),betaFits);
    csvwrite(strcat('data/',TEXT,'_ksFits.txt'),ksFits);
    
    csvwrite(strcat('data/',TEXT,'_betaFitsNext.txt'),betaFitsNext);
    csvwrite(strcat('data/',TEXT,'_ksFitsNext.txt'),ksFitsNext);
    
end




%% Bias check - The fitted \beta as a function of the true \beta is considered.

if PlotCheck == true

% Nearest Neighbour Check
    MidBeta = mean(betaFits,2);
    StdBeta = std(betaFits,0,2);
    
    % Asymmetric error. skew1 and skew2 are the points above and below the
    % mean respectively. To find the asymmetric error, we simply calculate
    % the standard deviation of each separately, but mirrored around the
    % mean. That is, for the upper uncertainty, the points below the mean
    % are removed, and replaced by a mirror image of the points above the
    % mean.
    skew1 = zeros(BetN,1);skew2 = zeros(BetN,1);
    for iskew=1:BetN
        temp1 = betaFits(iskew,betaFits(iskew,:)>=MidBeta(iskew));
        temp2 = betaFits(iskew,betaFits(iskew,:)<=MidBeta(iskew));
        skew1(iskew) = std([temp1,2*MidBeta(iskew) - temp1],0,2);
        skew2(iskew) = std([temp2,2*MidBeta(iskew) - temp2],0,2);
    end
    % If the standard deviation is calculated from < 2 points, Matlab
    % returns nan. To plot them, we set the uncertainty to 0 instead.
    skew1(isnan(skew1)) = 0;
    skew2(isnan(skew2)) = 0;


    figure()
    % Plot straight line for comparison.
    plot(BetMin:BetStep:BetMax,BetMin:BetStep:BetMax,'LineWidth',2)
    axis([0 1 0 1])
    hold on
    xlabel('True \beta')
    ylabel('Fitted \beta')
    plot(BetMin:BetStep:BetMax,MidBeta,'r-','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBeta+skew1,'r--','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBeta-skew2,'r--','LineWidth',2)
    title('Asymmetric Error - Nearest Neighbour')
    axis on
    box on
    grid on
    hold off


    figure()
    % Plot straight line for comparison.
    plot(BetMin:BetStep:BetMax,BetMin:BetStep:BetMax,'LineWidth',2)
    axis([0 1 0 1])
    hold on
    xlabel('True \beta')
    ylabel('Fitted \beta')
    plot(BetMin:BetStep:BetMax,MidBeta,'r-','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBeta+StdBeta,'r--','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBeta-StdBeta,'r--','LineWidth',2)
    title('Symmetric Error - Nearest Neighbour')
    axis on
    box on
    grid on
    hold off

    
% Next-to-Nearest Neighbour Check
    MidBetaNext = mean(betaFitsNext,2);
    StdBetaNext = std(betaFitsNext,0,2);

    % Asymmetric error
    skew1Next = zeros(BetN,1);skew2Next = zeros(BetN,1);
    for iskew=1:BetN
        temp1 = betaFitsNext(iskew,betaFitsNext(iskew,:)>=MidBetaNext(iskew));
        temp2 = betaFitsNext(iskew,betaFitsNext(iskew,:)<=MidBetaNext(iskew));
        skew1Next(iskew) = std([temp1,2*MidBetaNext(iskew) - temp1],0,2);
        skew2Next(iskew) = std([temp2,2*MidBetaNext(iskew) - temp2],0,2);
    end
    % If the standard deviation is calculated from < 2 points, Matlab
    % returns nan. To plot them, we set the uncertainty to 0 instead.
    skew1Next(isnan(skew1Next)) = 0;
    skew2Next(isnan(skew2Next)) = 0;

    figure()
    % Plot straight line for comparison.
    plot(BetMin:BetStep:BetMax,BetMin:BetStep:BetMax,'LineWidth',2)
    axis([0 1 0 1])
    hold on
    xlabel('True \beta')
    ylabel('Fitted \beta')
    plot(BetMin:BetStep:BetMax,MidBetaNext,'r-','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBetaNext+skew1Next,'r--','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBetaNext-skew2Next,'r--','LineWidth',2)
    title('Asymmetric Error - Next-to-Nearest Neighbour')
    axis on
    box on
    grid on
    hold off


    figure()
    % Plot straight line for comparison.
    plot(BetMin:BetStep:BetMax,BetMin:BetStep:BetMax,'LineWidth',2)
    axis([0 1 0 1])
    hold on
    xlabel('True \beta')
    ylabel('Fitted \beta')
    plot(BetMin:BetStep:BetMax,MidBetaNext,'r-','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBetaNext+StdBetaNext,'r--','LineWidth',2)
    plot(BetMin:BetStep:BetMax,MidBetaNext-StdBetaNext,'r--','LineWidth',2)
    title('Symmetric Error - Next-to-Nearest Neighbour')
    axis on
    box on
    grid on
    hold off

end