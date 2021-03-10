% This file treats a generated dataset similar to the positions of Buzzard
% nests. For conservation reasons, the real positions of the buzzards are
% not included in the library. Instead, a generated set with of the same
% form is included.
% Apart from this, the input is the Coulomb gas distances calculated with
% CoulombGasLoad.m and the uncertainties calculated with betaUncertainty.m.
% The Coulomb gas distances are necessary for the fit, but the
% uncertainties may be approximated with an included file in the library
% data.
% The output is the years grouped and fitted with NN and NNN spacing, both
% as a text file and in a figure with a linear fit.

clear;


% ---- Parameters for the Coulomb gas fit ----
% Number of particles
Nc = 200;
% Number of configurations
Nconf = 1e4;
% Number of iterations
Nsteps = 1e2;
% Step size
ep = Nc^(-1/3);

% \beta used in the fit.
BetMin = 0.0;
BetMax = 1.2;
BetStep = 0.1;


% The years looked at
YearsMin = 2000;
YearsMax = 2019;
Years = YearsMin:YearsMax ;
NYears = length(Years);

% The amount taken into the moving average.
Grouping = 5;

% Prelocations for grouped variables
BetaGroup = ones(NYears - Grouping+1,1)*NaN;
BetaSecondGroup = ones(NYears - Grouping+1,1)*NaN;
PopulationGroup = ones(NYears - Grouping+1,1)*NaN;
PopulationErrGroup = ones(NYears - Grouping+1,1)*NaN;
PopulationAll = zeros(NYears,1);
TEXTfull = ['Buzzard_Grouping',num2str(Grouping),'_Nconf',num2str(Nconf)];

% Targeted fit. This starts the fit at either 0.5 or the previous result
% and searches from there. It is faster than exploring the whole space, but
% may miss the global minimum if several local ones exist.
Targeted = false;

% Make the Coulomb gas fit. (Rather than simply loading the years and
% saving the distances. This may be done for speed.)
CoulombFit = true;

% Plot the spacing distributions of the individual groups.
PlotCoulomb = false;

% For the initial treatment of this data, a mask was applied to data set to
% cut away edge points. It is based on spheres drawn around each point. The
% points with spheres touching the edge are considered edge points.
% Implementation-wise, it is done on a grid rather than geometrically. The
% lattice sites are coloured by the spheres, which determines where the
% edge is. This seems to introduce a biased. The option is still included
% for completeness.
UseMask = false;

% For the unfolding, the width of the Gaussian has to be given. 4.5 is
% found by comparing to products of Ginibre.
a = 4.5;


%% Find the distances for all years
DistMinAllYears = cell(NYears , 1);
DistMinNextAllYears = cell(NYears , 1);
for iDist = 1:NYears
    
    % The first part is designed for the particular format of the buzzard
    % data. For other datasets, this part will need modification.
    disp(['Looking at the year ',num2str(Years(iDist))])
    DataCell = importdata(['data/BuzzardGenerated',num2str( Years(iDist) ),'.csv']);
    disp('Data imported');
    DataRaw = DataCell.data;

    % As the amount of nests is not the same for each year, we extract the
    % length.
    N = length(DataRaw(:,1));
    % As the buzzard nests are given in longitude and latitude, which are
    % very different this far north, we convert to kilometres in x- and
    % y-directions from the centre of the points.
    Middle = mean(DataRaw);
    
    % Prelocate for the positions in kilometres
    PosFull = zeros(size(DataRaw));

    % Convert from longitude and lattitude to kilometers from the centre.
    for iPoint = 1:N
        PosFull(iPoint,:) = [sign(DataRaw(iPoint,1) - Middle(1))*...
            deg2km(distance(Middle(1),Middle(2),DataRaw(iPoint,1),Middle(2)) ),...
                            sign(DataRaw(iPoint,2) - Middle(2))*...
            deg2km( distance(Middle(1),Middle(2),Middle(1),DataRaw(iPoint,2)) )];
    end

    disp('Data treated');

    % Find the non-unfolded distances to calculate the mean NN-spacing for
    % the Gaussian width in the unfolding
    disp('Calculating distances')
    DistancesNonUn = sqrt(abs((repmat(PosFull(:,1),1,length(PosFull(:,1)))...
                    -repmat(PosFull(:,1)',length(PosFull(:,1)),1)).^2 ...
                    + (repmat(PosFull(:,2),1,length(PosFull(:,2)))...
                    -repmat(PosFull(:,2)',length(PosFull(:,2)),1)).^2));

    % Find the smallest distances
    MinDistNonUn = NaN * ones(1,N);
    for iMin = 1:N
        % Cut out the diagonal before finding the minimum.
        DistTemp =  DistancesNonUn([1:(iMin-1),(iMin+1):N],iMin);
        MinDistNonUn(iMin) = min(DistTemp);
    end
    
    
    % It is possible to give each point an individual width, but this does
    % not seem to work as well as giving all points the same width.
    sigs = ones(N,1) * mean(MinDistNonUn) * a;
    
    % Find the unfolded distances
    Distances = UnfoldingGaussian(sigs,PosFull);


    % Find the minimum of the unfolded distances
    MinDist = NaN * ones(1,N);
    MinNextDist = NaN * ones(1,N);
    for iMin = 1:N
        % Cut out the diagonal before finding the minimum.
        DistTemp =  Distances([1:(iMin-1),(iMin+1):N],iMin);
        [minTemp,minTempPlace] = min(DistTemp);
        
        % Adjust for the cutout of the diagonal for the index of the
        % minimum.
        if minTempPlace >= iMin
            minTempPlace = minTempPlace + 1;
        end

        % To cut out any duplicates in the data, we set a threshold for
        % what is considered a duplicate. (Duplicates are present in the
        % buzzard data.)
        threshold = 1/N^2;

        % If the minimum is larger than the threshhold, simply accept it.
        % If not, check whether we have excluded the other point at the
        % location. If we have, simply find the next distance. Otherwise
        % exclude the point. (We do not want to exclude both points in case
        % of a duplicate.)
        if minTemp > threshold
            MinDist(iMin) = minTemp;
            MinNextDistTemp = mink(DistTemp,2);
            MinNextDist(iMin) = MinNextDistTemp(2);
        elseif isnan(MinDist(minTempPlace)) == true
            MinDist(iMin) = min(DistTemp(DistTemp>threshold));
            MinNextDistTemp = mink(DistTemp(DistTemp>threshold),2);
            MinNextDist(iMin) = MinNextDistTemp(2);
        end
    end

    disp('Distances calculated')

    % Here the mask is applied.
    if UseMask == true
        % Radius of spheres, chosen heuristically
        R = min( max(MinDistNonUn), 2.5*mean(MinDistNonUn));
        % Resolution of grid
        res = R/30;
        
        % CutOutBuzz takes the input Pos, R, res, PlotPos, and PlotThings.
        % The last two are for visualisation of the mask. (This is Figure
        % 5.4, left in Adam Mielke's PhD-thesis.)
        CutEntr = CutOutBuzz(PosFull,R,res,PosFull,false);
        
        % Choose only points in the bulk
        DistCut = MinDist(:,CutEntr);
        MinNextDistCut = MinNextDist(:,CutEntr);
    else
        DistCut = MinDist;
        MinNextDistCut = MinNextDist;
    end

    % Remove NaN. (These will be points that are too close to each other.)
    DistCut = (DistCut(~isnan(DistCut)));
    MinNextDistCut = (MinNextDistCut(~isnan(MinNextDistCut)));
    % The new amount of points.
    Ncut = length(DistCut);

    % Rescale the distances of the individual years such that the first
    % moment of the NN-distance is 1.
    DistCutScale = mean(DistCut);
    DistCut = sort(DistCut)/DistCutScale;
    MinNextDistCut = sort(MinNextDistCut)/DistCutScale;
    
    % Save the number of nests that are not duplicates. (And in the case
    % the mask was applied, it also excludes edge points.)
    PopulationAll(iDist) = Ncut;
    
    % Save the distances of the individual years in the cells, so they may
    % be more easily extracted for the grouping.
    DistMinAllYears{iDist} = DistCut;
    DistMinNextAllYears{iDist} = MinNextDistCut;
    
end


%% ----------------- Loop over the group starting points -----------------
for iStart = 0:(NYears - Grouping)

    disp(['Looking at the years ',num2str(YearsMin+iStart),'-',num2str(YearsMin+Grouping-1+iStart)]);
    
    % Text for loading the individual groups
    TEXT = ['Buzzard_start',num2str(iStart)];
    if UseMask == true
        TEXT = [TEXT,'_Mask'];
    end

    % Prelocation for the distances and populations in the group
    DistSave = [];
    DistNextSave = [];

    % Which years are in this group?
    YearsGroup = (YearsMin+iStart):(YearsMin+Grouping-1+iStart);

    % Loop over the years in the group.
    for iFiles = 1:length(YearsGroup)
        
        DistCut = DistMinAllYears{iStart + iFiles};
        MinNextDistCut = DistMinNextAllYears{iStart + iFiles};
        
        % Add distances to the group
        DistSave(end+1:end+length(DistCut)) = DistCut;
        DistNextSave(end+1:end+length(DistCut)) = MinNextDistCut;

    end
    
    % The population levels in this group
    Population = PopulationAll( (iStart + 1) : (iStart + Grouping) );
    
    % Scale the group distances.
    Scale = mean(DistSave);
    DistSave = sort(DistSave)/Scale;
    DistNextSave = sort(DistNextSave)/Scale;

    % -------------------------- Coulomb gas fit --------------------------
    if CoulombFit==true

        % For details on the Coulomb fit and the use of the parameters
        % passed to it, see the individual functions.
        if Targeted == true
            % Either choose 0.5 or the previously found \beta.
            if iStart == 0
                StartingPoint = 0.5;
                StartingPointNext = 0.5;
            else 
                StartingPoint = BetaFit;
                StartingPointNext = BetaFit2;
            end
            
            [BetaFit,ksDist] = CoulombGasFitTargeted(StartingPoint,Nc,...
                Nsteps,Nconf,ep,DistSave,BetMin,BetMax,BetStep);
            [BetaFit2,ksDist2] = CoulombGasFit2Targeted(StartingPointNext,...
                Nc,Nsteps,Nconf,ep,DistNextSave,BetMin,BetMax,BetStep);
        else
            [BetaFit,ksDist] = CoulombGasFit(TEXT,Nc,Nsteps,Nconf,ep,...
                DistSave,BetMin,BetMax,BetStep,PlotCoulomb,...
                ['Using the years: ',num2str(min(Years)),'-',num2str(max(Years))]);
            [BetaFit2,ksDist2] = CoulombGasFit2(TEXT,Nc,Nsteps,Nconf,ep,...
                DistNextSave,BetMin,BetMax,BetStep,PlotCoulomb,...
                ['Next-to-nearest, using the years: ',num2str(min(Years)),'-',num2str(max(Years))]);
        end


        BetaGroup(iStart+1) = BetaFit;
        BetaSecondGroup(iStart+1) = BetaFit2;

        %Save the fitted beta
        csvwrite(strcat('data/',TEXTfull,'_Beta.txt'),BetaGroup);
        csvwrite(strcat('data/',TEXTfull,'_BetaSecond.txt'),BetaSecondGroup);
    end

    % Save the mean population of the group.
    PopulationGroup(iStart+1) = mean(Population);
    % Propagation of Poisson error to mean. (No extra sqrt(N), because this
    % is just a generalisation of the error of an average.)
    PopulationErrGroup(iStart+1) = sqrt( sum(Population) ) / Grouping;
    % Save the population data
    csvwrite(strcat('data/',TEXTfull,'_Population.txt'),PopulationGroup);
end


%% ----------------------- Plot Population and Beta -----------------------
if length(BetaGroup)>2 % (The linear fit should only be done if we have more than 2 groups)
    
    % The middle years of the groups are calculated for the plot.
    MiddleYear = (YearsMin+(Grouping-1)/2):(YearsMax-(Grouping-1)/2);


    if CoulombFit==true
        % Try to load the error estimates for this particular configuration
        % of \beta and group size. If this does not exist, an approximation
        % is made with another file, where the error is rescaled by the
        % square root of the ratio between the chosen group size and the
        % group size of the backup file (5).
        % If the correct error is not used, the figure files will get the
        % label "_RoughError" to distinguish them.
        try
            BetaFitMat = importdata(['data/CoulombBetaFit_Group',num2str(Grouping),...
                '_N200_Nsteps100_Nconf',num2str(Nconf),'_BetMin',num2str(BetMin),...
                '_BetMax',num2str(BetMax),'_BetStep',num2str(BetStep),'_betaFits.txt']);
            BetaFitNextMat = importdata(['data/CoulombBetaFit_Group',num2str(Grouping),...
                '_N200_Nsteps100_Nconf',num2str(Nconf),'_BetMin',num2str(BetMin),...
                '_BetMax',num2str(BetMax),'_BetStep',num2str(BetStep),'_betaFitsNext.txt']);
            ErrorFound = true;
        catch
            disp('Correct uncertainty file not found, using rough estimate.')
            TEXTfull = [TEXTfull,'_RoughError'];
            BetaFitMat = importdata('data/CoulombBetaFit_Group5_N200_Nsteps100_Nconf1000_BetMin0_BetMax1.5_BetStep0.1_betaFits.txt')/sqrt(Grouping/5);
            BetaFitNextMat = importdata('data/CoulombBetaFit_Group5_N200_Nsteps100_Nconf1000_BetMin0_BetMax1.5_BetStep0.1_betaFitsNext.txt')/sqrt(Grouping/5);
            ErrorFound = false;
        end

    % Load the file made by betaUncertainty.m and find the uncertainty of
    % the points. First prelocate both the errors and the covariance.
    BetaFitErr = zeros(size(BetaGroup));
    BetaFitNextErr = zeros(size(BetaSecondGroup));
    BetaFitErrMat = zeros(length(BetaGroup));
    BetaFitNextErrMat = zeros(length(BetaSecondGroup));

    for iBeta = 1:length(BetaFitErr)
        % Calculate the errors. If the backup file is used, the index is
        % generated from the value of the fitted \beta. Otherwise we simply
        % find the correct index.
        if ErrorFound == true
            idBeta = find((BetMin:BetStep:BetMax)==BetaGroup(iBeta));
            idBetaNext = find((BetMin:BetStep:BetMax)==BetaSecondGroup(iBeta));
            BetaFitErr(iBeta) = std( BetaFitMat( idBeta ,: ) );
            BetaFitNextErr(iBeta) = std( BetaFitNextMat( idBetaNext ,: ) );

            % Set the diagonal of the covariance matrix to the error.
            BetaFitErrMat(iBeta,iBeta) = BetaFitErr(iBeta);
            BetaFitNextErrMat(iBeta,iBeta) = BetaFitNextErr(iBeta);
            % The off-diagonal entries are found by averaging over entries
            % with the given distance between them. This is done in the
            % function CovCons(). The iBeta2-index loops over the upper
            % triangle.
            for iBeta2 = (iBeta+1):length(BetaFitErr)
                idBeta2 = find((BetMin:BetStep:BetMax)==BetaGroup(iBeta2));
                CovVecThis = CovCons(BetaFitMat( idBeta,: ),iBeta2-1);
                CovVecThat = CovCons(BetaFitMat( idBeta2,: ),iBeta2-1);
                CovVecNextThis = CovCons(BetaFitMat( round(10*BetaSecondGroup(iBeta)+1),: ),iBeta2-1);
                CovVecNextThat = CovCons(BetaFitMat( round(10*BetaSecondGroup(iBeta2)+1),: ),iBeta2-1);
                % The combination of two different \beta is done with the
                % following heuristic meassure.
                BetaFitErrMat(iBeta,iBeta2) = sqrt((CovVecThis(end)^2 + CovVecThat(end)^2)/2);
                BetaFitNextErrMat(iBeta,iBeta2) = sqrt((CovVecNextThis(end)^2 + CovVecNextThat(end)^2)/2);
            end
        else
            BetaFitErr(iBeta) = std( BetaFitMat( round(10*BetaGroup(iBeta)+1),: ) );
            BetaFitNextErr(iBeta) = std( BetaFitNextMat( round(10*BetaSecondGroup(iBeta)+1),: ) );

            % Set the diagonal of the covariance matrix to the error.
            BetaFitErrMat(iBeta,iBeta) = BetaFitErr(iBeta);
            BetaFitNextErrMat(iBeta,iBeta) = BetaFitNextErr(iBeta);
            % The off-diagonal entries are found by averaging over entries
            % with the given distance between them. This is done in the
            % function CovCons(). The iBeta2-index loops over the upper
            % triangle.
            for iBeta2 = (iBeta+1):length(BetaFitErr)
                CovVecThis = CovCons(BetaFitMat( round(10*BetaGroup(iBeta)+1),: ),iBeta2-1);
                CovVecThat = CovCons(BetaFitMat( round(10*BetaGroup(iBeta2)+1),: ),iBeta2-1);
                CovVecNextThis = CovCons(BetaFitMat( round(10*BetaSecondGroup(iBeta)+1),: ),iBeta2-1);
                CovVecNextThat = CovCons(BetaFitMat( round(10*BetaSecondGroup(iBeta2)+1),: ),iBeta2-1);
                % The combination of two different \beta is done with the
                % following heuristic meassure.
                BetaFitErrMat(iBeta,iBeta2) = sqrt((CovVecThis(end)^2 + CovVecThat(end)^2)/2);
                BetaFitNextErrMat(iBeta,iBeta2) = sqrt((CovVecNextThis(end)^2 + CovVecNextThat(end)^2)/2);
            end
        end
        
        % Fill the lower triangle.
        BetaFitErrMat = (BetaFitErrMat + BetaFitErrMat')/2;
        BetaFitNextErrMat = (BetaFitNextErrMat + BetaFitNextErrMat')/2;
    end

    % \beta-fit with fitting method from Barlow, using fitting matrix.
    C = [ones(length(MiddleYear),1),MiddleYear'];
    ParamFit = inv(C' * inv(BetaFitErrMat) * C) * C' * inv(BetaFitErrMat) * BetaGroup;

    % Make the same fit for the NNN-\beta
    ParamFitNext = inv(C' * inv(BetaFitNextErrMat) * C) * C' * inv(BetaFitNextErrMat) * BetaSecondGroup;
    
    end


    % Make the linear fit for the population. The fit is done to the
    % individual years rather than the groups, so we don't have to deal
    % with the error propagation.
    Cpop = [ones(YearsMax-YearsMin+1,1),(YearsMin:YearsMax)'];
    Vpop = diag(sqrt(PopulationAll));
    ParamPop = inv(Cpop' * inv(Vpop) * Cpop) * Cpop' * inv(Vpop) * PopulationAll;

    
    %% Plot the fitted \beta
    if CoulombFit==true
        % Initialise figure and make plot with error bars of NN-spacing.
        BetaFig = figure();
        errorbar(MiddleYear,BetaGroup,BetaFitErr,'bx');

        % Plot NNN-spacing on top as well as the linear fits.
        hold on
        errorbar(MiddleYear,BetaSecondGroup,BetaFitNextErr,'r+');
        x = (1998:2021)';
        yBeta = (ParamFit(2)*x + ParamFit(1))';
        yBetaNext = (ParamFitNext(2)*x + ParamFitNext(1))';
        plot(x,yBeta,'b-','LineWidth',2);
        plot(x,yBetaNext,'r-','LineWidth',2);
        hold off

        % Set ticks
        xticks(2000:2019);
        yticks(BetMin:BetStep:(BetMax+1));

        axis([2000 2019 0 (max(BetaGroup)+max(BetaFitErr))*1.32]);

        % Remove every other tick on the axes to make the plot less
        % cluttered.
        ax = gca;
        xlabels = string(ax.XAxis.TickLabels); % extract
        xlabels(2:2:end) = nan; % remove every other one
        ax.XAxis.TickLabels = xlabels; % set
        ylabels = string(ax.YAxis.TickLabels); % extract
        ylabels(2:2:end) = nan; % remove every other one
        ax.YAxis.TickLabels = ylabels; % set

        axis([MiddleYear(1)-0.5 MiddleYear(end)+0.5 0 (max(BetaGroup)+max(BetaFitErr))*1.32]);

        % Make the figure look nice.
        box on
        grid on
        legend('Nearest Neighbour','Next-to-Nearest Neighbour','Location','northwest');
        xlabel('Middle Year')
        ylabel('Fitted \beta')
        set(gca,'fontsize', 18);

        % Save figure, both as .pdf and as .fig. The first three lines make
        % sure that the pdf is cut down to size rather than being save in
        % the middle of an A4 page.
        set(BetaFig,'Units','Inches');
        PlotPos = get(BetaFig,'Position');
        set(BetaFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[PlotPos(3), PlotPos(4)])
        savefig(['figure/',TEXTfull,'_Beta.fig']);
        print(BetaFig,['figure/',TEXTfull,'_Beta.pdf'],'-dpdf');
    end

    
    %% Plot population and fit
    % Initialise figure and make plot with error bars.
    PopFig = figure();
    errorbar(MiddleYear,PopulationGroup,PopulationErrGroup,'kx');

    hold on
    x = (1998:2021)';
    yPop = (ParamPop(2)*x + ParamPop(1))';
    plot(x,yPop,'k-','LineWidth',2);
    hold off

    xticks(2000:2019);
    yticks(0:20:280);

    axis([2000 2019 0 280]);

    % Remove every other tick on the axes to make the plot less
    % cluttered.
    ax = gca;
    xlabels = string(ax.XAxis.TickLabels); % extract
    xlabels(2:2:end) = nan; % remove every other one
    ax.XAxis.TickLabels = xlabels; % set
    ylabels = string(ax.YAxis.TickLabels); % extract
    ylabels(2:2:end) = nan; % remove every other one
    ax.YAxis.TickLabels = ylabels; % set

    axis([MiddleYear(1)-0.5 MiddleYear(end)+0.5 0 200]);

    % Make the figure look nice.
    box on
    grid on
    legend off;
    xlabel('Middle Year')
    ylabel('Average Population')
    set(gca,'fontsize', 18);

    % Save figure, both as .pdf and as .fig. The first three lines make
    % sure that the pdf is cut down to size rather than being save in
    % the middle of an A4 page.
    set(PopFig,'Units','Inches');
    PlotPos = get(PopFig,'Position');
    set(PopFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[PlotPos(3), PlotPos(4)])
    savefig(['figure/',TEXTfull,'_Population.fig']);
    print(PopFig,['figure/',TEXTfull,'_Population.pdf'],'-dpdf');

end
