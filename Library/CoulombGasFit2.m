% Next-to-nearest neighbour fit for Coulomb gas, where each value of \beta
% is tested.

function [BetaFit,ksDist] = CoulombGasFit2(TEXT,Nc,Nsteps,Nconf,ep,DistSave,BetMin,BetMax,BetStep,PlotThings,TitleText)
    
    % NNN-spacing Poisson
    PoissonNext1 = @(s) pi^2/8 .* s.^3 .* exp(-pi/4*s.^2);
    

    disp('Coulomb Gas Fit - Next to Nearest')

    % Starting point, Kolmogorov distance set to 1.
    BetaFit = 0;
    ksDist = 1;

    % Go through all \beta and find the one with smallest Kolmogorov
    % distance.
    for betaNum = fliplr(BetMin:BetStep:BetMax)

        disp(['Comparing to \beta=',num2str(betaNum)]);

        ColorVec = [ max(1 - (betaNum - BetMin) / (BetMax - BetMin) ,0) , 0 ,...
            max((betaNum - BetMin) / (BetMax - BetMin) ,0)];

        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(betaNum),'_ep',num2str(ep)];

            dists = importdata(['data/',TEXTc,'_NextDists.txt']);

            [~,~,ksDistTemp] = kstest2( DistSave , sort(dists(:)) );

        disp(['Kolmogorov distance is ',num2str(ksDistTemp)]);

        if ksDistTemp < ksDist
            BetaFit = betaNum;
            ksDist = ksDistTemp;
        end

    end

    disp(['Best fit is \beta=',num2str(BetaFit)])

    % Plot Poisson, Ginibre, and the Coulomb fit
    if PlotThings == true
        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(BetaFit),'_ep',num2str(ep)];

        if BetaFit ~= 0
            dists = importdata(['data/',TEXTc,'_NextDists.txt']);

            [Counts,Edges] = histcounts(dists,40,'Normalization','pdf');

            x = Edges(1:end-1)+diff(Edges)/2;
            y = Counts;
        end

        % Plot the Ginibre spacing distribution
        xGin = 0:0.01:3;
        yGin = xGin;

        for iy = 1:length(yGin)
            yGin(iy) = 1.1429*GinibreNNN(xGin(iy)*1.1429,50);
        end

        figBeta = figure();
        histogram(DistSave,30,'Normalization','pdf','FaceColor','cyan')
        hold on
        plot(xGin,yGin,'b--','LineWidth',2)
        if BetaFit ~= 0
            plot(x,y,'k','LineWidth',2)
        end
        plot(xGin,PoissonNext1(xGin),'r:','Linewidth',2);
        hold off
        if BetaFit == 0
            leg = legend('Data','Ginibre','Poisson');
        else
            leg = legend('Data','Ginibre',['\beta=',num2str(BetaFit)],'Poisson');
        end

        set(gca,'FontSize',20)
        axis([0 3 0 1.4])
        axis square
        grid on
        box on
        leg.FontSize = 18;
        title([TitleText,' with ks=',num2str(ksDist)])
        savefig(['figure/',TEXT,'_Next_beta.fig']);
        print(figBeta,['figure/',TEXT,'_Next_beta.pdf'],'-dpdf');
    end

end