% Nearest neighbour fit for Coulomb gas, where each value of \beta is
% tested.


function [BetaFit,ksDist] = CoulombGasFit(TEXT,Nc,Nsteps,Nconf,ep,DistSave,BetMin,BetMax,BetStep,PlotThings,TitleText)
    
    
    % Ginibre spacing function
    GinibreSpacing = @(s,cut) 2*s*sum( ((s.^2).^(1:cut)) ./ (gammainc( s.^2,1 + (1:cut),'upper').*factorial(1:cut)) )...
                * exp( sum(log(gammainc(s^2, 1 + (0:cut),'upper') ) ) );

    % Wigner Surmise (beta=1 is 2D Poisson)
    WS = @(beta,x) 2*gamma(beta/2+1).^(beta+1) ./ gamma((beta+1)/2).^(beta+2).*...
                     x.^beta .* exp(-(x.*gamma(1+beta/2)/gamma((beta+1)/2)).^2) ;

    disp('Coulomb Gas Fit')

    
    % Starting point, Kolmogorov distance set to 1.
    BetaFit = 0;
    ksDist = 1;

    % Go through all \beta and find the one with smallest Kolmogorov
    % distance.
    for betaNum = fliplr(BetMin:BetStep:BetMax)

        disp(['Comparing to \beta=',num2str(betaNum)]);

        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(betaNum),'_ep',num2str(ep)];

        dists = importdata(['data/',TEXTc,'_dists.txt']);

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

        % If the fit is Poisson, simply plot this.
        if BetaFit ~= 0
            dists = importdata(['data/',TEXTc,'_dists.txt']);

            [Counts,Edges] = histcounts(dists,40,'Normalization','pdf');

            x = Edges(1:end-1)+diff(Edges)/2;
            y = Counts;
        end

        % Plot the Ginibre spacing distribution
        xGin = 0:0.01:3;
        yGin = xGin;

        for iy = 1:length(yGin)
            yGin(iy) = 1.1429*GinibreSpacing(xGin(iy)*1.1429,50);
        end

        figBeta = figure();
        histogram(DistSave,30,'Normalization','pdf','FaceColor','cyan')
        hold on
        plot(xGin,yGin,'b--','LineWidth',2)
        if BetaFit ~= 0
            plot(x,y,'k','LineWidth',2)
        end
        plot(xGin,WS(1,xGin),'r:','Linewidth',2);
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
        savefig(['figure/',TEXT,'_beta.fig']);
        print(figBeta,['figure/',TEXT,'_beta.pdf'],'-dpdf');
    end

end