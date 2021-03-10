% Nearest neighbour spacing fit for Coulomb, where a starting point is used
% to increase the speed of the fit. This method may miss the global
% minimum if several local ones exist.

function [BetaFit,ksDist] = CoulombGasFitTargeted(StartingPoint,Nc,Nsteps,Nconf,ep,DistSave,BetMin,BetMax,BetStep)

    disp('Coulomb Gas Fit')
    
    
    BetVec = BetMin:BetStep:BetMax;
    
    BetMaxInd = length(BetVec);

    BetaFit = StartingPoint;
    

    % Load distances.
    disp(['Comparing to \beta=',num2str(BetaFit)]);
    
    TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
            '_Nconf',num2str(Nconf),'_beta',num2str(BetaFit),'_ep',num2str(ep)];

    dists = importdata(['data/',TEXTc,'_dists.txt']);

    [~,~,ksDist] = kstest2( DistSave , sort(dists(:)) );
    disp(['Kolmogorov distance is ',num2str(ksDist)]);

    % Index of starting point.
    iBet = find(BetaFit == BetVec);
    
    % Test the \beta below the start.
    BetMoveMinus = iBet-1;
    % If we are at the smallest index, set the Kolmogorov distance to
    % 1, so this direction is not picked.
    if BetMoveMinus >= 1
        
        BetaLook = BetVec(BetMoveMinus);
        disp(['Comparing to \beta=',num2str(BetaLook)]);
        
        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(BetaLook),'_ep',num2str(ep)];
        
        dists = importdata(['data/',TEXTc,'_dists.txt']);
        
        [~,~,ksDistMinus] = kstest2( DistSave , sort(dists(:)) );
        disp(['Kolmogorov distance is ',num2str(ksDistMinus)]);
    else
        ksDistMinus = 1;
    end

    % Test the \beta above
    BetMovePlus = iBet+1;
    % If we are at the highest index, set the Kolmogorov distance to 1,
    % so this direction is not picked.
    if BetMovePlus <= BetMaxInd
        BetaLook = BetVec(BetMovePlus);
        disp(['Comparing to \beta=',num2str(BetaLook)]);
        
        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(BetaLook),'_ep',num2str(ep)];
        
        dists = importdata(['data/',TEXTc,'_dists.txt']);
        
        [~,~,ksDistPlus] = kstest2( DistSave , sort(dists(:)) );
        disp(['Kolmogorov distance is ',num2str(ksDistPlus)]);
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
    if DirMove ~= 0 && iBetMove > 1 && iBetMove < BetMaxInd
        BetaLook = BetVec(iBetMove + DirMove);
        disp(['Comparing to \beta=',num2str(BetaLook)]);
        
        TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                '_Nconf',num2str(Nconf),'_beta',num2str(BetaLook),'_ep',num2str(ep)];
        
        dists = importdata(['data/',TEXTc,'_dists.txt']);
        
        [~,~,ksDistTemp] = kstest2( DistSave , sort(dists(:)) );
        disp(['Kolmogorov distance is ',num2str(ksDistTemp)]);
        
        while iBetMove >= 1 && iBetMove <= BetMaxInd && ksDistTemp < ksDist
            ksDist = ksDistTemp;
            iBetMove = iBetMove + DirMove;
            if iBetMove > 1 && iBetMove < BetMaxInd
                BetaLook = BetVec(iBetMove + DirMove);
                disp(['Comparing to \beta=',num2str(BetaLook)]);

                TEXTc = ['Coulomb_N',num2str(Nc),'_Nsteps',num2str(Nsteps),...
                        '_Nconf',num2str(Nconf),'_beta',num2str(BetaLook),'_ep',num2str(ep)];

                dists = importdata(['data/',TEXTc,'_dists.txt']);

                [~,~,ksDistTemp] = kstest2( DistSave , sort(dists(:)) );
                disp(['Kolmogorov distance is ',num2str(ksDistTemp)]);
            end
        end
    end
    
    BetaFit = BetVec(iBetMove);
    
    disp(['Best fit is \beta=',num2str(BetaFit)])

end