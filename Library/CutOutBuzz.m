% Mask for determining the bulk of a dataset with a non-trivial edge, such
% as holes or simply an unknown shape. For the buzzard nests, this seems to
% introduce a bias, but it may be useful for other datsets.

function CutEntr = CutOutBuzz(Pos,R,res,PlotPos,PlotThings)

    % Make a grid around the points.
    gridX = min(Pos(:,1))-R : res : max(Pos(:,1))+R;
    gridY = min(Pos(:,2))-R : res : max(Pos(:,2))+R;
    Lx = length(gridX);
    Ly = length(gridY);
    
    % All grid coordinates.
    grid = repmat(gridX',Ly,2);

    for iGrid = 1:Ly
        grid(1+(iGrid-1)*Lx:iGrid*Lx,2) = gridY(iGrid);
    end

    % Squared distances from data points to each grid site.
    distSq = abs( (repmat(Pos(:,1),1,length(grid(:,1)))-repmat(grid(:,1)',length(Pos(:,1)),1)).^2 ...
          + (repmat(Pos(:,2),1,length(grid(:,2)))-repmat(grid(:,2)',length(Pos(:,2)),1)).^2 );

    % Which lattice sites are closer than R to the data points. This is the
    % mask that estimates the support of the dataset.
    GridIsland = reshape(min(distSq)<R^2,Lx,Ly);
    
    % Which grid sites are on the edge of the support?
    GridBorder = zeros(size(GridIsland));
    % The outer part of the grid is at least a border.
    GridBorder(:,1) = 1; GridBorder(:,end) = 1;
    GridBorder(1,:) = 1; GridBorder(end,:) = 1;
    
    % Go through grid sites and see if it is on the border of the support.
    for ix = 2:Lx-1
    for iy = 2:Ly-1
        
        if GridIsland(ix,iy) && (GridIsland(ix-1,iy-1) == false ||...
           GridIsland(ix+1,iy-1) == false || GridIsland(ix-1,iy+1) == false...
           || GridIsland(ix+1,iy+1) == false )
       
                    GridBorder(ix,iy) = 1;
                    
        end
        
    end
    end
    
    % Accept only the data points that are more than R from the border.
    CutEntr = min(distSq(:,logical(reshape(GridBorder,1,Lx*Ly))),[],2) >= R^2;
    
    % Visualise the mask.
    if PlotThings == true
        figure();
        plot(grid(min(distSq)<R^2,1),grid(min(distSq)<R^2,2),'.','Color', [0.7 0.7 0.7]);
        hold on
        plot(Pos(:,1),Pos(:,2),'go','LineWidth',2);
        plot(Pos(CutEntr,1),Pos(CutEntr,2),'yx','LineWidth',2);
        plot(PlotPos(:,1),PlotPos(:,2),'bo','LineWidth',2);
        plot(PlotPos(CutEntr,1),PlotPos(CutEntr,2),'rx','LineWidth',2);
        hold off
        axis([-15 15 -10 10])
    end
    
end