function [CI] = CI(ND_PopObj,TSObj)
    % Non-dominated sorting of archived data
    [FrontNo,~] = NDSort(TSObj,inf);
    ND_TSObj = TSObj(FrontNo==1,:);
    
    % angle 
    angle  = acos(1-pdist2((max(ND_PopObj,0)),max(ND_TSObj,0),'cosine'));   % Calculate the angle between the two sets of vectors
    [Angle,~] = min(angle,[],2); % Calculate the minimum angle between each candidate solution and the ND sample

    % Objective Normalization
    ND_PopObj = (ND_PopObj - min(TSObj,[],1))./(max(TSObj,[],1)-min(TSObj,[],1));
    TSObj = (TSObj- min(TSObj,[],1))./(max(TSObj,[],1)-min(TSObj,[],1));
    Z = min(TSObj,[],1);

    % CN
    ddt  = pdist2(ND_PopObj, TSObj,'euclidean');
    ddt(ddt==0) = inf;
    CN = min(ddt,[],2);

    % DC
    DC  = pdist2(ND_PopObj, Z,'euclidean');

    % Indicators Normalization
    N_angle = (Angle - min(Angle,[],1))./(max(Angle,[],1)-min(Angle,[],1));
    N_CN = (CN - min(CN,[],1))./(max(CN,[],1)-min(CN,[],1));
    N_DC = (DC - min(DC,[],1))./(max(DC,[],1)-min(DC,[],1));
    
    % CI
    r1 = rand;
    r2 = rand;
    r3 = rand;
    CI = r1 * N_CN + r2 * N_angle - r3 * N_DC;
end

