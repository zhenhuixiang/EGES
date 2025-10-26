function Popreal = ES(PopDec,PopObj,TSDec,TSObj) % PopDec: candidate solution data; TSDec: archived data;
     % Non-dominated sorting of archived data
     [FrontNo,~] = NDSort(TSObj,inf);
     ND_TSObj = TSObj(FrontNo==1,:);
     ND_TSDec = TSDec(FrontNo==1,:);
     
     % Non-dominated sorting of candidate solution data
     [FrontNo1,~] = NDSort(PopObj,inf);
     ND_PopObj = PopObj(FrontNo1==1,:);
     ND_PopDec = PopDec(FrontNo1==1,:);   
     
     % parameters
     Popreal = [];
     rand1 = 1;

     %% CI sampling
     if rand < rand1
         if ~isempty(PopDec)
             CI_predicted = CI(ND_PopObj,TSObj);
             [~, index] = max(CI_predicted);
             Popreal = [Popreal; PopDec(index,:)];
         end
     end
end


