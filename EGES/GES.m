function [PopDec,PopObj] = GES(A1,Model,Problem)
    typeOfGroups = 2;
    numberOfGroups = 4;


    [V, Problem.N] = UniformPoint(Problem.N,Problem.M); % Reference Points
    
    Zmin   = min(A1.objs,[],1);
    if size(A1.decs,1) >= Problem.N
         Next = NSGAIIIEnvironmentalSelection(A1,Problem.N,V,Zmin); 
    end
    
    Population = Next;
    Zmin           = min(Population(all(Population.cons<=0,2)).objs,[],1); % Ideal point
    
    w = 1; wmax = 20;
    while w <= wmax
        MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));                        % Tournament Selection
        Offspring  = GESoperating(Problem, Population(MatingPool), numberOfGroups, typeOfGroups, Model, A1);    % GES
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);                             % Ideal point
        Population = NSGAIIIEnvironmentalSelection([Population,Offspring],Problem.N,V,Zmin);                % GLMO_NSGAIIIEnvironmentalSelection
        w = w+1;
    end
    PopDec = Population.decs;
    PopObj = Population.objs;
end