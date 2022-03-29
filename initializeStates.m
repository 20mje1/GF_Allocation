        function  [combs,sizeS] = initializeStates()
        %Function responsible for initializing the state space, which is
        %all possible states for the buffer statuses of all users. The size
        %of the state space = 3^(Number of users)
        vectors = { [1 2 3], [1 2 3], [1 2 3] , [1,2,3],[1 2 3],[1 2 3],[1 2 3],[1 2 3],[1 2 3],[1 2 3]}; 
        n = numel(vectors); 
        combs = cell(1,n); 
        [combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); 
        combs = cat(n+1, combs{:}); 
        combs = reshape(combs,[],n);
        sizeS = length(combs);
        
        end
