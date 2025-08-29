function NDDP = PariyangatS_Lab5(file,target,eps_s)
%input: name of file that contains data set, target point of interpolation,
%convergence criteria epsilon
%output: displays y-approximations, final divided differences table,
%returns a row vector with final approximation and level of convergence

    data = load(file);

    %i: index of the first data point
    i = 1;  
    %check if target is out of bounds of the data set
    found = false;
    if(data(1) > target)
        found = true;
    elseif(data(length(data)) < target)
        i = length(data)-1;
        found = true;
    end

    %find the first two data points and add to table
    count = 1;
    while(~found)
        if(data(count) > target)
            i = count-1;
            found = true;
        end
        count = count + 1; 
    end
    T = data(i:i+1, :);
    
    %calculate first divided difference
    T(1,3) = (T(2,2)-T(1,2)) / (T(2,1) - T(1,1));
    numrows = 2;
    %calculate first linear approximation
    yappx = T(1,2) + T(1,3)*(target-T(1,1)) 
   
    %indices of candidates 1 and 2
    cand_1 = i - 1; 
    cand_2 = i + 2; 
    %number of data points to the left and right of target
    numleft = 0; 
    numright = 0; 
    
    %keep track of convergence 
    converged = 0; 

    complete = 0; 
    while(~complete)
        %if-else tree to select next point from data set
        %pointer: index of new point 
        if(cand_1 <= 0) && (cand_2 > (length(data)))
            complete = 1;
        elseif(cand_2 > length(data))
            pointer = cand_1; 
            cand_1 = cand_1 - 1;
            numleft = numleft + 1;
        elseif(cand_1 <= 0)
            pointer = cand_2;
            cand_2 = cand_2 + 1;
            numright = numright + 1; 
        elseif(abs(target-data(cand_1)) < abs(target-data(cand_2)))
            pointer = cand_1;
            cand_1 = cand_1 - 1;
            numleft = numleft + 1; 
        elseif(abs(target-data(cand_1)) > abs(target-data(cand_2)))
            pointer = cand_2;
            cand_2 = cand_2 + 1;
            numright = numright + 1; 
        elseif(numleft < numright)
            pointer = cand_1;
            cand_1 = cand_1 - 1;
            numleft = numleft + 1; 
        else
            pointer = cand_2;
            cand_2 = cand_2 + 1;
            numright = numright + 1;
        end
    
        if(complete == 0)
            T(numrows+1,1:2) = data(pointer, :);
            numrows = numrows + 1;

            %calculate divided differences
            col = 3;
            for row = (numrows-1):-1:1
                T(row,col) = (T(row+1,col-1)-T(row,col-1)) / (T(numrows,1)-T(row,1));
                col = col + 1;
            end
        
            %calculate next term in NDDP
            E = T(1,numrows+1);
            for idx = 1:numrows-1
                E = E * (target-T(idx,1));
            end
            
            %calculate next approximation 
            yappx = yappx + E
                  
            %check for convergence
            eps = abs(E/yappx);
            if(eps < eps_s)
                fprintf('converged\n')
                converged = 1;
                complete = 1;
            end
        end    
    end
    
    NDDP = [converged,yappx,eps];
    T

end
