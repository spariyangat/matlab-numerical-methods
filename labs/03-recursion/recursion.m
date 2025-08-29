function x = PariyangatS_Lab6(A,b)
%input: non-singular square matrix of numbers A, and a column vector of
%numbers b
%output: displays A and b matrices recieved from calling program, returns
%solution x vector

    %display A and b matrices recieved from calling program
    fprintf('Initial Matrix A:\n');
    disp(A)
    fprintf('Initial Matrix B:\n');
    disp(b)
    
    %check for terminal instance
    if(size(A,1) <= 1)
        x = b/A;   

    %if the instance is non-terminal
    else
        m = cat(2,A,b);
        
        %find largest magnitude value in column 1
        [~,rowlargest] = max(abs(A(1:size(A,1),1)));

        %swap current row with row with largest first value
        if(rowlargest ~= 1)
            temp = m(1,:);
            m(1,:) = m(rowlargest,:);
            m(rowlargest,:) = temp;
        end
    
        %perform row operations
        for idx = 2:size(m,1)
            r = m(idx,1) / m(1,1);
            m(idx,:) = m(idx,:) - ((r)*(m(1,:)));
        end

        %identify the smaller system
        smallA = m(2:size(m,1),2:(size(m,2)-1));
        smallb = m(2:size(m,1),size(m,2));
        
        %spawn the child
        x = PariyangatS_Lab6(smallA, smallb);

        %back substition to solve for next x
        nextx = m(1,size(m,2));
        for idx = 1:(size(m,2)-2)
            count1 = size(m,2)-1-idx;
            nextx = nextx - (x(count1)*m(1,count1+1));
        end
        nextx = nextx / m(1,1);
        x = [nextx;x];

    end
end

