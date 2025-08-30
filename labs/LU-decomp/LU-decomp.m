function PariyangatS_Lab6a(a,name)
%input: square numeric matrix, name of the output file
%output: writes a concatenated matrix of L,U,and bgen to the output file

    %size of the matrix
    N = size(a,1);

    %initialize bgen vector and empty l and u matrices
    bgen = transpose(1:N);
    l = zeros(N,N);
    u = zeros(N,N);
                                                        
    %iterate through entire a matrix
    for k=1:N  
        %find all u(k,k) candidates 
        uvals = [];
        for idx=k:N
            uvals(idx-k+1) = a(idx,k) - dot(l(idx,1:(k-1)),u(1:(k-1),k));
        end
        
        %find max u(k,k) value
        [~,maxind] = max(abs(uvals));
        maxind = maxind + (k-1); %adjust for row index
        
        %perform partial pivoting if necessary
        if(maxind ~= k && k < N)
            %swap rows of bgen
            bgen([k,maxind],:) = bgen([maxind,k],:);
            
            %swap rows of a
            a([k,maxind],:) = a([maxind,k],:);
  
            %swap rows of l
            l([k,maxind],:) = l([maxind,k],:);
        end
        
        %calculate u values
        for idx = k:N
            u(k,idx) = a(k,idx)-dot(l(k,1:(k-1)),u(1:(k-1),idx));
        end

        %calculate l values
        for idx = (k+1):N
            l(idx,k) = (a(idx,k)-dot(l(idx,1:(k-1)),u(1:(k-1),k)))/u(k,k);
        end
        l(k,k) = 1;
    end

    %concatenate final L-U-b matrix and write to file
    LUb = cat(2,l,u,bgen);
    writematrix(LUb,name, 'WriteMode', 'append');
    
end
