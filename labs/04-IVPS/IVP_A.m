function PariyangatS_Lab7(g,s_a,h,numsteps,eps_s)
%input: a cell array with function handles, the initial augmented state
%array, step size, number of steps to take, convergence criteria
%output: print historical augmented state array to an output file, 
%plot each state variable vs independent variable

%% RK4 Code
% Use RK4 to seed 4 rows in the augmented state array, corresponding to
% our first four time points. 
    N = length(g);
    slopes = zeros(4,N);
    for iter = 1:3
        state = s_a(iter,:);
        for c = 1:N
         slopes(1,c) = g{c}(state);
        end
        state = s_a(iter,:) + (h/2)*[1, slopes(1,:)];
        for c = 1:N
            slopes(2,c) = g{c}(state);
        end
        state = s_a(iter,:) + (h/2)*[1, slopes(2,:)];
        for c = 1:N
            slopes(3,c) = g{c}(state);
        end
        state = s_a(iter,:) + h*[1, slopes(3,:)];
        for c = 1:N
            slopes(4,c) = g{c}(state);
        end
        RK = (slopes(1,:) + 2*slopes(2,:) + 2*slopes(3,:) + slopes(4,:)) / 6;   
        s_a(iter+1,:) = s_a(iter,:) + h*[1, RK];       
    end

%% ABM Method
    %print first 4 rows of the historical augmented state array
    filename = 'PariyangatS_timetrace.csv';
    for idx=1:size(s_a,1)
        writematrix(s_a(idx,:),filename, 'WriteMode', 'append');
    end
    
    %predictor and corrector coefficients 
    pbeta = [55/24, -59/24, 37/24, -3/8];
    pmod = 251/270;
    cbeta = [3/8, 19/24, -5/24, 1/24];
    cmod = -19/270; 
    
    N = length(g); %size of the system
    predictors = zeros(2,N); %matrix of predictor values

    maincounter = 1;
    while((size(s_a,1)-1) < numsteps)
        currentind = size(s_a,1); %current row number
        current = [s_a(currentind,:)]; %state variables at current time step
        
        %calculate predictors
        predictors(1,:) = predictors(2,:);
        for idx = 1:N
            predictors(2,idx) = s_a(currentind,idx+1) + h*(pbeta(1)*(g{idx}(current))...
            + pbeta(2)*g{idx}(s_a(currentind-1,:)) + pbeta(3)*g{idx}(s_a(currentind-2,:)) ...
            + pbeta(4)*g{idx}(s_a(currentind-3,:)));
        end

        %calculate predictor modifier
        modifiers = zeros(2,N); %matrix of modifier values 
        if(maincounter ~= 1)
            for idx = 1:N
                predictormodifiers(idx) = pmod*(correctors(idx) - predictors(1,idx));
            end       
        else
            predictormodifiers = zeros(1,N);
        end
        for idx = 1:N
            modifiers(2,idx) = predictors(2,idx) + predictormodifiers(idx);
        end
    
        %iterating corrector to convergence
        complete = 0; 
        while(~complete) 
            currenttime = s_a(currentind,1)+h;
            current = cat(2,currenttime,modifiers(2,:));
            
            %calculate corrector and corrector modifiers
            modifiers(1,:) = modifiers(2,:);
            for idx = 1:N
                correctors(idx) = s_a(currentind,idx+1) + h*(cbeta(1)*(g{idx}(current)) + ...
                    cbeta(2)*g{idx}(s_a(currentind,:)) + cbeta(3)*g{idx}(s_a(currentind-1,:)) + ...
                    cbeta(4)*g{idx}(s_a(currentind-2,:)));
                modifiers(2,idx) = correctors(idx) + cmod*(correctors(idx)-predictors(2,idx));
            end
            
            %check convergence
            for idx = 1:N
                eps(idx) = abs((modifiers(2,idx) - modifiers(1,idx))/(modifiers(2,idx)));
            end
            check = all(eps < eps_s);
            if(check == 1)
                complete = 1;
            end
        end

        %final values 
        s_a(currentind + 1,:) = cat(2,current(1),modifiers(2,:));
        writematrix(s_a(currentind + 1,:),filename,'WriteMode','append');
        
        maincounter = maincounter + 1; 
    end
    
    %plot of state variables vs. independent variable
    for idx = 1:N
        figure(idx)
        plot(s_a(:,1),s_a(:,(idx+1)));
        titletext = "Variable #" + idx;
        title(titletext);
    end
    
    
end