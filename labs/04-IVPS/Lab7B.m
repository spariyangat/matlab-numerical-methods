function Lab7B(g,s_a,h,numsteps,eps_s)

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

        maincounter = maincounter + 1; 
    end
    


    
    %finding step where time reaches 19 seconds
    dampedstep = 1;
    found = false;
    while(~found)
        if(s_a(dampedstep,1) >= 19)
            found = true; 
        end
        dampedstep = dampedstep + 1; 
    end

    %finding step where time reaches 20 seconds
    startstep = dampedstep;
    found = false;
    while(~found)
        if(s_a(startstep,1) >= 20)
            found = true; 
        end
        startstep = startstep + 1; 
    end
    
    %finding step where time reaches 20.3 seconds
    midstep = startstep;
    found = false;
    while(~found)
        if(s_a(midstep,1) >= 20.3)
            found = true; 
        end
        midstep = midstep + 1; 
    end

    %finding step where time reaches 21.3 seconds
    endstep = startstep;
    found = false;
    while(~found)
        if(s_a(endstep,1) >= 21.3)
            %this is where one full cycle of the
            %oscillation is completed
            found = true; 
        end
        endstep = endstep + 1; 
    end

    %finding max displacement for mass 1 after 19 seconds
    maxm1 = max(s_a(dampedstep:(dampedstep+100),2))

    %finding max displacement for mass 2 after 19 seconds
    maxm2 = max(s_a(dampedstep:(dampedstep+100),4))

    maxm2/maxm1
    
    figure(1)
    subplot(3,1,1)
    plot(s_a(dampedstep:(dampedstep+250),1),s_a(dampedstep:(dampedstep+250),2));
    title("x1 vs time");
    xlabel("time");
    ylabel("position of m1");

    %figure(2)
    subplot(3,1,2)
    plot(s_a(dampedstep:(dampedstep+250),1),s_a(dampedstep:(dampedstep+250),4));
    title("x2 vs time");
    xlabel("time");
    ylabel("position of m2");

    %plotting driving force using Fsin(omegat)    
    %figure(3)
    subplot(3,1,3)
    t = s_a(dampedstep:(dampedstep+250),1);
    drivingforce = 10*sin(5.101*t);
    plot(t,drivingforce);
    title("driving force");
    xlabel("time");
    ylabel("driving force")
    
    %find phi values
    [~,maxidx1] = max(s_a(startstep:endstep,2));
    timemax1 = s_a(startstep+maxidx1,1) %time of the max displacement for m1

    [~,maxidx2] = max(s_a(startstep:endstep,4));
    timemax2 = s_a(startstep+maxidx2,1) %time of the max displacement for m2
    [~,maxidxdriving] = max(drivingforce(1:(midstep-dampedstep)));
    timemaxdriving = s_a(dampedstep+maxidxdriving,1) %time of the max value of driving force

    %time lag for 1
    tl1 = timemax1 - timemaxdriving 

    %time lag for 2
    tl2 = timemax2 - timemaxdriving 


    %find step where driving force is 0
    firstzero = 0;
    counter = dampedstep + 1;
    found = 0;
    drivingforce(counter-dampedstep);
    while(~found)
        if(drivingforce(counter-dampedstep) <= 0)
            firstzero = counter;
            found = true;
        end
        counter = counter + 1;
    end
    
    %find the second time where driving force is 0 
    secondzero = firstzero;
    counter = firstzero + 1;
    found = 0;
    while(~found)
        if(drivingforce(counter-dampedstep) >= 0)
            secondzero = counter;
            found = true;
        end
        counter = counter + 1;
    end

    
    %calculate period for driving force

    period = 2*(s_a(secondzero,1) - s_a(firstzero,1))
    
    %calculate phase lag (radians) from time lag (seconds)
    phi1 = (tl1 * (2*pi)) / period
    phi2 = (tl2 * (2*pi)) / period


end