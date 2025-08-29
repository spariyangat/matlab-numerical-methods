function r = PariyangatS_Lab4(flag,handles,vector,eps_s)
%input: a method flag, handles needed for desired method, vector of 
% necessary starting values, convergence criteria epsilon as a decimal
%output: root of the passed in function (handle) to the user-specified
% level of convergence

    switch flag
        case 1 %False Position Method
            xl = vector(1);
            xu = vector(2);

            %run algorithm once before the while loop
            xr = (((-handles(xl))/(handles(xu)-handles(xl)))*(xu-xl))+xl;
            fxr = handles(xr);
            if((fxr * handles(xu)) < 0)
                xl = xr;
            else
                xu = xr;
            end
            
            %starter flag technique to enter the while loop
            first = 1;
            while (first == 1) || (abs((xr - xrprev)/(xr)) > eps_s)
                xrprev = xr;
                xr = (((-handles(xl))/(handles(xu)-handles(xl)))*(xu-xl))+xl;
                fxr = handles(xr);
                %replace lower or upper bound with xr
                if((fxr * handles(xu)) < 0)
                    xl = xr;
                else
                    xu = xr;
                end
                first = 0;
            end

        case 2 %Secant Method  
            xrprev = vector(1); 
            xr = vector(2);
            while (abs((xr - xrprev)/(xr)) > eps_s)
                temp = xr; 
                %applying the secant method equation to calculate new guess
                xr = (((-handles(xrprev))/(handles(xr)-handles(xrprev)))*(xr-xrprev))+xrprev;
                xrprev = temp; 
            end

        case 3 %Muller's Method
            xr2prev = vector(1);
            xr1prev = vector(2);
            xr = vector(3);
            while(abs((xr - xr1prev)/(xr)) > eps_s)
                %calculate a,b,c values using Muller's method equations 
                a = (((handles(xr)-handles(xr1prev))/(xr-xr1prev)) - ...
                ((handles(xr1prev)-handles(xr2prev))/(xr1prev-xr2prev))) / (xr - xr2prev);    
                b = (a*(xr-xr2prev)) + (handles(xr)-handles(xr2prev))/(xr-xr2prev);
                c = handles(xr);

                %calculate new guess using the polynomial
                temp = xr;
                xr = xr - ((2*c)/(b-sqrt((b^2)-(4*a*c))));   
                xr2prev = xr1prev;
                xr1prev = temp;
            end
            
        case 4 %Multi-dimensional Modified Secant Method
            delta = 0.01;
            roots = vector';

            %gateway flag technique
            converged = 0; 
            while(~converged)
                %create function vector f with handles
                for idx = 1:length(handles)
                    f(idx) = handles{idx}(roots); 
                end

                %calculate partials and create partials matrix
                for count1 = 1:length(handles)
                    for count2 = 1:length(roots)
                        modified_vals = roots;
                        modified_vals(count2) = modified_vals(count2)*(1+delta);
                        a = handles{count1}(modified_vals);
                        b = handles{count1}(roots);
                        c = (roots(count2))*(delta);
                        partials(count1,count2) = (a-b)/c; 
                    end
                end  

                %calculate guess and keep track of previous guesses
                prevroots = roots; 
                roots = (roots) - (partials)\(f');   
        
                %gateway flag technique
                counter = 1;
                flag2 = 1;                
                while(flag2 && counter <= length(roots))
                    eps_a = abs((roots(counter)-prevroots(counter))/roots(counter));
                    %check for convergence
                    if(eps_a > eps_s)
                        flag2 = 0; 
                    elseif((eps_a <= eps_s) && (counter == length(roots)))
                        flag2 = 0;
                        converged = 1;
                    end
                    counter = counter+1;       
                end
            end
            xr = roots;

    end
    
    r = xr;
   
end