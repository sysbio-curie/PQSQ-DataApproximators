function drawPotential(intervals, potential_function_handle, quadratic)

R = max(intervals);
mx = R*1.5;

xi = 0;
yi = 0;


for x=0:mx/100:mx
   
    k=1;
    
    a = 1;
    b = 0;
    
    %f(xi)
    
    for i=1:size(intervals,2)-1
        if((intervals(i)<=x)&(x<intervals(i+1)))
                xk = intervals(i);
                xk1 = intervals(i+1);
                a = (potential_function_handle(xk)-potential_function_handle(xk1))/(xk*xk-xk1*xk1);
                b = (potential_function_handle(xk1)*xk*xk-potential_function_handle(xk)*xk1*xk1)/(xk*xk-xk1*xk1);
         end
    end
    
    if(x>intervals(size(intervals,2)))
        a = 0;
        b = potential_function_handle(R);
    end
    
    xi1 = x;
    yi1 = b+a*x*x;
    
    fxi = potential_function_handle(xi);
    fxi1 = potential_function_handle(xi1);
    
    if(quadratic)
        plot([xi*xi xi1*xi1],[yi yi1],'b-','LineWidth',4); hold on;
    else
        plot([xi xi1],[yi yi1],'b-','LineWidth',4); hold on;
    end
    
    xi = xi1;
    yi = yi1;
    
end

xi = 0;
yi = 0;
    for i=1:size(intervals,2)-1
                xk = intervals(i);
                xk1 = intervals(i+1);
                a = (potential_function_handle(xk)-potential_function_handle(xk1))/(xk*xk-xk1*xk1);
                b = (potential_function_handle(xk1)*xk*xk-potential_function_handle(xk)*xk1*xk1)/(xk*xk-xk1*xk1);
                for x=0:mx/100:mx
                    xi1 = x;
                    yi1 = b+a*x*x;
                
                if(yi1<potential_function_handle(R*1.2))
                if(quadratic)
                    plot([xi*xi xi1*xi1],[yi yi1],'b.','LineWidth',1); hold on;
                else
                    plot([xi xi1],[yi yi1],'b.','LineWidth',1); hold on;
                end
                    xi = xi1;
                    yi = yi1;
                end
                end
    end

xi = 0;
fxi = 0;
for x=0:mx/10:mx
    fxi1 = potential_function_handle(x);
    xi1 = x;
if(xi1<R*1.3)
if(quadratic)    
    plot([xi*xi xi1*xi1],[fxi fxi1],'r--','LineWidth',4); hold on;
else
    plot([xi xi1],[fxi fxi1],'r--','LineWidth',4); hold on;
end
    xi = x;
    fxi = fxi1;
end 
end


    for i=1:size(intervals,2)-1
        if(quadratic)    
        plot([intervals(i+1)*intervals(i+1) intervals(i+1)*intervals(i+1)],[-potential_function_handle(R*1.2)/20 potential_function_handle(R*1.2)],'k-'); 
        else
        plot([intervals(i+1) intervals(i+1)],[-potential_function_handle(R*1.2)/20 potential_function_handle(R*1.2)],'k-');
        end
    end

axis off;

diff = 0;
for x=0:mx/100:mx
    if(x < R)
    delta = (PQSQ(x,intervals,potential_function_handle)-potential_function_handle(x))/(potential_function_handle(x)+0.001);
    %display(sprintf('Delta = %f',delta));
    diff = diff + abs(delta);
    end
end

display(sprintf('Mean difference = %f',diff/100));


end


