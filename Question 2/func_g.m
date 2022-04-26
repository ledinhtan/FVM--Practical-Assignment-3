function [g]= func_g(x,y,k)
if(k==2)
    g= 2*x*y^2 - x*(x^2 - y^2) + 3*x^2*y;
end

   if(k==3)
    g= 0;
end
 