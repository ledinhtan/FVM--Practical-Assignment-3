function [h]= func_h(x,y,k)
if(k==2)
    h= - (x*(x^2 - 6*y^2))/2 - y*(x^2 - y^2) - 2*x^2*y - x^3/2;
end
if(k==3)
    h=0;
end