function f1=f(x,y,k)
if(k==1)
   f1 = 4*pi^2*cos(2*pi*x)*cos(2*pi*y)+4*pi^2*cos(2*pi*x)*cos(2*pi*y);
end
if(k==2)
    f1 = 4*pi^2*cos(2*pi*x)*cos(5*pi*y)+25*pi^2*cos(2*pi*x)*cos(5*pi*y);
end
if(k==3)
    f1 = 16*pi^2*cos(4*x*pi)+16*pi^2*cos(4*y*pi);
end