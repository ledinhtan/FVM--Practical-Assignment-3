function f=f(x,y,k)
if(k==1)
   f=-2*x*(x-1)-2*y*(y-1);
end
if(k==2)
   f=-4;
end
if(k==3)
    f=4*pi^2*sin(2*pi*(x-1))+20;
end