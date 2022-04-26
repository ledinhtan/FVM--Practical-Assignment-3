
function uex=u_exact(x,y,k)
if(k==1)
    uex= cos(2*pi*x)*cos(2*pi*y);
end
if(k==2)
    uex=cos(2*pi*x)*cos(5*pi*y);
end
if(k==3)
    uex = cos(4*x*pi)+cos(4*y*pi);
end