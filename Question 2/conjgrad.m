function x = conjgrad(A, b, tol)
    dim = size( b); 
    if ( dim(2) > 1) 
       b = b'; 
    end
    if nargin < 3
        tol = 1e-10;
    end
    x = b;
    r = b - A * x;
    if norm(r) < tol
        return
    end
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
  
    for k = 1:numel(b);
       r = r - t*z;
       if( norm(r) < tol )
            return;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
    end
 end