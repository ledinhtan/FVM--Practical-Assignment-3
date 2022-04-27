function [A] = stiffnessMatrix(a,b,c,d,j,Nx,Ny)
A=zeros(Nx);
%Create matrix A1 
    if j == 1
        for i=1:Nx
           
            if i==1 
                A(1,1)= b(1)+d(j);
                A(1,2)= -b(1);
            elseif i== Nx
                A(Nx,Nx)  = a(Nx)+d(j);
                A(Nx,Nx-1)= -a(Nx);
            else
                A(i,i)  = b(i)+a(i)+d(j);
                A(i,i-1)= -a(i);
                A(i,i+1)= -b(i);
            end
        end
        
    elseif j ==  Ny
        for i=1:Nx
         
            if i==1 
                A(1,1)= b(1)+c(Ny);
                A(1,2)= -b(1);
            elseif i== Nx
                A(Nx,Nx)  = a(Nx)+c(Ny);
                A(Nx,Nx-1)= -a(Nx);
            else
                A(i,i)  = b(i)+c(j)+a(i);
                A(i,i-1)= -a(i);
                A(i,i+1)= -b(i);
            end
        end
    
    else
         for i=1:Nx
            
            if i==1 
                A(1,1)= b(1)+c(j)+d(j);
                A(1,2)= -b(1);
            elseif i== Nx
                A(Nx,Nx)  = a(Nx)+d(j)+c(j);
                A(Nx,Nx-1)= -a(Nx);
            else
                A(i,i)  = a(i)+b(i)+c(j)+d(j);
                A(i,i-1)= -a(i);
                A(i,i+1)= -b(i);
            end
        end
    end
end
