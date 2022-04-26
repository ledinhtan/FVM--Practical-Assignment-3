clear all
clc
close all
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
cases=1;
Nx=10;% Number of control volume
Ny= 8;
M=4;% number of iteration when refine mesh
norml2=zeros(M,1); % norm l2;
normh1=zeros(M,1); % norrm h1

ll=zeros(M,1);


for jj=1:M
    dx=(bx-ax)/Nx;
    dy=(by-ay)/Ny;
    
    % Create the mesh point
    x=zeros(Nx+1,1);
    for i=1:Nx+1
        x(i)=ax+(i-1)*dx;%% x_1/2..
    end
    
    y=zeros(Ny+1,1);
    for j=1:Ny+1
        y(j)=ay+(j-1)*dy;%% y_1/2..
    end
    
    % create control point
    x_cp=zeros(Nx+2,1);
    for i=1:Nx+2
        if(i==1)
            x_cp(i)=x(i);
        else
            if(i==Nx+2)
                x_cp(i)=x(i-1);
            else
                x_cp(i)=(x(i-1)+x(i))/2.0;
            end
        end
    end
    
     y_cp=zeros(Ny+2,1);
    for j=1:Ny+2
        if(j==1)
            y_cp(j)=y(j);
        else
            if(j==Ny+2)
                y_cp(j)=y(j-1);
            else
                y_cp(j)=(y(j-1)+y(j))/2.0;
            end
        end
    end
    
    cpoint=zeros((Nx+2)*(Ny+2),2);
    for j=1:Ny+2
        for i=1:Nx+2
            cpoint((j-1)*(Nx+2)+i,1)=x_cp(i);
            cpoint((j-1)*(Nx+2)+i,2)=y_cp(j);
        end
    end
    
    h=zeros(Nx,1);
    a=zeros(Nx,1);
    b=zeros(Nx,1);
    for i=1:Nx 
        h(i)=x(i+1)-x(i);
        a(i)= 1/(h(i)*(x_cp(i+1)-x_cp(i)));
        b(i)= 1/(h(i)*(x_cp(i+2)-x_cp(i+1)));

    end
   k=zeros(Ny,1);
    c=zeros(Ny,1);
    d=zeros(Ny,1);
    for j=1:Ny
        k(j)=y(j+1)-y(j);
        c(j)=1/(k(j)*(y_cp(j+1)-y_cp(j)));
        d(j)=1/(k(j)*(y_cp(j+2)-y_cp(j+1)));

   end
    % Creare the matrix A
  
   A=sparse(Nx*Ny); 
    for j=1:Ny
        [A1] = maTranLapGhep(a,b,c,d,j,Nx,Ny)
        if(j==1)    
            A((j-1)*Nx+1:j*Nx,(j-1)*Nx+1:j*Nx)=A1;
            A((j-1)*Nx+1:j*Nx,j*Nx+1:(j+1)*Nx)=-d(1)*eye(Nx);
        else
            if(j==Ny)
                A((j-1)*Nx+1:j*Nx,(j-1)*Nx+1:j*Nx)=A1;
                A((j-1)*Nx+1:j*Nx,(j-2)*Nx+1:(j-1)*Nx)=-c(Ny)*eye(Nx);
            else
                A((j-1)*Nx+1:j*Nx,(j-1)*Nx+1:j*Nx)=A1;
                A((j-1)*Nx+1:j*Nx,j*Nx+1:(j+1)*Nx)=-d(j)*eye(Nx);
                A((j-1)*Nx+1:j*Nx,(j-2)*Nx+1:(j-1)*Nx)=-c(j)*eye(Nx);
            end
        end
    end
  
    % Create vector b
  %% Tim vector F
    
      F=zeros(Nx*Ny,1);
    for j=1:Ny
        for i=1:Nx
            F((j-1)*Nx+i)=f(cpoint(j*(Nx+2)+i+1,1),cpoint(j*(Nx+2)+i+1,2),cases);
          
         end
    end
    B=zeros(Nx*Ny,1);
    sum=0;

    for j=1:Ny
        for i=1:Nx
            sum= sum+h(i)*k(j)*F((j-1)*Nx+i);            
        end
    end
   for j=1:Ny
       for i=1:Nx
         B((j-1)*Nx+i)= F((j-1)*Nx+i)-sum;  
       end
   end
    
      u = conjgrad(A,B,10^-6);
      
    % Create exact solution
    u_ex=zeros(Nx+2,Ny+2);
    for i=1:Nx+2
        for j=1:Ny+2
        u_ex(i,j)=u_exact(cpoint((j-1)*(Nx+2)+i,1),cpoint((j-1)*(Nx+2)+i,2),cases);
        end
    end
    
    %Create discrete solution

u_dis=u_ex;
    for j=1:Ny
        for i=1:Nx
            u_dis(i+1,j+1)=u((j-1)*Nx + i);
        end
    end
  
   figure
    subplot(121)
    surf(y_cp,x_cp,u_dis)
    title({['solution discrese khi Nx =',num2str(Nx), ' Ny =',num2str(Ny)]});
    subplot(122)
    surf(y_cp,x_cp,u_ex)
    title('solution exact')
        

%% Calculate the error on H1
    for i=1:Nx+1
        for j=1:Ny+1
            norml2(jj) = norml2(jj) + (u_dis(i,j) - u_ex(i,j))^2*(x_cp(i+1)-x(i))*(y_cp(j+1)-y(j));
        end
    end
    norml2(jj)=sqrt(norml2(jj));
    
%% Calculate the error on H1
    for i=1:Nx+1
        for j=1:Ny+1
             normh1(jj) = normh1(jj) + ((u_dis(i+1,j)-u_ex(i+1,j))-...
                           (u_dis(i,j)-u_ex(i,j)))^2/(x_cp(i+1)-x_cp(i))*(y_cp(j+1)-y_cp(j));
        end
    end
    normh1(jj)=sqrt(normh1(jj));
         %% Refine mesh
    
    Nx = 2*Nx;
    Ny =  2*Ny;
    ll(jj)= Ny;
end
figure
plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll),1.5*log(ll)+2, 'black', log(ll), 2*log(ll)+1,'g');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')

   