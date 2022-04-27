clear all
clc
close all
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
cases=2;
Nx=8;% Number of control volume
Ny= 10;
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
    
    
    a=zeros(Nx,1);
    b=zeros(Nx,1);
    for i=1:Nx 
        h=x(i+1)-x(i);
        a(i)= 1/(h*(x_cp(i+1)-x_cp(i)));
        b(i)= 1/(h*(x_cp(i+2)-x_cp(i+1)));

    end
   
    c=zeros(Ny,1);
    d=zeros(Ny,1);
    for j=1:Ny
        k=y(j+1)-y(j);
        c(j)=1/(k*(y_cp(j+1)-y_cp(j)));
        d(j)=1/(k*(y_cp(j+2)-y_cp(j+1)));

   end
    % Creare the matrix A
    A=sparse(Nx*Ny);
    
    A(1:Nx,1:Nx)=stiffnessMatrix(a,b,c,d,1,Nx);
    A(1:Nx,Nx+1:2*Nx)= -d(1)*eye(Nx);
  
    A((Ny-1)*Nx+1:Ny*Nx,(Ny-1)*Nx+1:Ny*Nx)=maTranLapGhep(a,b,c,d,Ny,Nx);
    A((Ny-1)*Nx+1:Ny*Nx,(Ny-2)*Nx+1:(Ny-1)*Nx)= -c(Ny)*eye(Nx);
    
    for j=2:Ny-1
        A((j-1)*Nx+1:j*Nx,(j-1)*Nx+1:j*Nx)=maTranLapGhep(a,b,c,d,j,Nx);
        A((j-1)*Nx+1:j*Nx,j*Nx+1:(j+1)*Nx)= -d(j)*eye(Nx);
        A((j-1)*Nx+1:j*Nx,(j-2)*Nx+1:(j-1)*Nx)= -c(j)*eye(Nx);
    end
    
    % Create vector b
  %% Tim vector F

     F=zeros(Nx*Ny,1);       
    for j=1:Ny
        for i=1:Nx
            s = 0;
            if (j==1)
                s = c(1)*u_exact(cpoint(i+1,1),cpoint(i+1,2),cases);
                if (i==1)
                    s = s + a(1)*u_exact(cpoint(Nx+3,1),cpoint(Nx+3,2),cases);
                end
                if (i==Nx)
                    s = s + b(Nx)*u_exact(cpoint(2*(Nx+2),1),cpoint(2*(Nx+2),2),cases);
                end
            end
            
            if (j==Ny)
                temp = (j+1)*(Nx+2)+1;
                s = d(Ny)*u_exact(cpoint(temp+i,1),cpoint(temp+i,2),cases);
                if (i==1)
                    s = s + a(1)*u_exact(cpoint(temp-Nx-2,1),cpoint(temp-Nx-2,2),cases);
                end
                if (i==Nx)
                    s = s + b(Nx)*u_exact(cpoint(temp-1,1),cpoint(temp-1,2),cases);
                end
            end
            
            if(j~=1 && j ~= Ny)
                if (i==1)
                    s = a(1)*u_exact(ax,y_cp(j+1),cases);
                end
                if (i==Nx)
                    s = b(Nx)*u_exact(bx,y_cp(j+1),cases);
                end
            end
            F((j-1)*Nx+i)=s + f(cpoint(j*(Nx+2)+i+1,1),cpoint(j*(Nx+2)+i+1,2),cases);
         end
    end
    u=A\F;
      
    % Create exact solution
    u_ex=zeros(Nx+2,Ny+2);
    for i=1:Nx+2
        for j=1:Ny+2
        u_ex(i,j)=u_exact(cpoint((j-1)*(Nx+2)+i,1),cpoint((j-1)*(Nx+2)+i,2),cases);
        end
    end
    
    %Create discrete solution
    u_dis=zeros(Nx+2,Ny+2);
    k=1;
       for j=1:Ny+2
           for i=1:Nx+2
               if (i==1||i==Nx+2)
                   u_dis(i,j)=u_ex(i,j);
             
               elseif(j==1||j==Ny+2)
                   u_dis(i,j)=u_ex(i,j);
               else
                   u_dis(i,j)=u(k);
                   k=k+1;
               end
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
plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll),1.5*log(ll)+3, 'black', log(ll), 2*log(ll)+3,'green');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')
