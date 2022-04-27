clear all
clc
close all
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;

N1=10;% Number of control volume x
N2=10;% Number of control volume y

M=4;% number of iteration when refine mesh
norml2=zeros(M,1); % norm l2
normh1=zeros(M,1); % norrm h1

ll=zeros(M,1);

for jj=1:M
    dx=(bx-ax)/N1;
    
    % Create the mesh point
    x=zeros(N1+1,1);
    for i_iter=1:N1+1
        x(i_iter)=ax+(i_iter-1)*dx;
    end
    
    % create control point
    
    x_cp=zeros(N1+2,1);
    for i_iter=1:N1+2
        if(i_iter==1)
            x_cp(i_iter)=x(i_iter);
        else
            if(i_iter==N1+2)
                x_cp(i_iter)=x(i_iter-1);
            else
                x_cp(i_iter)=(x(i_iter-1)+x(i_iter))/2.0;
            end
        end
    end
    
    dy=(by-ay)/N2;
 
    % Create the mesh point
    y=zeros(N2+1,1);
    for i_iter=1:N2+1
        y(i_iter)=ay+(i_iter-1)*dy;
    end
    
    % create control point
    y_cp=zeros(N2+2,1);
    for i_iter=1:N2+2
        if(i_iter==1)
            y_cp(i_iter)=y(i_iter);
        else
            if(i_iter==N2+2)
                y_cp(i_iter)=y(i_iter-1);
            else
                y_cp(i_iter)=(y(i_iter-1)+y(i_iter))/2.0;
            end
        end
    end
    
    a=zeros(N1,1);
    for i_iter=1:N1
        a(i_iter)=1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+1)-x_cp(i_iter)));
    end
    
    b=zeros(N1,1);
    for i_iter=1:N1
        b(i_iter)=1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+2)-x_cp(i_iter+1)));
    end
    
    c=zeros(N2,1);
    for i_iter=1:N2
        c(i_iter)=1/((y(i_iter+1)-y(i_iter))*(y_cp(i_iter+1)-y_cp(i_iter)));
    end
    
    d=zeros(N2,1);
    for i_iter=1:N2
        d(i_iter)=1/((y(i_iter+1)-y(i_iter))*(y_cp(i_iter+2)-y_cp(i_iter+1)));
    end
    
    % Creare the Matrix
    A=sparse(N1*N2);
    for i=1:N2
        [Ai,Ci,Di]=maTranLapGhep(a,b,c,d,i,N1);
        if i==1
            A((i-1)*N1+1:N1,i*1:N1)=Ai;
            A((i-1)*N1+1:N1,(N1 +1):2*N1)=Di;
        elseif i==N2
            A(((i-1)*N1+1):i*N1,((i-1)*N1 +1):N2*N1)=Ai;
            A(((i-1)*N1+1):i*N1,((i-2)*N1 +1):(i-1)*N1 )=Ci;
          
        else
            A((i-1)*N1+1:(i-1)*N1+ N1,(i-1)*N1-N1+1:(i-1)*N1)=Ci;
            A((i-1)*N1+1:(i-1)*N1+ N1,(i-1)*N1+1:i*N1)=Ai;
            A((i-1)*N1+1:(i-1)*N1+ N1,i*N1+1:(i+1)*N1)=Di;
            
            
        end
    end

    
    % Create vector b
    F=zeros(N1*N2,1);     
    for i=1:N2
        s1=zeros(N1,1);
        for j=1:N1
            s1(j)=f(x_cp(j+1),y_cp(i+1));
        end
        
        F((N1*i-N1 +1):N1*i)=s1;
    end
   
    u=zeros(N1*N2,1);
    u=A\F;
    
    u_ex=zeros(N2+2,N1+2);
    for i=1:N1+2
        for j=1:N2+2
            u_ex(j,i)=u_exact(x_cp(i),y_cp(j));
        end
    end
    
    u_dis=zeros(N2+2,N1+2);

    for i=2:N2+1
        for j=2:N1+1
            u_dis(i,j)=u((i-2)*N1+j-1);
        end
    end
    figure
    subplot(1,2,1)
    surf(x_cp,y_cp,u_dis)
    title('u_dis')
    subplot(1,2,2)
    surf(x_cp,y_cp,u_ex)
    title(['u_ex, with N1=','N2 =',num2str(N1,N2)])
    
    for i_iter=1:N2
        for j_jter=1:N1
            norml2(jj)=norml2(jj)+(u_dis(i_iter+1,j_jter+1)-...
                       u_ex(i_iter+1,j_jter+1))^2*...
                       (x(j_jter+1)-x(j_jter))*(y(i_iter+1)-y(i_iter));
        end
    end
    norml2(jj)=sqrt(norml2(jj));
    
    for i_iter=1:N2+1
        for j_jter=1:N1+1
            normh1(jj)=normh1(jj)+((u_dis(i_iter+1,j_jter)-...
                       u_ex(i_iter+1,j_jter))-(u_dis(i_iter,j_jter)-...
                       u_ex(i_iter,j_jter)))^2/(x_cp(j_jter+1)-...
                       x_cp(j_jter))*(y_cp(i_iter+1)-y_cp(i_iter));
        end
    end

    normh1(jj)=sqrt(normh1(jj)); 
    
    ll(jj)=N1;
    N1=2*N1;
    N2=2*N2;
     
end

figure
plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll),1.5*log(ll)+2, 'black', log(ll), 2*log(ll)+1.5,'green');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')
