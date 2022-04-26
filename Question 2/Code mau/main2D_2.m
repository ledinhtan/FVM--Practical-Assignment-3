% Solve 2D Laplace equation Neumann Boundary Condition
% luoi deu.
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
N=8; % Number of control volume
number_mesh=4; % number of iteration when refine mesh
number_mesh_point=zeros(number_mesh,1);
norm_l2=zeros( number_mesh,1); % norm l2;
norm_h1=zeros( number_mesh,1); % norrm h1

 %% Solve discrite solution and refine mesh
 for inumber_mesh=1:number_mesh
     number_mesh_point(inumber_mesh)=N;
     delta_x=(bx-ax)/N;
     delta_y=(by-ay)/N;
% Create the mesh point x(i+1/2)
    x=zeros(N+1,1);
    for i_iter=1:N+1
        x(i_iter)=ax+(i_iter-1)*delta_x; %luoi deu
    end
% Create the mesh point y(i+1/2)
    y=zeros(N+1,1);
    for j_iter=1:N+1
        y(j_iter)=ay+(j_iter-1)*delta_y; %luoi deu
    end    
% create control point x(i)
    x_cp=zeros(N+2,1);
    for i_iter=1:N+2
        if(i_iter==1)
            x_cp(i_iter)=x(i_iter);% voi  x_1/2=x_0
        else
            if(i_iter==N+2)
                x_cp(i_iter)=x(i_iter-1);% voi x_N1+1/2=x_N1+1
            else
                x_cp(i_iter)=(x(i_iter-1)+x(i_iter))/2.0; % x_cp(i_iter) la trung diem cua x(i_iter-1)va x(i_iter)
            end
        end
    end
    
% create control point y(i)
    y_cp=zeros(N+2,1);
    for j_iter=1:N+2
        if(j_iter==1)
            y_cp(j_iter)=y(j_iter);% voi  y_1/2=y_0
        else
            if(j_iter==N+2)
                y_cp(j_iter)=y(j_iter-1);% voi y_N1+1/2=y_N2+1
            else
                y_cp(j_iter)=(y(j_iter-1)+y(j_iter))/2.0; % y_cp(i_iter) la trung diem cua y(j_iter-1)va y(j_iter)
            end
        end
    end

 %% QUA TRINH TAO LAP MATRIX A.
   A=cell(N);
   A_1=sparse(N,N); % tao matran A_i voi i=1
   A_i=sparse(N,N); % tao matran A_i voi i chay tu 2 den (N_1-1)
   A_N1=sparse(N,N); % tao matran A_i voi i=N_1
    for i_iter=1:N
        for j_iter=1:N
            %% By setting.
            a_i=1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+1)-x_cp(i_iter)));
            b_i=1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+2)-x_cp(i_iter+1)));
            c_j=1/((y(j_iter+1)-y(j_iter))*(y_cp(j_iter+1)-y_cp(j_iter)));
            d_j=1/((y(j_iter+1)-y(j_iter))*(y_cp(j_iter+2)-y_cp(j_iter+1)));
            s_ij=a_i+b_i+c_j+d_j;
            
            %% Creare the Matrix C_i
            C_i=-eye(N)*a_i;  
            
            %% Creare the Matrix D_i
            D_i=-eye(N)*b_i;
            
            %% Creare the Matrix 
            if(j_iter==1)
                %% Creare the Matrix A_1
                A_1(j_iter,j_iter+1)=-d_j;
                A_1(j_iter,j_iter)= -a_i-c_j+s_ij;
                %% Creare the Matrix A_i
                A_i(j_iter,j_iter+1)=-d_j;
                A_i(j_iter,j_iter)= -c_j+s_ij;
                %% Creare the Matrix A_N1
                A_N1(j_iter,j_iter+1)=-d_j;
                A_N1(j_iter,j_iter)= -c_j-b_i+s_ij;
            else
                if(j_iter==N)
                    %% Creare the Matrix A_1
                    A_1(j_iter,j_iter-1)=-c_j;
                    A_1(j_iter,j_iter)= -a_i-d_j+s_ij;
                    %% Creare the Matrix A_i
                    A_i(j_iter,j_iter-1)=-c_j;
                    A_i(j_iter,j_iter)= -d_j+s_ij;
                    %% Creare the Matrix A_N1
                    A_N1(j_iter,j_iter-1)=-c_j;
                    A_N1(j_iter,j_iter)= -b_i-d_j+s_ij;
                else
                    %% Creare the Matrix A_1
                    A_1(j_iter,j_iter-1)=-c_j;
                    A_1(j_iter,j_iter+1)=-d_j;
                    A_1(j_iter,j_iter)= -a_i+s_ij;
                    %% Creare the Matrix A_i
                    A_i(j_iter,j_iter-1)=-c_j;
                    A_i(j_iter,j_iter+1)=-d_j;
                    A_i(j_iter,j_iter)= s_ij;
                    %% Creare the Matrix A_N1
                    A_N1(j_iter,j_iter-1)=-c_j;
                    A_N1(j_iter,j_iter+1)=-d_j;
                    A_N1(j_iter,j_iter)= -b_i+s_ij;
                end
            end
            
            %% Creare the Matrix A using cell2mat
            A{i_iter,j_iter}=sparse(N,N);
       
            if i_iter==1
                A{i_iter,i_iter}=A_1;
                A{i_iter,i_iter+1}=D_i;
            elseif i_iter==N
                A{i_iter,i_iter}=A_N1;
                A{i_iter,i_iter-1}=C_i;
            else
                A{i_iter,i_iter}=A_i;
                A{i_iter,i_iter-1}=C_i;
                A{i_iter,i_iter+1}=D_i;
            end
        end
    end
    %% dung cell2mat() de luu matran A duoi dang khoi.
    A=cell2mat(A);

 %% Create the Matrix b
    b=sparse(N);
    for i_iter=1:N 
        for j_iter=1:N
            b(i_iter,j_iter)=(functionf(x(i_iter),y(j_iter))+functionf(x(i_iter+1),y(j_iter))+functionf(x(i_iter),y(j_iter+1))...
                +functionf(x(i_iter+1),y(j_iter+1)))/4.0; % Trepozoidal rule 
        end
    end

 %% Create the Matrix F
    F=zeros((N)^2,1);
    r=1;
    t=1;
    while r<(N)^2
          F(r:r+N-1,1)=b(t,:);
          r=r+N;
          t=t+1;
    end
    
    
%% Solve discrete solution

  u=zeros(N*N,1);
  u=lsqr(A,F,10^-9,4*N);
  
  
%% exactly solution
    u_ex=sparse(N+2);
    for i_iter=1:N+2
        for j_iter=1:N+2
            u_ex(i_iter,j_iter)=exact_solution(x_cp(i_iter),y_cp(j_iter));
        end
    end
    
 %% Create discrete solution with boundary 
   u_dis=u_ex;
   %% cach chuyen vector sang matran.
   r_dis=2;
   t_dis=1;
    for i_iter=1:N+2
       for j_iter=1:N+2
           u_dis(i_iter,1)=u_dis(i_iter,2); % dieu kien bien: u_y(x,0)=0
           u_dis(i_iter,N+2)=u_dis(i_iter,N+1); % dieu kien bien: u_y(x,1)=0
           u_dis(1,j_iter)=u_dis(2,j_iter); % dieu kien bien : u_x(0,y)=0
           u_dis(N+2,j_iter)=u_dis(N+1,j_iter); % dieu kien bien : u_x(1,y)=0
       end
   end
   while r_dis < N+2
       u_dis(r_dis,2:N+1)=u(t_dis:t_dis+(N-1),1);
       r_dis=r_dis+1;
       t_dis=t_dis+N;
   end
 
  %%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i_iter=1:N
       for j_iter=1:N
            norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i_iter+1,j_iter+1)-u_ex(i_iter+1,j_iter+1))^2*(x(i_iter+1)-x(i_iter))*(y(j_iter+1)-y(j_iter));
       end
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    Norm_L2=norm_l2(inumber_mesh)
    
  %% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i_iter=1:N+1
       for j_iter=1:N+1
            norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+((u_dis(i_iter+1,j_iter+1)-u_ex(i_iter+1,j_iter+1))-(u_dis(i_iter,j_iter)-u_ex(i_iter,j_iter)))^2/...
                (x_cp(i_iter+1)-x_cp(i_iter))*(y_cp(j_iter+1)-y_cp(j_iter));
       end
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
    Norm_H1=norm_h1(inumber_mesh)
    
%% Figure exact and dicrete solutions 

     figure
    [X,Y]=meshgrid(x_cp,y_cp);
    subplot(1,2,1)
    Z_1=u_dis;
    surf(X,Y,Z_1)
    title('Discrete solution')
    subplot(1,2,2)
    Z_2=u_ex;
    surf(X,Y,Z_2)
    title('Exact solution')
 
%% Refine mesh (increse mesh point)   
    N=2*N;
 end
%% Figure for errors respect to number of mesh point

figure
plot(log(number_mesh_point),-log(norm_l2),'r', log(number_mesh_point), -log(norm_h1),'blue', log(number_mesh_point),1.5*log(number_mesh_point)+2, 'black', log(number_mesh_point), 2*log(number_mesh_point)+1.5,'green');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x','Location','NorthEastOutside');