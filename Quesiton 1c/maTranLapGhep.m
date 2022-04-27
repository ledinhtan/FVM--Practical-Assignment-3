function [Ai,Ci,Di]=stiffnessMatrix(a,b,c,d,j,N1)
Ai=zeros(N1);
Ci =zeros(N1);
Di =zeros(N1);
%---------------------------------------------------------
%Tao ma tran Ai
%--------------------------------------------------------- 
for i=1:N1
    if i==1
        Ai(1,1)=a(i)+b(i)+c(j)+d(j);
        Ai(1,2)=-b(1);
    else if i==N1
            Ai(N1,N1)=a(i)+b(i)+c(j)+d(j);
            Ai(N1,N1-1)=-a(N1);
        else
            Ai(i,i)=a(i)+b(i)+c(j)+d(j);
            Ai(i,i-1)=-a(i);
            Ai(i,i+1)=-b(i);
        end
    end
    %---------------------------------------------------------
    %Tao ma tran Bi va Ci
    %--------------------------------------------------------- 
Ci(i,i) = -c(j);
Di(i,i) = -d(j); 
end





