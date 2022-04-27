function [A]=stiffnessMatrix(a,b,c,d,j,N1)
%---------------------------------------------------------
%Create the A matrix
%--------------------------------------------------------- 
for i=1:N1
    if i==1
        A(1,1)=a(1)+b(1)+c(j)+d(j);
        A(1,2)=-b(1);
    else if i==N1
            A(N1,N1-1)=-a(N1);
            A(N1,N1)=a(N1)+b(N1)+c(j)+d(j);
        else
            A(i,i)=a(i)+b(i)+c(j)+d(j);
            A(i,i-1)=-a(i);
            A(i,i+1)=-b(i);
        end
    end
end
end
