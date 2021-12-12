function [B] = Bernstein(i,n,u)
%Nurbs book algorithm A1.2
%Compute the value of a Bernstein polynomial. *I
%Input: i,n,u
%Output: B
temp = zeros([n,1]);
for j=0:n %compute the columns of Table 1.1 in a temporary array
    temp(j+1) = 0.0;
end
temp(n-i+1) = 1.0;
u1 = 1.0-u;
for k=1:n
    for j=n:k
        temp(j+1) = u1*temp(j+1) + u*temp(j-1+1);
    end
end
B = temp(n+1);
end
