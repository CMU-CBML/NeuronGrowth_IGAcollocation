% %  This Function evaluates the basis functions and first derivatives 
% %               functions at a given parameter value u.
% % 
% %  Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: 
% %     Berlin 1995; pp. 72-73.
% % 
% %  June 17, 2003
% %  J. Austin Cottrell
% %  CES Graduate Student
% %  Texas Institute for Computational Engineering Science
% %  University of Texas at Austin
% %
% %  Modified to codes Matlab by :
% %  Chien Thai Hoang & Hung Nguyen Xuan
% %
% %   Faculty of Mathematics & Informatics, University of Natural Sciences
% %   Vietnam   National University?HCM
% 
% % 
function [ders]=dersbasisfuns3_modified(i,p,order_deriv,u,u_knot)
% 
% %     --------------variable declarations--------------------------------
% % Input: (modified by Kuanren - 02172021)
% %     i -> knot span
% %     p -> degree of curve
% %     order_deriv -> order of derivative
% %     u -> evaluation value
% %     u_knot -> knot vector
% 
% % Out:
% %     ders(1,:)  -> values of the basis function N_(i-p) ... N_(i) at the evaluation point u
% %     ders(2,:)  -> first derivatives
% %     ders(3,:)  -> second derivatives
% 
% 
ders=zeros(order_deriv+1,p+1);
N=zeros(p+1,p+1);
N(1,1) = 1;
a=zeros(order_deriv+1,p+1);
left=zeros(1,p+1);
right=zeros(1,p+1);

%     -------------------------------------------------------------------

for j = 1:p    
   left(j+1) = u - u_knot(i+1-j+1); %updated index - Kuanren - 02172021
   right(j+1) = u_knot(i+j+1) - u;
   saved = 0;
   for r = 0:j-1
       N(j+1,r+1) = right(r+2) + left(j-r+1);
       temp = N(r+1,j)/N(j+1,r+1);%%%%% error here -> check line 43&44
       N(r+1,j+1) = saved + right(r+2)*temp;
       saved = left(j-r+1)*temp;
   end 
   N(j+1,j+1) = saved;
end 

% load basis functions
for j = 0:p
    ders(1,j+1) = N(j+1,p+1);
end 

% compute derivatives
for r = 0:p % loop over function index
    s1 = 0;
    % alternate rows in array a
    s2 = 1;
    a(1,1) = 1;

    for k = 1:order_deriv % loop to compute kth derivative
        d = 0;
        rk = r-k;
        pk = p-k;
        
        if(r >= k)
            a(s2+1,1) = a(s1+1,1)./N(pk+2,rk+1);
            d = a(s2+1,1).*N(rk+1,pk+1);        
        end
        
        if(rk >= -1)
            j1 = 1;
        else
            j1 = fix(-rk);
        end
        
        if((r-1) <= pk)
            j2 = fix(k-1);
        else
            j2 = fix(p-r);
        end
        
        for j = j1:j2;
            a(s2+1,j+1) =(a(s1+1,j+1) - a(s1+1,j))./N(pk+2,rk+j+1);
            d = d + a(s2+1,j+1).*N(rk+j+1,pk+1);
        end % j =fix(j2+1);
        
        if(r <= pk)
            a(s2+1,k+1) = -a(s1+1,k)./N(pk+2,r+1);
            d = d + a(s2+1,k+1).*N(r+1,pk+1);
        end
        
        ders(k+1,r+1) = d;
        j = s1;
        s1 = s2;
        % switch rows
        s2 = j;
    end 
end    

%   multiply through by the correct factors;
r = p;
for k = 1:order_deriv
    for j = 0:p
        ders(k+1,j+1) = ders(k+1,j+1).*r;
    end
    r = fix(r.*(p-k));
end



