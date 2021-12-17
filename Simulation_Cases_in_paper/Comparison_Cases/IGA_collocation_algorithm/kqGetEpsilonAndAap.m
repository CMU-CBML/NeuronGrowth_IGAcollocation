function [epsilon, epsilon_deriv, aap, P_dy, P_dx] = kqGetEpsilonAndAap(epsilonb,delta,phi,xtheta,L_NuNv,U_NuNv,NuN1v,N1uNv)
% kqGet_epsilon_and_aap(phi,theta,NuN1v,N1uNv) 
% This function calculates epsilon and aap (a*a') based on phi, theta,
% NuN1v, and N1uNv.
% kqGet_epsilon_and_aap() output epsilon aap in the format of [epsilon,aap]

aniso = 6;

P_dy = NuN1v*phi;
P_dx = N1uNv*phi;

sz = sqrt(length(P_dx));

P_dy = reshape(P_dy,sz*sz,1);
P_dx = reshape(P_dx,sz*sz,1);

atheta =(full(atan2(P_dy,P_dx)));

epsilon = epsilonb.*(1.0+delta*cos(aniso*(atheta-xtheta)));
epsilon_deriv = -epsilonb.*(aniso*delta*sin(aniso.*(atheta-xtheta)));
aap = epsilon.*epsilon_deriv;

epsilon = U_NuNv\(L_NuNv\epsilon);
aap = U_NuNv\(L_NuNv\aap);

epsilon = reshape(epsilon,sz*sz,1);
aap = reshape(aap,sz*sz,1);