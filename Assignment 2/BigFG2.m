function [BigF,BigG] = BigFG2(h,tau,A,Q,B)
%H Compute vertices of polytopic overapproximation
%(sampling time,maximum delay,system mat.,eigenvector mat.)
%System parameters
l=eig(A);
l1=l(1);l2=l(2);%eigenvalues
Lh=[exp(l1*h),0;0,exp(l2*h)];
% Ass 2
alpha_1 = expm(l2*(h-tau));
alpha_2 = expm(l1*(h-tau));

Q_inv = inv(Q);

F0_part = Q_inv*[exp(h), 0; 0, 1/3.7*exp(4.3*h)]*Q*B;
F1_part = Q_inv*[-1, 0; 0, 0]*Q*B;
F2_part = Q_inv*[0, 0; 0, -1/3.7]*Q*B;

F0 = [expm(A*h), F0_part; zeros(1,3)];
F1 = [zeros(2), F1_part; zeros(1,3)];
F2 = [zeros(2), F2_part; zeros(1,3)];

BigF = F0 + alpha_1*F1 + alpha_2*F2;

G0 = [Q_inv*[-1, 0; 0, -1/3.7]*Q*B; 1];
G1 = [Q_inv*[1, 0; 0, 0]*Q*B; 0];
G2 = [Q_inv*[0, 0; 0, 1/3.7]*Q*B;0];

BigG = G0 + alpha_1*G1 + alpha_2*G2;
end

