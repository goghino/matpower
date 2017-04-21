function [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam)
%D2SBUS_DV2   Computes 2nd derivatives of power injection w.r.t. voltage.
%   [GAA, GAV, GVA, GVV] = D2SBUS_DV2(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex bus power injections. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Gaa = (d/dVa (dSbus_dVa.')) * lam
%       Gav = (d/dVm (dSbus_dVa.')) * lam
%       Gva = (d/dVa (dSbus_dVm.')) * lam
%       Gvv = (d/dVm (dSbus_dVm.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

n = length(V);
Ibus    = Ybus * V;
diaglam = sparse(1:n, 1:n, lam, n, n);
diagV   = sparse(1:n, 1:n, V, n, n);

A = sparse(1:n, 1:n, lam .* V, n, n);
B = Ybus * diagV;
C = A * conj(B);
D = Ybus' * diagV;
E = conj(diagV) * (D * diaglam - sparse(1:n, 1:n, D*lam, n, n));
F = C - A * sparse(1:n, 1:n, conj(Ibus), n, n);
G = sparse(1:n, 1:n, ones(n, 1)./abs(V), n, n);

Gaa = E + F;
Gva = 1j * G * (E - F);
Gav = Gva.';
Gvv = G * (C + C.') * G;

%% alternative implementation of Gaa
% Gaa1 = sparse(n,n);
% for s = 1:n %for each eq constraint/row of Ybus
% for i=1:n 
%     for j=1:n
%         
%         tmp = 0;
%         
%         if(s == i && i == j)
%             tmp = tmp - (Ibus(s)); %1j*1j
%         end
%         
%         if (s == i)
%             tmp = tmp + (Ybus(s,j)*V(j));
%         end
%         
%         if(s == j)
%             tmp = tmp + (Ybus(s,i)*V(i));
%         end
%         
%         if (i == j)
%             tmp = tmp - (Ybus(s,j)*V(j)); %1j*1j
%         end
%         
%         
%         Gaa1(i,j) = Gaa1(i,j) + V(s)*lam(s)*conj(tmp);
% 
%     end
% end
% end
% 
% Gav1 = sparse(n,n);
% for s = 1:n %for each eq constraint/row of Ybus
% for i=1:n 
%     for j=1:n
%         
%         tmp = 0;
%         
%         if(s == i && i == j)
%             tmp = tmp +  (1/abs(V(s)) * Ibus(s));
%         end
%         
%         if (s == i)
%             tmp = tmp + (Ybus(s,j)*V(j)/abs(V(j)));
%         end
%         
%         if(s == j)
%             tmp = tmp - (1/abs(V(s)) * Ybus(s,i)*V(i));
%         end
%         
%         if (i == j)
%             tmp = tmp - (Ybus(s,j)*V(j)/abs(V(j)));
%         end
%         
%         
%         Gav1(i,j) = Gav1(i,j) + 1j*V(s)*lam(s)*conj(tmp);
% 
%     end
% end
% end
% 
% Gvv1 = sparse(n,n);
% for s = 1:n %for each eq constraint/row of Ybus
% for i=1:n 
%     for j=1:n
%         
%         tmp = 0;
%         
%         if (s == i)
%             tmp = tmp + (Ybus(s,j)*V(j)/abs(V(j)));
%         end
%         
%         if(s == j)
%             tmp = tmp + ( Ybus(s,i)*V(i)/abs(V(i)));
%         end
% 
%         Gvv1(i,j) = Gvv1(i,j) + V(s)/abs(V(s))*lam(s)*conj(tmp);
% 
%     end
% end
% end
