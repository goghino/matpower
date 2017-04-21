function [Haa, Hav, Hva, Hvv] = d2Sbr_dV2(Cbr, Ybr, V, lam)
%D2SBR_DV2   Computes 2nd derivatives of complex power flow w.r.t. voltage.
%   [HAA, HAV, HVA, HVV] = D2SBR_DV2(CBR, YBR, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex branch power flows. Takes sparse connection matrix CBR, sparse
%   branch admittance matrix YBR, voltage vector V and nl x 1 vector of
%   multipliers LAM. Output matrices are sparse.
%
%   Example:
%       f = branch(:, F_BUS);
%       Cf =  sparse(1:nl, f, ones(nl, 1), nl, nb);
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       Cbr = Cf;
%       Ybr = Yf;
%       [Haa, Hav, Hva, Hvv] = d2Sbr_dV2(Cbr, Ybr, V, lam);
%
%   Here the output matrices correspond to:
%       Haa = (d/dVa (dSbr_dVa.')) * lam
%       Hav = (d/dVm (dSbr_dVa.')) * lam
%       Hva = (d/dVa (dSbr_dVm.')) * lam
%       Hvv = (d/dVm (dSbr_dVm.')) * lam
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
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

%% define
nl = length(lam);
nb = length(V);

diaglam = sparse(1:nl, 1:nl, lam, nl, nl);
diagV   = sparse(1:nb, 1:nb, V, nb, nb);

A = Ybr' * diaglam * Cbr;
B = conj(diagV) * A * diagV;
D = sparse(1:nb, 1:nb, (A*V) .* conj(V), nb, nb);
E = sparse(1:nb, 1:nb, (A.'*conj(V)) .* V, nb, nb);
F = B + B.';
G = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

Haa = F - D - E;
Hva = 1j * G * (B - B.' - D + E);
Hav = Hva.';
Hvv = G * F * G;

%% alternative representation of the hessians
% Vf = Cbr * V;
% If = Ybr * V;
% Haa1 = sparse(nb,nb);
% for row = 1:nl %1:nl
% for i = 1:nb
% for j = 1:nb
%     
%     tmp = 0;
%     
%     if (Cbr(row,i) == 1 && i == j)
%         %tmp = tmp + lam(i)*1j*1j*Vf(i)* conj(If(i));
%         tmp = tmp - (If(row));
%     end
%     
%     if (Cbr(row,i) == 1)
%         %tmp = tmp + lam(i)*1j*Vf(i)* conj(Ybr(i,k)*V(k)*1j);
%         tmp = tmp +  (Ybr(row,j)*V(j));
%     end
%     
%     if (Cbr(row,j) == 1)
%         %tmp = tmp + lam(i)*1j*Vf(i)* conj(Ybr(i,j)*V(j)*1j);
%         tmp = tmp +  (Ybr(row,i)*V(i));
%     end
%     
%     if (i == j)
%         %tmp = tmp + lam(i)*Vf(i)* conj(Ybr(i,k)*V(k)*1j*1j);
%         tmp = tmp -  (Ybr(row,j)*V(j));
%     end
% 
%     %Haa1(j,k) = Haa1(j,k) + tmp;
%     Haa1(i,j) = Haa1(i,j) + lam(row) * Vf(row) * conj(tmp);
% end
% end
% end
% 
% Hav1 = sparse(nb,nb);
% for row = 1:nl %1:nl
% for i = 1:nb
% for j = 1:nb
%     
%     tmp = 0;
%     
%     if (Cbr(row,i) == 1 && i == j)
%         tmp = tmp + (If(row)/abs(Vf(row)));
%     end
%     
%     if (Cbr(row,i) == 1)
%         tmp = tmp +  (Ybr(row,j)*V(j)/abs(V(j)));
%     end
%     
%     if (Cbr(row,j) == 1)
%         tmp = tmp - (Ybr(row,i)*V(i)/abs(Vf(row)));
%     end
%     
%     if (i == j)
%         tmp = tmp -  (Ybr(row,j)*V(j)/abs(V(j)));
%     end
% 
%     
%     Hav1(i,j) = Hav1(i,j) + Vf(row)*1j*lam(row)*conj(tmp);
% end
% end
% end
% 
% Hvv1 = sparse(nb,nb);
% for row = 1:nl %1:nl
% for i = 1:nb
% for j = 1:nb
%     
%     tmp = 0;
%     
%     if (Cbr(row,i) == 1)
%         tmp = tmp +  (Ybr(row,j)*V(j)/abs(V(j)));
%     end
%     
%     if (Cbr(row,j) == 1)
%         tmp = tmp +  (Ybr(row,i)*V(i)/abs(V(i)));
%     end
% 
%     
%     Hvv1(i,j) = Hvv1(i,j) + lam(row)*Vf(row)/abs(Vf(row))* conj(tmp);
% end
% end
% end
