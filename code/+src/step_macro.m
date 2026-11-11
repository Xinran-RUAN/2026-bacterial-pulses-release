function Unew = step_macro(U, grid, p)
%STEP_MACRO One implicit upwind step for the limiting macroscopic model.
rho = U.rho; S = U.S; N = U.N; dt = p.dt;
R = model.compute_response(rho, S, N, grid, p);
D = utils.vint(grid.v.^2, grid.w)/grid.V;
bS = R.BS/grid.V; bN = R.BN/grid.V;
uS = -bS; uN = -bN;
% The limiting scheme is rho_t = D rho_xx - d_x(uS rhoS) - d_x(uN rhoN)
% In the rho-matrix form use cS=bS=-uS and cN=bN=-uN.
cS = bS; cN = bN; cS([1,end]) = 0; cN([1,end]) = 0;
A = model.build_rho_matrix(grid, dt, D, cS, cN, uS, uN);
rhs = rho(:)/dt;
rho_new = (A\rhs).';
if isfield(p,'macro_limiter') && p.macro_limiter
    rho_new = max(rho_new,0);
    m0 = sum(rho)*grid.dx; m1 = sum(rho_new)*grid.dx;
    if m1 > 0, rho_new = rho_new*(m0/m1); end
end
[Snew,Nnew] = src.step_chemicals(S, N, rho_new, grid, p);
Unew = struct('rho',rho_new,'S',Snew,'N',Nnew,'g',[],'t',U.t+dt);
end
