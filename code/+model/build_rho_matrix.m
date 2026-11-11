function A = build_rho_matrix(grid, dt, a, cS, cN, uS, uN)
%BUILD_RHO_MATRIX Matrix for rho/dt - d_x(a d_x rho) - d_x(cS rhoS) - d_x(cN rhoN).
% cS,cN,uS,uN are 1-by-(Nx+1) interface arrays. Upwind is chosen by uS,uN.
N = grid.Nx; dx = grid.dx;
rows = []; cols = []; vals = [];
% identity
rows = [rows, 1:N]; cols = [cols, 1:N]; vals = [vals, (1/dt)*ones(1,N)];
% interior interfaces e=2,...,N between cells L=e-1 and R=e
for e = 2:N
    Lc = e-1; Rc = e;
    % diffusion flux a*(rho_R-rho_L)/dx
    coeff = a/dx^2;
    rows = [rows, Lc, Lc, Rc, Rc];
    cols = [cols, Lc, Rc, Lc, Rc];
    vals = [vals, coeff, -coeff, -coeff, coeff];
    % S drift flux cS*rho_up, upwind by uS
    [rr,cc,vv] = local_adv_entries(Lc,Rc,cS(e),uS(e),dx);
    rows = [rows, rr]; cols = [cols, cc]; vals = [vals, vv];
    % N drift flux
    [rr,cc,vv] = local_adv_entries(Lc,Rc,cN(e),uN(e),dx);
    rows = [rows, rr]; cols = [cols, cc]; vals = [vals, vv];
end
A = sparse(rows, cols, vals, N, N);
end

function [rows,cols,vals] = local_adv_entries(Lc,Rc,c,u,dx)
% contribution of -d_x(c*rho_up) across one interface.
if u >= 0
    up = Lc;
else
    up = Rc;
end
rows = [Lc, Rc];
cols = [up, up];
vals = [-c/dx, c/dx];
end
