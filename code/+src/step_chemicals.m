function [Snew,Nnew] = step_chemicals(S, N, rho, grid, p)
Ncells = grid.Nx; dt = p.dt;
L = utils.lap_neumann_matrix(Ncells, grid.dx); % -d_xx
I = speye(Ncells);
AS = I/dt + p.DS*L + p.alpha*I;
rhsS = S(:)/dt + p.beta*rho(:);
Snew = (AS\rhsS).';
AN = I/dt + p.DN*L + p.gamma*spdiags(rho(:),0,Ncells,Ncells);
rhsN = N(:)/dt;
Nnew = (AN\rhsN).';
end
