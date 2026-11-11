function R = compute_response(rho, S, N, grid, p)
%COMPUTE_RESPONSE Centered psi_S, psi_N and effective velocities at interfaces.
%
% The response is evaluated as
%
%     phi(eps_response*q_t + v*q_x),
% where eps_response is set by run_paper_case:
%     macro   : eps_response = 0,
%     kinetic : eps_response = p.eps.
% This keeps the macro solver consistent with the analytical epsilon->0
% travelling-wave formula, while the kinetic solver still uses the finite
% epsilon response.

v = grid.v;
w = grid.w;
V = grid.V;
Nx = grid.Nx;
Nv = grid.Nv;
dx = grid.dx;

% Cell-to-edge values and derivatives.
rho_e = utils.cell_to_edge(rho);
S_e   = utils.cell_to_edge(S);
N_e   = utils.cell_to_edge(N);
Sx_e  = utils.grad_cell_to_edge(S, dx);
Nx_e  = utils.grad_cell_to_edge(N, dx);

% Approximate time derivatives at cell centers from the PDE, then average to edges.
L = utils.lap_neumann_matrix(Nx, dx);   % -d_xx
Sxx = -(L*S(:)).';
Nxx = -(L*N(:)).';
St_cell = p.DS*Sxx - p.alpha*S + p.beta*rho;
Nt_cell = p.DN*Nxx - p.gamma*rho.*N;
St_e = utils.cell_to_edge(St_cell);
Nt_e = utils.cell_to_edge(Nt_cell);

% Epsilon multiplying the time derivative in the response.  If absent,
% use p.eps for backward compatibility; run_paper_case sets it explicitly.
if isfield(p, 'response_eps') && ~isempty(p.response_eps)
    eps_response = p.response_eps;
elseif isfield(p, 'eps') && ~isempty(p.eps)
    eps_response = p.eps;
else
    eps_response = 0;
end

phi_scale = model.response_scale(grid, p);

DtS = eps_response*St_e + v*Sx_e;   % Nv-by-(Nx+1), implicit expansion
DtN = eps_response*Nt_e + v*Nx_e;
rawS = p.chiS * phi_scale * model.phi_response(DtS, p);
rawN = p.chiN * phi_scale * model.phi_response(DtN, p);

% Centering in velocity keeps <psi>=0.  It does not change the first moment
% because the velocity grid is symmetric and <v>=0.
meanS = utils.vint(rawS, w)/V;
meanN = utils.vint(rawN, w)/V;
psiS = rawS - ones(Nv,1)*meanS;
psiN = rawN - ones(Nv,1)*meanN;
psi = psiS + psiN;
BS = utils.vint(v.*psiS, w);
BN = utils.vint(v.*psiN, w);
uS = -BS/V;
uN = -BN/V;

R = struct('rho_e',rho_e,'S_e',S_e,'N_e',N_e, ...
           'Sx_e',Sx_e,'Nx_e',Nx_e,'St_e',St_e,'Nt_e',Nt_e, ...
           'eps_response',eps_response,'phi_scale',phi_scale, ...
           'psiS',psiS,'psiN',psiN,'psi',psi, ...
           'BS',BS,'BN',BN,'uS',uS,'uN',uN);
end
