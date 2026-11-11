function [Unew, diagnostic] = step_kinetic_AP(U, grid, p)
%STEP_KINETIC_AP One step of the decoupled AP micro-macro scheme with exact decomposition of T=1+eps*psi.
rho = U.rho; S = U.S; N = U.N; g = U.g;
dt = p.dt; epsv = p.eps;
reflection = model.reflection_indices(grid);
g(:,[1,end]) = 0.5*(g(:,[1,end])+g(reflection,[1,end]));

R = model.compute_response(rho, S, N, grid, p);
Mres = model.upwind_g_transport(g, grid);
nonlin = model.project_zero_mean(R.psi .* g, grid);
Kold = model.K_equilibrium(rho, R, grid);
lam = dt/(epsv^2 + dt);

gtilde = (epsv^2*g - epsv*dt*Mres + dt*Kold - epsv*dt*nonlin)/(epsv^2 + dt);
gtilde(:,[1,end]) = 0.5*(gtilde(:,[1,end])+gtilde(reflection,[1,end]));

D = utils.vint(grid.v.^2, grid.w)/grid.V;
bS = R.BS/grid.V;
bN = R.BN/grid.V;
a = lam*D;
cS = lam*bS; cN = lam*bN;
% no-flux boundaries
cS([1,end]) = 0; cN([1,end]) = 0;

A = model.build_rho_matrix(grid, dt, a, cS, cN, R.uS, R.uN);
% known RHS
Jtilde_raw = utils.vint(grid.v.*gtilde,grid.w);
Jtilde = Jtilde_raw; Jtilde([1,end]) = 0;
v_dgtilde = diff(Jtilde)/grid.dx;
Fdiff_old = a * utils.grad_cell_to_edge(rho, grid.dx);
rhoS_old = model.upwind_rho_edges(rho, R.uS);
rhoN_old = model.upwind_rho_edges(rho, R.uN);
FoldS = cS .* rhoS_old; FoldN = cN .* rhoN_old;
rhs = rho(:)/dt - v_dgtilde(:) ...
    - model.div_flux_center(Fdiff_old, grid).' ...
    - model.div_flux_center(FoldS, grid).' ...
    - model.div_flux_center(FoldN, grid).';

rho_new = (A\rhs).';
Knew = model.K_equilibrium(rho_new, R, grid);
g_new = gtilde + lam*(Knew-Kold);
% Close the provisional boundary before testing the reconstructed distribution.
g_new(:,[1,end]) = 0.5*(g_new(:,[1,end])+g_new(reflection,[1,end]));
mass_before_limiter = sum(rho_new)*grid.dx;
rho_after_density = rho_new;

if p.use_limiter
    [rho_new, g_new, limiter_info] = model.positivity_limiter(rho_new, g_new, grid, p);
else
    limiter_info = struct('theta_cells',0,'theta_min',1,'negative_density_cells',0, ...
        'eta_triggered',false,'eta_edges',0,'eta_min',1,'min_f_before_eta',NaN, ...
        'early_return',false,'correction_applied',false,'mass_scale',1);
end

if nargout>1
    % These are the fluxes used in THIS density solve, before post-processing.
    Fnew = a*utils.grad_cell_to_edge(rho_after_density,grid.dx) ...
        + cS.*model.upwind_rho_edges(rho_after_density,R.uS) ...
        + cN.*model.upwind_rho_edges(rho_after_density,R.uN);
    Jused = Jtilde-Fnew+Fdiff_old+FoldS+FoldN;
    diagnostic = struct('mass_before',sum(rho)*grid.dx, ...
        'mass_after_density',mass_before_limiter,'mass_after',sum(rho_new)*grid.dx, ...
        'Jtilde_raw_boundary',Jtilde_raw([1,end]),'Jused_boundary',Jused([1,end]), ...
        'mass_residual',mass_before_limiter-sum(rho)*grid.dx-dt*(Jused(1)-Jused(end)), ...
        'limiter_mass_change',sum(rho_new)*grid.dx-mass_before_limiter,'limiter',limiter_info);
end

% Update chemicals. Use rho_new by default. This is slightly more stable for long runs.
rho_for_chem = rho;
if p.update_chem_with_new_rho
    rho_for_chem = rho_new;
end
[S_new,N_new] = src.step_chemicals(S, N, rho_for_chem, grid, p);

Unew = struct('rho',rho_new,'S',S_new,'N',N_new,'g',g_new,'t',U.t+dt);
end
