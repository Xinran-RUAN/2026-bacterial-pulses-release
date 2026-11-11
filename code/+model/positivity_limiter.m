function [rho_new, g_edge_new, info] = positivity_limiter(rho, g_edge, grid, p)
%POSITIVITY_LIMITER Zhang--Shu style limiter for f = rho/|V| + eps*g.
V = grid.V; w = grid.w; epsv = p.eps;
Nx = grid.Nx;
assert(isfinite(epsv) && epsv>0, 'Kinetic limiter requires epsilon>0.');
reflection = model.reflection_indices(grid);
mass0 = sum(rho)*grid.dx;
assert(isfinite(mass0) && mass0>0, 'Limiter requires positive total input mass.');
info = struct('theta_cells',0,'theta_min',1,'negative_density_cells',0, ...
    'eta_triggered',false,'eta_edges',0,'eta_min',1,'min_f_before_eta',NaN, ...
    'early_return',false,'correction_applied',false,'mass_scale',1);
assert(all(isfinite(rho(:))) && all(isfinite(g_edge(:))), ...
    'Limiter input must be finite.');
gc = 0.5*(g_edge(:,1:end-1) + g_edge(:,2:end));
F = ones(grid.Nv,1)*(rho/V) + epsv*gc;
assert(all(isfinite(F(:))), 'Reconstructed limiter input must be finite.');
if all(rho(:)>=0) && all(F(:)>=0)
    rho_new = rho;
    g_edge_new = g_edge;
    info.early_return = true;
    info.min_f_before_eta = min(F,[],'all');
    return;
end
info.correction_applied = true;
Fstar = F;
for j = 1:Nx
    rhoj = utils.vint(F(:,j), w);
    feq = rhoj/V;
    mj = min(F(:,j));
    if rhoj > 0 && mj < 0
        theta = min(1, feq/(feq - mj));
        info.theta_cells = info.theta_cells + (theta<1);
        info.theta_min = min(info.theta_min,theta);
        Fstar(:,j) = feq + theta*(F(:,j)-feq);
        Fstar(:,j) = max(Fstar(:,j),0);
    elseif rhoj <= 0
        info.negative_density_cells = info.negative_density_cells + 1;
        Fstar(:,j) = 0;
    end
end
rho_star = utils.vint(Fstar, w);
mass_star = sum(rho_star)*grid.dx;
assert(isfinite(mass_star) && mass_star>0, 'Positive mass correction failed.');
info.mass_scale = mass0/mass_star;
Fstar = info.mass_scale*Fstar;
rho_new = utils.vint(Fstar, w);
assert(all(isfinite(rho_new)) && all(rho_new>=0), 'Density correction failed.');
gc_new = (Fstar - ones(grid.Nv,1)*(rho_new/V))/epsv;
gc_ex = [gc_new(reflection,1), gc_new, gc_new(reflection,end)];
g_edge_new = 0.5*(gc_ex(:,1:end-1) + gc_ex(:,2:end));
% remove tiny residual velocity mean on interfaces
mean_g = utils.vint(g_edge_new, w)/V;
g_edge_new = g_edge_new - ones(grid.Nv,1)*mean_g;
% Keep the old output exactly when its actual reconstructed f is nonnegative.
F_return = ones(grid.Nv,1)*(rho_new/V) ...
    + 0.5*epsv*(g_edge_new(:,1:end-1)+g_edge_new(:,2:end));
info.min_f_before_eta = min(F_return,[],'all');
if any(F_return(:)<0)
    q_edge = min([rho_new(1),rho_new],[rho_new,rho_new(end)])/V;
    min_g = min(g_edge_new,[],1);
    eta_edge = ones(1,Nx+1);
    negative = min_g<0;
    eta_edge(negative) = min(1,q_edge(negative)./(-epsv*min_g(negative)));
    g_edge_new = g_edge_new .* eta_edge;
    info.eta_triggered = true;
    info.eta_edges = sum(eta_edge<1);
    info.eta_min = min(eta_edge);
end
end
