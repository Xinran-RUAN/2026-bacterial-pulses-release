function K = K_equilibrium(rho, R, grid)
%K_EQUILIBRIUM local equilibrium K on interfaces for current rho and response R.
drho = utils.grad_cell_to_edge(rho, grid.dx);
rhoS = model.upwind_rho_edges(rho, R.uS);
rhoN = model.upwind_rho_edges(rho, R.uN);
K = -(grid.v/grid.V).*drho ...
    - (ones(grid.Nv,1)*(rhoS/grid.V)).*R.psiS ...
    - (ones(grid.Nv,1)*(rhoN/grid.V)).*R.psiN;
end
