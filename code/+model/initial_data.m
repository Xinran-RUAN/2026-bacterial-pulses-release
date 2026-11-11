function U = initial_data(grid, p)
x = grid.x;
rho = p.rho0_amp * exp(-p.rho0_lambda*x);
S = p.S0 * ones(1, grid.Nx);
N = p.N0 * ones(1, grid.Nx);
g = zeros(grid.Nv, grid.Nx+1);
U = struct('rho',rho,'S',S,'N',N,'g',g,'t',0);
end
