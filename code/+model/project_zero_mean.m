function Z = project_zero_mean(F, grid)
%PROJECT_ZERO_MEAN Apply I-Pi to an Nv-by-Nx array, Pi h = <h>/|V|.
meanF = utils.vint(F, grid.w)/grid.V;
Z = F - ones(grid.Nv,1)*meanF;
end
