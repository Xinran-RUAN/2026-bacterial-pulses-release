function Mres = upwind_g_transport(g, grid)
%UPWIND_G_TRANSPORT Compute (I-Pi_h)(v d_x g) on the g interface grid.
dx = grid.dx; v = grid.v;
reflection = model.reflection_indices(grid);
gex = [g(reflection,2), g, g(reflection,end-1)];
gL = gex(:,1:end-2);
gC = gex(:,2:end-1);
gR = gex(:,3:end);
vp = max(v,0);
vm = max(-v,0);
M = vp.*(gC-gL)/dx - vm.*(gR-gC)/dx;
Mres = model.project_zero_mean(M, grid);
end
