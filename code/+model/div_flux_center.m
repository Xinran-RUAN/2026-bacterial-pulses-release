function dF = div_flux_center(F, grid)
%DIV_FLUX_CENTER divergence of interface flux into cell centers.
dF = (F(2:end)-F(1:end-1))/grid.dx;
end
