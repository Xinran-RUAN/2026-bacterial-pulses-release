function grid = make_grid(p)
Nx = round(p.L/p.dx);
x = ((1:Nx)-0.5)*p.dx;
xe = (0:Nx)*p.dx;
vel = utils.velocity_grid(p.vmin, p.vmax, p.dv, p.velocity_mode);
grid = struct('Nx',Nx,'x',x,'xe',xe,'dx',p.dx, ...
              'v',vel.v,'w',vel.w,'V',vel.V,'Nv',vel.Nv);
end
