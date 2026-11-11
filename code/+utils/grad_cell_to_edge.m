function gu = grad_cell_to_edge(u, dx)
%GRAD_CELL_TO_EDGE derivative at interfaces with homogeneous Neumann extension.
uex = [u(1), u(:).', u(end)];
gu = diff(uex)/dx;
end
