function L = lap_neumann_matrix(N, dx)
%LAP_NEUMANN_MATRIX returns -d_xx with homogeneous Neumann BC on cell centers.
e = ones(N,1);
L = spdiags([-e, 2*e, -e], -1:1, N, N)/dx^2;
L(1,1) = 1/dx^2;
L(end,end) = 1/dx^2;
end
