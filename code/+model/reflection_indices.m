function idx = reflection_indices(grid)
%REFLECTION_INDICES Validate quadrature and pair v with -v, including old grids.
persistent oldv oldw oldV oldidx
v = grid.v(:); w = grid.w(:);
if isequal(v,oldv) && isequal(w,oldw) && isequal(grid.V,oldV)
    idx = oldidx; return;
end
assert(numel(v)==numel(w) && all(isfinite(v)) && all(isfinite(w)) ...
    && all(w>0), 'Reflection requires finite velocities and positive weights.');
tolv = 64*eps(max(1,max(abs(v))));
tolw = 64*eps(max(1,max(abs(w))));
assert(abs(sum(w)-grid.V)<=64*eps(max(1,abs(grid.V))) && grid.V>0, ...
    'grid.V must equal the positive sum of quadrature weights.');
[~,order] = sort(v); idx = zeros(size(v)); idx(order) = flipud(order);
assert(max(abs(v+v(idx)))<=tolv && max(abs(w-w(idx)))<=tolw, ...
    'Specular reflection requires paired velocities and paired weights.');
assert(isequal(idx(idx),(1:numel(v)).'), 'Reflection map must be an involution.');
oldv=v; oldw=w; oldV=grid.V; oldidx=idx;
end
