function vel = velocity_grid(vmin, vmax, dv, mode)
%VELOCITY_GRID Build velocity nodes and quadrature weights.
% mode='continuous' uses a uniform trapezoid rule on [vmin,vmax].
% mode='four' uses {-1,-0.5,0.5,1} with equal weights summing to |V|=2.
if nargin < 4, mode = 'continuous'; end
switch lower(mode)
    case 'four'
        v = [-1; -0.5; 0.5; 1];
        w = 0.5*ones(4,1);      % sum(w)=2
    otherwise
        v = (vmin:dv:vmax).';
        w = dv*ones(size(v));
        w(1) = 0.5*dv; w(end) = 0.5*dv;
end
vel.v = v;
vel.w = w;
vel.V = sum(w);
vel.Nv = numel(v);
end
