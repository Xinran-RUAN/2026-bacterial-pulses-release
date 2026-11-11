function scale = response_scale(grid, p)
%RESPONSE_SCALE Multiplicative scale applied to phi in the response.
%
% In the sign-response macro limit,
%
%   u = - chi/|V| * int_V v * phi(v*q_x) dv.
%
% With phi(y)=-sign(y), this gives
%
%   u = chi * scale * int_V |v|dv / |V| * sign(q_x).
%
% Therefore scale=|V|/int_V |v|dv makes the macroscopic drift amplitude
% equal to the chi used in the analytical formulas.

scale = 1;
if isfield(p, 'phi_scale') && ~isempty(p.phi_scale)
    scale = p.phi_scale;
end

if isfield(p, 'normalize_phi_amplitude') && p.normalize_phi_amplitude
    abs_moment = utils.vint(abs(grid.v), grid.w);
    if abs_moment > 0
        scale = scale * grid.V / abs_moment;
    end
end
end
