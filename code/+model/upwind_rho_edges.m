function rho_up = upwind_rho_edges(rho, u)
%UPWIND_RHO_EDGES Upwind rho at interfaces for velocity u on interfaces.
rho_ex = [rho(1), rho(:).', rho(end)];
left = rho_ex(1:end-1);
right = rho_ex(2:end);
rho_up = left;
rho_up(u < 0) = right(u < 0);
% enforce no advective boundary flux by a harmless convention
rho_up(1) = rho(1);
rho_up(end) = rho(end);
end
