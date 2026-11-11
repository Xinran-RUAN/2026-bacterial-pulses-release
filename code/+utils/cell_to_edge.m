function ue = cell_to_edge(u)
%CELL_TO_EDGE Neumann extension followed by local averaging.
uex = [u(1), u(:).', u(end)];
ue = 0.5*(uex(1:end-1) + uex(2:end));
end
