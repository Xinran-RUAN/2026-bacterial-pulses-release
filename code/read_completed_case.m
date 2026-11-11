function d = read_completed_case(file)
%READ_COMPLETED_CASE Read a final MAT or a complete historical checkpoint.
vars = who('-file', file);
if any(strcmp(vars, 'out'))
    d = load(file, 'out', 'stats', 'settings', 'events', 'micro', 'history');
else
    d = load(file, 'U', 'stats', 'snap', 'settings', 'events', 'nevents', ...
        'micro', 'history', 'n0', 'elapsed_prior');
    assert(d.n0 == 50000 && abs(d.U.t - 50) < 1e-7 && numel(d.snap) == 501, ...
        'Only a complete T=50 checkpoint is admissible.');
    grid = model.make_grid(d.settings.p);
    d.stats.mass_final = sum(d.U.rho) * grid.dx;
    d.stats.seconds = d.elapsed_prior;
    d.events = d.events(1:d.nevents,:);
    d.out = struct('grid', grid, 'U', d.U, 'snap', d.snap, 'cm_t', [d.snap.t], ...
        'cm_x', arrayfun(@(s) sum(grid.x.*s.rho)/sum(s.rho), d.snap), ...
        'params', d.settings.p, 'kind', d.settings.case.kind);
end
end
