function root = startup_paper
%STARTUP_PAPER Add the release's production solver and return the release root.
root = fileparts(fileparts(mfilename('fullpath')));
code = fullfile(root, 'code');
addpath(code, '-begin');
assert(strcmp(which('model.positivity_limiter'), fullfile(code, '+model', 'positivity_limiter.m')) ...
    && strcmp(which('src.step_kinetic_AP'), fullfile(code, '+src', 'step_kinetic_AP.m')), ...
    'Another MATLAB package shadows the release solver.');
end
