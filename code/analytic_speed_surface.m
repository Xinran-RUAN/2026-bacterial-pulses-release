function speed = analytic_speed_surface(chiS, chiN)
%ANALYTIC_SPEED_SURFACE Solve Eq. (5.2) on the Figure 2(a) grid.
% Rows index chiS and columns index chiN. Use the original scalar root guess.
DS = 2;
alpha = 0.05;
speed = zeros(numel(chiS), numel(chiN));
for j = 1:numel(chiN)
    for k = 1:numel(chiS)
        equation = @(x) chiN(j) - x - chiS(k) * x / sqrt(4 * DS * alpha + x^2);
        speed(k,j) = fzero(equation, 1);
    end
end
end
