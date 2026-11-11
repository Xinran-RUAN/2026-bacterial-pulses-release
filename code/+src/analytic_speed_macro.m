function c = analytic_speed_macro(chiS, chiN, DS, alpha)
% Solve chiN - c = chiS*c/sqrt(4*DS*alpha+c^2). Return NaN if no positive root.
fun = @(c) chiN - c - chiS*c./sqrt(4*DS*alpha + c.^2);
if fun(0) <= 0
    c = NaN; return;
end
cmax = max(2*chiN+2, 2);
if fun(cmax) > 0
    c = NaN; return;
end
c = fzero(fun,[0,cmax]);
end
