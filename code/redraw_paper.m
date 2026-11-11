function figures = redraw_paper(selection, output_dir, recomputed_dir)
%REDRAW_PAPER Draw a figure (1-7), 'table', or 'all' from saved data.
% Output defaults to a new temporary directory. Supply a recomputed data
% root as the third argument to use rerun files for Figures 1 and 3-7
% and Table 1. Figure 2 uses Eq. (5.2) and the supplied centroid data.

root = startup_paper;
if nargin < 1 || isempty(selection), selection = 'all'; end
if nargin < 2 || isempty(output_dir), output_dir = tempname; end
if nargin < 3, recomputed_dir = []; end
data_dir = fullfile(root, 'data');
case_dir = data_dir;
if ~isempty(recomputed_dir), case_dir = recomputed_dir; end
assert(strcmp(selection, 'all') || strcmp(selection, 'table') || ...
    (isnumeric(selection) && isscalar(selection) && ismember(selection, 1:7)), ...
    'Select a figure number from 1 to 7, ''table'', or ''all''.');
if ~isfolder(output_dir), mkdir(output_dir); end

figures = gobjects(0);
indices = 1:7;
if isnumeric(selection), indices = selection; end
if strcmp(selection, 'table'), indices = []; end
cases = paper_cases();
blue = [0 0.447 0.741];
orange = [0.85 0.325 0.098];
gold = [0.929 0.694 0.125];
colors = {blue, orange, gold, [0 0 0]};
styles = {'-', '--', '-.', '-'};
epsilons = [0.1 0.05 0.01 0];
epsilon_labels = {'$\varepsilon=0.1$', '$\varepsilon=0.05$', '$\varepsilon=0.01$', 'macro'};
times = [1 10 20 30 40 50];

for index = indices
    switch index
        case 1
            s = load_case(case_dir, cases, 0, 0);
            [chiS, chiN] = src.effective_chi_macro(s.params, s.grid);
            speed = src.analytic_speed_macro(chiS, chiN, s.params.DS, s.params.alpha);
            D = utils.vint(s.grid.v.^2, s.grid.w)/s.grid.V;
            comparison_times = [30 50];
            limits = [14 33; 28 47];
            for panel = 1:2
                rho = s.rho(times == comparison_times(panel), :);
                [~, peak] = max(rho);
                z = s.grid.x - s.grid.x(peak);
                unit = exp((chiS+chiN-speed)/D*min(z,0) + (chiN-chiS-speed)/D*max(z,0));
                analytical = unit*sum(rho)/sum(unit);
                [fig, ax] = new_axes([760 560], [.13 .14 .775 .785], 28);
                draw_line(ax,s.grid.x,rho,'Numerical',[0 0 1],'-',2.2);
                draw_line(ax,s.grid.x,analytical,'Analytical',orange,'--',2.2);
                set(ax, 'XLim', limits(panel,:), 'YLim', [0 .115], 'YTick', 0:.02:.1);
                labels(ax, '$x$', '$\rho$');
                legend(ax, 'Location','northeast', 'Interpreter','latex', 'FontSize',28);
                figures(end+1) = output_figure(fig, output_dir, sprintf('figure01%c', 'a'+panel-1), nargout);
            end
        case 2
            centroids = load(fullfile(data_dir, 'figure02', 'historical_centroids.mat'), ...
                'chiS', 'chiN', 'macro_x', 'kinetic_x');
            chiS = centroids.chiS;
            chiN = centroids.chiN;
            macro = centroid_speeds(centroids.macro_x, 10);
            kinetic = centroid_speeds(centroids.kinetic_x, 6);
            analytic_grid = 0:.02:2;
            speeds = {analytic_speed_surface(analytic_grid, analytic_grid), macro, kinetic};
            xs = {analytic_grid, chiN, chiN};
            ys = {analytic_grid, chiS, chiS};
            titles = {'Predicted wave speed', 'Wave speed for the macro model', ...
                'Wave speed for the kinetic model'};
            for panel = 1:3
                [fig, ax] = new_axes([560 420], [.1207 .2019 .6995 .683], 30);
                pcolor(ax, xs{panel}, ys{panel}, speeds{panel});
                shading(ax, 'interp');
                colormap(ax, jet(256));
                clim(ax, [0 2]);
                colorbar(ax, 'Position', [.85 .2019 .036 .683]);
                set(ax,'XLim',[0 2],'YLim',[0 2],'XTick',0:2,'YTick',0:2);
                labels(ax, '$\chi_N$', '$\chi_S$');
                title(ax, titles{panel}, 'Interpreter','latex', 'FontSize',22);
                figures(end+1) = output_figure(fig, output_dir, sprintf('figure02%c', 'a'+panel-1), nargout);
            end
        case 3
            states = load_cases(case_dir, cases, zeros(1,4), epsilons);
            [fig, ax] = new_axes([1200 620], [.13 .2674 .775 .6707], 26);
            legend_lines = gobjects(1,4);
            for k = 1:numel(times)
                for j = 1:numel(states)
                    h = draw_line(ax,states(j).grid.x,states(j).rho(k,:),epsilon_labels{j}, ...
                        colors{j},styles{j},2.2 + .8*(j == 4));
                    if k == 1, legend_lines(j) = h; end
                end
            end
            evolution_axes(ax);
            label_times(ax,states(4),times);
            lg = legend(ax, legend_lines, 'Location','southoutside', 'NumColumns',4, ...
                'Interpreter','latex','FontSize',20);
            lg.Position(1:2) = [.3115 .0908];
            figures(end+1) = output_figure(fig, output_dir, 'figure03', nargout);
        case 4
            es = [.001 .002 .005 .01 .02 .05 .1];
            macro = load_case(case_dir, cases, 0, 0);
            errors = zeros(size(es));
            for j = 1:numel(es)
                s = load_case(case_dir, cases, 0, es(j));
                errors(j) = sqrt(s.grid.dx*sum((s.rho(end,:)-macro.rho(end,:)).^2));
            end
            [fig, ax] = new_axes([560 430], [.23 .2 .65 .67], 22);
            h = draw_line(ax,es,errors,'numerical error',blue,'-',1.9);
            set(h,'Marker','o','MarkerFaceColor','white','MarkerSize',7);
            draw_line(ax,es,errors(end)*(es/es(end)),'$O(\varepsilon)$','k','--',1.3);
            draw_line(ax,es,errors(end)*(es/es(end)).^2,'$O(\varepsilon^2)$','k','-.',1.3);
            set(ax, 'XScale','log','YScale','log','XLim',[0.0008317637711026709 0.12022644346174127], ...
                'YLim',[1e-5 1], 'XTick',[.001 .01 .1], 'YTick',[1e-5 1e-4 .001 .01 .1 1]);
            grid(ax,'on');
            grid(ax,'minor');
            labels(ax, '$\varepsilon$', '$\|\rho^\varepsilon(\cdot,t)-\rho(\cdot,t)\|_{L^2}$');
            title(ax,'Density error at $t=50$','Interpreter','latex');
            legend(ax,'Location','northwest','Interpreter','latex','FontSize',16);
            figures(end+1) = output_figure(fig, output_dir, 'figure04', nargout);
        case 5
            states = load_cases(case_dir, cases, zeros(1,4), epsilons);
            [fig, ax] = new_axes([720 460], [.13 .1357 .775 .7869], 22);
            h = gobjects(1,5);
            h(4) = draw_line(ax,states(4).cm_t,speed_curve(states(4)),epsilon_labels{4},'k','-',2.4);
            speed = src.analytic_speed_macro(.5,1.1,2,.05);
            h(5) = draw_line(ax,[4 50],[speed speed],'$\sigma_{\mathrm{pred}}=0.724$','k','--',1.6);
            for j = 1:3
                h(j) = draw_line(ax,states(j).cm_t,speed_curve(states(j)),epsilon_labels{j},colors{j},styles{j},1.8);
            end
            set(ax,'XLim',[0 50],'YLim',[.55 1],'XTick',0:5:50,'YTick',.55:.05:1);
            grid(ax,'on');
            labels(ax,'$t$','wave speed');
            title(ax,'Traveling speed comparison','Interpreter','latex');
            legend(ax,h([4 5 1 2 3]),'Location','northeast','Interpreter','latex','FontSize',16);
            figures(end+1) = output_figure(fig, output_dir, 'figure05', nargout);
        case 6
            s = load_case(case_dir, cases, .001, .01);
            [fig, ax] = new_axes([1200 620], [.13 .2674 .775 .6707], 26);
            for k = 1:numel(times)
                draw_line(ax,s.grid.x,s.rho(k,:),sprintf('$t=%d$',times(k)),'k','-',2.2);
            end
            evolution_axes(ax);
            label_times(ax,s,times);
            figures(end+1) = output_figure(fig, output_dir, 'figure06', nargout);
        case 7
            a_values = [.001 .0005 .0001 0];
            states = {load_cases(case_dir, cases, .001*ones(1,4), epsilons), ...
                load_cases(case_dir, cases, a_values, .01*ones(1,4))};
            legend_labels = {epsilon_labels, {'$a=0.001$', '$a=0.0005$', '$a=0.0001$', '$a=0$'}};
            positions = [.1 .23 .36 .7; .58 .23 .36 .7];
            [fig, ax] = new_axes([1200 540],positions(1,:),22);
            for panel = 1:2
                if panel == 2, ax = new_axes_on(fig,positions(2,:),22); end
                for j = 1:4
                    s = states{panel}(j);
                    draw_line(ax,s.grid.x,s.rho(end,:),legend_labels{panel}{j},colors{j},styles{j},2.2+.6*(j == 4));
                end
                set(ax,'XLim',[20 60],'YLim',[0 .12],'XTick',20:10:60,'YTick',0:.02:.12);
                labels(ax,'$x$','$\rho(x,50)$');
                grid(ax,'on');
                legend(ax,'Location','northeast','Interpreter','latex','FontSize',17);
                annotation(fig,'textbox',[positions(panel,1)+.135 .045 .09 .075], ...
                    'String',sprintf('(%c)','a'+panel-1),'LineStyle','none', ...
                    'HorizontalAlignment','center','FontSize',28.1);
            end
            figures(end+1) = output_figure(fig, output_dir, 'figure07', nargout);
    end
end

if strcmp(selection,'table') || strcmp(selection,'all')
    values = [.001 .0005 .0001 0];
    states = load_cases(case_dir, cases, values, .01*ones(1,4));
    file = fopen(fullfile(output_dir,'table_1_wave_speeds.csv'),'w');
    assert(file >= 0, 'Cannot write the Table 1 CSV.');
    close_file = onCleanup(@() fclose(file));
    fprintf(file,'a,numerical_propagation_speed,L2_density_to_sign\n');
    for j = 1:numel(states)
        s = states(j);
        fit = s.cm_t >= 30-1e-8 & s.cm_t <= 50+1e-8;
        coeff = polyfit(s.cm_t(fit),s.cm_x(fit),1);
        density_error = sqrt(s.grid.dx*sum((s.rho(end,:)-states(end).rho(end,:)).^2));
        fprintf(file,'%.17g,%.17g,%.17g\n',values(j),coeff(1),density_error);
    end
end
fprintf('Outputs: %s\n', output_dir);
end

function s = load_case(data_dir, cases, a, epsilon)
c = cases([cases.a] == a & [cases.epsilon] == epsilon);
assert(isscalar(c),'Unknown or ambiguous paper case.');
file = fullfile(data_dir, sprintf('figure%02d',c.figure), [c.id '.mat']);
assert(isfile(file),'Required case is missing: %s',file);
d = read_completed_case(file);
o = d.out;
assert(o.params.dx == .025 && o.params.dt == .001 && o.params.dv == .02 && ...
    abs(o.U.t-50) < 1e-7 && numel(o.snap) == 501 && o.params.response_eps == epsilon, ...
    'Grid, time or epsilon does not match the paper case: %s',file);
if a == 0
    assert(strcmp(o.params.phi_type,'sign'),'Wrong response function: %s',file);
else
    assert(strcmp(o.params.phi_type,'exp_sign') && o.params.phi_alpha == a, ...
        'Wrong exponential response parameter: %s',file);
end
snapshots = o.snap(10*[1 10 20 30 40 50]+1);
s = struct('grid',o.grid,'params',o.params,'rho',vertcat(snapshots.rho), ...
    'cm_t',o.cm_t,'cm_x',o.cm_x);
end

function states = load_cases(data_dir, cases, a_values, epsilons)
states = arrayfun(@(a,e) load_case(data_dir,cases,a,e),a_values,epsilons);
end

function result = centroid_speeds(centroids, tail)
result = zeros(size(centroids,1),size(centroids,2));
for j = 1:size(centroids,2)
    for k = 1:size(centroids,1)
        increments = diff(reshape(centroids(k,j,end-tail:end),1,[]));
        result(k,j) = mean(increments);
    end
end
end

function speed = speed_curve(s)
[t, ix] = unique(s.cm_t(:),'stable');
x = s.cm_x(:);
times = (0:.1:50).';
x = interp1(t,x(ix),times,'pchip',NaN);
speed = nan(size(times));
for k = 1:numel(times)
    keep = abs(times-times(k)) <= .5 & isfinite(x);
    if nnz(keep) >= 2
        fit = polyfit(times(keep),x(keep),1);
        speed(k) = fit(1);
    end
end
end

function [fig, ax] = new_axes(canvas, position, font_size)
fig = figure('Visible','off','Color','white','Position',[100 100 canvas]);
fig.Position = [100 100 canvas];
ax = new_axes_on(fig,position,font_size);
end

function ax = new_axes_on(fig, position, font_size)
ax = axes(fig,'Position',position,'FontSize',font_size,'Box','on','LineWidth',.8, ...
    'TickLabelInterpreter','latex');
hold(ax,'on');
end

function line_handle = draw_line(ax, x, y, name, color, style, width)
line_handle = plot(ax,x,y,'Color',color,'LineStyle',style,'LineWidth',width,'DisplayName',name);
end

function labels(ax, x, y)
xlabel(ax,x,'Interpreter','latex');
ylabel(ax,y,'Interpreter','latex');
end

function evolution_axes(ax)
set(ax,'XLim',[0 40],'YLim',[0 .2],'XTick',0:10:40,'YTick',0:.05:.2);
ytickformat(ax,'%.2f');
grid(ax,'on');
labels(ax,'$x$','$\rho(x,t)$');
end

function label_times(ax, s, times)
[peak_heights, peak_cells] = max(s.rho, [], 2);
for k = 1:numel(times)
    offset = -.4;
    align = 'right';
    if k == 1
        offset = 1.4;
        align = 'left';
    end
    text(ax,s.grid.x(peak_cells(k))+offset,peak_heights(k),sprintf('$t=%d$',times(k)), ...
        'HorizontalAlignment',align,'Interpreter','latex','FontSize',17);
end
end

function fig = output_figure(fig, output_dir, name, keep_open)
base = fullfile(output_dir,name);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
exportgraphics(fig,[base '.png'],'Resolution',150);
if ~keep_open
    close(fig);
    fig = gobjects(1);
end
end
