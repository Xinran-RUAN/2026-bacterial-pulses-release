function ok=run_paper_case(case_id,output_dir,smoke)
%RUN_PAPER_CASE Solve one fine-grid configuration and save MATLAB results.
% The optional smoke run uses the production mesh and stops at t=0.1.
root=startup_paper;
if nargin<3,smoke=false;end
assert(nargin>=2,'Specify an output directory outside the release data.');
cases=paper_cases;
c=cases(strcmp({cases.id},case_id));
assert(numel(c)==1,'Unknown production case: %s',case_id);

baseline=fullfile(root,'data','figure01','baseline.mat');
assert(isfile(baseline),'Missing data/figure01/baseline.mat. Install the Zenodo dataset as described in README.md.');
ref=load(baseline);p=ref.p;
p.Tfinal=50;
if smoke,p.Tfinal=0.1;end
p.output_times=0:.1:p.Tfinal;
p.verbose=false;p.report_dt=1;p.save_snapshots_to_disk=false;p.use_limiter=true;
p.phi_type='sign';if c.a>0,p.phi_type='exp_sign';p.phi_alpha=c.a;end
p.response_eps=c.epsilon;if strcmp(c.kind,'kinetic'),p.eps=c.epsilon;end
grid=model.make_grid(p);assert(isequal(grid,ref.g),'Grid differs from the saved baseline.');

case_dir=fullfile(output_dir,sprintf('figure%02d',c.figure));
published=fullfile(root,'data');
check_output_directory(case_dir,published);
if ~isfolder(case_dir),mkdir(case_dir);end
steps=round(p.Tfinal/p.dt);
settings=struct('case',c,'p',p,'version',version,'computer',computer,'Nx',grid.Nx,'Nv',grid.Nv);
finalfile=fullfile(case_dir,[c.id '.mat']);
checkpoint=fullfile(case_dir,[c.id '_checkpoint.mat']);
if isfile(finalfile)
    z=read_completed_case(finalfile);
    assert(isequal(z.settings.p,p)&&isequal(z.settings.case,c)&&isequal(z.out.grid,grid) ...
        &&abs(z.out.U.t-p.Tfinal)<1e-7, ...
        'Existing case settings differ. Choose a separate output directory.');
    ok=true;fprintf('Case %s is already complete.\n',case_id);return;
end
U=model.initial_data(grid,p);if strcmp(c.kind,'macro'),U.g=[];end;M0=sum(U.rho)*grid.dx;
stats=struct('mass_initial',M0,'mass_final',M0,'max_relative_mass_drift',0,'min_rho',min(U.rho),'min_final_f',min(U.rho)/grid.V,'min_scaled_final_f',1, ...
 'max_boundary_flux',0,'max_raw_Jtilde_boundary',0,'max_zero_mean',0,'max_scaled_zero_mean',0,'max_mass_residual',0,'max_limiter_mass_change',0,'max_used_boundary_flux',0, ...
 'early_return_steps',0,'correction_steps',0,'theta_steps',0,'theta_cells',0,'theta_min',1,'eta_steps',0,'eta_edges',0,'eta_min',1,'negative_density_cells',0,'seconds',0);
snap=struct('t',0,'rho',U.rho,'S',U.S,'N',U.N,'mass',M0);micro=struct('t',{},'g',{});if ~isempty(U.g),micro(1)=struct('t',0,'g',U.g);end
history=zeros(numel(p.output_times),13);history(1,:)=[0 M0 M0 M0 0 0 0 0 0 0 0 0 1];events=zeros(steps,11);nevents=0;n0=0;elapsed_prior=0;
if isfile(checkpoint)
 z=load(checkpoint);assert(isequal(z.settings,settings),'Checkpoint settings differ. Choose a separate output directory.');U=z.U;stats=z.stats;snap=z.snap;micro=z.micro;history=z.history;events=z.events;nevents=z.nevents;n0=z.n0;elapsed_prior=z.elapsed_prior;
 fprintf('RESUME %s t=%g\n',case_id,U.t);
end
timer=tic;
fprintf('START %s Nx=%d Nv=%d dx=%.8g dt=%.8g T=%g MATLAB=%s\n',case_id,grid.Nx,grid.Nv,p.dx,p.dt,p.Tfinal,version);
for n=n0+1:steps
  if strcmp(c.kind,'kinetic')
   [candidate,d]=src.step_kinetic_AP(U,grid,p);
   F=candidate.rho/grid.V+.5*p.eps*(candidate.g(:,1:end-1)+candidate.g(:,2:end));
   scale=abs(candidate.rho/grid.V)+.5*p.eps*(abs(candidate.g(:,1:end-1))+abs(candidate.g(:,2:end)));
   scaled=min(F./max(scale,realmin),[],'all');assert(all(isfinite(F(:)))&&scaled>=-512*eps,'Final reconstructed distribution failed positivity.');
   means=utils.vint(candidate.g,grid.w);mscale=utils.vint(abs(candidate.g),grid.w);rscale=max([candidate.rho(1),candidate.rho],[candidate.rho,candidate.rho(end)])/p.eps;
   scaledmean=max(abs(means)./max(mscale+rscale,realmin));assert(scaledmean<=512*eps,'Zero velocity mean failed.');
   J=utils.vint(grid.v.*candidate.g,grid.w);z=d.limiter;
   stats.min_final_f=min(stats.min_final_f,min(F,[],'all'));stats.min_scaled_final_f=min(stats.min_scaled_final_f,scaled);
   stats.max_zero_mean=max(stats.max_zero_mean,max(abs(means)));stats.max_scaled_zero_mean=max(stats.max_scaled_zero_mean,scaledmean);
   stats.max_boundary_flux=max(stats.max_boundary_flux,max(abs(J([1 end]))));stats.max_raw_Jtilde_boundary=max(stats.max_raw_Jtilde_boundary,max(abs(d.Jtilde_raw_boundary)));
   stats.max_used_boundary_flux=max(stats.max_used_boundary_flux,max(abs(d.Jused_boundary)));stats.max_limiter_mass_change=max(stats.max_limiter_mass_change,abs(d.limiter_mass_change));
   stats.early_return_steps=stats.early_return_steps+z.early_return;stats.correction_steps=stats.correction_steps+z.correction_applied;
   stats.theta_steps=stats.theta_steps+(z.theta_cells>0);stats.theta_cells=stats.theta_cells+z.theta_cells;stats.theta_min=min(stats.theta_min,z.theta_min);
   stats.eta_steps=stats.eta_steps+z.eta_triggered;stats.eta_edges=stats.eta_edges+z.eta_edges;stats.eta_min=min(stats.eta_min,z.eta_min);stats.negative_density_cells=stats.negative_density_cells+z.negative_density_cells;
   if z.correction_applied
    nevents=nevents+1;events(nevents,:)=[n*p.dt z.min_f_before_eta min(candidate.rho) z.theta_cells z.theta_min z.mass_scale z.eta_triggered z.eta_edges z.eta_min d.limiter_mass_change scaled];
   end
  else
   candidate=src.step_macro(U,grid,p);d=struct('mass_before',sum(U.rho)*grid.dx,'mass_after_density',sum(candidate.rho)*grid.dx,'mass_residual',sum(candidate.rho-U.rho)*grid.dx,'limiter_mass_change',0,'Jtilde_raw_boundary',[0 0]);
  end
  assert(all(isfinite(candidate.rho))&&all(candidate.rho>=0),'Final density failed positivity.');
  U=candidate;M=sum(U.rho)*grid.dx;stats.min_rho=min(stats.min_rho,min(U.rho));stats.max_relative_mass_drift=max(stats.max_relative_mass_drift,abs(M-M0)/M0);stats.max_mass_residual=max(stats.max_mass_residual,abs(d.mass_residual));
  assert(stats.max_relative_mass_drift<1e-9,'Mass drift threshold exceeded.');
  if mod(n,100)==0
   k=n/100+1;snap(k)=struct('t',n*p.dt,'rho',U.rho,'S',U.S,'N',U.N,'mass',M);
   history(k,:)=[n*p.dt d.mass_before d.mass_after_density M d.mass_residual d.limiter_mass_change max(abs(d.Jtilde_raw_boundary)) stats.max_boundary_flux stats.max_zero_mean stats.early_return_steps stats.correction_steps stats.eta_steps stats.eta_min];
  end
  if any(n==[1000 10000 20000 30000 40000 50000])&&~isempty(U.g),micro(end+1)=struct('t',n*p.dt,'g',U.g);end
  if mod(n,1000)==0
   elapsed=elapsed_prior+toc(timer);
   fprintf('PROGRESS %s t=%g seconds=%.1f drift=%.3g bypass=%d correction=%d minF=%.3g\n',case_id,n*p.dt,elapsed,stats.max_relative_mass_drift,stats.early_return_steps,stats.correction_steps,stats.min_final_f);
  end
  if mod(n,5000)==0
   saved=struct('U',U,'stats',stats,'snap',snap,'micro',micro,'history',history, ...
       'events',events,'nevents',nevents,'n0',n,'elapsed_prior',elapsed_prior+toc(timer), ...
       'settings',settings);
   save_case(checkpoint,saved);
   clear saved
  end
 end
 stats.mass_final=sum(U.rho)*grid.dx;stats.seconds=elapsed_prior+toc(timer);if strcmp(c.kind,'macro'),stats.min_final_f=NaN;stats.min_scaled_final_f=NaN;end
 out=struct('grid',grid,'U',U,'snap',snap,'cm_t',[snap.t],'cm_x',arrayfun(@(s)sum(grid.x.*s.rho)/sum(s.rho),snap),'params',p,'kind',c.kind);
 events=events(1:nevents,:);
 save_case(finalfile,struct('out',out,'stats',stats,'settings',settings, ...
     'micro',micro,'history',history,'events',events));
 delete_if_present(checkpoint);
 fprintf('COMPLETE %s seconds=%.3f mass=%.16g drift=%.3g\n',case_id,stats.seconds,stats.mass_final,stats.max_relative_mass_drift);
 ok=true;
end

function check_output_directory(folder,published)
% Resolve the nearest existing ancestor before creating the destination.
parent=folder;
while ~isfolder(parent)
    parent=fileparts(parent);
    if isempty(parent),parent=pwd;end
end
[found,info]=fileattrib(parent);
assert(found && ~strcmp(info.Name,published) && ~startsWith(info.Name,[published filesep]), ...
    'Choose an output directory outside the release data.');
end

function save_case(destination,variables)
temporary=[tempname(fileparts(destination)) '.mat'];
cleanup=onCleanup(@() delete_if_present(temporary));
save(temporary,'-struct','variables','-v7.3');
[ok,message]=movefile(temporary,destination,'f');
assert(ok,'%s',message);
end

function delete_if_present(file)
if isfile(file),delete(file);end
end

