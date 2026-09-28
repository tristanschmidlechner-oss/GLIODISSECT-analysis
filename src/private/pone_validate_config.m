function cfg=pone_validate_config(cfg,N)
assert(isstruct(cfg)&&isscalar(cfg),'PONE:Config','Configuration must be a scalar struct.');
names={'controls','full_required','training_required','q2_required'};
assert(all(isfield(cfg,names))&&numel(fieldnames(cfg))==4,'PONE:Config','Use the four explicit configuration fields.');
for k=1:4,validateattributes(cfg.(names{k}),{'numeric'},{'scalar','integer','finite','positive'});end
assert(cfg.controls==N,'PONE:ControlCount','Input and configured control counts differ.');
assert(N>=4 && cfg.full_required>=3 && cfg.full_required<=N && ...
 cfg.training_required>=3 && cfg.training_required<=N-1 && cfg.q2_required>=2 && cfg.q2_required<=N, ...
 'PONE:Config','Invalid coverage counts.');
end
