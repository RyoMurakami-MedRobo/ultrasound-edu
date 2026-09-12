function validate_execution_trace()
% Verify recorded channel increments against the actual signed DAS sum.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
cfg.probe=struct('Nelements',32,'pitch',.3e-3,'fc',5e6,'bandwidth',75);
cfg.medium=struct('c',1540); cfg.acq=struct('fs_factor',4);
cfg.tx=struct('scheme','plane','angle_deg',0,'focus_mm',18);
cfg.scat=struct('x',0,'z',.018,'rc',1);
cfg.recon=struct('zmax',.028,'fnumber',1.5,'rx_apod','rect');
cfg.options=struct('force_mock',true);
S=sim_engine(cfg); gz=linspace(.010,.026,121);
fns={@das_reference,@das_custom_template};
for j=1:numel(fns)
    fh=fns{j};
    baseline=fh(S.RF(:,:,1),S.tx(1),S.rx_pos,0,gz,S.c,S.fs);
    [img,delays,tr]=fh(S.RF(:,:,1),S.tx(1),S.rx_pos,0,gz,S.c,S.fs);
    assert(isequal(img,baseline),'Tracing changed the image.');
    assert(max(abs(sum(tr.contributions,2)-tr.coherent(:)))<1e-12,'Trace sum differs.');
    assert(all(tr.contributions(~isfinite(delays))==0),'Excluded channel contributed.');
end
fprintf('MATLAB execution trace tests OK\n');
end
