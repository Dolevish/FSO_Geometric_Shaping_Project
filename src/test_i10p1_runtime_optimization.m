function report = test_i10p1_runtime_optimization()
%TEST_I10P1_RUNTIME_OPTIMIZATION  Numerical regression for performance refactor.
%
% Commit I.10.1 must change runtime only. These tests compare the optimized
% cached implementation against a private copy of the pre-I.10.1 log-h
% evaluator/grid path and verify that cache reuse is actually exercised.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.1 - log-h runtime optimization regression suite\n');
    fprintf('============================================================\n');

    tests={@test_cached_path_matches_legacy,@test_cache_hit_observable,@test_no_turbulence_identity};
    names={'Cached adaptive evaluator matches legacy numerical path', ...
        'Repeated channel case reuses cached log-h quadrature', ...
        'No-turbulence path remains numerically unchanged'};

    passed=false(numel(tests),1); elapsed=nan(numel(tests),1); messages=strings(numel(tests),1);
    for k=1:numel(tests)
        fprintf('[%02d/%02d] %s ... ',k,numel(tests),names{k}); t=tic;
        try
            tests{k}(); elapsed(k)=toc(t); passed(k)=true; messages(k)="PASS";
            fprintf('PASS (%.2fs)\n',elapsed(k));
        catch ME
            elapsed(k)=toc(t); messages(k)=string(ME.message);
            fprintf('FAIL (%.2fs)\n  %s\n',elapsed(k),ME.message);
        end
    end

    report=struct('names',{names},'passed',passed,'elapsedSeconds',elapsed, ...
        'messages',messages,'nPassed',sum(passed),'nTotal',numel(tests));
    fprintf('\nPassed %d/%d tests in %.2fs.\n',report.nPassed,report.nTotal,sum(elapsed,'omitnan'));
    if ~all(passed)
        error('test_i10p1_runtime_optimization:Failure','Commit I.10.1 regression suite failed.');
    end
end


function test_cached_path_matches_legacy()
    cases=[4 20 0.1;16 30 0.2;32 30 0.3];
    for k=1:size(cases,1)
        M=cases(k,1); snr=cases(k,2); sig=cases(k,3);
        [dt,~]=select_logh_dt(M,snr,sig);
        cfg=build_fso_config(M,1,snr,sig,'minGap',0.01,'pinZero',true, ...
            'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive', ...
            'YGridScale',1.1,'loghDt',dt,'loghSpanSigma',8,'loghYBlockSize',512);

        X={define_constellation(M,1,true,"mean"), wide_geometry(M,1,0.01)};
        for j=1:numel(X)
            x=X{j}(:);
            AMI_functions.assert_constellation_feasible(x,cfg,1e-10);
            legacy=legacy_adaptive_eval(x,cfg.px,cfg,dt,8,1.1,512);
            optimized=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,1.1,512);
            err=abs(optimized-legacy);
            assert(err<=5e-12*max(1,abs(legacy)), ...
                'Numerical mismatch M=%d SNR=%.1f sigma=%.3f: %.3e',M,snr,sig,err);
        end
    end
end


function test_cache_hit_observable()
    M=8; snr=20; sig=0.1; dt=0.01;
    cfg=build_fso_config(M,1,snr,sig,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive');
    x=define_constellation(M,1,true,"mean");
    y=AMI_functions.build_noCSI_y_grid(cfg,1.1*max(abs(x)));

    clear AMI_noCSI_fast_logh_grid
    [a,d1]=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,8,y,512);
    [b,d2]=AMI_noCSI_fast_logh_grid(x,cfg.px,cfg,dt,8,y,512);
    assert(abs(a-b)<1e-14,'Repeated evaluation must be deterministic.');
    assert(~d1.quadratureCacheHit,'First evaluation after clear should miss the cache.');
    assert(d2.quadratureCacheHit,'Second identical evaluation should hit the cache.');
end


function test_no_turbulence_identity()
    M=8; cfg=build_fso_config(M,1,25,0,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive');
    x=wide_geometry(M,1,0.01);
    legacy=legacy_adaptive_eval(x,cfg.px,cfg,0.01,8,1.1,512);
    optimized=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,0.01,8,1.1,512);
    assert(abs(optimized-legacy)<1e-13,'No-turbulence result changed.');
end


function mi=legacy_adaptive_eval(x,px,params,dt,spanSigma,scale,block)
    y=AMI_functions.build_noCSI_y_grid(params,scale*max(abs(x)));
    mi=legacy_logh_grid(x,px,params,dt,spanSigma,y,block);
end


function mi_bits=legacy_logh_grid(x,px,params,dt,spanSigma,y_grid,yBlockSize)
    x=real(x(:).'); px=real(px(:).'); px=px/sum(px); y_grid=real(y_grid(:).');
    M=numel(x); Ny=numel(y_grid);

    if params.sig_t<1e-12
        means=params.R*x(:);
        d2=bsxfun(@minus,y_grid,means).^2;
        logPyx=-0.5*log(2*pi*params.sigma_n_sq)-d2/(2*params.sigma_n_sq);
        mi_bits=legacy_mi(logPyx,px,y_grid); return;
    end

    tMin=params.mu_t-spanSigma*params.sig_t; tMax=params.mu_t+spanSigma*params.sig_t;
    tGrid=tMin:dt:tMax;
    if isempty(tGrid)||tGrid(end)<tMax-100*eps(max(1,abs(tMax)))
        tGrid=[tGrid tMax]; %#ok<AGROW>
    else
        tGrid(end)=tMax;
    end
    tGrid=unique(tGrid,'stable');
    qWeights=legacy_trap_weights(tGrid);
    fT=exp(-0.5*((tGrid-params.mu_t)/params.sig_t).^2)/(sqrt(2*pi)*params.sig_t);
    q=qWeights.*fT; keep=q>0&isfinite(q); tGrid=tGrid(keep); q=q(keep);
    logQ=log(q(:)); hVals=exp(tGrid(:));
    inv2s2=1/(2*params.sigma_n_sq); logNorm=-0.5*log(2*pi*params.sigma_n_sq);
    logPyx=zeros(M,Ny);
    for j=1:M
        means=params.R*hVals*x(j);
        for b0=1:yBlockSize:Ny
            b1=min(Ny,b0+yBlockSize-1); yb=y_grid(b0:b1);
            d2=bsxfun(@minus,yb,means).^2;
            logTerms=bsxfun(@plus,logQ+logNorm,-inv2s2*d2);
            logPyx(j,b0:b1)=legacy_lse(logTerms,1);
        end
    end
    mi_bits=legacy_mi(logPyx,px,y_grid);
end


function mi=legacy_mi(logPyx,px,y)
    logPxPyx=bsxfun(@plus,log(px(:)),logPyx);
    logPy=legacy_lse(logPxPyx,1); M=size(logPyx,1); MI=0;
    for j=1:M
        Pyxj=exp(logPyx(j,:)); lr=(logPyx(j,:)-logPy)/log(2);
        z=Pyxj.*lr; z(~isfinite(z))=0; MI=MI+px(j)*trapz(y,z);
    end
    mi=max(0,MI);
end

function w=legacy_trap_weights(x)
    x=x(:).'; n=numel(x); if n==1,w=1;return;end
    dx=diff(x); w=zeros(1,n); w(1)=dx(1)/2; w(end)=dx(end)/2;
    if n>2,w(2:end-1)=(dx(1:end-1)+dx(2:end))/2;end
end
function s=legacy_lse(A,dim)
    m=max(A,[],dim); s=m+log(sum(exp(bsxfun(@minus,A,m)),dim));
end


function x=wide_geometry(M,Pavg,dmin)
    K=max(1,ceil(M/8)); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus<0,error('test_i10p1_runtime_optimization:Infeasible','Infeasible wide geometry.');end
    q=zeros(M,1); q(end-K+1:end)=surplus*M/K; x=b+q;
end
