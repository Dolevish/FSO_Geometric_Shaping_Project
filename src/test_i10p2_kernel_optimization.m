function report = test_i10p2_kernel_optimization()
%TEST_I10P2_KERNEL_OPTIMIZATION  Regression checks for Commit I.10.2.

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.2 - likelihood-kernel runtime regression suite\n');
    fprintf('============================================================\n');

    tests={@test_runtime_tuning_map,@test_i10p2_matches_i10p1, ...
        @test_zero_fast_path_agreement,@test_symbol_batch_agreement};
    names={'Deterministic runtime tuning map is exact', ...
        'I.10.2 auto path matches I.10.1 numerical path', ...
        'x=0 fast path preserves AMI', ...
        'Symbol batch sizes preserve AMI'};

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
        error('test_i10p2_kernel_optimization:Failure','Commit I.10.2 regression suite failed.');
    end
end


function test_runtime_tuning_map()
    [b4,s4,m4]=select_logh_runtime_tuning(4);
    [b8,s8]=select_logh_runtime_tuning(8);
    [b16,s16]=select_logh_runtime_tuning(16);
    [b32,s32,m32]=select_logh_runtime_tuning(32);
    assert(b4==1024 && s4==2);
    assert(b8==1024 && s8==2);
    assert(b16==2048 && s16==1);
    assert(b32==256 && s32==2);
    assert(strcmp(m4.ruleVersion,'I10.2-v1') && strcmp(m32.ruleVersion,'I10.2-v1'));
end


function test_i10p2_matches_i10p1()
    cases=[4 20 0.1;16 30 0.2;32 30 0.3];
    for k=1:size(cases,1)
        M=cases(k,1); snr=cases(k,2); sig=cases(k,3); [dt,~]=select_logh_dt(M,snr,sig);
        cfg=build_fso_config(M,1,snr,sig,'minGap',0.01,'pinZero',true, ...
            'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive', ...
            'YGridScale',1.1,'loghDt',dt,'loghSpanSigma',8,'loghYBlockSize',512);
        X={define_constellation(M,1,true,"mean"),wide_geometry(M,1,0.01)};
        for j=1:numel(X)
            x=X{j}(:);
            ref=cfg; ref.loghRuntimeKernel='i10p1'; ref.loghRuntimeBlockPolicy='fixed';
            ref.loghSymbolBatchSize=1; ref.loghZeroFastPath=false;
            vRef=AMI_noCSI_fast_logh_adaptive_grid(x,ref.px,ref,dt,8,1.1,512);
            vNew=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,1.1,512);
            err=abs(vNew-vRef);
            assert(err<=5e-12*max(1,abs(vRef)), ...
                'I.10.2 mismatch M=%d SNR=%.1f sigma=%.3f geometry=%d: %.3e',M,snr,sig,j,err);
        end
    end
end


function test_zero_fast_path_agreement()
    M=16; snr=30; sig=0.2; dt=0.01;
    cfg=build_fso_config(M,1,snr,sig,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive');
    cfg.loghRuntimeKernel='i10p2'; cfg.loghRuntimeBlockPolicy='fixed'; cfg.loghSymbolBatchSize=1;
    x=wide_geometry(M,1,0.01);
    a=cfg; a.loghZeroFastPath=false;
    b=cfg; b.loghZeroFastPath=true;
    va=AMI_noCSI_fast_logh_adaptive_grid(x,a.px,a,dt,8,1.1,512);
    [vb,d]=AMI_noCSI_fast_logh_adaptive_grid(x,b.px,b,dt,8,1.1,512);
    assert(abs(va-vb)<=5e-12*max(1,abs(va)),'Zero fast path changed AMI.');
    assert(d.nZeroFastSymbols>=1 && d.zeroFastPath,'Expected pinned-zero fast path to be exercised.');
end


function test_symbol_batch_agreement()
    M=32; snr=30; sig=0.3; dt=0.005;
    cfg=build_fso_config(M,1,snr,sig,'minGap',0.01,'pinZero',true, ...
        'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive');
    cfg.loghRuntimeKernel='i10p2'; cfg.loghRuntimeBlockPolicy='fixed'; cfg.loghZeroFastPath=true;
    x=wide_geometry(M,1,0.01);
    vals=nan(3,1); batches=[1 2 4];
    for k=1:numel(batches)
        c=cfg; c.loghSymbolBatchSize=batches(k);
        vals(k)=AMI_noCSI_fast_logh_adaptive_grid(x,c.px,c,dt,8,1.1,256);
    end
    assert(max(abs(vals-vals(1)))<=5e-12*max(1,abs(vals(1))), ...
        'Symbol batching changed AMI beyond tolerance.');
end


function x=wide_geometry(M,Pavg,dmin)
    K=max(1,ceil(M/8)); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus<0,error('test_i10p2_kernel_optimization:Infeasible','Infeasible wide geometry.');end
    q=zeros(M,1); q(end-K+1:end)=surplus*M/K; x=b+q;
end
