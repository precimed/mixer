figure(1); clf; hold on;
for hvec_model = 1:2
    if libisloaded('bgmg'), unloadlibrary('bgmg'); end;

    load('test.frq.mat')
    load('test.ld.mat'); index_A = index_A + 1; index_B = index_B + 1;
    load('test_ldsum.l2.ldscore.mat'); ldscore_sum_r2 = annomat;
    load('test_ldsum.l4.ldscore.mat'); ldscore_sum_r4 = annomat;
    load('test_ldsum_per_allele.l2.ldscore.mat'); ldscore_sum_r2_per_allele = annomat;
    load('test_ldsum_per_allele.l4.ldscore.mat'); ldscore_sum_r4_per_allele = annomat;
    defvec=false(size(mafvec)); defvec(1:3:end)=true;
    %defvec=true(size(mafvec));

    bgmg_shared_library = 'H:\GitHub\BGMG\src\build_win\bin\debug\bgmg.dll';
    bgmg_shared_library_header = 'H:\GitHub\BGMG\src\bgmg_matlab.h';

    if ~libisloaded('bgmg'), fprintf('Loading bgmg library: %s, %s... ', bgmg_shared_library, bgmg_shared_library_header); loadlibrary(bgmg_shared_library, bgmg_shared_library_header, 'alias', 'bgmg');  fprintf('OK.\n'); end;

    calllib('bgmg', 'bgmg_init_log', ['test.bgmglib.log']);

    m2c = @(x)(x-1);
    check = @()fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', 0));
    calllib('bgmg', 'bgmg_set_tag_indices', 0, length(defvec), sum(defvec), m2c(find(defvec)));  check();
    calllib('bgmg', 'bgmg_set_option', 0,  'r2min', 0.01); check();
    calllib('bgmg', 'bgmg_set_option', 0,  'kmax', 10); check();
    calllib('bgmg', 'bgmg_set_option', 0,  'max_causals', length(defvec)); check();  
    calllib('bgmg', 'bgmg_set_option', 0,  'num_components', 1); check();
    calllib('bgmg', 'bgmg_set_option', 0,  'cache_tag_r2sum', true); check();

    calllib('bgmg', 'bgmg_set_ld_r2_coo', 0, length(r2), m2c(index_A), m2c(index_B), r2); check();
    calllib('bgmg', 'bgmg_set_ld_r2_csr', 0);  check();

    calllib('bgmg', 'bgmg_set_weights_randprune', 0, 10, 0.8);  check();
    if hvec_model == 1
        hvec = mafvec .* (1-mafvec) * 2;
    else
        hvec = ones(size(mafvec));
    end
    calllib('bgmg', 'bgmg_set_hvec', 0, length(hvec), hvec);  check();
    calllib('bgmg', 'bgmg_set_option', 0,  'diag', 0); check();

    pBuffer = libpointer('singlePtr', zeros(sum(defvec), 1, 'single'));
    calllib('bgmg', 'bgmg_retrieve_weights', 0, sum(defvec), pBuffer);  check(); bgmg_weights = pBuffer.Value;
    calllib('bgmg', 'bgmg_retrieve_ld_tag_r2_sum', 0, sum(defvec), pBuffer);  check(); bgmg_sum_r2 = pBuffer.Value;
    calllib('bgmg', 'bgmg_retrieve_ld_tag_r4_sum', 0, sum(defvec), pBuffer);  check(); bgmg_sum_r4 = pBuffer.Value;
    clear pBuffer

    if hvec_model == 1
        plot(bgmg_sum_r2, bgmg_sum_r2 - ldscore_sum_r2_per_allele(defvec), '.')
        plot(bgmg_sum_r4, bgmg_sum_r4 - ldscore_sum_r4_per_allele(defvec), '.')
    else
        plot(bgmg_sum_r2, bgmg_sum_r2 - ldscore_sum_r2(defvec), '.')
        plot(bgmg_sum_r4, bgmg_sum_r4 - ldscore_sum_r4(defvec), '.')
    end
end