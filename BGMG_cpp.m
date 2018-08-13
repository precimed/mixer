classdef BGMG_cpp
  properties
    Context
    verbose
  end
  properties (Dependent)
    defvec              % instead of tag_indices
    num_snp
    num_tag
    chrnumvec
    mafvec
    weights
    zvec1
    zvec2
    nvec1
    nvec2
    ld_tag_r2_sum
    ld_tag_r4_sum
  end
  
  methods
    function obj = BGMG_cpp(context)
        if ~libisloaded('bgmg'), error('call BGMG_cpp.load(...) before creating BGMG_cpp objects'); end;
        if ~exist('context', 'var'), context = 0; end
        obj.Context = context;
        obj.verbose = true;
    end

    function check(obj)
        if obj.verbose, fprintf('RESULT: %s; STATUS: %s\n', calllib('bgmg', 'bgmg_get_last_error'), calllib('bgmg', 'bgmg_status', obj.Context)); end;
    end
    
    function dispose(obj)
        calllib('bgmg', 'bgmg_dispose', obj.Context); obj.check();
    end
    
    % set/get defvec
    function obj = set.defvec(obj, val)
        calllib('bgmg', 'bgmg_set_tag_indices', obj.Context, length(val), sum(val), find(val)-1); obj.check();
    end
    function val = get.defvec(obj)
        pBuffer = libpointer('int32Ptr', zeros(obj.num_tag, 1, 'int32'));
        calllib('bgmg', 'bgmg_retrieve_tag_indices', obj.Context, obj.num_tag, pBuffer); obj.check();
        tag_indices = pBuffer.Value; clear pBuffer
        tag_indices = tag_indices + 1;  % convert to matlab unity-based indices
        val = false(obj.num_snp, 1);    % generate defvec
        val(tag_indices) = true;
    end
    function val = get.num_tag(obj)
        val = calllib('bgmg', 'bgmg_get_num_tag', obj.Context); obj.check();
    end
    function val = get.num_snp(obj)
        val = calllib('bgmg', 'bgmg_get_num_snp', obj.Context); obj.check();
    end

    % get/set chrnumvec
    function obj = set.chrnumvec(obj, val)
        calllib('bgmg', 'bgmg_set_chrnumvec', obj.Context, length(val), val);  obj.check();
    end
    function val = get.chrnumvec(obj)
        pBuffer = libpointer('int32Ptr', zeros(obj.num_snp, 1, 'int32'));
        calllib('bgmg', 'bgmg_retrieve_chrnumvec', obj.Context, obj.num_snp, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    
    % set/get mafvec
    function obj = set.mafvec(obj, val)
        calllib('bgmg', 'bgmg_set_mafvec', obj.Context, length(val), val);  obj.check();
    end
    function val = get.mafvec(obj)
        pBuffer = libpointer('singlePtr', zeros(obj.num_snp, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_mafvec', obj.Context, obj.num_snp, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    
    % set/get zvec
    function obj = set_zvec(obj, trait_index, val)
        calllib('bgmg', 'bgmg_set_zvec', obj.Context, trait_index, length(val), val);  obj.check();
    end
    function val = get_zvec(obj, trait_index)
        pBuffer = libpointer('singlePtr', zeros(obj.num_tag, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_zvec', obj.Context, trait_index, obj.num_tag, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    function val = get.zvec1(obj)
        trait_index = 1; val = obj.get_zvec(trait_index);
    end
    function obj = set.zvec1(obj, val)
        trait_index = 1; obj.set_zvec(trait_index, val);
    end
    function val = get.zvec2(obj)
        trait_index = 2; val = obj.get_zvec(trait_index);
    end
    function obj = set.zvec2(obj, val)
        trait_index = 2; obj.set_zvec(trait_index, val);
    end

    % set/get nvec
    function obj = set_nvec(obj, trait_index, val)
        calllib('bgmg', 'bgmg_set_nvec', obj.Context, trait_index, length(val), val);  obj.check();
    end
    function val = get_nvec(obj, trait_index)
        pBuffer = libpointer('singlePtr', zeros(obj.num_tag, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_nvec', obj.Context, trait_index, obj.num_tag, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    function val = get.nvec1(obj)
        trait_index = 1; val = obj.get_nvec(trait_index);
    end
    function obj = set.nvec1(obj, val)
        trait_index = 1; obj.set_nvec(trait_index, val);
    end
    function val = get.nvec2(obj)
        trait_index = 2; val = obj.get_nvec(trait_index);
    end
    function obj = set.nvec2(obj, val)
        trait_index = 2; obj.set_nvec(trait_index, val);
    end
    
    % set/get weights
    function set_weights_randprune(obj, randprune_n, randprune_r2)
        calllib('bgmg', 'bgmg_set_weights_randprune', obj.Context, randprune_n, randprune_r2); obj.check();
    end
    function obj = set.weights(obj, val)
        calllib('bgmg', 'bgmg_set_weights', obj.Context, length(val), val);  obj.check();
    end
    function val = get.weights(obj)
        pBuffer = libpointer('singlePtr', zeros(obj.num_tag, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_weights', obj.Context, obj.num_tag, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end

    % get ld scores (sum of r2 and r4)
    function val = get.ld_tag_r2_sum(obj)
        pBuffer = libpointer('singlePtr', zeros(obj.num_tag, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_ld_tag_r2_sum', obj.Context, obj.num_tag, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    function val = get.ld_tag_r4_sum(obj)
        pBuffer = libpointer('singlePtr', zeros(obj.num_tag, 1, 'single'));
        calllib('bgmg', 'bgmg_retrieve_ld_tag_r4_sum', obj.Context, obj.num_tag, pBuffer); obj.check();
        val = double(pBuffer.Value); clear pBuffer
    end
    
    % set LD structure
    function set_ld_r2_coo_from_matlab_indices(obj, index_A, index_B, r2)
        % convert matlab to c++ indices
        index_A = index_A - 1;
        index_B = index_B - 1;
        calllib('bgmg', 'bgmg_set_ld_r2_coo', obj.Context, length(r2), index_A, index_B, r2); obj.check(); 
    end
    function set_ld_r2_coo_from_file(obj, filename)
        fprintf('Loading %s...', filename);
        calllib('bgmg', 'bgmg_set_ld_r2_coo_from_file', obj.Context, filename); obj.check(); 
        fprintf('OK.\n'); 
    end
    function set_ld_r2_csr(obj, chr_label)
        if ~exist('chr_label', 'var'), chr_label = -1; end;  % -1 means to finalize all chromosomes
        calllib('bgmg', 'bgmg_set_ld_r2_csr', obj.Context, chr_label); obj.check();
    end
    
    % get LD structure
    function [tag, r2] = get_ld_r2_snp(obj, snp_index)
      num_ld_r2 = calllib('bgmg', 'bgmg_num_ld_r2_snp', obj.Context, snp_index); obj.check();
      pBuffer_tag = libpointer('int32Ptr', zeros(num_ld_r2, 1, 'int32'));
      pBuffer_r2 = libpointer('singlePtr', zeros(num_ld_r2, 1, 'single'));
      calllib('bgmg', 'bgmg_retrieve_ld_r2_snp', obj.Context, snp_index, num_ld_r2, pBuffer_tag, pBuffer_r2); obj.check();
      tag = pBuffer_tag.Value; clear pBuffer_tag
      r2 = pBuffer_r2.Value; clear pBuffer_r2
    end

    function [snp, tag, r2] = get_ld_r2_chr(obj, chr_label)
      num_ld_r2 = calllib('bgmg', 'bgmg_num_ld_r2_chr', obj.Context, chr_label); obj.check();
      pBuffer_snp = libpointer('int32Ptr', zeros(num_ld_r2, 1, 'int32'));
      pBuffer_tag = libpointer('int32Ptr', zeros(num_ld_r2, 1, 'int32'));
      pBuffer_r2 = libpointer('singlePtr', zeros(num_ld_r2, 1, 'single'));
      calllib('bgmg', 'bgmg_retrieve_ld_r2_chr', obj.Context, chr_label, num_ld_r2, pBuffer_snp, pBuffer_tag, pBuffer_r2); obj.check();
      snp = pBuffer_snp.Value; clear pBuffer_snp
      tag = pBuffer_tag.Value; clear pBuffer_tag
      r2 = pBuffer_r2.Value; clear pBuffer_r2
    end

    % set option
    function set_option(obj, option, value)
        calllib('bgmg', 'bgmg_set_option', obj.Context, option, value); obj.check();
    end
    
    % log likelihood cache
    function clear_loglike_cache(obj)
        calllib('bgmg', 'bgmg_clear_loglike_cache', obj.Context); obj.check();
    end
    function loglike_trajectory = extract_univariate_loglike_trajectory(obj) 
        num_loglike_entries = calllib('bgmg', 'bgmg_get_loglike_cache_size', obj.Context);
        pBuffer_pivec = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2_zero = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2_beta = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_cost = libpointer('doublePtr', zeros(1, 1, 'double'));
        loglike_trajectory=[];
        loglike_trajectory.pivec = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2_zero = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2_beta = zeros(num_loglike_entries, 1);
        loglike_trajectory.cost = zeros(num_loglike_entries, 1);
        for i=1:num_loglike_entries
            calllib('bgmg', 'bgmg_get_loglike_cache_univariate_entry', obj.Context, i-1, pBuffer_pivec, pBuffer_sig2_zero, pBuffer_sig2_beta, pBuffer_cost);  obj.check(); 
            loglike_trajectory.pivec(i) = pBuffer_pivec.Value;
            loglike_trajectory.sig2_zero(i) = pBuffer_sig2_zero.Value;
            loglike_trajectory.sig2_beta(i) = pBuffer_sig2_beta.Value;
            loglike_trajectory.cost(i) = pBuffer_cost.Value;
        end
        loglike_trajectory.cost(loglike_trajectory.cost > 1e99) = nan;
        clear pBuffer_pivec pBuffer_sig2_zero pBuffer_sig2_beta pBuffer_cost
    end
    function loglike_trajectory = extract_bivariate_loglike_trajectory(obj) 
        num_loglike_entries = calllib('bgmg', 'bgmg_get_loglike_cache_size', obj.Context);
        pBuffer_pivec = libpointer('singlePtr', zeros(1, 3, 'single'));
        pBuffer_sig2_beta = libpointer('singlePtr', zeros(1, 2, 'single'));
        pBuffer_rho_beta = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2_zero = libpointer('singlePtr', zeros(1, 2, 'single'));
        pBuffer_rho_zero = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_cost = libpointer('doublePtr', zeros(1, 1, 'double'));
        loglike_trajectory=[];
        loglike_trajectory.pivec = zeros(num_loglike_entries, 3);
        loglike_trajectory.sig2_beta = zeros(num_loglike_entries, 2);
        loglike_trajectory.rho_beta = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2_zero = zeros(num_loglike_entries, 2);
        loglike_trajectory.rho_zero = zeros(num_loglike_entries, 1);
        loglike_trajectory.cost = zeros(num_loglike_entries, 1);
        for i=1:num_loglike_entries
            calllib('bgmg', 'bgmg_get_loglike_cache_bivariate_entry', obj.Context, i-1, 3, pBuffer_pivec, 2, pBuffer_sig2_beta, pBuffer_rho_beta, 2, pBuffer_sig2_zero, pBuffer_rho_zero, pBuffer_cost);  obj.check(); 
            loglike_trajectory.pivec(i, :) = pBuffer_pivec.Value;
            loglike_trajectory.sig2_beta(i, :) = pBuffer_sig2_beta.Value;
            loglike_trajectory.rho_beta(i) = pBuffer_rho_beta.Value;
            loglike_trajectory.sig2_zero(i, :) = pBuffer_sig2_zero.Value;
            loglike_trajectory.rho_zero(i) = pBuffer_rho_zero.Value;
            loglike_trajectory.cost(i) = pBuffer_cost.Value;
        end
        loglike_trajectory.cost(loglike_trajectory.cost > 1e99) = nan;
        clear pBuffer_pivec pBuffer_sig2_zero pBuffer_sig2_beta pBuffer_cost pBuffer_rho_beta pBuffer_rho_zero
    end
    
    % cost functions
    function cost = calc_univariate_cost(obj, trait_index, pi_vec, sig2_zero, sig2_beta)
        cost = calllib('bgmg', 'bgmg_calc_univariate_cost', obj.Context, trait_index, pi_vec, sig2_zero, sig2_beta); obj.check();
    end
    function pdf  = calc_univariate_pdf(obj, trait_index, pi_vec, sig2_zero, sig2_beta, zgrid)
        pBuffer = libpointer('singlePtr', zeros(length(zgrid), 1, 'single'));
        calllib('bgmg', 'bgmg_calc_univariate_pdf', obj.Context, trait_index, pi_vec, sig2_zero, sig2_beta, length(zgrid), zgrid, pBuffer);  obj.check(); 
        pdf = pBuffer.Value'; clear pBuffer
    end
    function cost = calc_bivariate_cost(obj, pi_vec, sig2_beta, rho_beta, sig2_zero, rho_zero)
        cost = calllib('bgmg', 'bgmg_calc_bivariate_cost', obj.Context, 3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero); obj.check();
    end
    function pdf  = calc_bivariate_pdf(obj, pi_vec, sig2_beta, rho_beta, sig2_zero, rho_zero, zgrid1, zgrid2)
        pBuffer = libpointer('singlePtr', zeros(length(zgrid1), 1, 'single'));
        calllib('bgmg', 'bgmg_calc_bivariate_pdf', obj.Context, length(pi_vec), pi_vec, length(sig2_beta), sig2_beta, rho_beta, length(sig2_zero), sig2_zero, rho_zero, length(zgrid1), zgrid1, zgrid2, pBuffer); obj.check(); 
        pdf = reshape(pBuffer.Value', size(zgrid1)); clear pBuffer
    end
  end
  
  methods(Static)
    function load(bgmg_shared_library, bgmg_shared_library_header)
        if ~libisloaded('bgmg'), 
            fprintf('Loading bgmg library: %s, %s... ', bgmg_shared_library, bgmg_shared_library_header); 
            loadlibrary(bgmg_shared_library, bgmg_shared_library_header, 'alias', 'bgmg');  
            fprintf('OK.\n');
        end;
    end
    
    function unload
        if libisloaded('bgmg'), unloadlibrary('bgmg'); end;
    end
    
    function init_log(log_file)
        calllib('bgmg', 'bgmg_init_log', log_file);
    end
    
    function dispose_all_context_ids
        calllib('bgmg', 'bgmg_dispose', -1);
    end
    
    function log(varargin)
        fprintf(varargin{:});
        if libisloaded('bgmg'), calllib('bgmg', 'bgmg_log_message', regexprep(sprintf(varargin{:}),'[\n\r]+', '\t')); end;
    end
  end
end

% ToDo:
% DLL_PUBLIC int64_t bgmg_retrieve_tag_r2_sum(int context_id, int component_id, float num_causal, int length, float* buffer);
% DLL_PUBLIC int64_t bgmg_retrieve_ld_tag_r2_sum(int context_id, int length, float* buffer);  // LD scores (r2 and r4)
% DLL_PUBLIC int64_t bgmg_retrieve_ld_tag_r4_sum(int context_id, int length, float* buffer);
% DLL_PUBLIC int64_t bgmg_retrieve_weighted_causal_r2(int context_id, int length, float* buffer);
% DLL_PUBLIC double bgmg_calc_univariate_cost_with_deriv(int context_id, int trait_index, double pi_vec, double sig2_zero, double sig2_beta, int deriv_length, double* deriv);
