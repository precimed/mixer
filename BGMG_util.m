classdef BGMG_util
  methods(Static)
    function [pi0, pi_vec] = normalize_pi_vec(pi_vec)
        num_mix  = size(pi_vec, 2);
        if (num_mix) == 1
            pi0 = 1-pi_vec;
        else
            pi0 = prod(1-pi_vec')';
            pi_vec = pi_vec .* (1-repmat(pi0, [1, num_mix])) ./ repmat(sum(pi_vec, 2), [1, num_mix]);
        end
    end

    function ok = validate_params(params)
        ok = true;

        % Validate required fields
        if ~isfield(params, 'pi_vec'), error('pi_vec field is required'); end;
        if ~isfield(params, 'sig2_zero'), error('sig2_zero field is required'); end;
        if ~isfield(params, 'sig2_beta'), error('sig2_beta field is required'); end;

        % Validate orientation and report size mismatches
        if size(params.pi_vec, 1) ~= 1, error('pi_vec must be row vector, one value per mixture component'); end;
        if size(params.sig2_zero, 2) ~= 1, error('sig2_zero must be column vector, one value per trait'); end;
        if size(params.pi_vec, 2) ~= size(params.sig2_beta, 2), error('size mismatch between pi_vec and sig2_beta'); end;
        if size(params.sig2_zero, 1) ~= size(params.sig2_beta, 1), error('size mismatch between sig2_zero and sig2_beta'); end;
        if isfield(params, 'rho_beta') && size(params.rho_beta, 1) ~= 1, error('rho_beta must be row vector, one value per mixture component'); end;
        if isfield(params, 'rho_beta') && size(params.rho_beta, 2) ~= size(params.pi_vec,2), error('size mismatch between pi_vec and rho_beta'); end;

        % Validate parameter values
        if any(params.pi_vec < 0), warning('pi_vec < 0'); ok = false; end;
        if any(params.pi_vec > 1), warning('pi_vec > 1'); ok = false; end;
        if any(params.sig2_zero(:) <= 0), warning('sig2_zero <= 0'); ok = false; end;
        if any(params.sig2_beta(:) < 0), warning('sig2_beta < 0'); ok = false; end;
        if isfield(params, 'rho_zero') && ((params.rho_zero < -1) || (params.rho_zero > 1)), warning('rho_zero must be in [-1, 1]'); ok = false; end;
        if isfield(params, 'rho_beta') && any((params.rho_beta < -1) | (params.rho_beta > 1)), warning('rho_beta must be in [-1, 1]'); ok = false; end;
    end

    function values = op_power(values, operation, power)
        % Operation power.
        % Calculates (operation^power)(values).
        %
        % Examples:
        % op_power(2,          @(x,y)(x*y), 10) = 2^10 = 1024
        % op_power(3,          @(x,y)(x*y), 4)  = 3^4  = 81
        % op_power(5,          @(x,y)(x+y), 7)  = 5*7  = 35
        % op_power([1 2; 2 1], @(A,B)(A*B), 2)  = [5 4; 4 5] - matrix power
        if (power == 0), error('not supported'); end;
        if (power ~= floor(power)), error('power must be integer'); end;
        bi = de2bi(power - 1);
        powN = values;
        for i=1:length(bi)
            if bi(i), values = operation(values, powN); end;
            powN = operation(powN, powN);
        end
    end

    function y = exp_amd(x,s)
        if ~exist('s', 'var'), s = 1; end

        minval = eps;
        if s == 0
          x = minval+x;
          y = log(x);
        else
          y = exp(x);
          y = (y-exp(log(minval)));
        end
    end

    function y = sigmf_of(x, s)
        if ~exist('s', 'var'), s = 0; end
        if s == 0
            y = BGMG_util.logit_amd((x + 1)/2, s);
        else
            y = BGMG_util.logit_amd(x, s) * 2 - 1;
        end
    end

    function y = logit_amd(x,s)
        if ~exist('s', 'var'), s = 0; end

        minval = eps;
        if s == 0
          x = minval+x*(1-2*minval);
          y = BGMG_util.logit(x,s);
        else
          y = BGMG_util.logit(x,s);
          y = (y-BGMG_util.logit(BGMG_util.logit(minval,0),1))/(1-2*minval);
        end
    end

    function y = logit(x,invflag)
        % logit transformation
        %
        % Input:
        % ------
        % x,            input value
        % invflag,      logical, inverse?
        %
        % Return:
        % ------
        % y,            transformed value
        if exist('invflag','var') & invflag
          y = exp(x)./(1+exp(x));
        else
          y = log(x./(1-x));
        end
    end

    function y = erf_amd(x,s)
        if ~exist('s', 'var'), s = 1; end

        minval = eps;
        if s == 0
          x = x*(1-minval);
          y = erfinv(x);
        else
          y = erf(x);
          y = y/(1-minval);
        end
    end

    function cost = UGMG_CPP_fminsearch_cost(iv)
        ov = BGMG_util.UGMG_mapparams1(iv);
        cost = calllib('bgmg', 'bgmg_calc_univariate_cost', 0, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
        fprintf('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
    end

    function ov = UGMG_mapparams1(iv, options)
        % mapparams for univariate mixture with a single causal component
        if ~exist('options', 'var'), options=[]; end;
        if ~isfield(options, 'pi_vec'), options.pi_vec = nan; end;
        if ~isfield(options, 'sig2_zero'), options.sig2_zero = nan; end;
        if ~isfield(options, 'sig2_beta'), options.sig2_beta = nan; end;

        is_packing = isstruct(iv); cnti = 1;
        if is_packing, ov = []; else ov = struct(); end;

        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.logit_amd, 'pi_vec');
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_zero');
        [ov, ~] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
    end

    function [ov, cnti] = mapparams(iv, ov, cnti, options, transform, field)
        idx_pack = isnan(options.(field));
        transform_forward = 0;
        transform_backward = 1;
        if isstruct(iv)
            % transform from struct to vector
            if any(idx_pack(:))
                ov = cat(2,ov,BGMG_util.rowvec(transform(iv.(field)(idx_pack),transform_forward)));
            end
        else
            % transform from vector to struct
            ov.(field) = options.(field);
            if any(idx_pack(:))
                ov.(field)(idx_pack) = transform(iv(cnti : (cnti + sum(idx_pack(:)) - 1)), transform_backward);
                cnti=cnti+sum(idx_pack(:));
            end
        end
    end

    function v = colvec(M)
        v = reshape(M,[numel(M) 1]);
    end

    function v = rowvec(M)
        v = reshape(M,[1 numel(M)]);
    end
    
    function [lamGC] = lamGCfromQQ(logqvec,hv_logp)

      qMedian=0.5;
      xMedian = -log10(qMedian);                % == 0.3;  At q==0.5 (i.e., x=0.3), 50% of SNPs lie abiove and 50% lie below.
      xd = BGMG_util.colvec(logqvec);                     % -log10(q) for Data.
      
      ixd = find(xd>xMedian,1);                 % q==0.5 is the median by definition -- with 50% proportion above and below.
      
      if ixd<=10 | isempty(ixd),  lamGC=nan;   return;  end
      if ixd>length(xd)-10,       lamGC=nan;   return;  end
      
      %xd = xd(ixd-1:ixd);                      % xd minimally brackets the median.
      %yd = hv_logp(ixd-1:ixd);                 % The y-axis values corresponding to xd.
      xd = xd(ixd-10:ixd+10);
      
      if length(find(~isfinite(xd))),  lamGC=nan;   return;  end
      
      yd = hv_logp(ixd-10:ixd+10);
      XD = [ones(length(xd),1) xd] ;
      b1 = XD\yd'; %'
      xx = xd;
      yy = XD*b1;  %   yy = b1(2)*xd + b1(1)
      %figure(99); plot(xx,yy,'*');
      ixd = find(xx>xMedian,1);                 % q==0.5 is the median by definition -- with 50% proportion above and below.
      xxd = xx(ixd-1:ixd);                      % xxd minimally brackets the median.
      yyd = yy(ixd-1:ixd);                      % The y-axis values corresponding to xxd.
      
      % (xd,yd) forms a jaggidy staircase whereas (xxd,yyd) forms a straight line.
      %ydMedian = interp1(xd,yd,xMedian);       % The y-axis value corresponding to xMedian.
      ydMedian = interp1(xxd,yyd,xMedian);      % The y-axis value corresponding to xMedian.
      pdMedian = 10^-ydMedian;                  % The data  p-value corresponding to qMedian;
      z2_medianFromQQ_data = chi2inv(1-pdMedian,1);  % Note, hv_logp = -log10(2*normcdf(-abs(hv_z))); z2_data_median is the z^2 value at ydMedian.
      z2_medianNull        = chi2inv(   qMedian,1);     % 0.454
      % Equivalently:
      %z2_medianFromQQ_data = norminv(pdMedian/2)^2;
      %z2_medianNull        = norminv( qMedian/2)^2;
      
      lamGC = z2_medianFromQQ_data/z2_medianNull;
    end

    function [univariate_ci_funcs, bivariate_ci_funcs, all_funcs] = find_extract_funcs(options)
        if ~exist('options', 'var'), options  = []; end
        if ~isfield(options, 'total_het'), options.total_het = nan; end;
        univariate_ci_funcs.sig2_zero         = @(params)(params.sig2_zero);
        univariate_ci_funcs.sig2_zero_minus1  = @(params)(params.sig2_zero - 1);
        univariate_ci_funcs.sig2_beta         = @(params)(params.sig2_beta);
        univariate_ci_funcs.pi_vec            = @(params)(params.pi_vec);
        univariate_ci_funcs.h2                = @(params)((params.sig2_beta*params.pi_vec')*options.total_het);

        bivariate_ci_funcs.sig2_zero_T1       = @(params)(params.sig2_zero(1));
        bivariate_ci_funcs.sig2_zero_T1_minus1= @(params)(params.sig2_zero(1) - 1);
        bivariate_ci_funcs.sig2_zero_T2       = @(params)(params.sig2_zero(2));
        bivariate_ci_funcs.sig2_zero_T2_minus1= @(params)(params.sig2_zero(2) - 1);
        bivariate_ci_funcs.pi_vec_C1          = @(params)(params.pi_vec(1));
        bivariate_ci_funcs.pi_vec_C2          = @(params)(params.pi_vec(2));
        bivariate_ci_funcs.pi_vec_C3          = @(params)(params.pi_vec(3));
        bivariate_ci_funcs.h2_T1              = @(params)((params.sig2_beta(1, :)*params.pi_vec')*options.total_het);
        bivariate_ci_funcs.h2_T2              = @(params)((params.sig2_beta(2, :)*params.pi_vec')*options.total_het);
        bivariate_ci_funcs.rho_zero           = @(params)(params.rho_zero);
        bivariate_ci_funcs.sig2_beta_T1       = @(params)(params.sig2_beta(1,1));
        bivariate_ci_funcs.sig2_beta_T2       = @(params)(params.sig2_beta(2,2));
        bivariate_ci_funcs.rho_beta           = @(params)(params.rho_beta(3));
        bivariate_ci_funcs.pi1u               = @(params)(sum(params.pi_vec([1 3])));
        bivariate_ci_funcs.pi2u               = @(params)(sum(params.pi_vec([2 3])));
        bivariate_ci_funcs.rg                 = @(params)(params.rho_beta(3) * params.pi_vec(3) / sqrt(sum(params.pi_vec([1 3])) * sum(params.pi_vec([2 3]))));
        bivariate_ci_funcs.pi12_minus_pi1u_times_pi2u = @(params)(params.pi_vec(3) - sum(params.pi_vec([1 3])) * sum(params.pi_vec([2 3])));
        bivariate_ci_funcs.pi12_over_pi1u     = @(params)(params.pi_vec(3) / sum(params.pi_vec([1 3])));
        bivariate_ci_funcs.pi12_over_pi2u     = @(params)(params.pi_vec(3) / sum(params.pi_vec([2 3])));
        bivariate_ci_funcs.pi1_over_pi1u      = @(params)(params.pi_vec(1) / sum(params.pi_vec([1 3])));
        bivariate_ci_funcs.pi2_over_pi2u      = @(params)(params.pi_vec(2) / sum(params.pi_vec([2 3])));
        bivariate_ci_funcs.pi1u_over_pi2u     = @(params)(sum(params.pi_vec([1 3])) / sum(params.pi_vec([2 3])));
        bivariate_ci_funcs.pi2u_over_pi1u     = @(params)(sum(params.pi_vec([2 3])) / sum(params.pi_vec([1 3])));
        bivariate_ci_funcs.h2pleio_T1         = @(params)(params.sig2_beta(1, 3)*params.pi_vec(3)*options.total_het);
        bivariate_ci_funcs.h2pleio_over_h2_T1 = @(params)((params.sig2_beta(1, 3)*params.pi_vec(3)) ./ (params.sig2_beta(1, :)*params.pi_vec'));
        bivariate_ci_funcs.h2pleio_T2         = @(params)(params.sig2_beta(2, 3)*params.pi_vec(3)*options.total_het);
        bivariate_ci_funcs.h2pleio_over_h2_T2 = @(params)((params.sig2_beta(2, 3)*params.pi_vec(3)) ./ (params.sig2_beta(2, :)*params.pi_vec'));
        
        fns = fieldnames(univariate_ci_funcs);
        for i=1:length(fns)
            all_funcs.(['UGMG_T1_', fns{i}]) = @(result)(univariate_ci_funcs.(fns{i})(result.univariate{1}.params));
            all_funcs.(['UGMG_T2_', fns{i}]) = @(result)(univariate_ci_funcs.(fns{i})(result.univariate{2}.params));
        end
        
        fns = fieldnames(bivariate_ci_funcs);
        for i=1:length(fns)
            all_funcs.(['BGMG_', fns{i}]) = @(result)(bivariate_ci_funcs.(fns{i})(result.bivariate.params));
        end
    end
    
    function ci = extract_ci_funcs(ci_params, ci_funcs, params, ci_alpha)
        ci_func_names = fieldnames(ci_funcs);
        for i=1:length(ci_func_names)
            ci_func_name = ci_func_names{i};
            ci_func      = ci_funcs.(ci_func_name);

            pe = ci_func(params);  % pe = point estimate

            if ~isempty(ci_params)
                dist = nan(length(ci_params), numel(pe));
                for j=1:length(ci_params), dist(j, :) = BGMG_util.rowvec(ci_func(ci_params{j})); end
            else
                dist = nan(2, numel(pe));
            end

            ci_result.point_estimate = pe;
            ci_result.mean = reshape(mean(dist), size(pe));
            ci_result.median = reshape(median(dist), size(pe));
            ci_result.lower = reshape(quantile(dist,     ci_alpha/2), size(pe));
            ci_result.upper = reshape(quantile(dist, 1 - ci_alpha/2), size(pe));
            ci_result.se = reshape(std(dist), size(pe));
            ci_result.pval = reshape(2*normcdf(-abs(ci_result.mean ./ ci_result.se)), size(pe));

            ci.(ci_func_name) = ci_result;
        end
    end

    % extract point estimates
    function [header, data] = result2str_point_estimates(result, options)
        [~, ~, funcs] = BGMG_util.find_extract_funcs(options);

        data = ''; header = '';
        ci_func_names = fieldnames(funcs);
        for j=1:length(ci_func_names)
            func = funcs.(ci_func_names{j});
            try
                val = mat2str(func(result));
            catch
                val = '';
            end
            data = sprintf('%s%s\t', data, val);
            header = sprintf('%s%s\t', header, ci_func_names{j});
        end
    end

    function result2str(fileID, result)
        cis = {};
        if isfield(result, 'univariate')
            for i=1:length(result.univariate)
                if isfield(result.univariate{i}, 'ci')
                    cis{end+1, 1} = result.univariate{i}.ci;
                end
            end
        end
        if isfield(result, 'bivariate') && isfield(result.bivariate, 'ci')
            cis{end+1, 1} = result.bivariate.ci;
        end

        for t=1:length(cis)
            fs=fieldnames(cis{t});
            for i=1:length(fs)
                values = cis{t}.(fs{i});

                element = values.mean;
                if (size(element, 1) > 1) && (size(element, 2) > 1), continue; end;

                for k=1:length(element)
                    cs = fieldnames(values);

                    header = ''; data = '';
                    header = [header 'measure\t'];
                    data = [data fs{i}];
                    if length(element) > 1, data = [data sprintf('(%i)', k)]; end;
                    data = [data '\t'];

                    for j=1:length(cs)
                        header = [header cs{j} '\t'];
                        value = values.(cs{j});
                        value = value(k);
                        data = [data mat2str(value) '\t'];
                    end

                    if i==1 && k==1, fprintf(fileID, [header '\n']); end
                    fprintf(fileID, [data '\n']);
                end
            end
            fprintf(fileID, '\n');
        end
    end
    
    function s = vec2str(vec, digits)
        if ~exist('digits', 'var'), digits = 6; end;
        format = ['%.', int2str(digits), 'f'];
        if length(vec) == 1
            s=sprintf(format, vec);
        else
            s=['[', sprintf(format, vec(1))];
            for i=2:length(vec)
                s=[s, ', ', sprintf(format, vec(i))];
            end
            s = [s ']'];
        end
    end

    function so = struct_to_display(si)
        so = struct();
        si_fields = fieldnames(si);
        for field_index=1:length(si_fields)
            field = si_fields{field_index};
            if size(si.(field), 1) > 1
                so.(field) = mat2str(si.(field), 3);
            else
                so.(field) = si.(field);
            end
        end
    end
    
    function loglike_trajectory = extract_univariate_loglike_trajectory() 
        num_loglike_entries = calllib('bgmg', 'bgmg_get_loglike_cache_size', 0);
        pBuffer_pivec = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2zero = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2beta = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_cost = libpointer('doublePtr', zeros(1, 1, 'double'));
        loglike_trajectory=[];
        loglike_trajectory.pivec = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2zero = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2beta = zeros(num_loglike_entries, 1);
        loglike_trajectory.cost = zeros(num_loglike_entries, 1);
        for i=1:num_loglike_entries
            calllib('bgmg', 'bgmg_get_loglike_cache_univariate_entry', 0, i-1, pBuffer_pivec, pBuffer_sig2zero, pBuffer_sig2beta, pBuffer_cost);  check(); 
            loglike_trajectory.pivec(i) = pBuffer_pivec.Value;
            loglike_trajectory.sig2zero(i) = pBuffer_sig2zero.Value;
            loglike_trajectory.sig2beta(i) = pBuffer_sig2beta.Value;
            loglike_trajectory.cost(i) = pBuffer_cost.Value;
        end
        clear pBuffer_pivec pBuffer_sig2zero pBuffer_sig2beta pBuffer_cost
    end
    
    function loglike_trajectory = extract_bivariate_loglike_trajectory() 
        num_loglike_entries = calllib('bgmg', 'bgmg_get_loglike_cache_size', 0);
        pBuffer_pivec = libpointer('singlePtr', zeros(1, 3, 'single'));
        pBuffer_sig2beta = libpointer('singlePtr', zeros(1, 2, 'single'));
        pBuffer_rho_beta = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_sig2zero = libpointer('singlePtr', zeros(1, 2, 'single'));
        pBuffer_rho_zero = libpointer('singlePtr', zeros(1, 1, 'single'));
        pBuffer_cost = libpointer('doublePtr', zeros(1, 1, 'double'));
        loglike_trajectory=[];
        loglike_trajectory.pivec = zeros(num_loglike_entries, 3);
        loglike_trajectory.sig2beta = zeros(num_loglike_entries, 2);
        loglike_trajectory.rho_beta = zeros(num_loglike_entries, 1);
        loglike_trajectory.sig2zero = zeros(num_loglike_entries, 2);
        loglike_trajectory.rho_zero = zeros(num_loglike_entries, 1);
        loglike_trajectory.cost = zeros(num_loglike_entries, 1);
        for i=1:num_loglike_entries
            calllib('bgmg', 'bgmg_get_loglike_cache_bivariate_entry', 0, i-1, 3, pBuffer_pivec, 2, pBuffer_sig2beta, pBuffer_rho_beta, 2, pBuffer_sig2zero, pBuffer_rho_zero, pBuffer_cost);  check(); 
            loglike_trajectory.pivec(i, :) = pBuffer_pivec.Value;
            loglike_trajectory.sig2beta(i, :) = pBuffer_sig2beta.Value;
            loglike_trajectory.rho_beta(i) = pBuffer_rho_beta.Value;
            loglike_trajectory.sig2zero(i, :) = pBuffer_sig2zero.Value;
            loglike_trajectory.rho_zero(i) = pBuffer_rho_zero.Value;
            loglike_trajectory.cost(i) = pBuffer_cost.Value;
        end
        result.loglike_fit_trajectory = loglike_trajectory;
        clear pBuffer_pivec pBuffer_sig2zero pBuffer_sig2beta pBuffer_cost pBuffer_rho_beta pBuffer_rho_zero

    end
    
  end
end

