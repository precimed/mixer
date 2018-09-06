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

    % divide or multiply argument by a factor of 3
    % to avoid 0 == log(1). Specialy suited for sig2_zero.
    function y = exp3_amd(x,s)
        if ~exist('s', 'var'), s = 1; end

        minval = eps;
        if s == 0
          x = minval+3*x;
          y = log(x);
        else
          y = exp(x);
          y = (y-exp(log(minval)))/3;
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
    
    function ov = UGMG_mapparams1_decorrelated_parametrization(iv, options)
        % mapparams for univariate mixture with a single causal component
        if ~exist('options', 'var'), options=[]; end;
        if ~isfield(options, 'sig2_zero'), options.sig2_zero = nan; end;
        if ~isfield(options, 'pi_vec'), options.pi_vec = nan; end;
        if ~isfield(options, 'sig2_beta'), options.sig2_beta = nan; end;
        if isnan(options.pi_vec) ~= isnan(options.sig2_beta), error('isnan(options.pi_vec) ~= isnan(options.sig2_beta)'); end;
        
        is_packing = isstruct(iv); cnti = 1;
        if is_packing, ov = []; else ov = struct(); end;

        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp3_amd, 'sig2_zero');
        if is_packing
            ov = cat(2, ov, log(atanh(iv.pi_vec)) + log(iv.sig2_beta));
            ov = cat(2, ov, log(atanh(iv.pi_vec)) - log(iv.sig2_beta));
        else
            ov.pi_vec = tanh(exp((iv(cnti) + iv(cnti+1))/2));
            ov.sig2_beta = exp((iv(cnti) - iv(cnti+1))/2);
        end
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
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp3_amd, 'sig2_zero');
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
        bivariate_ci_funcs.pi12_over_min_piXu = @(params)(params.pi_vec(3) / min(sum(params.pi_vec([1 3])), sum(params.pi_vec([2 3]))));
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
    
    function cost = UGMG_fminsearch_cost(ov, trait_index)
        cost = BGMG_cpp().calc_univariate_cost(trait_index, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
        BGMG_cpp.log('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
    end
    
    function cost = UGMG_CPP_fminsearch_cost(iv, trait_index)
        ov = BGMG_util.UGMG_mapparams1(iv);
        cost = BGMG_cpp().calc_univariate_cost(trait_index, ov.pi_vec, ov.sig2_zero, ov.sig2_beta);
        BGMG_cpp.log('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost);
    end
    
    function [cost, gradient] = UGMG_fminsearch_cost_with_gradient(ov, trait_index)
        if nargout == 1
          cost = BGMG_util.UGMG_fminsearch_cost(ov, trait_index);
          return
        end
        
        pBuffer = libpointer('doublePtr', zeros(3, 1, 'double'));
        cost = calllib('bgmg', 'bgmg_calc_univariate_cost_with_deriv', 0, trait_index, ov.pi_vec, ov.sig2_zero, ov.sig2_beta, 3, pBuffer);
        gradient = pBuffer.value;
        clear pBuffer
        BGMG_cpp.log('pi_vec=%.5e, sig2_zero=%.3f, sig2_beta=%.5e, cost=%.3f, deriv=%s\n', ov.pi_vec, ov.sig2_zero, ov.sig2_beta, cost, mat2str(gradient));
    end
    
    function cost = BGMG_fminsearch_cost(ov)
        if (length(ov.pi_vec) == 1)
            % pleiotropic model (one component with shared effects)
            cost = BGMG_cpp().calc_bivariate_cost([0 0 ov.pi_vec], ov.sig2_beta, ov.rho_beta, ov.sig2_zero, ov.rho_zero);
        elseif (length(ov.pi_vec) == 3)
            % full model
            cost = BGMG_cpp().calc_bivariate_cost(ov.pi_vec, ov.sig2_beta(:, 3), ov.rho_beta(3), ov.sig2_zero, ov.rho_zero);
        else
            error('not implemented');
        end

        filt = @(x)unique(x(x~=0));
        BGMG_cpp.log('Bivariate : pi_vec=[%s], rho_beta=[%s], sig2_beta1=[%s], sig2_beta2=[%s], rho_zero=%.3f, sig2_zero=[%s], cost=%.3e\n', ...
            sprintf('%.3e ', ov.pi_vec), ...
            sprintf('%.3f ', ov.rho_beta(end)), ...
            sprintf('%.2e ', filt(ov.sig2_beta(1, :))), ...
            sprintf('%.2e ', filt(ov.sig2_beta(2, :))), ...
            ov.rho_zero, ...
            sprintf('%.3f ', ov.sig2_zero), cost);
        if ~isfinite(cost), cost=1e99; end;
    end

    function ov = BGMG_mapparams1_rho(iv, options)
        % mapparams for bivaraite mixture with a one causal pleiotropic component
        % all params except rho_zero and rho_beta must be fixed by options
        if ~isfield(options, 'rho_zero'), options.rho_zero = nan; end;
        if ~isfield(options, 'rho_beta'), options.rho_beta = nan; end;

        is_packing = isstruct(iv); cnti = 1;
        if is_packing, ov = []; else ov = options; end;

        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_zero');
        [ov, ~] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_beta');
    end

    function ov = BGMG_mapparams3_decorrelated_parametrization(iv)
        % mapparams for BGMG model with 9 free parameters
        % - 3 parameters from UGMG_mapparams1_decorrelated_parametrization(trait1)
        % - 3 parameters from UGMG_mapparams1_decorrelated_parametrization(trait2)
        % - 3 parameters from BGMG_mapparams3_rho_and_pifrac

        is_packing = isstruct(iv);
        if is_packing, ov = []; else ov = struct(); end;
        
        if is_packing
            ov = cat(2, ov, BGMG_util.UGMG_mapparams1_decorrelated_parametrization(struct('pi_vec', sum(iv.pi_vec([1,3])), 'sig2_beta', iv.sig2_beta(1, 3), 'sig2_zero', iv.sig2_zero(1))));
            ov = cat(2, ov, BGMG_util.UGMG_mapparams1_decorrelated_parametrization(struct('pi_vec', sum(iv.pi_vec([2,3])), 'sig2_beta', iv.sig2_beta(2, 3), 'sig2_zero', iv.sig2_zero(2))));
            ov = cat(2, ov, BGMG_util.BGMG_mapparams3_rho_and_pifrac(iv));
        else
        	p1 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(iv(1:3));
            p2 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(iv(4:6));
            ov = BGMG_util.BGMG_mapparams3_rho_and_pifrac(iv(7:9), struct('pi_vec', [p1.pi_vec, p2.pi_vec], 'sig2_beta', [p1.sig2_beta, p2.sig2_beta], 'sig2_zero', [p1.sig2_zero, p2.sig2_zero]));
        end
    end

    function ov = BGMG_mapparams3_rho_and_pifrac(iv, options)
        % mapparams for BGMG model with 3 free parameters:
        % - pi12frac = pi12/min(pi1u, pi2u)
        % - rho_zero
        % - rho_beta
        %
        % options must contain univariate results: 
        % - pi_vec, vector 2x1, univariate polygenicitiyes pi1u, pi2u
        % - sig2_beta, vector 2x1
        % - sig2_zero, vector 2x1
        %
        % test
        % params  = struct('pi_vec', [0.1 0.2 0.3], 'rho_zero', 0.1, 'rho_beta', 0.2); 
        % options = struct('pi_vec', [sum(params.pi_vec([1 3])), sum(params.pi_vec([2 3]))], 'sig2_beta', [1e-2 1e-3], 'sig2_zero', [1.5 1.8]); 
        % x = BGMG_util.BGMG_mapparams3_rho_and_pifrac(params, options)
        % p = BGMG_util.BGMG_mapparams3_rho_and_pifrac(x, options)
        % options.rho_zero = 0.1; options.rho_beta = 0.2;

        if ~exist('options', 'var'), options=[]; end;
        if ~isfield(options, 'rho_zero'), options.rho_zero = nan; end;
        if ~isfield(options, 'rho_beta'), options.rho_beta = [0 0 nan]; end;

        transform_forward = 0;
        transform_backward = 1;
        is_packing = isstruct(iv); cnti = 1;
        if is_packing, ov = []; else ov = struct(); end;

        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_zero');
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_beta');
        
        if is_packing
            if ~isfield(options, 'pi_vec'), options.pi_vec = [sum(iv.pi_vec([1,3])), sum(iv.pi_vec([2,3]))]; end;
            assert(sum(iv.pi_vec([1,3])) == options.pi_vec(1));
            assert(sum(iv.pi_vec([2,3])) == options.pi_vec(2));
            if isfield(iv, 'sig2_beta') && isfield(options, 'sig2_beta'), assert(all(iv.sig2_beta(:, end) == BGMG_util.colvec(options.sig2_beta))); end;
            if isfield(iv, 'sig2_zero') && isfield(options, 'sig2_zero'), assert(all(iv.sig2_zero == BGMG_util.colvec(options.sig2_zero))); end;
            ov = cat(2, ov, BGMG_util.logit_amd(iv.pi_vec(end) / min(options.pi_vec), transform_forward));
        else
            pi12frac = BGMG_util.logit_amd(iv(cnti), transform_backward); cnti=cnti+1;
            pi12 = min(options.pi_vec) * pi12frac;
            ov.pi_vec = [options.pi_vec(1)-pi12, options.pi_vec(2)-pi12, pi12];
            ov.sig2_beta = [options.sig2_beta(1), 0, options.sig2_beta(1); ...
                            0, options.sig2_beta(2), options.sig2_beta(2)];
            ov.sig2_zero = BGMG_util.colvec(options.sig2_zero);
        end
    end
    
    function ov = BGMG_mapparams3_decorrelated_parametrization_9arguments(iv)
        % for ci computation (used only to unpack 9-component vector into param struct)
        p1 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(iv(1:3));
        p2 = BGMG_util.UGMG_mapparams1_decorrelated_parametrization(iv(4:6));
        ov = BGMG_util.BGMG_mapparams3_rho_and_pifrac(iv(7:9), struct(...
            'pi_vec', [p1.pi_vec, p2.pi_vec], ...
            'sig2_beta', [p1.sig2_beta, p2.sig2_beta], ...
            'sig2_zero', [p1.sig2_zero, p2.sig2_zero]));
    end
    

    function ov = BGMG_mapparams3(iv, options)
        % mapparams for saturated bivaraite mixture with a three causal component
        % (trait1-specific, trait2-specific, and pleiotropic components)
        % pleiotropic component re-use the same sig2_beta as trait-specific
        % components.

        if ~exist('options', 'var'), options=[]; end;
        if ~isfield(options, 'pi_vec'), options.pi_vec = [nan nan nan]; end;
        if ~isfield(options, 'sig2_zero'), options.sig2_zero = [nan; nan]; end;
        if ~isfield(options, 'sig2_beta'), options.sig2_beta = [nan; nan]; end;
        if ~isfield(options, 'rho_zero'), options.rho_zero = nan; end;
        if ~isfield(options, 'rho_beta'), options.rho_beta = [0 0 nan]; end;
        if all(size(options.sig2_beta) == [2 3]), options.sig2_beta = options.sig2_beta(:, 3); end;
        
        is_packing = isstruct(iv); cnti = 1;
        if is_packing, ov = []; else ov = struct(); end;

        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.logit_amd, 'pi_vec');
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp3_amd, 'sig2_zero');
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_zero');
        [ov, cnti] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.sigmf_of, 'rho_beta');

        % Tricks to map sig2_beta from two-vector into 2x3 matrix with two zero
        % elements and equal variances in trait-specific and pleiotropic components.
        if is_packing
            if isfield(iv, 'sig2_beta'), iv.sig2_beta = iv.sig2_beta(:,3); end;
            [ov, ~] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
        else
            [ov, ~] = BGMG_util.mapparams(iv, ov, cnti, options, @BGMG_util.exp_amd, 'sig2_beta');
            s1 = ov.sig2_beta(1); s2 = ov.sig2_beta(2);
            ov.sig2_beta = [s1 0 s1; 0 s2 s2];
        end
    end
    
    function defvec_tmp = find_hardprune_indices(defvec_tmp, hardprune_r2, mafvec, chrnumvec, hardprune_plink_ld_bin, chr_labels)
        % Use hard threshold to exlude sinonimous SNPs from fit. Just one
        % iteration of random pruning with very high r2 threshold. Non-selected
        % SNPs are excluded.
        BGMG_cpp.log('Excluding variants based on random pruning at %.3f threshold...\n', hardprune_r2);
        tag_indices_tmp = find(defvec_tmp);
        bgmglib=BGMG_cpp(1);
        bgmglib.dispose();
        bgmglib.defvec = defvec_tmp;
        bgmglib.mafvec = mafvec;
        bgmglib.chrnumvec = chrnumvec;
        for chr_index=1:length(chr_labels),
            bgmglib.set_ld_r2_coo_from_file(strrep(hardprune_plink_ld_bin,'@', sprintf('%i', chr_labels(chr_index))));
            bgmglib.set_ld_r2_csr(chr_labels(chr_index));
        end;
        hardprune_n = 1;
        bgmglib.set_weights_randprune(hardprune_n, hardprune_r2);
        weights_bgmg = bgmglib.weights;
        bgmglib.dispose();
        defvec_tmp(tag_indices_tmp(weights_bgmg==0)) = false;
        BGMG_cpp.log('Exclude %i variants after hard pruning at %.3f threshold (%i variants remain)\n', sum(weights_bgmg == 0), hardprune_r2, sum(defvec_tmp));
    end
  end
end

