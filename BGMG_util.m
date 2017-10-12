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

    function v = colvec(M)
        v = reshape(M,[numel(M) 1]);
    end

    function v = rowvec(M)
        v = reshape(M,[1 numel(M)]);
    end
  end
end
