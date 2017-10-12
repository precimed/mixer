classdef BGMG_util
methods(Static)
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
          y = logit(x,s);
        else
          y = logit(x,s);
          y = (y-logit(logit(minval,0),1))/(1-2*minval);
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
end
end
