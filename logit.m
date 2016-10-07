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
