'''
(c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''
# Utility classes for univariate and bivariate fit
# Contains
# _log_exp_converter, _logit_logistic_converter, _arctanh_tanh_converter - converters to map bounded parameters into -inf, +inf range
# UnivariateParams, BivariateParams - represent parameters, with basic functionality like "calculate cost"
# Several univariate and bivariate parametrizations - suitable for fitting and (some specific parametrization) for uncertainty calculation
# _calculate_univariate_uncertainty, _calculate_bivariate_uncertainty - estimates confidence intervals for parameters and their aggregates (like h2 or rg) 
 
import numpy as np
import numdifftools as nd

epsval = np.finfo(float).eps
minval = np.finfo(float).min
maxval = np.finfo(float).max

def _log_bounded(x): # transforms, roughtly speaking, [0.00001, 100000] to [-10, 10]
    if not np.isfinite(x): return x
    if x<epsval: x=epsval
    if x>maxval: x=maxval
    y = np.log(x)
    return y

def _exp_bounded(x): # transforms, roughtly speaking, [-10, 10] to [0.0001, 100000]
    if not np.isfinite(x): return x
    y = np.exp(x)
    if y<epsval: y=epsval
    if y>maxval: y=maxval
    return y

def _logit_bounded(x): # transforms, roughly speaking, [0.00001, 0.99999] to [-10, 10]
    if not np.isfinite(x): return x
    if x<epsval: x=epsval
    if x>(1-epsval): x=1-epsval
    y = _log_bounded(x / (1-x))
    return y
    
def _logistic_bounded(x): # transforms, roughly speaking, [-10, 10] to [0.00001, 0.99999]
    if not np.isfinite(x): return x
    y = _exp_bounded(x) / (1 + _exp_bounded(x))
    if y<epsval: y=epsval
    if y>(1-epsval): y=1-epsval
    return y
    
# converter that maps [0, +inf] domain to [-inf, inf], and back if infvlag=True
def _log_exp_converter(x, invflag=False):
    return _exp_bounded(x) if invflag else _log_bounded(x)

# converter that maps [0, 1] domain to [-inf, inf], and back if infvlag=True
def _logit_logistic_converter(x, invflag=False):
    return _logistic_bounded(x) if invflag else _logit_bounded(x)

# converter that maps [-1, 1] domain to [-inf, inf], and back if infvlag=True
# this has an additional property that arctanh(X) ~ X for small values of X.
def _arctanh_tanh_converter(x, invflag=False):
    return (2.0 * _logistic_bounded(2.0 * x) - 1.0) if invflag else 0.5*_logit_bounded(0.5*x + 0.5)

class UnivariateParams(object):
    def __init__(self, pi, sig2_beta, sig2_zero):
        self._pi = pi
        self._sig2_beta = sig2_beta
        self._sig2_zero = sig2_zero
        self._validate()
        
    def _validate(self):
        for val in [self._pi, self._sig2_beta, self._sig2_zero]:
            assert np.isscalar(val)
            assert np.greater_equal(val, 0)
        assert np.less_equal(self._pi, 1)
    
    def __str__(self):
        description = []
        for attr_name in '_pi', '_sig2_beta', '_sig2_zero':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'UnivariateParams({})'.format(', '.join(description))
    __repr__ = __str__

    def as_dict(self):
        return {'pi': self._pi, 'sig2_beta': self._sig2_beta, 'sig2_zero': self._sig2_zero}
    
    def cost(self, lib, trait):
        value = lib.calc_univariate_cost(trait, self._pi, self._sig2_zero, self._sig2_beta)
        return value if np.isfinite(value) else 1e100

class BivariateParams(object):
    def __init__(self, pi=None, sig2_beta=None, rho_beta=None, sig2_zero=None, rho_zero=None, params1=None, params2=None, pi12=None):
        if (params1 is not None) and (params2 is not None) and (pi12 is not None):
            self._pi = [params1._pi - pi12, params2._pi - pi12, pi12]
            self._sig2_beta = [params1._sig2_beta, params2._sig2_beta]
            self._sig2_zero = [params1._sig2_zero, params2._sig2_zero]
        else:
            self._pi = pi
            self._sig2_beta = sig2_beta
            self._sig2_zero = sig2_zero
        self._rho_beta = rho_beta
        self._rho_zero = rho_zero
        self._validate()

    def _rg(self):
        pi1u = self._pi[0] + self._pi[2]
        pi2u = self._pi[1] + self._pi[2]
        return self._rho_beta * self._pi[2] / np.sqrt(pi1u * pi2u)

    def _validate(self):
        assert len(self._pi) == 3
        assert len(self._sig2_beta) == 2
        assert len(self._sig2_zero) == 2
        assert np.isscalar(self._rho_zero)
        assert np.isscalar(self._rho_beta)
        
        assert np.greater_equal(self._pi, 0).all()
        assert np.less_equal(np.sum(self._pi), 1.0)

        assert np.greater_equal(self._sig2_beta, 0).all()
        assert np.greater_equal(self._sig2_zero, 0).all()
        
        assert np.greater_equal([self._rho_zero, self._rho_beta], -1).all()
        assert np.less_equal([self._rho_zero, self._rho_beta], 1).all()

    def __str__(self):
        description = []
        for attr_name in '_pi', '_sig2_beta', '_rho_beta', '_sig2_zero', '_rho_zero':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        description.append('rg: {}'.format(self._rg()))
        return 'BivariateParams({})'.format(', '.join(description))
    __repr__ = __str__
    
    def as_dict(self):
        return {'pi': self._pi, 'sig2_beta': self._sig2_beta, 'sig2_zero': self._sig2_zero,
                'rho_zero': self._rho_zero, 'rho_beta': self._rho_beta}

    def cost(self, lib):
        value = lib.calc_bivariate_cost(self._pi, self._sig2_beta, self._rho_beta, self._sig2_zero, self._rho_zero)
        return value if np.isfinite(value) else 1e100

# TBD: there is some inconsistency across "Parametrization" classes.
# The following classes are OK:
#   - UnivariateParametrization_natural_axis
#   - BivariateParametrization_constUNIVARIATE_natural_axis
#   - BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO
# All other classes should be corrected by removing the "fit" function, and moving this logic into mixer.py apply_fit_sequence

# Parametrization that constrains polygenicity
# Typical usage is as follows:
#   optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
#   params0, details = UnivariateParametrization_constPI(1.0, 1.5, 1e-4, libbgmg, 1).fit(optimizer)
class UnivariateParametrization_constPI(object):
    def __init__(self, const_pi, init_sig2_zero, init_sig2_beta, lib, trait):
        self._init_vec = [_log_exp_converter(init_sig2_beta, invflag=False),
                          _log_exp_converter(init_sig2_zero, invflag=False)] 
        self._const_pi = const_pi
        self._lib = lib
        self._trait = trait

    def _vec_to_params(self, vec):
        return UnivariateParams(pi=self._const_pi,
                                sig2_beta=_log_exp_converter(vec[0], invflag=True),
                                sig2_zero=_log_exp_converter(vec[1], invflag=True))
    
    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib, self._trait)
    
    # optimizer should return a result which contains an "x" field
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

# Parametrization that constraints pi and sig2_beta parameters
class UnivariateParametrization_constPI_constSIG2BETA(object):
    def __init__(self, init_sig2_zero, const_params, lib, trait):
        # const_params is params to derive constraints
        self._init_vec = [_log_exp_converter(init_sig2_zero, invflag=False)]
        self._const_params = const_params
        self._lib = lib
        self._trait = trait

    def _vec_to_params(self, vec):
        return UnivariateParams(pi=self._const_params._pi,
                                sig2_beta=self._const_params._sig2_beta,
                                sig2_zero=_log_exp_converter(vec[0], invflag=True))
    
    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib, self._trait)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

# Parametrization that constraints sig2_zero parameter, and a product of sig2_beta and pi
class UnivariateParametrization_constH2_constSIG2ZERO(object):
    def __init__(self, init_pi, const_params, lib, trait):
        # const_params is params to derive constraints
        self._init_vec = [_logit_logistic_converter(init_pi, invflag=False)]
        self._const_params = const_params
        self._lib = lib
        self._trait = trait

    def _vec_to_params(self, vec):
        pi = _logit_logistic_converter(vec[0], invflag=True)
        sig2_beta = (self._const_params._sig2_beta * self._const_params._pi) / pi
        return UnivariateParams(pi=pi,
                                sig2_beta=sig2_beta,
                                sig2_zero=self._const_params._sig2_zero)
    
    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib, self._trait)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

# Parametrization that constraints sig2_zero parameter, and a product of sig2_beta and pi
class UnivariateParametrization_constH2_constSIG2ZERO_boundedPI(object):
    def __init__(self, const_params, max_pi, lib, trait):
        # const_params is params to derive constraints
        self._const_params = const_params
        self._lib = lib
        self._trait = trait
        self._max_pi = max_pi

    def _vec_to_params(self, pi):
        if pi<epsval: pi=epsval
        if pi>self._max_pi: pi=self._max_pi
        sig2_beta = (self._const_params._sig2_beta * self._const_params._pi) / pi
        return UnivariateParams(pi=pi,
                                sig2_beta=sig2_beta,
                                sig2_zero=self._const_params._sig2_zero)
    
    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib, self._trait)
    
    def fit(self, scalar_optimizer):
        result = scalar_optimizer(self._calc_cost)
        return self._vec_to_params(result.x), result

# Unconstrained parametrization with "independent axis", i.e.
#   x1 = log(sig2_zero)
#   x2 = log(atanh(pi)) + log(sig2_beta)
#   x3 = log(atanh(pi)) - log(sig2_beta)
# The reason for atanh(pi) in the formulas is to make sure that the inverse transform always give a valid pi (e.g. between zero to one)
# atanh is particularly helpful because atanh(x)~x for small x, and we expect pi to be small.
class UnivariateParametrization(object):
    def __init__(self, init_params, lib, trait):
        self._init_vec = self._params_to_vec(init_params)
        self._lib = lib
        self._trait = trait

    def _params_to_vec(self, params):
        arctanh_pi = _arctanh_tanh_converter(params._pi, invflag=False)
        return [_log_exp_converter(params._sig2_zero, invflag=False),
                _log_exp_converter(arctanh_pi, invflag=False) + _log_exp_converter(params._sig2_beta, invflag=False),
                _log_exp_converter(arctanh_pi, invflag=False) - _log_exp_converter(params._sig2_beta, invflag=False)]
        
    def _vec_to_params(self, vec):
        sig2_zero=_log_exp_converter(vec[0], invflag=True)
        sig2_beta=_log_exp_converter((vec[1] - vec[2]) / 2.0, invflag=True)
        pi= _arctanh_tanh_converter(_log_exp_converter((vec[1] + vec[2]) / 2.0, invflag=True), invflag=True)
        return UnivariateParams(pi=pi, sig2_beta=sig2_beta, sig2_zero=sig2_zero)
    
    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib, self._trait)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

# Standard univariate parametrization
#   x1 = log(sig2_zero)
#   x2 = log(sig2_beta)
#   x3 = logit(pi)
class UnivariateParametrization_natural_axis(object):
    def __init__(self, lib, trait):
        self._lib = lib
        self._trait = trait

    def params_to_vec(self, params):
        return [_log_exp_converter(params._sig2_zero, invflag=False),
                _log_exp_converter(params._sig2_beta, invflag=False),
                _logit_logistic_converter(params._pi, invflag=False)]
        
    def vec_to_params(self, vec):
        sig2_zero=_log_exp_converter(vec[0], invflag=True)
        sig2_beta=_log_exp_converter(vec[1], invflag=True)
        pi = _logit_logistic_converter(vec[2], invflag=True)
        return UnivariateParams(pi=pi, sig2_beta=sig2_beta, sig2_zero=sig2_zero)
    
    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)
    
def _max_rg(pi1u, pi2u):
    return min(pi1u, pi2u) / np.sqrt(pi1u * pi2u)

# BGMG_cpp_fit_bivariate_fast (fits rho_zero, rho_beta, with constrain imposed by maxRG)
class BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG(object):
    def __init__(self, const_sig2_beta, const_sig2_zero, max_rg, init_rho_beta, init_rho_zero, lib):
        assert abs(init_rho_beta) <= abs(max_rg)
        assert abs(init_rho_zero) <= 1
        self._init_vec = [_arctanh_tanh_converter(init_rho_beta / abs(max_rg), invflag=False),
                          _arctanh_tanh_converter(init_rho_zero, invflag=False)]
        self._const_sig2_beta = const_sig2_beta
        self._const_sig2_zero = const_sig2_zero
        self._max_rg = abs(max_rg)
        self._lib = lib

    def _vec_to_params(self, vec):
        return BivariateParams(pi=[0,0,1], 
                               sig2_beta=self._const_sig2_beta,
                               sig2_zero=self._const_sig2_zero,
                               rho_beta=_arctanh_tanh_converter(vec[0], invflag=True) * self._max_rg,
                               rho_zero=_arctanh_tanh_converter(vec[1], invflag=True))

    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

# BGMG_cpp_fit_bivariate_fast (fits pi12 and rho_beta constrained on rg and all univariate params)
class BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO_boundedPI(object):
    def __init__(self, const_params1, const_params2, const_rg, const_rho_zero, lib):
        self._max_pi12 = min([const_params1._pi, const_params2._pi])
        self._min_pi12 = abs(const_rg) * np.sqrt(const_params1._pi * const_params2._pi)
        assert self._min_pi12 < self._max_pi12
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._const_rg = const_rg
        self._const_rho_zero = const_rho_zero
        self._lib = lib
        
    def _vec_to_params(self, pi12):
        # The following assertion doesn't work with Brent method, which may evaluate outside of the bracket range
        #assert (self._min_pi12 <= pi12) and (pi12 <= self._max_pi12)
        if pi12 < self._min_pi12: pi12 = self._min_pi12
        if pi12 > self._max_pi12: pi12 = self._max_pi12
        rho_beta = self._const_rg * np.sqrt(self._const_params1._pi * self._const_params2._pi) / pi12
        assert abs(rho_beta) <= 1
        return BivariateParams(pi=[self._const_params1._pi - pi12, self._const_params2._pi - pi12, pi12], 
                               sig2_beta=[self._const_params1._sig2_beta, self._const_params2._sig2_beta],
                               sig2_zero=[self._const_params1._sig2_zero, self._const_params2._sig2_zero],
                               rho_beta=rho_beta,
                               rho_zero=self._const_rho_zero)

    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib)
    
    # optimizer can be, for example
    # scalar_optimizer = lambda func, xLeft, xRight: scipy.optimize.minimize_scalar(func,  method='Brent', bracket=[xLeft, xRight])
    def fit(self, scalar_optimizer):
        result = scalar_optimizer(self._calc_cost)
        return self._vec_to_params(result.x), result

# BGMG_cpp_fit_bivariate_fast (fits pi12 and rho_beta constrained on rg and all univariate params)
# The difference from the above function (....boundedPI) is that here we map params into [-inf,+inf],
# while ...boundedPI is parametrized directly with PI.
class BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(object):
    def __init__(self, const_params1, const_params2, const_rg, const_rho_zero, lib):
        self._max_pi12 = min([const_params1._pi, const_params2._pi])
        self._min_pi12 = abs(const_rg) * np.sqrt(const_params1._pi * const_params2._pi)
        assert self._min_pi12 <= self._max_pi12
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._const_rg = const_rg
        self._const_rho_zero = const_rho_zero
        self._lib = lib

    def params_to_vec(self, params):
        assert(params._pi[2] >= self._min_pi12)
        assert(params._pi[2] <= self._max_pi12)
        return [_logit_logistic_converter((params._pi[2] - self._min_pi12) / (self._max_pi12 - self._min_pi12), invflag=False)]

    def vec_to_params(self, vec):
        if not hasattr(vec, "__len__"): vec = [vec]
        pi12 = self._min_pi12 + _logit_logistic_converter(vec[0], invflag=True) * (self._max_pi12 - self._min_pi12) 
        rho_beta = self._const_rg * np.sqrt(self._const_params1._pi * self._const_params2._pi) / pi12
        assert abs(rho_beta) <= 1
        return BivariateParams(pi=[self._const_params1._pi - pi12, self._const_params2._pi - pi12, pi12], 
                               sig2_beta=[self._const_params1._sig2_beta, self._const_params2._sig2_beta],
                               sig2_zero=[self._const_params1._sig2_zero, self._const_params2._sig2_zero],
                               rho_beta=rho_beta,
                               rho_zero=self._const_rho_zero)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)

# Bivariate parametrization with "independent axis", i.e.
#   x1 =     atanh(pi12/pi12max) *     atanh(rho_beta) 
#   x2 = log(atanh(pi12/pi12max) / abs(atanh(rho_beta)))
#   x3 = atanh(rho_zero)
# The first parameter, x1, defines defines the overal rg, 
# at least when both pi12 and rho_beta are sufficiently small
# The second parameter, x2, defines how rg is docomposed into pi12 and rho_beta,
# again when both pi12 and rho_beta are small enough.
# The reason for atanh(pi12/pi12max)) and atanh(rho_beta)) is to make sure that the invorse transform always give valid parameters.
# The reason fro log in x2 is to transform the ratio from [0, +inf] into [-inf, +inf] domain
# BGMG_cpp_fit_bivariate      (fits pi12, rho_beta, rho_zero)
class BivariateParametrization_constUNIVARIATE(object):
    def __init__(self, const_params1, const_params2, init_pi12, init_rho_beta, init_rho_zero, lib):
        max_pi12 = min(const_params1._pi, const_params2._pi)
        assert((init_pi12 >= 0) and (init_pi12 <= max_pi12))
        assert(abs(init_rho_beta) <= 1.0)
        assert(abs(init_rho_zero) <= 1.0)
        
        atanh_rho_beta =_arctanh_tanh_converter(init_rho_beta, invflag=False)
        atanh_pi12_frac = _arctanh_tanh_converter(init_pi12 / max_pi12, invflag=False)
        self._init_vec = [atanh_pi12_frac * atanh_rho_beta,
                          _log_exp_converter(atanh_pi12_frac / abs(atanh_rho_beta), invflag=False),
                          _arctanh_tanh_converter(init_rho_zero, invflag=False)]
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def _vec_to_params(self, vec, params1=None, params2=None):
        _params1 = params1 if params1 is not None else self._const_params1
        _params2 = params2 if params2 is not None else self._const_params2
        max_pi12 = min(_params1._pi, _params2._pi)
        atanh_pi12_frac = np.sqrt(np.abs(vec[0] * _log_exp_converter(vec[1], invflag=True)))
        atanh_rho_beta = vec[0] / atanh_pi12_frac
        pi12 = max_pi12 * _arctanh_tanh_converter(atanh_pi12_frac, invflag=True)
        rho_beta = _arctanh_tanh_converter(atanh_rho_beta, invflag=True)
        rho_zero = _arctanh_tanh_converter(vec[2], invflag=True)
        return BivariateParams(pi=[_params1._pi - pi12, _params2._pi - pi12, pi12], 
                               sig2_beta=[_params1._sig2_beta, _params2._sig2_beta],
                               sig2_zero=[_params1._sig2_zero, _params2._sig2_zero],
                               rho_beta=rho_beta,
                               rho_zero=rho_zero)

    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

class BivariateParametrization_constUNIVARIATE_natural_axis(object):
    def __init__(self, const_params1, const_params2, lib):
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def params_to_vec(self, params):
        max_pi12 = min(self._const_params1._pi, self._const_params2._pi)
        return [_arctanh_tanh_converter(params._rho_beta, invflag=False),
                _arctanh_tanh_converter(params._rho_zero, invflag=False),
                _logit_logistic_converter(params._pi[2] / max_pi12, invflag=False)]

    def vec_to_params(self, vec):
        max_pi12 = min(self._const_params1._pi, self._const_params2._pi)
        rho_beta = _arctanh_tanh_converter(vec[0], invflag=True)
        rho_zero = _arctanh_tanh_converter(vec[1], invflag=True)
        pi12 = max_pi12 * _logit_logistic_converter(vec[2], invflag=True)
        return BivariateParams(pi=[self._const_params1._pi - pi12, self._const_params2._pi - pi12, pi12], 
                               sig2_beta=[self._const_params1._sig2_beta, self._const_params2._sig2_beta],
                               sig2_zero=[self._const_params1._sig2_zero, self._const_params2._sig2_zero],
                               rho_beta=rho_beta, rho_zero=rho_zero)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)
    
# BGMG_cpp_fit_bivariate_fast_constrained (fits rho_zero to adapt to another reference)
class BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI(object):
    def __init__(self, const_params1, const_params2, const_pi12, const_rho_beta, init_rho_zero, lib):
        assert abs(init_rho_zero) <= 1
        assert (const_pi12 >= 0) and (const_pi12 <= min(const_params1._pi, const_params2._pi))
        self._init_vec = [_arctanh_tanh_converter(init_rho_zero, invflag=False)]
        self._const_pi12 = const_pi12
        self._const_rho_beta = const_rho_beta
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def _vec_to_params(self, vec):
        return BivariateParams(pi=[self._const_params1._pi - self._const_pi12,
                                   self._const_params2._pi - self._const_pi12,
                                   self._const_pi12], 
                               sig2_beta=[self._const_params1._sig2_beta, self._const_params2._sig2_beta],
                               sig2_zero=[self._const_params1._sig2_zero, self._const_params2._sig2_zero],
                               rho_beta=self._const_rho_beta,
                               rho_zero=_arctanh_tanh_converter(vec[0], invflag=True))

    def _calc_cost(self, vec):
        return self._vec_to_params(vec).cost(self._lib)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

def _hessian_robust(hessian, hessdiag):
    # for noisy functions hessian might be badly estimated
    # if we detect a problem with hessian fall back to hess diagonal
    hessdiag[hessdiag < 0] = 1e15
    hessdiag = np.diag(hessdiag)
    if not np.isfinite(hessian).all(): return hessdiag
    try:
        hessinv = np.linalg.inv(hessian)
    except np.linalg.LinAlgError as err:
        return hessdiag
    if not np.isfinite(hessinv).all(): return hessdiag
    if np.less_equal(np.linalg.eigvals(hessinv), 0).any(): return hessdiag
    return hessian

def _calculate_univariate_uncertainty(parametrization, alpha, totalhet, num_snps, num_samples):
    NCKoef = 0.319 # this koef gives proportion of causal variants that explain 90% of heritability. 
                   # it is specific to BGMG with single gaussian, with MAF specific model
    funcs = [('pi', lambda x: x._pi),
             ('nc', lambda x: x._pi * num_snps),
             ('nc@p9', lambda x: x._pi * num_snps * NCKoef),
             ('sig2_beta', lambda x: x._sig2_beta),
             ('sig2_zero', lambda x: x._sig2_zero),
             ('h2', lambda x: x._sig2_beta * x._pi * totalhet)]
    stats = [('mean', lambda x: np.mean(x)),
             ('median', lambda x: np.median(x)),
             ('std', lambda x: np.std(x)),
             ('lower', lambda x: np.percentile(x, 100.0 * (  alpha/2))),
             ('upper', lambda x: np.percentile(x, 100.0 * (1-alpha/2)))]
    hessian = _hessian_robust(nd.Hessian(parametrization._calc_cost)(parametrization._init_vec), 
                              nd.Hessdiag(parametrization._calc_cost)(parametrization._init_vec))
    x_sample = np.random.multivariate_normal(parametrization._init_vec, np.linalg.inv(hessian), num_samples)
    sample = [parametrization._vec_to_params(x) for x in x_sample]
    result = {}
    for func_name, func in funcs:
        result[func_name] = {'point_estimate': func(parametrization._vec_to_params(parametrization._init_vec))}
        param_vector = [func(s) for s in sample]
        for stat_name, stat in stats:
            result[func_name][stat_name] = stat(param_vector)
    return result, sample

def _calculate_bivariate_uncertainty(parametrization, ci_samples, alpha, totalhet, num_snps, num_samples):
    NCKoef = 0.319 # this koef gives proportion of causal variants that explain 90% of heritability. 
                   # it is specific to BGMG with single gaussian, with MAF specific model
    
    funcs = [('sig2_zero_T1', lambda x: x._sig2_zero[0]),
             ('sig2_zero_T2', lambda x: x._sig2_zero[1]),
             ('sig2_beta_T1', lambda x: x._sig2_beta[0]),
             ('sig2_beta_T2', lambda x: x._sig2_beta[1]),
             ('h2_T1', lambda x: x._sig2_beta[0] * (x._pi[0] + x._pi[2]) * totalhet),
             ('h2_T2', lambda x: x._sig2_beta[1] * (x._pi[1] + x._pi[2]) * totalhet),
             ('rho_zero', lambda x: x._rho_zero),
             ('rho_beta', lambda x: x._rho_beta),
             ('rg', lambda x: x._rho_beta * x._pi[2] / np.sqrt((x._pi[0] + x._pi[2]) * (x._pi[1] + x._pi[2]))),
             ('pi1', lambda x: x._pi[0]),
             ('pi2', lambda x: x._pi[1]),
             ('pi12', lambda x: x._pi[2]),
             ('pi1u', lambda x: x._pi[0] + x._pi[2]),
             ('pi2u', lambda x: x._pi[1] + x._pi[2]),
             ('nc1', lambda x: num_snps * x._pi[0]),
             ('nc2', lambda x: num_snps * x._pi[1]),
             ('nc12', lambda x: num_snps * x._pi[2]),
             ('nc1u', lambda x: num_snps * (x._pi[0] + x._pi[2])),
             ('nc2u', lambda x: num_snps * (x._pi[1] + x._pi[2])),
             ('nc1@p9', lambda x: NCKoef * num_snps * x._pi[0]),
             ('nc2@p9', lambda x: NCKoef * num_snps * x._pi[1]),
             ('nc12@p9', lambda x: NCKoef * num_snps * x._pi[2]),
             ('nc1u@p9', lambda x: NCKoef * num_snps * (x._pi[0] + x._pi[2])),
             ('nc2u@p9', lambda x: NCKoef * num_snps * (x._pi[1] + x._pi[2])),
             ('totalpi', lambda x: np.sum(x._pi)),
             ('totalnc', lambda x: num_snps * np.sum(x._pi)),
             ('totalnc@p9', lambda x: NCKoef * num_snps * np.sum(x._pi)),             
             ('pi1_over_totalpi', lambda x: x._pi[0] / np.sum(x._pi)),
             ('pi2_over_totalpi', lambda x: x._pi[1] / np.sum(x._pi)),
             ('pi12_over_totalpi', lambda x: x._pi[2] / np.sum(x._pi)),
             ('pi1_over_pi1u', lambda x: x._pi[0] / (x._pi[0] + x._pi[2])),
             ('pi2_over_pi2u', lambda x: x._pi[1] / (x._pi[1] + x._pi[2])),
             ('pi12_over_pi1u', lambda x: x._pi[2] / (x._pi[0] + x._pi[2])),
             ('pi12_over_pi2u', lambda x: x._pi[2] / (x._pi[1] + x._pi[2])),
             ('pi1u_over_pi2u', lambda x: (x._pi[0] + x._pi[2]) / (x._pi[1] + x._pi[2])),
             ('pi2u_over_pi1u', lambda x: (x._pi[1]  + x._pi[2]) / (x._pi[0] + x._pi[2]))]
              
    stats = [('mean', lambda x: np.mean(x)),
             ('median', lambda x: np.median(x)),
             ('std', lambda x: np.std(x)),
             ('lower', lambda x: np.percentile(x, 100.0 * (  alpha/2))),
             ('upper', lambda x: np.percentile(x, 100.0 * (1-alpha/2)))]
    
    hessian = _hessian_robust(nd.Hessian(parametrization._calc_cost)(parametrization._init_vec), 
                              nd.Hessdiag(parametrization._calc_cost)(parametrization._init_vec))
    x_sample = np.random.multivariate_normal(parametrization._init_vec, np.linalg.inv(hessian), num_samples)
    sample = [parametrization._vec_to_params(x, params1=ci_s1, params2=ci_s2) for ci_s1, ci_s2, x in zip(ci_samples[0], ci_samples[1], x_sample)]
    result = {}
    for func_name, func in funcs:
        result[func_name] = {'point_estimate': func(parametrization._vec_to_params(parametrization._init_vec))}
        param_vector = [func(s) for s in sample]
        for stat_name, stat in stats:
            result[func_name][stat_name] = stat(param_vector)
    return result, sample


if False:
    '''
    # test that converters effectively map the range of parameters into unbounded domain that will work with fminsearch
    def isclose(a, b, rel_tol=1e-09, abs_tol=2*epsval):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
    tests = [(precimed.mixer.utils._log_exp_converter,        0.0, maxval, [0.00001, 0.001, 1.2, 1000, 10000]), 
            (precimed.mixer.utils._logit_logistic_converter, 0.0, 1.0,   [0.00001, 0.001, 0.1, 0.9, 0.999, 0.99999]),
            (precimed.mixer.utils._arctanh_tanh_converter,    -1.0, 1.0,  [-0.9999, -0.991, -0.9, 0.9, 0.999, 0.99999])]
    for func, limLow, limHigh, values in tests:
        for index, val in enumerate([limLow] + values + [limHigh]):
            #print('{}: {} -> {} -> {}'.format(func, val, func(val, False), func(func(val, False), True)))
            assert np.isfinite(func(val, False)), 'Error in {}({})'.format(func, val)
            assert isclose(val, func(func(val, False), True)), 'Error in {}({})'.format(func, val)
            isLimit = (index==0) or (index==(len(values)+ 1))
            assert abs(func(val, False)) < (1000 if isLimit else 20), '{}({}) is too large ({})'.format(func, val, func(val, False))
            assert abs(func(val, False)) > (0.001 if isLimit else 0.1), '{}({}) is too small ({})'.format(func, val, func(val, False))
        assert isclose(func(-10000, True), limLow), '{} vs {}'.format(func(-10000, True), limLow)
        assert isclose(func( 10000, True), limHigh), '{} vs {}'.format(func(10000, True), limHigh)


    # validate that _vec_to_params and _params_to_vec complement each other in UnivariateParametrization 
    pp=UnivariateParametrization(UnivariateParams(0.001, 1e-4, 1.23), None, 1)
    UnivariateParametrization(pp._vec_to_params([2, -50, 4.5]), None, 1)._init_vec
    # >>> [2.0, -50.00000075435396, 4.499999245646041]
    '''
