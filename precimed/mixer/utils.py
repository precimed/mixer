'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''
# Utility classes for univariate and bivariate fit
# Contains
# _log_exp_converter, _logit_logistic_converter, _arctanh_tanh_converter - converters to map bounded parameters into -inf, +inf range
# UnivariateParams, BivariateParams - represent parameters, with basic functionality like "calculate cost"
# Several univariate and bivariate parametrizations - suitable for fitting and (some specific parametrization) for uncertainty calculation
# _calculate_univariate_uncertainty, _calculate_bivariate_uncertainty - estimates confidence intervals for parameters and their aggregates (like h2 or rg) 
 
import numpy as np
import numpy.matlib
import numdifftools as nd
import scipy.optimize
import scipy.stats
from scipy.interpolate import interp1d
import json

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

    def find_pi_mat(self, num_snp):
        return self._pi * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def find_sig2_mat(self, num_snp):
        return self._sig2_beta * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def cost(self, lib, trait):
        value = lib.calc_unified_univariate_cost(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp), 
                                                 sig2_zeroA=self._sig2_zero, sig2_zeroC=1, sig2_zeroL=0)
        return value if np.isfinite(value) else 1e100

    def aux(self, lib, trait):
        return lib.calc_unified_univariate_aux(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp), 
                                                 sig2_zeroA=self._sig2_zero, sig2_zeroC=1, sig2_zeroL=0)

    def pdf(self, lib, trait, zgrid):
        return lib.calc_unified_univariate_pdf(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp),
                                               sig2_zeroA=self._sig2_zero, sig2_zeroC=1, sig2_zeroL=0, zgrid=zgrid)

    def power(self, lib, trait, ngrid, zthresh=5.45):
        return lib.calc_unified_univariate_power(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp),
                                                 sig2_zeroA=self._sig2_zero, sig2_zeroC=1, sig2_zeroL=0, zthresh=zthresh, ngrid=ngrid)

# params for MAF-, LD-, and annotation-dependent architectures
# this also supports mixture of small and large effects (pass vector pi and sig2_beta)
# NB! Trick #1 np.cumsum(sig2_beta) is what gives variance per SNP
# NB! Trick #2 pi[0]==1 indicates that this component is responsible for "everything else"
class AnnotUnivariateParams(object):
    def __init__(self, pi=[None], sig2_beta=[None], sig2_annot=None, s=None, l=None, sig2_zeroA=None, sig2_zeroL=None, annomat=None, annonames=None, mafvec=None, tldvec=None):
        if annomat is not None: assert(annomat.ndim == 2) # 1D arrays don't work in np.dot as we want matrix multiplication
        self._pi = [pi] if ((pi is None) or isinstance(pi, (int, float))) else pi
        self._sig2_beta = [sig2_beta] if ((sig2_beta is None) or isinstance(sig2_beta, (int, float))) else sig2_beta
        self._sig2_annot = sig2_annot
        self._s = s
        self._l = l
        self._sig2_zeroA = sig2_zeroA
        self._sig2_zeroL = sig2_zeroL

        self._annomat = annomat
        self._annonames = annonames
        self._mafvec = mafvec
        self._tldvec = tldvec

    def as_string(self, attrs_list=['_pi', '_sig2_beta', '_sig2_annot', '_s', '_l', '_sig2_zeroA', '_sig2_zeroL']):
        description = []
        for attr_name in attrs_list:
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'AnnotUnivariateParams({})'.format(', '.join(description))

    def __str__(self):
        return self.as_string()

    __repr__ = __str__

    def as_dict(self):
        return {'pi': self._pi, 'sig2_beta': self._sig2_beta,
                'sig2_zeroA': self._sig2_zeroA, 'sig2_zeroL': self._sig2_zeroL,
                's': self._s, 'l': self._l,
                'sig2_annot': self._sig2_annot, 'annonames': self._annonames}

    def fit_sig2_annot(self, lib, trait):
        nnls_mat = self.nnls_mat(lib, trait=trait)

        defvec = np.isfinite(lib.get_zvec(trait)) & (lib.weights > 0)
        nnls_mat = nnls_mat[defvec, :]
        nnls_mat1 = np.concatenate([np.ones(shape=(nnls_mat.shape[0], 1), dtype=np.float32), nnls_mat], 1)
        z2 = np.square(lib.get_zvec(trait)[defvec])
        w = lib.weights[defvec]

        z2w = np.multiply(z2, np.sqrt(w))
        nnls_mat1w = np.multiply(nnls_mat1, np.matlib.repmat(np.sqrt(w.reshape(-1, 1)), 1, nnls_mat1.shape[1]))
        # note that in the above formulas the intercept is also weighted - that's correct 
        # weighted least squares is equivalent to pre-conditioning by sqrt(w), as one can validate as follows:
        # mod_wls = sm.WLS(z2, nnls_mat1, weights=w); res_wls = mod_wls.fit(); res_wls.params
        # x, res, rank, s = np.linalg.lstsq(nnls_mat1w, z2w, rcond=None)

        # https://en.wikipedia.org/wiki/Non-negative_least_squares
        #sig2_annot, cost=scipy.optimize.nnls(nnls_mat1, z2)
        sig2_annot, _=scipy.optimize.nnls(nnls_mat1w, z2w)
        self._sig2_zeroA=sig2_annot[0]  # save intercept
        self._sig2_annot=sig2_annot[1:]

    def drop_zero_annot(self):
        self._annomat=self._annomat[:, self._sig2_annot>0]
        self._annonames=[a for a, s in zip(self._annonames, self._sig2_annot) if s>0]
        self._sig2_annot=self._sig2_annot[self._sig2_annot > 0]
    
    def find_pi_mat(self, num_snp):
        pi = ([np.max([0, 1-np.sum(self._pi[1:])])] + self._pi[1:]) if (self._pi[0]==1) else self._pi  # Trick #2
        if (len(self._pi) > 1): raise(NotImplementedError("I suspect the matrix is wrongly oriented here"))
        return np.matmul(np.ones(shape=(num_snp, 1), dtype=np.float32), np.array(pi, dtype=np.float32).reshape(1, -1))

    def find_sig2_mat(self):
        sig2_beta = np.cumsum(self._sig2_beta) # Trick #1 (all subsequent components have larger variance)
        return np.matmul(
                np.multiply(np.dot(self._annomat, np.array(self._sig2_annot).astype(np.float32)),
                            np.multiply(np.power(np.float32(2.0) * self._mafvec * (1-self._mafvec), np.float32(self._s)),
                                        np.power(self._tldvec, np.float32(self._l)))).reshape(-1, 1),
                np.array(sig2_beta, dtype=np.float32).reshape(1, -1))

    def find_annot_h2(self, annomat):
        sig2_mat = self.find_sig2_mat()
        hetvec = np.float32(2.0) * self._mafvec * (1-self._mafvec)
        h2_vec = np.multiply(hetvec, np.matmul(sig2_mat, self.find_pi_mat(num_snp=1).reshape([-1, 1])).flatten())
        return np.matmul(h2_vec.reshape((1, len(h2_vec))), annomat)

    def find_annot_enrich(self, annomat):
        h2_annot = self.find_annot_h2(annomat)
        h2_total = h2_annot[0][0]
        snps_annot = np.sum(annomat, 0)
        snps_total = snps_annot[0]
        return np.divide(np.divide(h2_annot, h2_total), np.divide(snps_annot, snps_total))

    def nnls_mat(self, lib, trait):
        # NB! This function assumes an infinitesimal model
        assert(len(self._sig2_beta) == 1)

        num_annot = self._annomat.shape[1]
        
        pi_vec = np.ones(shape=(lib.num_snp, 1), dtype=np.float32)
        sig2_vec = np.multiply(np.power(np.float32(2.0) * self._mafvec * (1-self._mafvec), np.float32(self._s)),
                               np.power(self._tldvec, np.float32(self._l))) * self._sig2_beta[0]
        sig2_zeroL = 0 # force sig2_zeroL to zero (this may seem a bit tricky or contriversial, but later we add self._sig2_zeroL via infinitesimal model => it shouldn't contribute to annotation-specific sigma2_beta)
        sig2_zeroA = 0 # force sig2_zeroA to zero (this need not to contirbute to each TLD score across all annotation categories;
                       # the intercept term is added (and computed) later by fit_sig2_annot()
        sig2_zeroC = 1
        retval=np.zeros(shape=(lib.num_tag, num_annot), dtype=np.float32)
        lib.set_option('aux_option', 1)  # AuxOption_Ezvec2
        for annot_index in range(0, num_annot):
            retval[:, annot_index] = lib.calc_unified_univariate_aux(1, pi_vec, np.multiply(self._annomat[:, annot_index], sig2_vec), sig2_zeroA, sig2_zeroC, sig2_zeroL)
        lib.set_option('aux_option', 0)  # AuxOption_None
        return retval

    def cost(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        value = lib.calc_unified_univariate_cost(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        return value if np.isfinite(value) else 1e100

    def pdf(self, lib, trait, zgrid):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        return lib.calc_unified_univariate_pdf(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL, zgrid=zgrid)

    def tag_pdf(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', 2)  # AuxOption_TagPdf
        retval = lib.calc_unified_univariate_aux(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        lib.set_option('aux_option', 0)  # AuxOption_None
        return retval

    def tag_pdf_err(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', 3)  # AuxOption_TagPdfErr
        retval = lib.calc_unified_univariate_aux(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        lib.set_option('aux_option', 0)  # AuxOption_None
        return retval

    def power(self, lib, trait, ngrid, zthresh=5.45):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        return lib.calc_unified_univariate_power(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL, zthresh=zthresh, ngrid=ngrid)

class AnnotBivariateParams(object):
    def __init__(self):
        self._params1 = None # univariate architecture for each trait
        self._params2 = None # both params must be either single-component causal mixture, or an infinitesimal model
                             # both params are automatically alighed to share the list of annotation categories, and annotatinos not present in _params1 or _params2 receive sig2_annot=0.
        self._pi12 = None    # polygenic overlap parameter (must not exceed _params1._pi and _params2._pi)
        self._rho_annot = None # vector of effect size correlations for each annotation category
        self._rho_zeroA = None # correlation of sig2_zeroA
        self._rho_zeroL = None # correlation of sig2_zeroL

class AnnotUnivariateParametrization(object):
    def __init__(self, lib, trait, constraint):
        self._lib = lib
        self._trait = trait
        self._constraint = constraint # of type AnnotUnivariateParams, None indicate files that must be searched

    def params_to_vec(self, params):
        vec = []
        for constraint_pi, params_pi in zip(self._constraint._pi, params._pi):
            if constraint_pi is None: vec.append(_logit_logistic_converter(params_pi, invflag=False))
        for constraint_s2b, params_s2b in zip(self._constraint._sig2_beta, params._sig2_beta):
            if constraint_s2b is None: vec.append(_log_exp_converter(params_s2b, invflag=False))
        if self._constraint._s is None: vec.append(params._s)
        if self._constraint._l is None: vec.append(params._l)
        if self._constraint._sig2_zeroA is None: vec.append(_log_exp_converter(params._sig2_zeroA, invflag=False))
        if self._constraint._sig2_zeroL is None: vec.append(_log_exp_converter(params._sig2_zeroL, invflag=False))
        return vec

    def vec_to_params(self, vec):
        vec = list(vec)
        return AnnotUnivariateParams(
            pi=[(constraint_pi if (constraint_pi is not None) else _logit_logistic_converter(vec.pop(0), invflag=True)) for constraint_pi in self._constraint._pi],
            sig2_beta=[(constraint_s2b if (constraint_s2b is not None) else _log_exp_converter(vec.pop(0), invflag=True)) for constraint_s2b in self._constraint._sig2_beta],
            sig2_annot=self._constraint._sig2_annot,
            s=self._constraint._s if (self._constraint._s is not None) else vec.pop(0),
            l=self._constraint._l if (self._constraint._l is not None) else vec.pop(0),
            sig2_zeroA=self._constraint._sig2_zeroA if (self._constraint._sig2_zeroA is not None) else _log_exp_converter(vec.pop(0), invflag=True),
            sig2_zeroL=self._constraint._sig2_zeroL if (self._constraint._sig2_zeroL is not None) else _log_exp_converter(vec.pop(0), invflag=True),
            annomat=self._constraint._annomat,
            annonames=self._constraint._annonames,
            mafvec=self._constraint._mafvec,
            tldvec=self._constraint._tldvec)

    def calc_cost(self, vec):
        params = self.vec_to_params(vec)
        self._lib.log_message(params.as_string(attrs_list=['_pi', '_sig2_beta', '_s', '_l', '_sig2_zeroA', '_sig2_zeroL']))
        return params.cost(self._lib, self._trait)

# a specific approximation that constrains pi=1 and sig2_zeroL = sig2_beta.
class AnnotUnivariateParametrizationInf(object):
    def __init__(self, lib, trait, constraint):
        self._lib = lib
        self._trait = trait
        self._constraint = constraint # of type AnnotUnivariateParams

    def params_to_vec(self, params):
        return [_log_exp_converter(params._sig2_beta[0], invflag=False),
                _log_exp_converter(params._sig2_zeroA, invflag=False)]

    def vec_to_params(self, vec):
        vec = list(vec)
        sig2_beta = _log_exp_converter(vec.pop(0), invflag=True)
        sig2_zeroA = _log_exp_converter(vec.pop(0), invflag=True)
        return AnnotUnivariateParams(
            pi=1, sig2_beta=[sig2_beta], sig2_annot=[1], s=0, l=0, sig2_zeroA=sig2_zeroA, sig2_zeroL=sig2_beta,
            annomat=self._constraint._annomat,
            annonames=self._constraint._annonames,
            mafvec=self._constraint._mafvec,
            tldvec=self._constraint._tldvec)

    def calc_cost(self, vec):
        params = self.vec_to_params(vec)
        self._lib.log_message(params.as_string(attrs_list=['_pi', '_sig2_beta', '_s', '_l', '_sig2_zeroA', '_sig2_zeroL']))
        return params.cost(self._lib, self._trait)

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

    def _params1(self):
        return UnivariateParams(pi=self._pi[0] + self._pi[2], sig2_beta=self._sig2_beta[0], sig2_zero=self._sig2_zero[0])

    def _params2(self):
        return UnivariateParams(pi=self._pi[1] + self._pi[2], sig2_beta=self._sig2_beta[1], sig2_zero=self._sig2_zero[1])

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

    def find_pi_mat(self, num_snp):
        return np.matmul(np.array(self._pi, dtype=np.float32).reshape(-1, 1), np.ones(shape=(1, num_snp), dtype=np.float32))

    def find_sig2_mat(self, num_snp):
        return np.matmul(np.array(self._sig2_beta, dtype=np.float32).reshape(-1, 1), np.ones(shape=(1, num_snp), dtype=np.float32))

    def find_rho_vec(self, num_snp):
        return self._rho_beta * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def cost(self, lib):
        num_snp = lib.num_snp
        value = lib.calc_unified_bivariate_cost(self.find_pi_mat(num_snp), self.find_sig2_mat(num_snp), self.find_rho_vec(num_snp),
                                                sig2_zeroA=self._sig2_zero, sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=self._rho_zero, rho_zeroL=0)
        return value if np.isfinite(value) else 1e100

    def aux(self, lib):
        num_snp = lib.num_snp
        return lib.calc_unified_bivariate_aux(self.find_pi_mat(num_snp), self.find_sig2_mat(num_snp), self.find_rho_vec(num_snp),
                                              sig2_zeroA=self._sig2_zero, sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=self._rho_zero, rho_zeroL=0)

    def pdf(self, lib, zgrid):
        num_snp = lib.num_snp
        [zgrid1, zgrid2] = np.meshgrid(zgrid, zgrid)
        zgrid1=zgrid1[zgrid>=0, :]; zgrid2=zgrid2[zgrid>=0, :]

        pdf = lib.calc_unified_bivariate_pdf(self.find_pi_mat(num_snp), self.find_sig2_mat(num_snp), self.find_rho_vec(num_snp),
                                             sig2_zeroA=self._sig2_zero, sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=self._rho_zero, rho_zeroL=0,
                                             zvec1=zgrid1.flatten(), zvec2=zgrid2.flatten())

        pdf = pdf.reshape(zgrid1.shape)
        pdf = np.concatenate((np.fliplr(np.flipud(pdf[1:, :])), pdf))
        return pdf

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

def _calculate_univariate_uncertainty_funcs(alpha, totalhet, num_snps):
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
    return funcs, stats

def _calculate_univariate_uncertainty(parametrization, alpha, totalhet, num_snps, num_samples):
    funcs, stats = _calculate_univariate_uncertainty_funcs(alpha, totalhet, num_snps)
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

def _calculate_bivariate_uncertainty_funcs(alpha, totalhet, num_snps):
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
             ('dice', lambda x: (2 * x._pi[2]) / (x._pi[0] + x._pi[1] + 2*x._pi[2])),
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

    return funcs, stats

def _calculate_bivariate_uncertainty(parametrization, ci_samples, alpha, totalhet, num_snps, num_samples):
    funcs, stats = _calculate_bivariate_uncertainty_funcs(alpha, totalhet, num_snps)
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

def calc_qq_data(zvec, weights, hv_logp):
    # Step 0. calculate weights for all data poitns
    # Step 1. convert zvec to -log10(pvec)
    # Step 2. sort pval from large (1.0) to small (0.0), and take empirical cdf of this distribution
    # Step 3. interpolate (data_x, data_y) curve, as we don't need 10M data points on QQ plots
    data_weights = weights / np.sum(weights)                            # step 0
    data_y = -np.log10(2*scipy.stats.norm.cdf(-np.abs(zvec)))           # step 1
    si = np.argsort(data_y); data_y = data_y[si]                        # step 2
    data_x=-np.log10(np.flip(np.cumsum(np.flip(data_weights[si]))))     # step 2
    data_idx = np.not_equal(data_y, np.concatenate((data_y[1:], [np.Inf])))
    data_logpvec = interp1d(data_y[data_idx], data_x[data_idx],         # step 3
                            bounds_error=False, fill_value=np.nan)(hv_logp)
    return data_logpvec

def calc_qq_model(zgrid, pdf, hv_z):
    model_cdf = np.cumsum(pdf) * (zgrid[1] - zgrid[0])
    model_cdf = 0.5 * (np.concatenate(([0.0], model_cdf[:-1])) + np.concatenate((model_cdf[:-1], [1.0])))
    model_logpvec = -np.log10(2*interp1d(-zgrid[zgrid<=0], model_cdf[zgrid<=0],
                                         bounds_error=False, fill_value=np.nan)(hv_z))
    return model_logpvec

def calc_qq_plot(libbgmg, params, trait_index, downsample, mask=None, title=''):
    # mask can subset SNPs that are going into QQ curve, for example LDxMAF bin.
    if mask is None:
        mask = np.ones((libbgmg.num_tag, ), dtype=bool)

    zvec = libbgmg.get_zvec(trait_index)[mask]

    # Regular grid (vertical axis of the QQ plots)
    hv_z = np.linspace(0, np.min([np.max(np.abs(zvec)), 38.0]), 1000)
    hv_logp = -np.log10(2*scipy.stats.norm.cdf(-hv_z))

    # Empirical (data) QQ plot
    data_logpvec = calc_qq_data(zvec, libbgmg.weights[mask], hv_logp)

    # Estimated (model) QQ plots
    model_weights = libbgmg.weights
    mask_indices = np.nonzero(mask)[0]
    model_mask = np.zeros((len(model_weights), ), dtype=bool)
    model_mask[mask_indices[range(0, len(mask_indices), downsample)]] = 1
    model_weights[~model_mask] = 0
    model_weights = model_weights/np.sum(model_weights)

    zgrid = np.arange(0, 38.0, 0.05, np.float32)

    original_weights = libbgmg.weights
    libbgmg.weights = model_weights
    pdf = params.pdf(libbgmg, trait_index, zgrid)
    libbgmg.weights = original_weights

    pdf = pdf / np.sum(model_weights)
    zgrid = np.concatenate((np.flip(-zgrid[1:]), zgrid))  # extend [0, 38] to [-38, 38]
    pdf = np.concatenate((np.flip(pdf[1:]), pdf))
    model_logpvec = calc_qq_model(zgrid, pdf, hv_z)
    #plt.plot(model_logpvec, hv_logp, '-')
    #plt.plot(data_logpvec, hv_logp, '-')
    #plt.plot(hv_logp, hv_logp, '--k')

    return {'hv_logp': hv_logp, 'data_logpvec': data_logpvec, 'model_logpvec': model_logpvec,
            'n_snps': int(np.sum(mask)), 'sum_data_weights': float(np.sum(libbgmg.weights[mask])), 'title' : title}

def calc_power_curve(libbgmg, params, trait_index, downsample, nvec=None):
    power_nvec = np.power(10, np.arange(3, 8, 0.1)) if (nvec is None) else nvec

    original_weights = libbgmg.weights
    model_weights = libbgmg.weights

    mask = np.zeros((len(model_weights), ), dtype=bool)
    mask[range(0, len(model_weights), downsample)] = 1
    model_weights[~mask] = 0
    model_weights = model_weights/np.sum(model_weights)

    libbgmg.weights = model_weights     # temporary set downsampled weights
    try:
        power_svec = params.power(libbgmg, trait_index, power_nvec)
        libbgmg.weights = original_weights  # restore original weights
    except:
        libbgmg.weights = original_weights  # restore original weights
        return {}

    return {'nvec': power_nvec, 'svec': power_svec}

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if callable(obj):
            return str(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.float32):
            return np.float64(obj)
        return json.JSONEncoder.default(self, obj)

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
