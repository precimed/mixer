import os
import sys
import ctypes
import six
import logging
import numpy as np

def _p2n(arg_value):  # python2native
    if isinstance(arg_value, str):
        return ctypes.create_string_buffer(arg_value.encode('utf-8'))
    return arg_value

def _n2p(res_value):  # native2python
    if isinstance(res_value, bytes):
        if six.PY3:
            return res_value.decode('utf-8')
    return res_value

# Not implemented:
# - clear_loglike_cache
# - extract_univariate_loglike_trajectory
# - extract_bivariate_loglike_trajectory

class LibBgmg(object):
    def __init__(self, lib_name=None, context_id=0):
        self._context_id = context_id
        self.cdll, self._lib_name = self._load_cdll(lib_name)
        logging.info('__init__(lib_name={}, context_id={})'.format(self._lib_name, context_id))

        float32_pointer_type = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        float64_pointer_type = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        int32_pointer_type = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')

        # set function signatures ('restype' and 'argtype') for all functions that involve non-integer types
        # (pointers, floats, doubles, etc - either as input or as output)
        self.cdll.bgmg_get_last_error.restype = ctypes.c_char_p
        self.cdll.bgmg_status.restype = ctypes.c_char_p
        self.cdll.bgmg_set_tag_indices.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_retrieve_tag_indices.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_set_mafvec.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_mafvec.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_weights.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_weights.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_chrnumvec.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_retrieve_chrnumvec.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_set_zvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_zvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_nvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_nvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_causalbetavec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_causalbetavec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_snp_order.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type]
        self.cdll.bgmg_retrieve_snp_order.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type]
        self.cdll.bgmg_retrieve_k_pdf.argtypes = [ctypes.c_int, ctypes.c_int, float64_pointer_type]
        self.cdll.bgmg_set_option.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_double]
        self.cdll.bgmg_set_ld_r2_coo_from_file.argtypes = [ctypes.c_int, ctypes.c_char_p]
        self.cdll.bgmg_set_ld_r2_csr.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_set_weights_randprune.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float]
        self.cdll.bgmg_retrieve_ld_tag_r2_sum.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_ld_tag_r4_sum.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_num_ld_r2_snp.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r2_snp.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_num_ld_r2_chr.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r2_chr.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_num_ld_r2_snp_range.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r2_snp_range.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_univariate_cost.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        self.cdll.bgmg_calc_univariate_cost.restype = ctypes.c_double
        self.cdll.bgmg_calc_univariate_pdf.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_univariate_power.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_univariate_delta_posterior.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_bivariate_cost.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type, ctypes.c_int, float32_pointer_type, ctypes.c_float, ctypes.c_int, float32_pointer_type, ctypes.c_float]
        self.cdll.bgmg_calc_bivariate_cost.restype = ctypes.c_double
        self.cdll.bgmg_calc_bivariate_pdf.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type, ctypes.c_int, float32_pointer_type, ctypes.c_float, ctypes.c_int, float32_pointer_type, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type]

    def get_last_error(self):
        return _n2p(self.cdll.bgmg_get_last_error())

    def init_log(self, file):
        logging.info('init_log({})'.format(file)) 
        self.cdll.bgmg_init_log(_p2n(file))

    def log_message(self, message):
        logging.info('log_message({})'.format(message)) 
        self.cdll.bgmg_log_message(_p2n(message))

    def dispose(self):
        return self._check_error(self.cdll.bgmg_dispose(self._context_id))

    @property
    def status(self):
        return _n2p(self.cdll.bgmg_status(self._context_id))

    def init(self, bim_file, frq_file, chr_labels, trait1_file, trait2_file, exclude, extract):
        chr_labels_val = chr_labels if isinstance(chr_labels, str) else ' '.join([str(x) for x in chr_labels])
        return self._check_error(self.cdll.bgmg_init(
            self._context_id, _p2n(bim_file), _p2n(frq_file), _p2n(chr_labels_val), _p2n(trait1_file), _p2n(trait2_file), _p2n(exclude), _p2n(extract)))

    def set_option(self, option, value):
        if value is None: return None
        return self._check_error(self.cdll.bgmg_set_option(self._context_id, _p2n(option), value))

    def convert_plink_ld(self, plink_ld_gz, plink_ld_bin):
        return self._check_error(self.cdll.bgmg_convert_plink_ld(self._context_id, _p2n(plink_ld_gz), _p2n(plink_ld_bin)))

    def set_ld_r2_coo_from_file(self, filename):
        return self._check_error(self.cdll.bgmg_set_ld_r2_coo_from_file(self._context_id, _p2n(filename)))

    def set_ld_r2_csr(self, chr_label=-1):  # -1 means to finalize all chromosomes
        return self._check_error(self.cdll.bgmg_set_ld_r2_csr(self._context_id, chr_label))

    def set_weights_randprune(self, n, r2):
        return self._check_error(self.cdll.bgmg_set_weights_randprune(self._context_id, n, r2))

    @property
    def num_tag(self):
        return self._check_error(self.cdll.bgmg_get_num_tag(self._context_id))
   
    @property
    def num_snp(self):
        return self._check_error(self.cdll.bgmg_get_num_snp(self._context_id))

    @property
    def max_causals(self):
        return self._check_error(self.cdll.bgmg_get_max_causals(self._context_id))
   
    @property
    def k_max(self):
        return self._check_error(self.cdll.bgmg_get_k_max(self._context_id))

    @property
    def defvec(self):
        numpy_ndarray = np.zeros(shape=(self.num_tag,), dtype=np.int32)
        self._check_error(self.cdll.bgmg_retrieve_tag_indices(self._context_id, np.size(numpy_ndarray), numpy_ndarray))
        return numpy_ndarray

    @defvec.setter
    def defvec(self, val):
        indices = np.flatnonzero(np.array(val, dtype=np.bool)).astype(np.int32)
        self._check_error(self.cdll.bgmg_set_tag_indices(self._context_id, np.size(val), np.size(indices), indices))

    @property
    def mafvec(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_mafvec, np.float32, self.num_snp, trait=None)

    @mafvec.setter
    def mafvec(self, val):
        self._set_vec_impl(self.cdll.bgmg_set_mafvec, np.float32, val, trait=None)

    @property
    def weights(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_weights, np.float32, self.num_tag, trait=None)

    @weights.setter
    def weights(self, val):
        self._set_vec_impl(self.cdll.bgmg_set_weights, np.float32, val, trait=None)

    @property
    def chrnumvec(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_chrnumvec, np.int32, self.num_snp, trait=None)

    @chrnumvec.setter
    def chrnumvec(self, val):  
        self._set_vec_impl(self.cdll.bgmg_set_chrnumvec, np.int32, val, trait=None)

    @property
    def k_pdf(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_k_pdf, np.float64, self.k_max, trait=None)

    def set_zvec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_zvec, np.float32, val, trait=trait)

    def get_zvec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_zvec, np.float32, self.num_tag, trait=trait)

    def set_nvec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_nvec, np.float32, val, trait=trait)

    def get_nvec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_nvec, np.float32, self.num_tag, trait=trait)

    def set_causalbetavec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_causalbetavec, np.float32, val, trait=trait)

    def get_causalbetavec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_causalbetavec, np.float32, self.num_snp, trait=trait)

    def get_snp_order(self, component_id):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_snp_order, np.int32, self.k_max * self.max_causals, trait=component_id).reshape((self.k_max, self.max_causals))

    def set_snp_order(self, component_id, val):
        self._set_vec_impl(self.cdll.bgmg_set_snp_order, np.int32, val, trait=component_id)

    @property
    def zvec1(self):
        return self.get_zvec(trait=1)

    @zvec1.setter
    def zvec1(self, val):
        self.set_zvec(val, trait=1)

    @property
    def zvec2(self):
        return self.get_zvec(trait=2)

    @zvec2.setter
    def zvec2(self, val):
        self.set_zvec(val, trait=2)

    @property
    def nvec1(self):
        return self.get_nvec(trait=1)

    @nvec1.setter
    def nvec1(self, val):
        self.set_nvec(val, trait=1)

    @property
    def nvec2(self):
        return self.get_nvec(trait=2)

    @nvec2.setter
    def nvec2(self, val):
        self.set_nvec(val, trait=2)

    @property
    def causalbetavec1(self):
        return self.get_causalbetavec(trait=1)

    @causalbetavec1.setter
    def causalbetavec1(self, val):
        self.set_causalbetavec(val, trait=1)

    @property
    def causalbetavec2(self):
        return self.get_causalbetavec(trait=2)

    @causalbetavec2.setter
    def causalbetavec2(self, val):
        self.set_causalbetavec(val, trait=2)

    @property
    def ld_tag_r2_sum(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_ld_tag_r2_sum, np.float32, self.num_tag, trait=None)

    @property
    def ld_tag_r4_sum(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_ld_tag_r4_sum, np.float32, self.num_tag, trait=None)

    # return (tag, r2) tuple, representing LD structure of a given reference snp
    # tag and r2 are vectors of the same length, tag gives an index of a tag snp, r2 gives corresponding LD r2 correlation
    def get_ld_r2_snp(self, snp_index):
        num_ld_r2 = self._check_error(self.cdll.bgmg_num_ld_r2_snp(self._context_id, snp_index))
        tag_array = np.zeros(shape=(num_ld_r2,), dtype=np.int32)
        r2_array = np.zeros(shape=(num_ld_r2,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r2_snp(self._context_id, snp_index, num_ld_r2, tag_array, r2_array))
        return (tag_array, r2_array)

    # return (snp, tag, r2) tuple, representing LD structure of a given chromosome
    # snp, tag and r2 are vectors of the same length
    # snp gives an index of a reference snp
    # tag gives an index of a tag snp
    # r2 gives corresponding LD r2 correlation
    def get_ld_r2_chr(self, chr_label):
        num_ld_r2 = self._check_error(self.cdll.bgmg_num_ld_r2_chr(self._context_id, chr_label))
        snp_array = np.zeros(shape=(num_ld_r2,), dtype=np.int32)
        tag_array = np.zeros(shape=(num_ld_r2,), dtype=np.int32)
        r2_array = np.zeros(shape=(num_ld_r2,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r2_chr(self._context_id, chr_label, num_ld_r2, snp_array, tag_array, r2_array))
        return (snp_array, tag_array, r2_array)

    # return (snp, tag, r2) tuple, representing LD structure of a given range of snps
    # [from_snp, to_snp) - left snp inclusive, right snp exclusive, values between 0 and num_snp.
    # snp, tag and r2 are vectors of the same length
    # snp gives an index of a reference snp
    # tag gives an index of a tag snp
    # r2 gives corresponding LD r2 correlation
    def get_ld_r2_snp_range(self, from_snp, to_snp):
        num_ld_r2 = self._check_error(self.cdll.bgmg_num_ld_r2_snp_range(self._context_id, from_snp, to_snp))
        snp_array = np.zeros(shape=(num_ld_r2,), dtype=np.int32)
        tag_array = np.zeros(shape=(num_ld_r2,), dtype=np.int32)
        r2_array = np.zeros(shape=(num_ld_r2,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r2_snp_range(self._context_id, from_snp, to_snp, num_ld_r2, snp_array, tag_array, r2_array))
        return (snp_array, tag_array, r2_array)

    def calc_univariate_cost(self, trait, pi_vec, sig2_zero, sig2_beta):
        cost = self.cdll.bgmg_calc_univariate_cost(self._context_id, trait, pi_vec, sig2_zero, sig2_beta)
        self._check_error()
        return cost

    def calc_univariate_pdf(self, trait, pi_vec, sig2_zero, sig2_beta, zgrid):
        zgrid_data = (zgrid if isinstance(zgrid, np.ndarray) else np.array(zgrid)).astype(np.float32)
        pdf = np.zeros(shape=(np.size(zgrid),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_univariate_pdf(self._context_id, trait, pi_vec, sig2_zero, sig2_beta, np.size(zgrid), zgrid_data, pdf))
        return pdf

    def calc_univariate_power(self, trait, pi_vec, sig2_zero, sig2_beta, zthresh, ngrid):
        ngrid_data = (ngrid if isinstance(ngrid, np.ndarray) else np.array(ngrid)).astype(np.float32)
        svec = np.zeros(shape=(np.size(ngrid),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_univariate_power(self._context_id, trait, pi_vec, sig2_zero, sig2_beta, zthresh, np.size(ngrid), ngrid_data, svec))
        return svec    

    def calc_univariate_delta_posterior(self, trait, pi_vec, sig2_zero, sig2_beta):
        raise(RuntimeError('Disable calc_univariate_delta_posterior - for some reason it crashes in native c++ plugin'))
        c0 = np.zeros(shape=(np.size(self.num_tag),), dtype=np.float32)
        c1 = np.zeros(shape=(np.size(self.num_tag),), dtype=np.float32)
        c2 = np.zeros(shape=(np.size(self.num_tag),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_univariate_delta_posterior(self._context_id, trait, pi_vec, sig2_zero, sig2_beta, self.num_tag, c0, c1, c2))
        return (c0, c1, c2)

    def calc_bivariate_cost(self, pi_vec, sig2_beta, rho_beta, sig2_zero, rho_zero):
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        sig2_zero = (sig2_zero if isinstance(sig2_zero, np.ndarray) else np.array(sig2_zero)).astype(np.float32)
        cost = self.cdll.bgmg_calc_bivariate_cost(self._context_id, np.size(pi_vec), pi_vec, np.size(sig2_beta), sig2_beta, rho_beta, np.size(sig2_zero), sig2_zero, rho_zero)
        self._check_error()
        return cost

    def calc_bivariate_pdf(self, pi_vec, sig2_beta, rho_beta, sig2_zero, rho_zero, zvec1, zvec2):
        #, int length, float* zvec1, float* zvec2, float* pdf
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        sig2_zero = (sig2_zero if isinstance(sig2_zero, np.ndarray) else np.array(sig2_zero)).astype(np.float32)
        zvec1 = (zvec1 if isinstance(zvec1, np.ndarray) else np.array(zvec1)).astype(np.float32)
        zvec2 = (zvec2 if isinstance(zvec2, np.ndarray) else np.array(zvec2)).astype(np.float32)
        if np.size(zvec1) != np.size(zvec2): raise(RuntimeError("len(zvec1) != len(zvec2)"))
        pdf = np.zeros(shape=(np.size(zvec1),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_bivariate_pdf(self._context_id, np.size(pi_vec), pi_vec, np.size(sig2_beta), sig2_beta, rho_beta, np.size(sig2_zero), sig2_zero, rho_zero, np.size(zvec1), zvec1, zvec2, pdf))
        return pdf

    def __str__(self):
        description = []
        for attr_name in '_lib_name', '_context_id', 'num_snp', 'num_tag':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'LibBgmg({})'.format(', '.join(description))
    __repr__ = __str__

    def _set_vec_impl(self, func, arg_type, arg_value, trait=None):
        data = (arg_value if isinstance(arg_value, np.ndarray) else np.array(arg_value)).astype(arg_type)
        args = [self._context_id] + ([trait] if trait != None else []) + [np.size(data), data]
        self._check_error(func(*args))

    def _get_vec_impl(self, func, arg_type, arg_size, trait=None):
        numpy_ndarray = np.zeros(shape=(arg_size,), dtype=arg_type)
        args = [self._context_id] + ([trait] if trait != None else []) + [np.size(numpy_ndarray), numpy_ndarray]
        self._check_error(func(*args))
        return numpy_ndarray

    def _check_error(self, error_code=None):
        logging.debug('status:{}'.format(self.status))
        
        # if there is an error code, check that it is not negative
        if (error_code is not None) and (error_code < 0):
            raise RuntimeError(self.get_last_error())
        
        # if there is no error code, check that last error message is empty
        if (error_code is None) and (self.get_last_error()):
            raise RuntimeError(self.get_last_error())

        return error_code
  
    def _load_cdll(self, lib_name):
        # choose default library name
        default_lib_name = 'libbgmg.so'
        if sys.platform.startswith('win'):
            default_lib_name = 'bgmg.dll'
        if sys.platform.startswith('darwin'):
            default_lib_name = 'libbgmg.dylib'

        lib_names = []
        
        if lib_name is not None:
            lib_names.append(lib_name)
            
        env_lib_name = os.environ.get('BGMG_SHARED_LIBRARY')
        if env_lib_name is not None:
            lib_names.append(env_lib_name)
        
        lib_names.append(default_lib_name)
        
        # We look into 4 places: lib_name, BGMG_SHARED_LIBRARY, packaged default_lib_name
        # and then default_lib_name
        cdll = None
        for ln in lib_names:
            try:
                cdll = ctypes.CDLL(ln)
                break
            except OSError as e:
                if ln == default_lib_name:
                    exc = e
                continue
        if cdll is None:
            exception_message = (
                '{exc}\n'
                'Failed to load BGMG shared library from `{lib_names}`. '
                'Try to add the location of `{default_lib_name}` file into your PATH '
                'system variable, or to set BGMG_SHARED_LIBRARY - the specific system variable '
                'which may point to `{default_lib_name}` file, including the full path.'
            ).format(**locals())
            raise OSError(exception_message)

        return (cdll, ln)
