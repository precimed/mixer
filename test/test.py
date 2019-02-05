import cmmutils
import numpy as np

het = np.array([0.1, 0.25, 0.3],dtype='f8')
annot = np.array([1,0,0], dtype='u1')
is2 = np.array([0,2,5,6], dtype='u8')

r2 = np.array([0.4, 1, 0.2, 0.8, 1, 1], dtype='f4')
i21 = np.array([1, 0, 0, 2, 1, 2], dtype='u4')

ind_2use = np.array([0,2], dtype='u4')
ssize = np.array([4,7.1], dtype='f8')
is2_2use = np.array([0,2,3], dtype='u8')

s2_2use = np.empty(3, dtype='f4')
annot_s2_2use = np.empty(3, dtype='u1')

cmmutils.fill_s2(het, annot, is2, r2, i21, ind_2use, ssize, is2_2use, s2_2use, annot_s2_2use)

s2_2use_exact = np.array([0.4, 0.4, 2.13], dtype='f4')
annot_s2_2use_exact = np.array([0, 1, 0], dtype='u1')

assert (s2_2use_exact == s2_2use).all()
assert (annot_s2_2use_exact == annot_s2_2use).all()

print("Test pass")

# print(f"s2_2use: {s2_2use}")
# print("")
# print(f"annot_s2_2use: {annot_s2_2use}")
