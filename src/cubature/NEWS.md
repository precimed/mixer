# Release notes for cubature code

The following are the main changes in each subsequent tagged release
of the [cubature code by Steven G. Johnson](README.md).

## Version 1.0.3

* Transferred files and documentation to [Github](https://github.com/stevengj/cubature).

## Version 1.0.2

* Fix memory leak in `hcubature` on BSD (and MacOS) systems.  Thanks to
  Nicolas Tessore (@ntessore) for the bug report.

## Version 1.0.1

* `cubature.h` header now includes `<stdlib.h>`, to make sure `size_t` is defined

## Version 1.0

* Many API changes compared to pre-1.0 versions:
    - rename `adapt_integrate` -> `hcubature`
    - integrand now returns int to signal errors
    - `error_norm` argument for vector-valued integrands
    - `maxeval` and `npt` args are now `size_t`, not `int`

* New `pcubature` routines for p-adaptive (Clenshaw-Curtis) integration.

* Split `test.c` and some `#include` stuff from `hcubature.c` (to share
  with `pcubature.c`).
