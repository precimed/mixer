 SEMT - Compile-time symbolic differentiation via C++ templates
============================================================================

The SEMT library provides an easy way to define arbitrary functions
and obtain their derivatives via C++ expression templates.

This is an extension of Zvi Guttermans original work
accompanying his master thesis
["Symbolic Pre-computation for Numerical Applications"](http://www.cs.technion.ac.il/~ssdl/thesis/finished/2004/ZviGutterman/).

Additional features include:
* partial differentiation
* easy-to-use objects for multi-valued functions
* basic expression substitution that will simplify generated expressions
* advanced facilities to create expressions

Because of this origin, all code is free, see the LICENSE
file for details.

Simple example
----------------

```C++
    DVAR(x, 0);
    auto f = pow(x, INT(3)) + sin(x);
    cout << f << endl;
    cout << deriv_t(f, x) << endl;
```
will print

    ((x0)^(3) + sin(x0))
    ((3 * (x0)^(2)) + cos(x0))

to your console.

Usage
---------
Run `make doc` to create the documentation with doxygen.
You can view the documentation [here](http://www.numasoft.de/semt-doc/) or download it [here](https://github.com/downloads/st-gille/semt/semt-0-3-doc.tar.gz).
There you will find an introductory example and the comprehensive API documentation.

There are several example programs:

* `semt_examples` - These are used in the documentation.
* `semt_speed` - A test run on a mapping R^10 -> R^10.
* `semt_newton` - Implementation of Newton's method.
* `semt_check` - Some crude tests.

You can run
    `make [release | gprof | gcov] <binary>`
to compile a corresponding binary.

Using *release* mode enables compiler optimizations and places the binary in the `Release` folder.

Using *gprof* creates a call graph and with timings in `gprof/<binary>_<arguments>`, where you can supply arguments via the `${ARGS}` make variable.

Using *gcov* outputs coverage results into `gcov`.

Contact me at *gille@numasoft.de* for any questions.
