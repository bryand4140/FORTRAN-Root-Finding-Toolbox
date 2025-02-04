![FORTRAN-roots](media/logo.png)

## Description
A modern FORTRAN library for finding the roots of single variable functions and nth-order polynomials.

## Compiling and Executing
A Makefile is provided to build the code and has been written to be compatible with Windows and Linux OS. The Makefile has been configured to use main or test as the default executable. To build the code, use
```
make main.exe
```
By default, the code is built with double precision (`real64`). This can be changed within the module `MOD_Select_Real_Kind.f90`. After the code is built, module and executable files will be placed in a build directory. Use
````
./main
````
to run the code with your main program. Additionally, a `fpm.toml` file has been provided for users of the FORTRAN Package Manager (FPM) [Fortran Package Manager](https://github.com/fortran-lang/fpm).

## Methods
The code currently contains the following codes:
- `fun_zero` - Finding the root of a single-valued function within a bounded interval
- `poly_roots`- Finds the real and complex roots of an nth-degree polynomial

Detailed documentation on how to use the codes is provided in the PDF "Roots User Guid.pdf." 



