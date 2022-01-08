Affine Particle in Cell in 2D
================
![Screenshot](http://www.cs.columbia.edu/cg/raymond/apic2d.jpg)
APIC2D is an educational project to illustrate the affine-particle-in-cell algorithm in 2D, for water simulation.

The papers implemented here include:

Jiang, Chenfanfu, et al. "The affine particle-in-cell method." ACM Transactions on Graphics (TOG) 34.4 (2015): 51.

Batty, Christopher, Florence Bertails, and Robert Bridson. "A fast variational framework for accurate solid-fluid coupling." ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM, 2007.

Ando, Ryoichi, Nils Thurey, and Reiji Tsuruno. "Preserving fluid sheets with adaptively sampled anisotropic particles." IEEE transactions on visualization and computer graphics 18.8 (2012): 1202-1214.

Brackbill, Jeremiah U., and Hans M. Ruppel. "FLIP: A method for adaptively zoned, particle-in-cell calculations of fluid flows in two dimensions." Journal of Computational physics 65.2 (1986): 314-343.

Fei, Yun (Raymond), et al. "Revisiting Integration in the Material Point Method: A Scheme for Easier Separation and Less Dissipation." ACM Transactions on Graphics (TOG) 40.4 (2021): 109.

It contains multiple integrators that you may switch and compare through changing the `integration_scheme` variable in the code. Its value can be one of the following:
```
IT_PIC: original particle-in-cell (PIC)
IT_FLIP: original fluid-implicit-particle (FLIP)
IT_RPIC: rotational particle-in-cell (RPIC)
IT_APIC: affine particle-in-cell (APIC)
IT_AFLIP: affine fluid-implicit-particle (AFLIP)
IT_ASFLIP: affine separable fluid-implicit-particle (ASFLIP)
```

It also supports using different orders for velocity evaluation, where one may change the `velocity_order` variable in the code. Its value can be one of the following:
```
VO_EULER: first order evaluation
VO_RA2: Ralston's second order evaluation
VO_RK3: Runge Kutta's 3rd-order method
VO_RK4: Runge Kutta's 4rd-order method
```

Dependencies
--------------------
APIC2D depends on the Eigen libraries (included), TBB for parallelization, as well as FreeGLUT for simple visualization. You may use the Homebrew (on Mac) or APT package handling utility (on Ubuntu Linux) to install these dependencies. For example, with Homebrew on Mac OS X, these external dependencies can be installed through
```
brew install tbb freeglut
```

For Windows, CMake should automatically grab the dependencies from its submodule.

Compilation
-----------------
To compile APIC2D, you'll need CMake or CMake-GUI (https://cmake.org).

CMake:
1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..* (or *cmake -G Xcode ..* to generate Xcode project files on Mac; for Windows please use *cmake -G <generator> ..* to use the specific generators, or simply type *cmake -G* to list all the available generators)
2. Optionally you can adjust the options with *ccmake .*
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

CMake-GUI:
1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose the generator (for Windows, select the installed version of the Microsoft Visual Studio).
2. CMake should automatically find all the dependencies. If not, check the *Advanced* box and locate those missing libraries manually. On Windows, please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables.
4. open the solution (for Visual Studio or Xcode) and compile the code.
