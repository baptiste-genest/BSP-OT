## Official source code of the "BSP-OT: Sparse transport plans between discrete measures in loglinear time" paper (SIGGRAPH Asia 2025). 
Baptiste Genest, Nicolas Bonneel, Vincent Nivoliers, David Coeurjolly.

Standalone python binding and package : [kenshi84/bspot-python](https://github.com/kenshi84/bspot-python)

Python Optimal Transport binding : https://pythonot.github.io/master/auto_examples/plot_bsp_ot.html

CUDA version ***(10x speedup!!)*** : [baptiste-genest/BSP-OT-CUDA](https://github.com/baptiste-genest/BSP-OT-CUDA)

[![CMake on multiple platforms](https://github.com/baptiste-genest/BSP-OT/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/baptiste-genest/BSP-OT/actions/workflows/cmake-multi-platform.yml) 

[![](https://www.replicabilitystamp.org/logo/Reproducibility-small.png)](http://www.replicabilitystamp.org#https-github-com-baptiste-genest-bsp-ot) 

![teaser](https://github.com/baptiste-genest/BSP-OT/blob/main/teaser.jpg)

# BSP-OT Compilation Guide

This guide explains how to compile the BSPOT project using CMake, focusing on the top-level `CMakeLists.txt`. It also lists all dependencies and how to obtain them.

## Requirements

- **CMake** >= 3.12
- **C++20** compatible compiler (e.g., GCC 10+, Clang 10+, MSVC 2019+)
- **git** (for fetching some dependencies)
- **Internet connection** (for fetching dependencies using CPM)
- **Imagick** (only for color transfer, for image resizing and format conversion (not fetched))

### Dependencies

The following libraries are required and are automatically handled by the CMake build system via CPM or included CMake scripts:

- [Eigen3](https://gitlab.com/libeigen/eigen) (version 3.4.0, downloaded automatically)
- [OpenMP](https://www.openmp.org/) (for parallelization, usually provided by your compiler)
- [Polyscope](https://github.com/nmwsharp/polyscope)
- [geometry-central](https://github.com/nmwsharp/geometry-central)
- [spdlog](https://github.com/gabime/spdlog)
- [Spectra](https://github.com/yixuan/spectra)

All dependencies except Eigen3 are also included via CMake include scripts (see the `cmake/` directory).

## Step-by-Step Compilation

1. **Go in code folder**
   ```bash
   cd BSP-OT
   ```

2. **Create a build directory**
   ```bash
   mkdir build
   cd build
   ```

3. **Configure the project with CMake**
   ```bash
   cmake ..
   ```
   - This will download and configure all dependencies using CPM and the scripts in `cmake/`.

4. **Build the project**
   ```bash
   cmake --build .
   ```
   - This will build all executables defined in the main `CMakeLists.txt`.

## Available Executables

After compilation, the following programs will be built (if their sources are present):

- `bijections`
- `manifold_sampling`
- `persistance_diagrams_matching`
- `barycenters`
- `color_transfer`
- `stippling`
- `scale_rigid_registration`

Each corresponds to a source file in the `apps/` directory.

## Examples

Note that all parameters can be described with the help command on any the executable:
```bash
./any_exe --help
```
to reproduce the figure 8, you can execute
```bash
./bijections --mu_file ../data/point_clouds/armadillo.pts --nb_trees 64 --viz
```

the --viz parameter allows to see the results with polyscope. You should be able to do this:
![example](https://github.com/baptiste-genest/BSP-OT/blob/main/armadillo_example.gif)

For stippling 
```bash
./stippling --size_mu 10000 --nu_file ../data/images/fruits.png --res_grid 250 --output rslt.pts --viz
```

And color transfer
```bash
./color_transfer --target_image ../data/images/mountain.png --colors ../data/images/painting.jpg --iter 16 --output rslt.png
```

## Static parameters

To optimize performances, the code has some static parameters. For bijective applications (bijections, barycenters, persistance_diagrams_matching, color_transfer) you can compile with floats to get a speed-up without changing the quality. The other applications must use doubles. this is set by the type *scalar* defined in common/types.h. Double by default. Each main file in apps is compiled with a static dimension, if you want to try 2D examples, please set "static_dim = 2".

## Header only

If you want to easily import BSP-OT into your project, feel free to use the header only file `BSP-OT_header_only.h`. 
Note that eigen is still a dependancy. It can then easily be used via something like:

```cpp
#include "BSP-OT_header_only.h"


int main() {
	using namespace BSPOT;
	Points<2> A,B;
	A = Points<2>::Random(2,1000);
	B = Points<2>::Random(2,1000);


	auto cost = [&] (int i,int j) {
	return (A.col(i) - B.col(j)).squaredNorm();
	};

	// source points, target points, number of trees to compute and merge, cost between points
	auto T = computeGaussianBSPOT(A,B,64,cost);

	std::cout << "matching cost: " << T.evalMatching(cost) << std::endl;
	for (auto i : range(A.cols())) {
		std::cout << A.col(i).transpose() << " matched with " << B.col(T[i]).transpose() << std::endl;
	}

	return 0;
}
```

# Python binding

Python binding (maintained separately at [kenshi84/bspot-python](https://github.com/kenshi84/bspot-python)) can be installed via `pip`:

```bash
pip install bspot
```

Example usage:

```python
>>> import bspot
>>> import numpy as np
>>> A = np.random.randn(3,10000)
>>> B = np.random.randn(3,10000)
>>> bspot.set_num_threads(8)        # By default, use all available threads
>>> bspot.compute_matching(A, B, gaussian=True)
array([1586, 7207,  330, ..., 3329, 4056, 3637],
      shape=(10000,), dtype=int32)
```
