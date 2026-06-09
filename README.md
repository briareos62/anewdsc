# anewdsc

A C library for operations on polynomials and finding real roots.

## Dependencies

This library depends on the following:
- [GMP](https://gmplib.org/) (GNU Multiple Precision Arithmetic Library)
- [MPFR](https://www.mpfr.org/) (GNU Multiple Precision Floating-Point Reliably)
- [MPFI](http://mpfi.gforge.inria.fr/) (Multiple Precision Floating-point Interval library)

## Building

This project uses CMake for building. By default, both a shared library (`.so` on Linux, `.dll` on Windows) and a static library are built.

### Linux / macOS

By default, CMake configures a **Release** build (with optimizations, `-O3`).

```bash
mkdir build
cd build
cmake ..
make
```

To configure a **Debug** build (with debug symbols, `-g`), set `CMAKE_BUILD_TYPE`:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

### Windows (Shared Library / DLL)

To build a DLL on Windows, you can use MinGW/MSYS2 (defaults to Release):

```bash
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

For a Debug build with MinGW, append `-DCMAKE_BUILD_TYPE=Debug` to the `cmake` command.

Alternatively, you can build using Microsoft Visual Studio (MSVC):

```cmd
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

*(To build a Debug version with MSVC, simply use `--config Debug` instead).*

## Testing

To run the tests:

```bash
cd build
make test
```

## License

This project is available under a dual-licensing model:

**Open Source & Evaluation:** For open-source projects, academic research, teaching, and personal use, this library is licensed under the GNU GPLv3. Please note the strict copyleft conditions of this license. See the [LICENSE](LICENSE) file for details.

**Commercial Use:** If you wish to use this library in a commercial, proprietary product without having to release your own source code under the GPL, you require a commercial license. Please contact us at math[at]leipelt-hamburg.de.

> **Note for commercial users:** While a commercial license allows you to use `anewdsc` in proprietary products, please be aware that this library depends on GMP, MPFR, and MPFI, which are licensed under the LGPL. To easily comply with the LGPL requirements for these dependencies, we strongly recommend dynamically linking `anewdsc` and its dependencies (using shared libraries / DLLs).
