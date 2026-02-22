# Tests

## Prerequisites
- CMake ≥ 3.10
- Intel Fortran Compiler (e.g., oneAPI)

> **Important**: The starting point of the following commands must be the repository root directory (the directory containing the top-level `CMakeLists.txt`).

## 1. Create and enter a build directory

```
mkdir build
cd build
```

## 2. Configure the project with CMake

```
cmake ..
```

## 3. Build 

```
cmake --build . --config Debug
```

## 4. Run tests

```
ctest -C Debug --output-on-failure --verbose
```
