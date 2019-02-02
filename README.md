# Inverse Iterator

#### Description

InverseIterator is intended to perform inverse iteration algorithm to find the smallest eigenvalue of a given matrix, using AMGX library (https://github.com/NVIDIA/AMGX) to solve the system of linear equations. You can find more specific information in the documentation (folder bin/).

#### Build Inverse Iterator

Set up AMGX build path (variable AMGX_BUILD_PATH in the Makefile) and type:

```sh
$ make
```

#### Build the hamiltonian example

This example perform inverse iteration algorithm on Hamilton operator.
```sh
$ cd examples/
```
Set up build paths (variables AMGX_BUILD_PATH, II_PATH and CUDA_PATH in the Makefile). Set the architecture of your GPU (variable CFLAGS in the Makefile) and type:

```sh
$ make
```

#### Running the hamiltonian example

```sh
$ cd examples/
$ ./hamiltonian <size_of_the_hamiltonian> <value_substracted_from_diagonal> <accurancy_of_the_calculations>
```

#### License

MIT