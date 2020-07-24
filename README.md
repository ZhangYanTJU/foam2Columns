# foam2Columns

[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-7-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-7)
[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-8-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-8)

## Usage

It supports any number of fields and transforms them to this format (in each time directory):
```
L1: x y z var1 var2 ...
L2: x y z var1 var2 ...
.
.
.
Ln: x y z var1 var2 ...
```
where n is the cell number, and var1(p), var2(T) are the input variables:

```c++
foam2Columns -fields "(p T)"
```
