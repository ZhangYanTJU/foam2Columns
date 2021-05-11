# foam2Columns

[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-7-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-7)
[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-8-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-8)

## Usage

It supports any number of fields and transforms them to this format:
```
L1: x y z var1 var2 var3_x var3_y var3_z ...
L2: x y z var1 var2 var3_x var3_y var3_z ...
.
.
.
Ln: x y z var1 var2 var3_x var3_y var3_z ...
```
where `n` is the cell number, and `var1(p)`, `var2(T)` `var3(U)`are the input variables:

```c++
foam2Columns -fields "(p T)"
```

Yes, it supports volVectorField:
```c++
foam2Columns -vectorFields "(U)"
```

Or together:
```c++
foam2Columns -fields "(p T)" -vectorFields "(U)"
```

Thanks for ZmengXu.
With his help, it supports lagrangian fields:
```c++
foam2Columns -lagrangianFields "(d T)"
```

Or together:
```c++
foam2Columns -fields "(p T)" -vectorFields "(U)" -lagrangianFields "(d T)"
```

The results are save in `$FOAM_CASE/postProcessing/foam2Columns/$time/Eulerian_*`, and `$FOAM_CASE/postProcessing/foam2Columns/$time/Lagrangian_*` for lagrangian fields.