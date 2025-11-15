# PIMDkk

経路積分分子動力学プログラム

## Simulation modes

Specify the simulation engine through the `$Isimulation` keyword in `input.inp`.

- `0`: Path-integral molecular dynamics (PIMD)
- `1`: Ring-polymer molecular dynamics (RPMD)
- `2`: Centroid molecular dynamics (CMD)
- `3`: Path-integral hybrid Monte Carlo (PIHMC)
- `10`: Classical MD (if you want the classical limit just set `Nbead=1`)



