# LJ-PBC Pressure Output Check

Local, ignored smoke test for the LJ-PBC pressure-output path.

It runs a two-step classical `Iforce=32` argon box and checks that
`Result/pressure.dat` is created with finite pressure columns, the expected
`Nbead=1` equality `P_cv == P_prim`, and the expected fixed-cell volume.

Run:

```sh
make test
```

Expected result: `PASS`.
