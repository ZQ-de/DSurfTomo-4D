# Real-data input template

This directory provides an example configuration for applying DSurfTomo-4D to real observational data.

No real observational dataset is distributed with this repository. Users should prepare:

- a T1 reference model named `MOD`;
- matched common-path dispersion datasets for T1 and T2;
- a `DSurfTomo.in` file with `synthetic flag = 0`.

The filenames of the T1 and T2 dispersion datasets are user-defined and should match the first two lines of `DSurfTomo.in`.

The T1 reference model can be obtained from a conventional DSurfTomo inversion of the T1 data. The resulting model can be converted to the full `MOD` grid format using `make_4d_MOD_from_T1.py`.

After preparing the required files, run:

```bash
./DSurfTomo-4D
```

The main output is:

```
DSurfTomo.deltaVs.dat
```

which contains the recovered shear-wave velocity change between T1 and T2.
