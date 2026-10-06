# Synthetic example: −5% time-lapse velocity change

This example provides the ideal synthetic test used to demonstrate the basic performance of DSurfTomo-4D.

The T1 true model is kept unchanged, while a prescribed −5% shear-wave velocity anomaly is introduced into the T2 true model within a 3 × 3 × 3 grid region at depths of 0.6, 0.8, and 1.1 km.

The two epochs use the same observation geometry, station-pair configuration, and surface-wave period coverage. Therefore, differences between T1 and T2 arise only from the prescribed velocity change in the true model.

The same T1 and T2 synthetic observation setup is used for the Independent inversion, common-path comparison, and DSurfTomo-4D inversion.

## Directory structure

```text
m05_synthetic/
├── MOD.true_T1
├── MOD.true_T2
├── true_deltaVs.dat
├── true_deltaVs_percent.dat
├── indT1/
├── indT2/
└── 4D/
```

- `MOD.true_T1`: true velocity model at T1.
- `MOD.true_T2`: true velocity model at T2 with the prescribed −5% anomaly.
- `true_deltaVs.dat`: true absolute velocity change.
- `true_deltaVs_percent.dat`: true relative velocity change in percent.
- `indT1/`: conventional DSurfTomo inversion for T1.
- `indT2/`: conventional DSurfTomo inversion for T2.
- `4D/`: direct differential inversion using DSurfTomo-4D.

The T1 inversion result is converted to the full-grid reference model used by DSurfTomo-4D. The 4D inversion then uses the T1 and T2 observations together with this reference model to recover the time-lapse shear-wave velocity change.
