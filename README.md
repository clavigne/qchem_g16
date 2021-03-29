# qchem_g16: a Gaussian 16 External interface to QChem

Run QChem under the control of Gaussian 16 with `qchem_g16`.

For example, you can use `qchem_g16` to perform optimizations and transition
state searches using the algorithms in Gaussian but the quantum chemical
methods in QChem.

## Installation

`qchem_g16` is available as a [statically-linked linux binary.](https://github.com/clavigne/qchem_g16/releases/download/0.2/qchem_g16) Just download it
and put it somewhere on `$PATH`.

On other platforms, you can build it using [Cargo and rustc,](https://doc.rust-lang.org/cargo/getting-started/installation.html)
```bash
git clone https://github.com/clavigne/qchem_g16
cd qchem_g16
cargo --build --release
```

## Usage

The `qchem_g16` binary is a well-behaved [Gaussian external
program.](https://gaussian.com/external/) It requires only one supplemental
argument `--rem FILE`, where `FILE` is a QChem input file containing a `$rem`
group. For example, the Gaussian input would look like this,

```
%NProcShared=4
# external="qchem_g16 --rem params.rem"
  opt=(calcfc, noraman)

H2 optimization

0 1
C          1.38274       -0.22147        0.00555
C          0.50658       -1.30672       -0.00829
C         -0.87138       -1.09065       -0.01412
C         -1.37328        0.21076       -0.00468
C         -0.49712        1.29590        0.00988
C          0.88087        1.07986        0.01426
H          2.45631       -0.38993        0.00964
H          0.89762       -2.32073       -0.01448
H         -1.55394       -1.93613       -0.02612
H         -2.44679        0.37914       -0.00858
H         -0.88812        2.30994        0.01784
H          1.56346        1.92540        0.02438

```

and the QChem input (`params.rem`) like this,

```
$rem
method hf
basis 6-31g
$end
```

Provided that `qchem` and `qchem_g16` are both on `$PATH`, running `g16` on
this file will optimize benzene but with energy, gradient and Hessian
evaluations done by QChem.

Q-Chem will automatically use as many threads as `%NProcShared`. 

## Features

- [X] Run Q-Chem as an external program from Gaussian.
- [X] Correct number of threads for Q-Chem.
- [X] Provides energies, gradients and Hessians.
- [ ] Provides dipole moments and polarizability.
- [X] Works with implicit solvation models.

## License

`qchem_g16` is [free and unencumbered software released in the public
domain.](./LICENSE) 


## Credits

Made with ♥ using [Rust](https://doc.rust-lang.org/book/) 🦀 by me, Cyrille
Lavigne.
