// use clap::{crate_version, App, Arg};
use std::env;
use std::fs::{read_to_string, File};
use std::io::{Error, ErrorKind, Result, Write};
use std::path::Path;

const CONV_FACTOR: f64 = 4.46552493159e-4;

fn parse_floats_from_file(n: u8, p: &Path) -> Result<Vec<f64>> {
    // Parse a vector of floats from a file.
    let nums: std::result::Result<Vec<_>, _> = read_to_string(p)?
        .split_whitespace()
        .map(|x| x.parse::<f64>())
        .collect();

    match nums {
        Ok(i) => {
            if i.len() == n.into() {
                Ok(i)
            } else {
                Err(Error::new(ErrorKind::Other, "fewer floats than expected"))
            }
        }
        Err(_) => Err(Error::new(ErrorKind::Other, "failed to parse floats")),
    }
}

fn write_gaussian_output(natoms: u8, nder: u8, qchem_loc: &str, output_file: &str) -> Result<()> {
    let mut outfile = File::create(output_file)?;

    // energy
    outfile.write(format!("{:+20.12}", 0.0).as_bytes())?;

    // dipole
    outfile.write(format!("{:+20.12}{:+20.12}{:+20.12}\n", 0.0, 0.0, 0.0).as_bytes())?;

    // derivatives
    if nder > 0 {
        let mut data =
            parse_floats_from_file(3 * natoms, &Path::new(&qchem_loc).join("efield.dat"))?;
        for _ in 0..natoms {
            for el in data.drain(..3) {
                outfile.write(format!("{:+20.12}", el).as_bytes())?;
            }
            outfile.write("\n".as_bytes())?;
        }
        // polarizability + dip derivative (6 + 9 * Natoms)
        for _ in 0..(2 + 3 * natoms) {
            outfile.write(format!("{:+20.12}{:+20.12}{:+20.12}\n", 0.0, 0.0, 0.0).as_bytes())?;
        }
    }

    // hessian
    if nder > 1 {
        let n_hessian = (3 * natoms) * (3 * natoms + 1) / 2;
        let mut data =
            parse_floats_from_file(n_hessian, &Path::new(&qchem_loc).join("hessian.dat"))?;
        for _ in 0..(n_hessian / 3) {
            for el in data.drain(..3) {
                outfile.write(format!("{:+20.12}", el * CONV_FACTOR).as_bytes())?;
            }
            outfile.write("\n".as_bytes())?;
        }
    }
    Ok(())
}

fn main() {
    let nder = 2;
    let qchem_loc = env::var("QCHEM_RUNDIR").unwrap_or(".".to_string());
    let output = "gamout";

    write_gaussian_output(3, nder, &qchem_loc, output);

    // let matches = App::new("qchem-gaussian-collab")
    //     .author("Cyrille Lavigne <cyrille.lavigne@mail.utoronto.ca>")
    //     .about("bla")
    //     .version(crate_version!())
    //     .arg(
    //         Arg::new("nel")
    //             .value_name("NEL")
    //             .about("Total number of electrons.")
    //             .required(true)
    //             .takes_value(true),
    //     )
    //     .arg(
    //         Arg::new("homo")
    //             .long("homo")
    //             .value_name("IORBS")
    //             .multiple(true)
    //             .about("Indices of HOMO orbitals.")
    //             .required(true)
    //             .takes_value(true),
    //     )
    //     .arg(
    //         Arg::new("lumo")
    //             .long("lumo")
    //             .value_name("IORBS")
    //             .multiple(true)
    //             .about("Indices of LUMO orbitals.")
    //             .required(true)
    //             .takes_value(true),
    //     )
    //     .arg(
    //         Arg::new("norb")
    //             .long("norb")
    //             .value_name("NORBS")
    //             .about("Total number of orbitals.")
    //             .takes_value(true),
    //     )
    //     .get_matches();

    // // Get the orbital indices
    // let homos = BTreeSet::from_iter(
    //     matches
    //         .values_of_t::<usize>("homo")
    //         .unwrap_or_else(|e| e.exit()),
    // );
    // let lumos = BTreeSet::from_iter(
    //     matches
    //         .values_of_t::<usize>("lumo")
    //         .unwrap_or_else(|e| e.exit()),
    // );

    // // number of total electrons, active orbitals and active electrons
    // let nel = match matches.value_of_t::<usize>("nel") {
    //     Ok(i) => i,
    //     Err(e) => e.exit(),
    // };
}
