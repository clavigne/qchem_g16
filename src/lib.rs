use std::fs::{read_to_string, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result, Write};
use std::path::Path;

const CONV_FACTOR: f64 = 4.46552493159e-4;

fn parse_energy(p: &Path) -> Result<f64> {
    let f = File::open(p)?;
    let f = BufReader::new(f);
    let eline = f
        .lines()
        .filter_map(|x| match x {
            Ok(val) => val
                .starts_with(" The QM part of the energy is")
                .then(|| val),
            Err(_) => None,
        })
        .next();
    // now an Option(str)
    let out = match eline {
        Some(val) => Ok(val),
        None => Err(Error::new(
            ErrorKind::Other,
            "no energy line found in output file",
        )),
    };

    match out {
        Ok(val) => val
            .strip_prefix(" The QM part of the energy is")
            .unwrap()
            .trim()
            .parse::<f64>()
            .map_err(|_| Error::new(ErrorKind::Other, "failed to parse floats")),
        Err(e) => Err(e),
    }
}

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

pub fn qchem_translate_to_gaussian(
    natoms: u8,
    nder: u8,
    qchem_loc: &str,
    output_file: &str,
) -> Result<()> {
    let mut outfile = File::create(output_file)?;

    // energy
    let energy = parse_energy(&Path::new(&qchem_loc).join("qchem.out"))?;
    outfile.write(format!("{:+20.12}", energy).as_bytes())?;

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
