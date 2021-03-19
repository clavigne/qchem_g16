use std::convert::TryInto;
use std::fs::{read_to_string, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result, Write};
use std::path::Path;
use std::str::FromStr;

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

fn parse_nums_from_str<T: FromStr>(n: u8, data: String) -> Result<Vec<T>> {
    // Parse a vector of floats from a file.
    let nums: std::result::Result<Vec<_>, _> =
        data.split_whitespace().map(|x| x.parse::<T>()).collect();

    match nums {
        Ok(i) => {
            if i.len() == n.into() {
                Ok(i)
            } else {
                Err(Error::new(
                    ErrorKind::Other,
                    format!("expected {} values, got {}", n, i.len()),
                ))
            }
        }
        Err(_) => Err(Error::new(ErrorKind::Other, "failed to parse values")),
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
        let mut data = parse_nums_from_str::<f64>(
            3 * natoms,
            read_to_string(Path::new(&qchem_loc).join("efield.dat"))?,
        )?;
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
        let mut data = parse_nums_from_str::<f64>(
            n_hessian,
            read_to_string(Path::new(&qchem_loc).join("hessian.dat"))?,
        )?;
        for _ in 0..(n_hessian / 3) {
            for el in data.drain(..3) {
                outfile.write(format!("{:+20.12}", el * CONV_FACTOR).as_bytes())?;
            }
            outfile.write("\n".as_bytes())?;
        }
    }
    Ok(())
}

#[derive(Debug)]
pub struct Calculation {
    pub natoms: usize,
    pub nder: usize,
    pub charge: i8,
    pub spin: i8,
    pub z: Vec<u8>,
    pub coords: Vec<[f64; 3]>,
}

pub fn parse_gau_ein(infile: &str) -> Result<Calculation> {
    let gaussfile = read_to_string(infile)?;
    let mut gauss = gaussfile.lines();
    if let Some(header) = gauss.next() {
        // Parse
        let entries = parse_nums_from_str::<i8>(4, header.to_string())?;
        let natoms: usize = entries[0].try_into().unwrap();
        let nder: usize = entries[1].try_into().unwrap();
        let charge: i8 = entries[2];
        let spin: i8 = entries[3];
        let mut coords = Vec::new();
        let mut zvals = Vec::<u8>::new();

        for _ in 0..natoms {
            if let Some(line) = gauss.next() {
                let (start, end) = line.split_at(11);
                let atom = parse_nums_from_str::<u8>(1, start.to_string())?[0];
                let vals = parse_nums_from_str::<f64>(4, end.to_string())?;
                coords.push([vals[0], vals[1], vals[2]]);
                zvals.push(atom);
            } else {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Gaussian input file is truncated",
                ));
            }
        }
        Ok(Calculation {
            natoms: natoms,
            nder: nder,
            charge: charge,
            spin: spin,
            z: zvals,
            coords: coords,
        })
    } else {
        Err(Error::new(ErrorKind::Other, "Gaussian input file is empty"))
    }
}
