#[macro_use]
extern crate anyhow;
use anyhow::{Context, Result};
// std library
use std::convert::TryInto;
use std::fmt::Debug;
use std::str::FromStr;

fn parse_nums_from_str<T: FromStr + Debug, const N: usize>(data: String) -> Result<[T; N]> {
    // Parse a vector of floats from a file.
    let nums: std::result::Result<Vec<_>, _> =
        data.split_whitespace().map(|x| x.parse::<T>()).collect();

    match nums {
        Ok(vec) => {
            if vec.len() == N.into() {
                Ok({
                    let out: [T; N] = vec.try_into().unwrap();
                    out
                })
            } else {
                Err(anyhow!("expected {} values, got {}", N, vec.len()))
            }
        }
        Err(_) => Err(anyhow!("failed to parse values")),
    }
}

// QChem parsers
fn parse_energy(qchem_out: &str) -> Result<f64> {
    let tag = "Total energy in the final basis set =";
    let out = qchem_out
        .lines()
        .filter_map(|x| x.trim_start().strip_prefix(tag))
        .next()
        .context("no energy line found in qchem output")?
        .trim()
        .parse::<f64>()
        .context("failed to parse energy into float")?;
    Ok(out)
}

fn parse_gradient(natom: usize, qchem_chk: &str) -> Result<Vec<f64>> {
    let nlines = 3 * natom / 5 + 2;
    let grads: Result<Vec<_>, _> = qchem_chk
        .lines()
        .skip_while(|x| !x.trim_start().starts_with("Cartesian Forces"))
        .take(nlines)
        .flat_map(|x| x.split_whitespace())
        .filter(|x| x.contains("."))
        .map(|x| x.parse::<f64>())
        .collect();
    let grads = grads?;
    let mut vals = grads.iter();
    let mut out = vec![];
    for _ in 0..3 * natom {
        let k = vals.next().unwrap();
        out.push(*k);
    }
    Ok(out)
}

fn parse_hessian(natom: usize, qchem_chk: &str) -> Result<Vec<f64>> {
    let ncoord = natom * 3;
    let nel = ncoord * (ncoord + 1) / 2;
    let nlines = nel / 5 + 2;
    let hess: Result<Vec<_>, _> = qchem_chk
        .lines()
        .skip_while(|x| !x.trim_start().starts_with("Cartesian Force Constants"))
        .take(nlines)
        .flat_map(|x| x.split_whitespace())
        .filter(|x| x.contains("."))
        .map(|x| x.parse::<f64>())
        .collect();
    let hess = hess?;
    let mut vals = hess.iter();
    let mut out = vec![];
    for _ in 0..nel {
        let k = vals.next().unwrap();
        out.push(*k);
    }
    Ok(out)
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

impl Calculation {
    pub fn geometry(&self) -> String {
        let mut output = String::new();
        for i in 0..self.natoms {
            output.push_str(&format!(
                "{}   {}   {}   {}\n",
                self.z[i], self.coords[i][0], self.coords[i][1], self.coords[i][2]
            ));
        }
        output.trim().to_string()
    }

    pub fn qchem_molecule(&self) -> String {
        format!(
            "$molecule\n\
         {} {}\n\
         {}\n\
         $end\n",
            self.charge,
            self.spin,
            self.geometry()
        )
    }

    pub fn from_ext(gaussfile: &str) -> Result<Calculation> {
        let mut gauss = gaussfile.lines();
        let header = gauss.next().context("QChem output is truncated")?;

        // Parse
        let entries = parse_nums_from_str::<i8, 4>(header.to_string())?;
        let natoms: usize = entries[0].try_into().unwrap();
        let nder: usize = entries[1].try_into().unwrap();
        let charge: i8 = entries[2];
        let spin: i8 = entries[3];
        let mut coords = Vec::new();
        let mut zvals = Vec::<u8>::new();

        for _ in 0..natoms {
            if let Some(line) = gauss.next() {
                let (start, end) = line.split_at(11);
                let atom = parse_nums_from_str::<u8, 1>(start.to_string())?[0];
                let vals = parse_nums_from_str::<f64, 4>(end.to_string())?;
                coords.push([vals[0], vals[1], vals[2]]);
                zvals.push(atom);
            } else {
                return Err(anyhow!("Gaussian input file is truncated"));
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
    }

    pub fn translate_qchem(&self, qchem_out: &str, qchem_fchk: &str) -> Result<String> {
        let mut output = String::new();
        let nder = self.nder;
        let natoms = self.natoms;

        // energy
        eprintln!("\tparsing energy");
        output.push_str(&format!("{:+20.12}", parse_energy(qchem_out)?));
        eprintln!("\t\tdone");

        // dipole
        output.push_str(&format!("{:+20.12}{:+20.12}{:+20.12}\n", 0.0, 0.0, 0.0));

        // derivatives
        if nder > 0 {
            eprintln!("\tparsing gradient");
            let grads = parse_gradient(natoms, qchem_fchk)?;
            eprintln!("\t\tdone");
            let mut count = 0;
            for el in grads {
                output.push_str(&format!("{:+20.12}", el));
                count += 1;
                if count == 3 {
                    output.push('\n');
                    count = 0
                }
            }
            // polarizability + dip derivative (6 + 9 * Natoms)
            for _ in 0..(2 + 3 * natoms) {
                output.push_str(&format!("{:+20.12}{:+20.12}{:+20.12}\n", 0.0, 0.0, 0.0));
            }
        }

        // hessian
        if nder > 1 {
            eprintln!("\tparsing hessian");
            let hess = parse_hessian(natoms, qchem_fchk)?;
            eprintln!("\t\tdone");
            let mut count = 0;
            for el in hess {
                output.push_str(&format!("{:+20.12}", el));
                count += 1;
                if count == 3 {
                    output.push('\n');
                    count = 0
                }
            }
        }
        Ok(output)
    }
}
