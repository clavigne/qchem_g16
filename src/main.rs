use clap::{crate_version, App, Arg};
use extgauss::qchem_translate_to_gaussian;
use num_cpus;
use std::env;
use std::fs::{read_to_string, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result, Write};
use std::path::Path;

fn main() -> Result<()> {
    let matches = App::new("extgaussian-rs")
        .author("Cyrille Lavigne <cyrille.lavigne@mail.utoronto.ca>")
        .about("TODO")
        .version(crate_version!())
        .arg("<Layer>              'Layer of an ONIOM calculation.'")
        .arg("<InputFile>          'Input to external program.'")
        .arg("<OutputFile>         'Output from external program.'")
        .arg("<MsgFile>            'Messages for Gaussian.'")
        .arg("[FChkFile]           'Formatted checkpoint file.'")
        .arg("[MatElFile]          'Matrix elements.'")
        .arg(
            Arg::new("rem")
                .long("rem")
                .value_name("REMFILE")
                .about("File containing $rem options for QChem.")
                .required(true)
                .takes_value(true),
        )
        .get_matches();

    // Get InputFile
    let infile = matches.value_of("InputFile").unwrap();
    let outfile = matches.value_of("OutputFile").unwrap();
    let msgfile = matches.value_of("MsgFile").unwrap();
    let remfile = matches.value_of("rem").unwrap();
    let qchem_loc = env::var("EXTGAUSS_QCHEM_RUNDIR").unwrap_or(".".to_string());
    let num_threads = env::var("OMP_NUM_THREADS").unwrap_or(num_cpus::get().to_string());
    let qchem_call =
        env::var("EXTGAUSS_QCHEM_CALL").unwrap_or(format!("qchem -nt {}", num_threads));

    let mut msgs = File::create(msgfile)?;
    msgs.write(
        format!(
            "-+- extgaussian-rs v{} ----------------------------\n",
            crate_version!()
        )
        .as_bytes(),
    )?;
    msgs.write(format!(" |  input:     {}\n", infile).as_bytes())?;
    msgs.write(format!(" |  output:    {}\n", outfile).as_bytes())?;
    msgs.write(format!(" |  $rem file: {}\n", remfile).as_bytes())?;
    msgs.write(format!(" |  rundir:    {}\n", qchem_loc).as_bytes())?;
    msgs.write(format!(" |  calling:   {}\n", qchem_call).as_bytes())?;

    // Load rem lines

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

    // let nder = 2;
    // let output = "gamout";
    // qchem_translate_to_gaussian(3, nder, &qchem_loc, output)?;
    Ok(())
}
