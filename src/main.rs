use clap::{crate_version, App, Arg};
use extgauss::{parse_gau_ein, qchem_translate_to_gaussian};
use num_cpus;
use std::env;
use std::fs::{read_to_string, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result, Write};
use std::path::Path;
use std::process::Command;

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

    // Get argument values
    let infile = matches.value_of("InputFile").unwrap();
    let outfile = matches.value_of("OutputFile").unwrap();
    let msgfile = matches.value_of("MsgFile").unwrap();
    let remfile = matches.value_of("rem").unwrap();

    // get some env variables
    let qchem_loc = env::var("EXTGAUSS_QCHEM_RUNDIR").unwrap_or(".".to_string());
    let qchem_exe = env::var("EXTGAUSS_QCHEM_EXE").unwrap_or("qchem".to_string());
    let num_threads = env::var("OMP_NUM_THREADS").unwrap_or(num_cpus::get().to_string());
    let qchem_args = ["-nt", &num_threads];

    let mut msgs = File::create(msgfile)?;
    msgs.write(
        format!(
            "-+---------------------------------------------- extgaussian-rs v{} \n",
            crate_version!()
        )
        .as_bytes(),
    )?;
    msgs.write(format!(" |  input:     {}\n", infile).as_bytes())?;
    msgs.write(format!(" |  output:    {}\n", outfile).as_bytes())?;
    msgs.write(format!(" |  $rem file: {}\n", remfile).as_bytes())?;
    msgs.write(format!(" |  rundir:    {}\n", qchem_loc).as_bytes())?;
    msgs.write(format!(" |  calling:   {}\n", qchem_exe).as_bytes())?;
    msgs.write(format!(" |  args:      {:?}\n", qchem_args).as_bytes())?;

    // Load calculation details
    let calc = parse_gau_ein(infile)?;
    let jobtype = match calc.nder {
        0 => "sp",
        1 => "force",
        2 => "freq",
        _ => "",
    };

    let scf_guess = match Path::new(&qchem_loc).join("qchem.scratch").is_dir() {
        true => "scf_guess read\n",
        false => "",
    };

    // Make molecule data
    let mol = format!(
        "$molecule\n\
         {} {}\n\
         {}\n\
         $end\n",
        calc.charge,
        calc.spin,
        calc.get_geometry()
    );

    // Make qchem input
    let rem = format!(
        "{}\n\
         $rem\n{}\n\
         {}\
         jobtype {}\n\
         qm_mm true\n\
         qmmm_print true\n\
         hess_and_grad true\n\
         $end\n",
        mol,
        read_to_string(remfile)?.trim(),
        scf_guess,
        jobtype
    );

    msgs.write(
        "-+----------------------------------------------------- input to qchem\n\
         "
        .as_bytes(),
    )?;
    for line in rem.lines() {
        msgs.write(format!(" | {}\n", line).as_bytes())?;
    }
    msgs.write(
        "=+======================================================= qchem output\n\
         "
        .as_bytes(),
    )?;

    let qchem = Command::new(qchem_exe)
        .args(&qchem_args)
        .current_dir(&qchem_loc)
        .output();

    let qchem_out = match qchem {
        Ok(val) => std::str::from_utf8(&val.stdout).unwrap().to_string(),
        Err(e) => format!("Calling QChem failed\n {:?}\n", e),
    };
    msgs.write(&qchem_out.as_bytes());
    qchem_translate_to_gaussian(outfile, &calc, &qchem_loc, &qchem_out)?;

    Ok(())
}
