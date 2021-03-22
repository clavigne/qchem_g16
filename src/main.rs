use clap::{crate_version, App, Arg};
use num_cpus;
use qchem_g16::{parse_gau_ein, qchem_translate_to_gaussian};
use std::env;
use std::fs::{read_to_string, File};
use std::io::{Result, Write};
use std::path::Path;
use std::process::Command;

fn extract_rem(input: &str) -> (String, String) {
    let (before, rem1) = input.split_at(input.find("$rem").expect("no $rem input"));
    let (rem, after) = rem1.split_at(rem1.find("$end").expect("no $end to $rem input"));
    (rem.to_string(), [after, "\n", before].concat())
}

fn main() -> Result<()> {
    let matches = App::new("qchem_g16")
        .author("Cyrille Lavigne <cyrille.lavigne@mail.utoronto.ca>")
        .about(
            "\
\n\
This is a Gaussian 16 interface to use Q-Chem as an external calculator. It is
meant to be called from within Gaussian 16 using the External keyword and the
following invocation (where qchem_g16 is on $PATH),

# external=\"qchem_g16 --rem params.rem\"

Here, params.rem is a QChem parameter file. It should contain at least a $rem
section with method and basis arguments. It should *not* include jobtype, nor
should it include a $molecule section; these will be filled in by this script.
",
        )
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

    // paths that we will need
    let qchem_dir = Path::new(&qchem_loc);
    let qchem_inp = qchem_dir.join("qchem.inp");
    let qchem_out = qchem_dir.join("qchem.out");
    let qchem_scratch = qchem_dir.join("qchem.scratch");

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
    let (jobtype, hess_and_grad) = match calc.nder {
        0 => ("sp", ""),
        1 => ("force", ""),
        2 => ("freq", "hess_and_grad true\n"),
        _ => ("", ""),
    };

    let scf_guess = match qchem_scratch.is_dir() {
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
    let parameters = read_to_string(remfile)?;
    let (rem, extras) = extract_rem(parameters.trim());

    // Make qchem input
    let rem = format!(
        "{}\n\
         {}{}\
         jobtype {}\n\
         qm_mm true\n\
         qmmm_print true\n\
         input_bohr true\n\
         {}{}",
        mol, rem, scf_guess, jobtype, hess_and_grad, extras
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

    File::create(&qchem_inp)?.write(rem.as_bytes())?;

    let qchem = Command::new(qchem_exe)
        .args(&qchem_args)
        .arg(&qchem_inp)
        .arg(&qchem_out)
        .arg(&qchem_scratch)
        .current_dir(&qchem_dir)
        .output();

    let qchem_stdout = match qchem {
        Ok(val) => std::str::from_utf8(&val.stdout).unwrap().to_string(),
        Err(e) => format!("Calling QChem failed\n {:?}\n", e),
    };
    msgs.write(&qchem_stdout.as_bytes())?;
    qchem_translate_to_gaussian(outfile, &calc, &qchem_dir, &qchem_out)?;

    // delete dat files
    let _ = Command::new("rm")
        .arg("efield.dat")
        .arg("hessian.dat")
        .current_dir(&qchem_dir)
        .output();

    Ok(())
}
