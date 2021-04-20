use anyhow::{Context, Result};
use clap::{crate_version, App};
use qchem_g16::Calculation;
use std::env;
use std::fs::{read_to_string, File};
use std::io::Write;
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
        .arg("-r, --rem=<REM>      'File containing $rem options for QChem'")
        .arg("-d, --dir=[DIR]      'QChem run directory [default=.]'")
        .arg("-e, --exe=[EXE]      'QChem executable invocation [default=\"qchem\"]'")
        .arg("<Layer>              'Layer of an ONIOM calculation'")
        .arg("<InputFile>          'Input to external program'")
        .arg("<OutputFile>         'Output from external program'")
        .arg("<MsgFile>            'Messages for Gaussian'")
        .arg("[FChkFile]           'Formatted checkpoint file'")
        .arg("[MatElFile]          'Matrix elements'")
        .arg("--dry                'Do not actually run Q-Chem, but do parse any results present.'")
        .get_matches();

    // Get argument values
    let infile = matches.value_of("InputFile").unwrap();
    let outfile = matches.value_of("OutputFile").unwrap();
    let msgfile = matches.value_of("MsgFile").unwrap();
    let remfile = matches.value_of("rem").unwrap();
    let dry = matches.is_present("dry");

    // get some env variables
    let qchem_loc = matches.value_of("dir").unwrap_or(".");
    let qchem_exe = matches.value_of("exe").unwrap_or("qchem");
    let num_threads = env::var("OMP_NUM_THREADS").unwrap_or("1".to_string());
    let qchem_args = ["-nt", &num_threads];

    // paths that we will need
    let qchem_dir = Path::new(&qchem_loc);
    let qchem_inp = qchem_dir.join("qchem.inp");
    let qchem_out = qchem_dir.join("qchem.out");
    let qchem_scratch = qchem_dir.join("qchem.scratch");
    eprintln!("ready.");

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
    eprintln!("loading gaussian input");
    let gaussfile = read_to_string(infile).context(format!("No EIn file at {}", infile))?;
    let calc = Calculation::from_ext(&gaussfile)?;
    let (jobtype, hess_and_grad) = match calc.nder {
        0 => ("sp", ""),
        1 => ("force", ""),
        2 => ("freq", "hess_and_grad true\nvibman_print 6\n"),
        _ => ("", ""),
    };
    eprintln!("\tdone. calculation: {}", jobtype);

    eprintln!("looking for restart data");
    let scf_guess = match qchem_scratch.is_dir() {
        true => {
            eprintln!("\trestart data found");
            "scf_guess read\n"
        }
        false => {
            eprintln!("\tno restart data");
            ""
        }
    };

    eprintln!("reading & parsing rem file");
    let mol = calc.qchem_molecule();
    let parameters = read_to_string(remfile)?;
    let (rem, extras) = extract_rem(parameters.trim());
    eprintln!("\tdone");

    // Make qchem input
    eprintln!("building new rem group");
    let rem = format!(
        "{}\n\
         {}{}\
         jobtype {}\n\
         input_bohr true\n\
         sym_ignore true\n\
         gui 2\n\
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
    eprintln!("\tdone");

    let mut qchem_cli = Command::new(qchem_exe);
    let qchem_cli = qchem_cli
        .args(&qchem_args)
        // these are not relative paths because we are calling qchem from $run_dir
        .arg("qchem.inp")
        .arg("qchem.out")
        .arg("qchem.scratch")
        .current_dir(&qchem_dir);

    let (qchem_stdout, do_panic) = match dry {
        false => {
            eprintln!("qchem: {:#?}", qchem_cli);
            eprintln!("running...");
            let qchem = qchem_cli.output();
            eprintln!("\tdone");
            match qchem {
                Ok(val) => (std::str::from_utf8(&val.stdout).unwrap().to_string(), false),
                Err(e) => (format!("calling QChem failed\n {:?}\n", e), true),
            }
        }
        true => {
            eprintln!("qchem: {:#?}", qchem_cli);
            eprintln!("dry run!");
            eprintln!("\tdone");
            (
                "DRY RUN\n QChem will not be executed.\n\n".to_string(),
                false,
            )
        }
    };

    msgs.write(&qchem_stdout.as_bytes())?;
    if do_panic {
        panic!("failed to launch qchem");
    }

    eprintln!("loading qchem output");
    let qchem_output = read_to_string(&qchem_out)?;
    eprintln!("\tdone");

    eprintln!("translating qchem output");
    let translation = calc.translate_qchem(&qchem_output)?;
    eprintln!("\tdone");

    eprintln!("writing for gaussian and terminating");
    File::create(&outfile)?.write(translation.as_bytes())?;
    eprintln!("\tALL DONE");
    Ok(())
}
