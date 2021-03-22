#!/bin/sh
cargo build --target x86_64-unknown-linux-musl --release 
cp target/x86_64-unknown-linux-musl/release/qchem_g16 .
strip qchem_g16
