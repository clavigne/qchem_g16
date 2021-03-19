#!/bin/sh
cargo build --target x86_64-unknown-linux-musl --release 
cp target/x86_64-unknown-linux-musl/release/extgauss-rs .
strip extgauss-rs
