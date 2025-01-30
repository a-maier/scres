use std::env;
use std::path::PathBuf;

use cbindgen::Language;

fn main() {
    if cfg!(target_family = "unix") {
        write_c_header()
    }
}

fn write_c_header() {
    let out: PathBuf = [env::var("OUT_DIR").unwrap().as_str(), "scres.h"]
        .iter()
        .collect();

    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    let mut config = cbindgen::Config::default();
    config.cpp_compat = true;
    config.function.must_use =
        Some("__attribute__((warn_unused_result))".to_string());
    cbindgen::Builder::new()
        .with_config(config)
        .with_header(
            "/** C API for scres
 *
 * License: GPL 3.0 or later
 * Author: Andreas Maier <amaier@ifae.es>
*/",
        )
        .with_crate(crate_dir)
        .with_language(Language::C)
        .with_include_guard("SCRES_H")
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(out);
}
