scres
=====

This is a wrapper crate around [cres](https://crates.io/crates/cres)
for easier access to single-step cell resampling with an extended C
API.

Installation
------------

### Rust

Add this to your Cargo.toml:

```toml
[dependencies]
scres = "0.0.2"
```

### C

To use the C API, the easiest way is to download the header and
libraries from the [github
releases](https://github.com/a-maier/scres/releases). Copy all `lib*`
files to a directory where your C compiler looks for libraries and
`cres.h` and `scres.h` to a directory for headers.

The documentation consists of comments in `scres.h`. For an example,
see the
[examples](https://github.com/a-maier/cres/tree/master/examples)
subdirectory.

### C (manual installation)

Instead of downloading the precompiled libraries, you can build them
by running

```
git clone https://github.com/a-maier/scres
cd scres
cargo build --release
```

You can find the compiled libraries in `target/release` and the
headers `cres.h` and `scres.h` in subdirectories of
`target/release/build`.
