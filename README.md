# TaxI2
Calculates genetic differences between DNA sequences

This version (v2.1.0) includes the legacy application, plus the new backend.
For the most recent GUI, please visit: https://github.com/iTaxoTools/TaxIGui


## Scores for alignment
The file `TaxI/scores.json` (found in `%APPDATA%\iTaxoTools` or in `$XDG_CONFIG_HOME$/`)
contains the scores used in the sequence alignment.

Each line has the format:
```
score_identifier<Tab>value
```

The possible scores are:
* `gap penalty`: Score to open a gap in the middle of a sequence
* `gap extend penalty`: Score to extend an existing gap in the middle of a sequence
* `end gap penalty`: Score to create a gap at an end of a sequence.
* `end gap extend penalty`: Score to extend a gap at an end of a sequence.
* `match score`: Score for matching nucleotides
* `mismatch score`: Score for non-matching nucleotides

## Choosing the backend for calculating distances between sequences
If the file `TaxI/options.json` (found in `%APPDATA%\iTaxoTools` or in `$XDG_CONFIG_HOME$/`) exists,
it should contains the line:
```
distance_calculation<Tab>BACKEND
```
* If `BACKEND` is `Rust`, the distance calculation uses functions written in Rust for better performance. Requires that the Rust module is compiled (see below).
* If `BACKEND` is `Python`, the distance calculation uses functions written in Python. Require Biopython.

## Compiling the Rust module

TaxI2 includes a module for calculating distances between sequences, which is written in Rust. It might need to be compiled before it can be used.

### Rust installation
[Instructions on the Rust site](https://www.rust-lang.org/tools/install).

Detailed instruction below.

#### Linux installation

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup update
```

#### Windows installation
Download and run [rustup-init.exe](https://static.rust-lang.org/rustup/dist/i686-pc-windows-gnu/rustup-init.exe).

You might need to install [Microsoft Build Tools for C++](https://visualstudio.microsoft.com/de/visual-cpp-build-tools), which should include the c++ compiler, Windows 10 SDK and English language files.

Find out whether your Python interpreter is 32- or 64-bit. For this run:
```
python
>>> import platform
>>> platform.architecture
(bits, linkage)
```
The `bits` will contain the required information.

If Python is 32-bit, do:
```
rustup default stable-i686-pc-windows-msvc
```

If Python is 64-bit, do:
```
rustup default stable-x86_64-pc-windows-msvc
```

### Compilation

Rust module compilation is now automated via setuptools.
