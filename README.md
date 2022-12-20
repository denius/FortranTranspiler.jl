# Fortran to Julia transpiler

Still not yet production ready translator. There's still a lot to come, but the most of job it can do.
However, the script can convert BLAS and LAPACK.

Things that still need to be completed:
* Formatted I/O. Now simple formats are implemented through `println` and `@printf`.
Cyclic formats, formats in variables are not implemented. The FLOATS are output as C, not Fortran.
* Annotating arrays as function arguments, namely using `@view` and `OffsetArrays` where appropriate.

# Usage

To try use:
```sh
julia FortranTranspiler.jl someFortranFile.f
```

or use
```sh
julia FortranTranspiler.jl *.f
```

or even
```sh
julia FortranTranspiler.jl .
```
to proceed all Fortran files in current folder.

In Linux environment you can use it by `FortranTranspiler.jl *.f` if you drop
`FortranTranspiler.jl` as executable script somewhere in PATH.

Some packages will be required to run the script:
```julia
julia> ] add Compat DataStructures DocStringExtensions Espresso JuliaFormatter Printf Setfield Tokenize
```

Usage options `FortranTranspiler.jl --help`:
```
Transpiler¹ converts the FORTRAN77 and partly FORTRAN90 source code into Julia.
It uses the pipe of both parsing and naive Regex replacements to do as much as possible,
but the output may need to further refinement.

Usage:
    FortranTranspiler.jl [--lowercase | --uppercase] ... [--] <filename1.f> <filename2.f> <somedir> ...
    FortranTranspiler.jl [--lowercase | --uppercase] ... [--] .
    FortranTranspiler.jl -h | --help

Samples:
    FortranTranspiler.jl .
    FortranTranspiler.jl somefile.f somedir
    FortranTranspiler.jl *.f* *.F*

Options:
    -h, --help         Show this screen.
    --version          Show version.
    -q, --quiet        Suppress all console output.
    -v, --verbose      Be verbose.
    -vv, --verbose     Be more verbose.
    -vvv, --verbose    Be yet more verbose.
    --preserveext      Preserve files extensions (suffix), append .jl: SOMEFILE.f90.jl
    --uppercase        Convert all identifiers to upper case.
    --lowercase        Convert all identifiers to lower case.
    --greeks           Replace the greek letter names thats starts the var names
                       with the corresponding unicode symbol, like DELTA1 -> δ1.
    --subscripts       Replace tail suffixes in vars names with unicode subscripts
                       like SOMEVAR_1 => SOMEVAR₁, if exist. Can be applied few times.
    --greeksubscripts  SOMEVAR_gamma => SOMEVARᵧ
    --                 The rest of the args is only the filenames and dirs.
    --formatting       Try to format with JuliaFormatter package.
    --dontfixcontinue  Do not try to insert ommited CONTINUE in the ancient fortran DO LABEL loops.
    --packarrays       Insert also Arrays in returned values.
    --double           Evaluate 1.0E0 as Float64, despite in fortran 1.0E0 is Float32.
    -n, --dry-run      Make the processing but don't write output ".jl" files.

[^1]: [Source-to-source compiler](https://en.wikipedia.org/wiki/Source-to-source_compiler).
```

The options `--greeks`, `--subscripts` and `--greeksubscripts` are useful for
some converting of Fortran ASCII letters variables into pretty Unicode variable names.

The `-v` option (and so on) is used to show the internal working machinery and possible errors
in sourced Fortran files.

After processing, each resulting file will be tried to parse with the Julia translator,
and if the `--quiet` option is not specified, you can see the unsuccessfull result of this.

There is especial option `--formatting`, it uses
[JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to try to pretty format
resulted Julia code file.
This package has a really powerful syntax parser, and if it failed to process the file,
then most likely the `FortranTranspiler.jl` was not very successful.
Although most of the work has already been done and some will remain to handmade.

For developers there is also available DEBUG mode to discover internals:
```sh
JULIA_DEBUG="FortranTranspiler" FortranTranspiler.jl program.f90
```
or
```sh
JULIA_DEBUG="FortranTranspiler" FortranTranspiler.jl -vv program.f90
```

# Work In Progess

It is also necessary to work out the function arguments, or rather array arguments.
This will require a double pass of source.

For formatted I/O, you need the `FortranFormattedIO.jl` package, alike
[FortranFiles.jl](https://github.com/traktofon/FortranFiles.jl). There is no such package yet,
it needs to be developed. For this possibility there is a branch `v2`.


And Fortran90 needs a lot of work.

# Acknowledgment

Inspired by https://gist.github.com/rafaqz/fede683a3e853f36c9b367471fde2f56
