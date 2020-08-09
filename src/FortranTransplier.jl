#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

#=

# KNOWN ISSUES
* julia's `for` loop counter is local scope variable, unlike fortran.
  Thus the codes that use the value of counter after the loop will be broken
  and should be fixed via set loop counter before `for` to afterlast value. May be fixed!
* due to julia's GC, some implicit variables initialized inside the loop/ blocks
  may disappear upon they exit these blocks. Consider initilizing them with
  the appropriate value in the initial part of the function code.
* julia can't propagate back changed values of scalars via functions args, unlike fortran.
  Thus such changed scalars should be returned via functions `return` statement.
* all strings in fortran are `Vector{UInt8}` and should stay same type to preserve
  assignment statements. (Use `join()` and `collect()` to convert to strings and back)
* whitespaces matter despite of the old fortran, which can ignore them all:
    * in julia whitespaces between name of functions/arrays or braces are not allowed.
    * whitespaces in brackets is unwanted (due to https://github.com/JuliaLang/julia/issues/14853)
* some unserviceable comments are cutted off (eg. in expanded lines continuations).
* this script is mostly for the fixed-form fortran and in the free-form is not tested yet.
* in the `DATA` statement can occur uncatched repetitions like `DATA SOMEARRAY/8*0,1,2/`
* 'FORMAT' conversion is unaccomplished
* implied do-loops are not always caught
* there is no OPENMP in Julia

# TODO
    ~~custom precision like `dzero=0.0_psb_dpk_, done=1.0_psb_dpk_` in psb_const_mod.F90~~
        and `INTEGER,PARAMETER :: RP = SELECTEDREALKIND(15)`
    few PROGRAMs in one file
    use SOMEMODULE -> use .SOMEMODULE
    use SOMEMODULE, only: SOMEFUN -> import .SOMEMODULE: SOMEFUN
    RECURSIVE
    `FUNCTION f(x) RESULT(v)` can use as `f` as `v` as the returned value?
    parse DO...
    parse IF... ???
    insert running BLOCK DATA at startup in main()


    return intent(out):
    LSAME
    XERBLA
    P_VALUE
    RWK/IWK
    M_GETMEM

    see also https://www.math.fsu.edu/~dmandel/Fortran/Chap4.pdf
=#


module FortranTransplier

export convert_fortran

using Compat
using DataStructures
using DocStringExtensions
using Espresso
using JuliaFormatter
using Printf
using Tokenize
import Tokenize.Tokens.Token


# module for CLI running
module CLI

Base.Experimental.@optlevel 0

using ..FortranTransplier
using Compat
using DataStructures
using Glob
using Printf

const usage =
"""
    Transpiler¹ converts the fortran77 and partly fortran90 code into julia.
    It uses both parsing and naive regex replacements to do as much as possible,
    but the output may need further refinement.

    Usage:
        $(basename(@__FILE__)) [--lowercase | --uppercase] ... [--] <filename1.f> <filename2.f> ...
        $(basename(@__FILE__)) [--lowercase | --uppercase] ... [--] .
        $(basename(@__FILE__)) -h | --help

    Samples:
        $(basename(@__FILE__)) .
        $(basename(@__FILE__)) somefile.f somedir
        $(basename(@__FILE__)) *.f* *.F*

    Options:
        -h, --help         Show this screen.
        --version          Show version.
        -q, --quiet        Suppress all console output.
        -v, --verbose      Be verbose.
        -v, --verbose      Be more verbose.
        -v, --verbose      Be yet more verbose.
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
        --packarrays
        --strings
        --dropdeclarations
        --single           Evaluate 1.0E0 as Float32. Double precision in fortran is 1.0D0
        --omitimplicit     Omit implicit scalars initialization.
        -n, --dry-run      Make the processing but don't write output \".jl\" files.

    Inspired by https://gist.github.com/rafaqz/fede683a3e853f36c9b367471fde2f56

    [^1]: [Source-to-source compiler](https://en.wikipedia.org/wiki/Source-to-source_compiler).
"""

function main(args = String[])
    opts, fnames = parseargs(args)
    case = opts["--uppercase"] ? uppercase :
           opts["--lowercase"] ? lowercase : identity

    for fname in fnames
        opts["--veryverbose"] && println()
        opts["--verbose"] && println("$(fname):")
        code = open(fname) |> read |> String
        result = convert_fortran(code, casetransform=case, quiet=opts["--quiet"],
                                 verbose=opts["--verbose"], veryverbose=opts["--veryverbose"],
                                 veryveryverbose=opts["--veryveryverbose"],
                                 greeks=opts["--greeks"], subscripts=opts["--subscripts"],
                                 greeksubscripts=opts["--greeksubscripts"])
        opts["--dry-run"] || write(splitext(fname)[1] * ".jl", result)
        if opts["--formatting"]
            try
                result = JuliaFormatter.format_text(result, margin = 192)
                opts["--dry-run"] || write(splitext(fname)[1] * ".jl", result)
            catch
                @info("$(splitext(fname)[1] * ".jl")\nFORMATTING IS NOT SUCCESSFULL\n")
            end
        end
    end
    return nothing
end

function parseargs(args)

    function therecanbeonlyone!(opts, v1, v2)
        opts[v1] = opts[v1] || opts[v2]
        delete!(opts, v2)
    end

    fnames = Vector{String}()
    opts = OrderedDict{String, Any}()
    opts["--quiet"]           = opts["-q"]   = false
    opts["--verbose"]         = opts["-v"]   = false
    opts["--veryverbose"]     = opts["-vv"]  = false
    opts["--veryveryverbose"] = opts["-vvv"] = false
    opts["--dry-run"]         = opts["-n"]   = false
    opts["--uppercase"]       = false
    opts["--lowercase"]       = false
    opts["--greeks"]          = false
    opts["--greeksubscripts"] = false
    opts["--subscripts"]      = 0
    opts["--formatting"]      = false
    #opts["--returnnothing"]   = false
    opts["--packarrays"]      = false
    opts["--strings"]         = false
    #opts["--"]  = false

    while length(args) > 0
        if first(args) == "--help" || first(args) == "-h"
            println("$(basename(@__FILE__)):\n$usage\n")
            exit(0)
        elseif haskey(opts, first(args))
            opts[first(args)] += 1
            popfirst!(args)
        elseif first(opts) == "--" # end of options
            popfirst!(args)
            break
        elseif occursin(r"^--\w\w+$", first(args)) || occursin(r"^-\w$", first(args))
            @error("unrecognized option in opts: '$(first(args))'")
            exit(1)
        else
            push!(fnames, first(args))
            popfirst!(args)
        end
    end
    while length(args) > 0
        push!(fnames, first(args))
        popfirst!(args)
    end

    rxfortranfiles = (r"\.for$"i, r"\.f$"i, r"\.f77$"i, r"\.f90$"i)

    flist = Vector{String}()
    for fname in fnames
        if isfile(fname)
            push!(flist, normpath(fname))
        elseif isdir(fname)
            for (root, dirs, files) in walkdir(fname), f in files
                any(map(r->occursin(r,f), rxfortranfiles)) && append!(flist, (joinpath(root, f),))
            end
        else
            @warn "No such file or directory: $(fname)"
        end
    end
    unique!(flist)

    if opts["--verbose"] == 2 || opts["-v"] == 2 || opts["-vv"] == 1
        opts["--veryverbose"] = opts["--verbose"] = opts["-vv"] = opts["-v"] = true
    elseif opts["--verbose"] == 3 || opts["-v"] == 3 || opts["-vvv"] == 1
        opts["--veryveryverbose"] = opts["--veryverbose"] = opts["--verbose"] =
            opts["-vvv"] = opts["-vv"] = opts["-v"] = true
    end
    for (k,v) in opts
        !isa(opts[k], Bool) && opts[k] == 1 && (opts[k] = true)
    end

    therecanbeonlyone!(opts, "--veryveryverbose", "-vvv")
    therecanbeonlyone!(opts, "--veryverbose"    , "-vv")
    therecanbeonlyone!(opts, "--verbose"        , "-v" )
    therecanbeonlyone!(opts, "--quiet"          , "-q" )
    therecanbeonlyone!(opts, "--dry-run"        , "-n" )

    if opts["--veryveryverbose"]
        for (key,val) in opts @printf("  %-13s  =>  %-5s\n", key, repr(val)); end
        println("  fnames: $flist")
    end

    return opts, flist
end

# Auxiliary utility for testing purpose only.
function trytoinclude(mask::AbstractString = "*.jl")
    files = glob(mask)
    for f in files
        try
            include(f)
        catch e
            @error("in \"$f\":\n$e\n")
        end
    end
    return nothing
end

end # of module CLI

function convert_fortran(code; casetransform=identity, quiet=true,
                         verbose=false, veryverbose=false, veryveryverbose=false,
                         greeks=false, subscripts=0, greeksubscripts=false)

    # ORDER OF PROCESSING IS FRAGILE

    # Preprocess and cleanup fortran source

    # should be '\n' newlines only
    for rx in (r"\r\n" => @s_str("\n"), r"\r" => @s_str("\n"))
        code = replace(code, rx)
    end

    # expand tabs
    code = replace(code, r"^[ ]{0,5}\t"m => "      ")
    code = replace(code, r"\t"m => "  ")

    isfixedformfortran = !occursin(r"^[^cC*!#\n][^\n]{0,3}[^\d\h\n]"m, code)
    @debug "isfixedformfortran = $isfixedformfortran"

    # replace some symbols with the watermarks
    code = savespecialsymbols(code)

    # convert fortran comments to #
    code = convertfortrancomments(code, isfixedformfortran)

    commentstrings = OrderedDict{String,String}()
    code, commentstrings = savecomments(code, commentstrings, casetransform)

    # convert lines continuations marks to "\t\r\t"
    code = replacecontinuationmarks(code, isfixedformfortran)
    #write("test1.jl", code)

    code, formatstrings = collectformatstrings(code)
    veryveryverbose && printformatstrings(formatstrings)

    code, strings = savestrings(code, stringtype = "F8")

    # some simple fixes and conversions
    isfixedformfortran && (code = masklabels(code))
    write("test1.jl", code)
    code = foldl(replace, collect(repairreplacements), init=code)
    write("test2.jl", code)
    isfixedformfortran && (code = unmasklabels(code))

    # tokenize test
    errs = filter(x->x.token_error!=Tokens.TokenError(0), collect(tokenize(code)))
    length(errs) > 0 && @show errs

    # mark lines of code occupied by each module, subroutine and so on
    blocks = splitbyblocks(code, veryveryverbose)

    for i = 1:length(blocks)

        b = blocks[i]
        code = tostring(b[codelines])

        veryveryverbose && println()
        veryverbose && print("[$(b[startline]):$(b[lastline])] $(b[blocktype]): $(b[blockname])\n")

        # extract necessary information
        lines = splitonlines(concatcontinuedlines(stripcomments(code)))
        scalars, arrays, stringvars, vars = collectvars(lines)
        commons = collectcommon(lines)
        dolabels, gotolabels = collectlabels(lines)

        #veryveryverbose && length(gotolabels) > 0&& println("     labels : $gotolabels")
        #veryveryverbose && length(dolabels) > 0&& println("   dolabels : $dolabels")
        veryveryverbose && length(scalars) > 0 && println("    scalars : $scalars")
        veryveryverbose && length(arrays) > 0  && println("    arrays  : $arrays")
        veryveryverbose && length(commons) > 0 && println("    COMMONs : $commons")

        # replace array's parentheses with square brackets
        #write("test1.jl", code)
        code = replacearraysbrackets(code, arrays)
        #write("test2.jl", code)

        # READ and WRITE statements
        code = processiostatements(code, formatstrings)
        code = foldl(replace, reverse(collect(formatstrings)), init=code)

        #already done: code, strings = savestrings(code, stringtype = "F8")

        #write("test1.jl", code)
        code = insertabsentreturn(code)
        #write("test2.jl", code)

        code = processcommon(code, commons, arrays)

        code, commentstrings = commentoutdeclarations(code, commentstrings)

        code = processdostatements(code)

        #write("test1.jl", code)
        code = replacedocontinue(code, dolabels, gotolabels)
        #write("test2.jl", code)

        code = processconditionalgotos(code)

        #write("test1.jl", code)
        code = processifstatements(code)
        #write("test2.jl", code)
        code = processarithmif(code)
        #write("test2.jl", code)

        code = shapinglabels(code)

        # process replacements
        code = foldl(replace, collect(trivialreplacements), init=code)

        code = processimplieddoloops(code)

        code = processselectcase(code)

        code = processparameters(code)

        code = processdatastatement(code, vcat(arrays, "view"))

        code = splatprintviews(code)

        #write("test1.jl", code)
        code = processlinescontinuation(code)
        #write("test2.jl", code)

        lines, comments = splitoncomment(code)

        #write("test1.jl", tostring(lines))
        # remains straight syntax conversions
        for rx in replacements
            lines = map(a->replace(a, rx), lines)
        end
        #write("test2.jl", tostring(lines))

        # concat code lines back together and restore saved comments
        code = tostring(map(*, lines, comments))

        # strip all whitespaces inside brackets due to
        # https://github.com/JuliaLang/julia/issues/14853
        code = stripwhitesinbrackets(code)

        code = restorecomments(code, commentstrings)

        # removing unnecessaries
        for rx in removal
            code = replace(code, rx => "")
        end

        # create functions declaration with custom header
        for rx in headersprocessing
            code = replace(code, rx)
        end

        # restore saved formats, strings and symbols
        code = foldl(replace, reverse(collect(strings)), init=code)
        #code = foldl(replace, reverse(collect(formatstrings)), init=code)
        code = restorespecialsymbols(code)
        code = replace(code, r"\t?(\n|\r)\t|\t(\n|\r)\t?" => "\n")

        try
            Meta.parse(code, 1)
        catch e
            quiet || @error("$(b[1]) $(b[2])\n$e\n")
        end

        b[codelines] = code

    end # of loop over subroutines

    # transpiled code will be stored in the string `result`
    result = "using FortranFiles\nusing OffsetArrays\nusing Parameters\nusing Printf"
    #result = "using FortranFiles\nusing FortranStrings\nusing OffsetArrays\nusing Parameters\nusing Printf"

    # concat blocks of code and substitute code in place pointed by references
    i = 1
    while i <= length(blocks)
        b = blocks[i]
        result *= '\n' * b[codelines]
        for m in collect(eachmatch(r"[\n]?\h*#=(\w+):(\w+):(\d+):(\d+)=#\h*[\n]?"mi, result))
            k = findfirst(blocks) do a
                a[blocktype]      == m.captures[1] &&
                a[blockname]      == m.captures[2] &&
                "$(a[startline])" == m.captures[3] &&
                "$(a[lastline])"  == m.captures[4]
            end
            result = replace(result, m.match => '\n' * blocks[k][codelines] * '\n')
            splice!(blocks, k)
        end
        i += 1
    end

    return result
end

const SS = SubstitutionString

function mask(smbl)
    suffix = "iZjAcpPokM"
    if smbl == '@' || smbl == "@"
        return "at_$suffix"
    elseif smbl == '!' || smbl == "!"
        return "bang_$suffix"
    elseif smbl == '$' || smbl == "\$"
        return "dollar_$suffix"
    elseif smbl == '"' || smbl == "\""
        return "quote_$suffix"
    elseif smbl == '#' || smbl == "#"
        return "sha_$suffix"
    elseif smbl == '\\' || smbl == "\\"
        return "slashsymbol_$suffix"
    elseif smbl == "''"
        return "GhUtwoap_$suffix"
    elseif smbl == "end"
        return "lastpos_$suffix"
    elseif smbl == "return"
        return "ret_$suffix"
    else
        return smbl * "$suffix"
    end
end
mask() = mask("")

function savespecialsymbols(code::AbstractString)
    # save '#', '@', '!', "''"

    # save "''" inside strings 'SOME''STRING'
    # https://regex101.com/r/YrkC9C/1
    rx = r"('((?>[^']*|(?1))*)'){2,}"m
    #rx = r"('((?>[^'\n]*|(?1))*)'){2,}"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = m.match
        len = ncodeunits(code)
        while true
            str = replace(str, r"('((?>[^']*|(?1))*)')('((?>[^']*|(?1))*)')" => SS("'\\2$(mask("''"))\\4'"))
            len == (len = ncodeunits(str)) && break
        end
        code = replace(code, m.match => str)
    end

    # stitch F90 "&\n  &" continued lines marks in strings
    rx = r"('((?>[^']*|(?1))*)')"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"&\h*\n\h*&"m => "")
        code = replace(code, m.match => str)
    end

    # save '!', '"', '\\', '$' inside strings 'SOME!"\\$STRING'
    rx = r"('((?>[^'\n]*|(?1))*)')"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"!"    => mask('!')  )
        str = replace(str,     r"\""   => mask('"')  )
        str = replace(str,     r"[\\]" => mask('\\') )
        str = replace(str,     r"\$"   => mask('$')  )
        code = replace(code, m.match => str)
    end

    # replace '@' and '#' everywhere
    code  = replace(code, r"#" => mask('#'))
    code  = replace(code, r"@" => mask('@'))
    return code
end
function restorespecialsymbols(code::AbstractString)
    # restore saved symbols
    code = replace(code, mask('#')   => "#"    )
    code = replace(code, mask('@')   => "@"    )
    code = replace(code, mask('!')   => "!"    )
    code = replace(code, mask('$')   => "\\\$" )
    code = replace(code, mask('"')   => "\\\"" )
    code = replace(code, mask("''")  => "'"    )
    code = replace(code, mask("end") => "end"  )
    code = replace(code, mask('\\')  => "\\\\" )
    return code
end

function savecomments(code, commentstrings, casetransform=identity)

    rx = r"[ ]*#(?!CMMNT\d{10}).*(?:\t|)$" # spaces goes to comments
    lines = splitonlines(code)

    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            str = m.match
            key = @sprintf "CMMNT%05d%05d" i mod(hash(str), 2^16)
            commentstrings[key] = String(str)
            # case conversion while comments are detached
            lines[i] = casetransform(replace(lines[i], m.match => s"")) * " #$key"
        else
            # cut out all whitespaces not masked by the comments
            lines[i] = casetransform(replace(lines[i], r"[ ]*(\t|)$" => s"\1"))
        end
    end
    return sameasarg(code, lines), commentstrings
end
function restorecomments(code, commentstrings)
    for (k,v) in commentstrings
        code = replace(code, Regex("[ ]?#?$k") => v)
    end
    return code
end

"""
Convert lines continuations marks to "\t\r\t"
"""
function replacecontinuationmarks(code::AbstractString, fortranfixedform)

    code = replace(code, r"(?:&)([ ]*(?:#?CMMNT\d{10}|))(?:\n|\t\r\t)"m => SS("\\1\t\r\t"))
    code = replace(code, r"(?:\t\r\t|\n)([ ]*)&"m => SS("\t\r\t\\1 "))
    fortranfixedform && (code = replace(code, r"[\n\r][ ]{5}\S"m => SS("\t\r\t      ")))

    # also include empty and commented lines inside the block of continued lines in:
    # free-form
    rx = r"\t\r\t[ ]*(#?CMMNT\d{10}|)(\n[ ]*(#?CMMNT\d{10}|))*\n[ ]*(#?CMMNT\d{10}|)\t\r\t"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"\n" => SS("\t\r\t"))
        code = replace(code, m.match => str)
    end
    # fixed form
    rx = r"\n[ ]*(#?CMMNT\d{10}|)\t\r\t"m
    while true
        matches = collect(eachmatch(rx, code))
        length(matches) > 0 || break
        for m in reverse(matches)
            str = replace(m.match, r"\n" => SS("\t\r\t"))
            code = replace(code, m.match => str)
        end
    end
    return code
end

function collectformatstrings(code::AbstractString)
    formatstrings = OrderedDict{String,String}()

    # collect FORMAT(some,format,args)
    rx = r"^\h*\d+\h+format\h*\("mi
    rxfmt = r"^(\h*(\d+)\h+format\h*)\(([^\n]*)\)(\h*#?CMMNT\d+|\h*)$"mi
    for mx in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, mx.offset)
        str = stripcommentsbutlast(rstrip(concatcontinuedlines(code[r])))
        #str = replace(str, mask("''") => "''")
        m = match(rxfmt, str); fmt = m.captures[3]
        key = @sprintf "FMT%06d" mod(hash(fmt), 2^19)
        formatstrings[key] = fmt
        code = code[1:prevind(code,r[1])] * "$(m.captures[1])($key)$(m.captures[4])\n" *
               code[nextind(code,r[end]):end]
    end

    # collect remaining format strings: '(someformatstring)'
    # https://regex101.com/r/uiP6Is/5
    rx = r"('\(((?>(?!'\(|\)')*[^\n]|(?1))*?)\)')"mi
    #rx = r"('\(((?>[\s\S](?!'\(|\)')*|(?1))*?)\)')"mi
    for m in reverse(collect(eachmatch(rx, code)))
        str = String(m.match)
        #str = replace(m.match, mask("''") => "''")
        key = @sprintf "FMT%06d" mod(hash(str), 2^19)
        formatstrings[key] = str
        code = replace(code, m.match => key)
    end
    return code, formatstrings
end

function collectformats(lines::AbstractVector)
    # https://regex101.com/r/Lopg3O/1
    rx = r"^\h*(\d+)\h+format\h*(\(((?>[^\(\)]++|(?2))*)\))"i
    formats = Dict{String,String}()
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            formats[m.captures[1]] = m.captures[3]
        end
    end
    return formats
end

function printformatstrings(formatstrings)
    if length(formatstrings) > 0
        println("IO FORMAT strings:")
        foreach(reverse(collect(values(formatstrings)))) do s
            println(replace(restorespecialsymbols(s), r"\t\r?\t"=>"\n"))
        end
    end
    return nothing
end

function masklabels(code::AbstractString)
    # mask F77 labels: "1000" => "L1000:" to be like in fortran90
    # https://regex101.com/r/FG2iyI/4
    rx = r"(?=\n[ ]{0,4}\d[ ]{0,4})([\d ]{0,5})"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"(\d+)"m => s"L\1:")
        code = replace(code, m.match => str)
    end
    return code
end

unmasklabels(code::AbstractString) = replace(code, r"(\n\h*)L(\d+):(\h+)"m => s"\1\2\3")

"""
Simple names for Blocks structure fields indices
"""
@enum Blocks blocktype=1 blockname=2 startline=3 lastline=4 codelines=5 parent=6 firstincluded=7
Base.to_index(a::Blocks) = Int(a)
Base.promote_rule(T::Type, ::Type{Blocks}) = T
Base.convert(T::Type, a::Blocks) = T(Int(a))
for op in (:*, :/, :+, :-, :(:))
    @eval begin
        Base.$op(a::Number, b::Blocks) = $op(promote(a, b)...)
        Base.$op(a::Blocks, b::Number) = $op(promote(a, b)...)
        #Base.$op(a::Blocks, b::Blocks) = Blocks($op(Int(a), Int(b)))
    end
end

"""
Split code text into the blocks of code corresponding each module, subroutine and so on.
The resulted tree is flatted into Vector of the blocks and stored in the order of appearance.
Blocks of code: ["BlockType", "BlockName", firstline, lastline, ["Code Text"], Ref(ParentBlock)]
"""
function splitbyblocks(code, veryveryverbose=false)

    lines = splitonlines(code)
    lines .= rstrip.(stripcomments.(lines))

    # Initial blocks tree:
    # [blocktype, blockname, blockbegin, blockend, TEXT, parent,
    #     [blocktype, blockname, blockbegin, blockend, TEXT, parent, [...], [...], ], [...], ]
    sourcetree = ["FILE", "", 1, length(lines), String[], nothing]
    i = 1
    while (b = markbyblock(lines, i)) !== nothing
        push!(sourcetree, b)
        i = sourcetree[end][lastline] + 1
    end

    if veryveryverbose
        for b in blocksflattening(sourcetree)
            print("[$(b[startline]):$(b[lastline])] $(b[blocktype]): $(b[blockname])\n")
        end
    end

    sourcetree[codelines] = splitonlines(code)
    blocks = splitbyblock!(sourcetree)
    flatted = blocksflattening(blocks)

    if (i = findfirst(a -> a[blockname]==mask("NAMELESSPROGRAM"), flatted)) !== nothing
        #flatted[i][startline:lastline] .+= 1
        #prepend!(flatted[i][codelines], ("      PROGRAM NAMELESSPROGRAM",))
        #flatted[i][blockname] = "NAMELESSPROGRAM"
        flatted[i][codelines][1] = "      PROGRAM NAMELESSPROGRAM\n" * flatted[i][codelines][1]
    end

    return flatted
end

function markbyblock(lines, blockbegin)
    # split on: module | program | subroutine | function | blockdata

    # https://regex101.com/r/aWXhP4/5
    rxbegin = r"^(?:\h*|)((?:(?:\w+(?:\*\w+|)\h+|)(?:recursive\h+|)(?!\bend)function|(?:recursive\h+|)(?!\bend)subroutine|program|block\h*data|module(?!\h*procedure))\b)(?:\h+(\w+\b|)).*$"mi
    rxend = r"^(?:\h*|)(\bend(?:|function|subroutine|program|module)\b)(\h+\w+\b|).*$"mi

    # ["BlockType", "BlockName", firstline, lastline, String[], [], [...], ]
    blocks = Vector{Any}()
    headercatched = false
    lastnonempty = 1
    i = blockbegin
    while i <= length(lines)
        l, n = concatcontinuedlines(lines, i)
        i += n - 1 # skip multi-line
        if (m = match(rxbegin, l)) !== nothing
            if !headercatched
                blocktype = uppercase(replace(m.captures[1], r"^.*(\b\w+)$"=>s"\1"))
                blockname = uppercase(m.captures[2])
                append!(blocks, [blocktype, blockname, blockbegin, 0, String[], []])
                headercatched = true
            else
                nextblockbegin = min(i, lastnonempty+1)
                push!(blocks, markbyblock(lines, nextblockbegin))
                i = blocks[end][lastline]
            end

        elseif (m = match(rxend, l)) !== nothing
            if length(blocks) == 0
                @warn "Find only \"END\" without start of the program or subroutine"
                append!(blocks, ["PROGRAM", mask("NAMELESSPROGRAM"), blockbegin, 0, String[], []])
            end
            blocks[lastline] = i
            return blocks
        end

        isempty(lines[i]) || (lastnonempty = i)
        i += 1
    end

    if length(blocks) > 0
        @warn "Can't find final \"END\" of the $(blocks[blocktype]) $(blocks[blockname])"
        blocks[lastline] = length(lines)
        return blocks
    else
        return nothing
    end

end

function splitbyblock!(blocks)
    for i in length(blocks):-1:Int(firstincluded)
        b = blocks[i]
        ind = blocks[startline]
        b[codelines] = splice!(blocks[codelines], b[startline]-ind+1:b[lastline]-ind+1,
                       ["#=$(b[blocktype]):$(b[blockname]):$(b[startline]):$(b[lastline])=#"])
        b[parent] = Ref(blocks)
        blocks[i] = splitbyblock!(blocks[i])
    end
    return blocks
end

function blocksflattening(blocks)
    result = Vector{Any}()
    push!(result, blocks[1:firstincluded-1])
    for i in firstincluded:length(blocks)
        append!(result, blocksflattening(blocks[i]))
    end
    return result
end

function convertfortrancomments(code, isfixedformfortran = false)
    code1 = tostring(code)
    # '!' => '&' in 6-th col.
    isfixedformfortran && (code1 = replace(code1, r"(\n[ ]{5})!"m => s"\1&"))
    rx = r"(^[cC*]|!)([^\n]*)"m
    code1 = replace(code1, rx => s"#\2")
    code1 = replace(code1, r"#=" => "# ") # accident multiline comment
    return sameasarg(code, code1)
end

mutable struct Var
    decl::String
    name::String
    dim::String
    val::String
end
function Base.show(io::IO, v::Var)
    val = isempty(v.val) ? "" : "=$(v.val)"
    dim = isempty(v.dim) ? "" : "($(v.dim))"
    print(io, "$(v.decl)::$(v.name)$dim$val")
    return
end


function collectvars(lines::AbstractVector)

    lines1 = typeof(lines)()
    savedstrings = OrderedDict{String, String}()

    # save initial values strings
    for l in lines
        l1, savedstrings = savestrings(l, savedstrings)
        append!(lines1, (l1,))
    end

    declkeywords = ("DIMENSION", "INTEGER", "LOGICAL", "CHARACTER", "REAL", "COMPLEX",
                    "DOUBLEPRECISION", "INTENT", "TYPE", "EXTERNAL", "PARAMETER", )

    matched90 = Dict{Int,Any}()
    matched   = Dict{Int,Any}()
    for (i,l) in enumerate(lines1)
        ts = collect(tokenize(l))
        ts = filter(a->a.kind!=Tokens.WHITESPACE, ts)
        if any(a->a.kind == Tokens.DECLARATION, ts)
            matched90[i] = ts
        elseif ts[1].kind == Tokens.IDENTIFIER && ts[1].val in declkeywords
            matched[i] = ts
        end
    end

    vars    = Dict{String, Var}()
    #vars    = Dict{String, NamedTuple{(:decl,:name,:dim,:val),Tuple{String,String,String,String}}}()
    for (i,ts) in matched
        vars = parsevardeclstatement(lines1[i], ts, vars)
    end
    for (i,ts) in matched90
        vars = parsevardecl90statement(lines1[i], ts, vars)
    end
    #@debug for (k,v) in vars println("$k: \"$v\"") end
    #@debug "vars = \"$vars\""

    arrays = Vector{String}()
    scalars = Vector{String}()
    strings = Vector{String}()
    # catch CHARACTER???? https://regex101.com/r/1fZkiz/3
    rxch = r"^CHARACTER(?:\*?(\(((?>[^()]++|(?1))*)\))|\*(\d+)|(?=\h))"i

    # separate vars on types
    for (k,v) in vars
        occursin(r"\bEXTERNAL\b", v.decl) && (pop!(vars, k); continue) # exclude EXTERNAL
        length(v.dim) == 0 ? push!(scalars, k) : push!(arrays, k)
        occursin(r"\bCHARACTER\b", v.decl) && push!(strings, k)
        #occursin(rxch, v.decl) && push!(strings, k)
    end
    # in FORTRAN even single CHARACTER can be accessed by index:
    # `CH(:1)` or `CH(1:1)` thus any CHARACTER variable is an array
    append!(arrays, strings)

    @debug "scalars = \"$scalars\""
    @debug "arrays  = \"$arrays\""
    @debug "strings = \"$strings\""
    @debug "vars = \"$vars\""

    # restore initial strings values
    for s in strings
        vars[s].val = restorestrings(vars[s].val, savedstrings)
    end

    return scalars, arrays, strings, vars
end


function parsevardeclstatement(line, ts, vars)
    # parse strings like this:
    # INTEGER X(10),Y(*),Z(2,*),I
    # CHARACTER*(*) S(2)
    # CHARACTER*8 S/'12345678'/,R
    # CHARACTER JBCMPZ*2/'AB'/,R2*(*)

    startvarsind = 1

    @debug "parcedline = \"$line\""
    @debug for (i,z) in enumerate(ts) print("$i: "); show(z); println() end

    # split onto decl and vars
    if ts[1].val == "PARAMETER"
        @assert ts[2].kind == Tokens.LPAREN && ts[end-1].kind == Tokens.RPAREN
        splice!(ts, 2:2)
        splice!(ts, length(ts)-1:length(ts)-1)
        startvarsind = 2
    else
        for i = 1:length(ts)-1
            if ts[i+1].kind == Tokens.IDENTIFIER &&
                ts[i].kind in (Tokens.IDENTIFIER, Tokens.INTEGER, Tokens.RPAREN)
                startvarsind = i + 1; break
            end
        end
    end

    decltype0 = join(untokenize.(ts[1:startvarsind-1]))

    i = startvarsind
    while i <= length(ts)
        tt = ts[i]
        #@debug "tt.kind=\"$(tt.kind)\""
        if tt.kind == Tokens.ENDMARKER || tt.kind == Tokens.COMMENT
            break
        elseif tt.kind == Tokens.COMMA || tt.kind == Tokens.WHITESPACE
            i += 1
        elseif tt.kind == Tokens.IDENTIFIER
            varname = tt.val
            varname == "FUNCTION" && break
            t = ts[i+=1]
            if t.kind == Tokens.STAR #'*' # CHARACTER S*2
                any(a->a.val == "CHARACTER", ts[1:startvarsind-1]) || @error "unexpected token '*' in \"$(restore(ts))\""
                t = ts[i+=1]
                if t.kind == Tokens.LPAREN # '(' # CHARACTER S*(SOMEEXPR)
                    i, i0 = findmatchedbrace(ts, i), i
                    decltype = "CHARACTER" * restore(ts[i0:i]) # type overwritten
                else
                    t.kind == Tokens.INTEGER || @error "unexpected token \"$(t.val)\" in \"$(restore(ts))\""
                    decltype = "CHARACTER*$(t.val)"
                end
                t = ts[i+=1]
            else
                decltype = decltype0
            end
            if t.kind == Tokens.LPAREN #'(' # ARRAY(SOMEEXPR)
                i, i0 = skipbraces(ts, i), i
                dim = restore(ts[i0+1:i-2])
                t = ts[i]
            else
                dim = ""
            end
            if t.kind == Tokens.FWD_SLASH # '/' # initial value /42/
                i, i0 = findnext(a->a.kind==Tokens.FWD_SLASH, ts, i+1), i
                val = restore(ts[i0+1:i-1])
                i += 1
            elseif t.kind == Tokens.EQ # '=' # PARAMETER value
                i += 1; i0 = i
                # catch ',' or EOL
                i = something(findnexttoken(Tokens.COMMA, ts, i), length(ts)+1)
                val = restore(ts[i0:i-1])
            else
                val = ""
            end
            k = uppercase(varname)
            if haskey(vars, k)
                vars[k].decl *= (isempty(vars[k].decl) ? "" : ",") * decltype
                isempty(vars[k].dim) && (vars[k].dim = dim)
                isempty(vars[k].val) && (vars[k].val = val)
            else
                push!(vars, k => Var(decltype,varname,dim,val))
            end
        else
            @error "unexpected kind \"$(tt.kind)\" of token \"$tt\""
            exit(1)
        end
    end

    return vars
end


function parsevardecl90statement(line, ts, vars)
    # parse strings like this:
    # character*(*),parameter::NAME='FUNNAM'
    # Integer,Dimension(3)::K(2)=(/12,13/),I
    # Real(kind=kdp),Dimension(:),Intent(InOut)::XDONT

    @debug "parcedline = \"$line\""
    @debug for (i,z) in enumerate(ts) print("$i: "); show(z); println() end

    startvarsind = findfirst(t->t.kind==Tokens.DECLARATION, ts) + 1 # "::"
    decltype = restore(@view ts[1:startvarsind-2])
    i = findfirst(t->t.val=="DIMENSION", @view ts[1:startvarsind-2])
    dim0 = i!==nothing && ts[i+1].kind == LPAREN ? restore(@view ts[i+1:findmatchedbrace(ts,i)-1]) : ""
    # catch braces and its contents https://regex101.com/r/inyyeW/2
    #rxdim = r"DIMENSION(\(((?>[^()]++|(?1))*)\))"i
    #dim0 = (m = match(rxdim, decltype)) !== nothing ? m.captures[2] : ""

    i = startvarsind
    while i <= length(ts)
        tt = ts[i]
        #@debug "tt.kind=\"$(tt.kind)\""
        if tt.kind == Tokens.ENDMARKER || tt.kind == Tokens.COMMENT
            break
        elseif tt.kind == Tokens.COMMA || tt.kind == Tokens.WHITESPACE
            i += 1
        elseif tt.kind == Tokens.PAIR_ARROW # "=>"
            break
        elseif tt.kind == Tokens.IDENTIFIER #
            varname = tt.val
            t = ts[i+=1]
            if t.kind == Tokens.LPAREN #'(' # arrays
                # TODO?: catch character(len=10)
                i, i0 = findmatchedbrace(ts, i), i
                dim = restore(ts[i0+1:i-1])
                t = ts[i+=1]
            elseif dim0 != ""
                dim = dim0
            else
                dim = ""
            end
            if t.kind == Tokens.EQ # '=' # initial value
                i += 1; i0 = i
                # catch ',' or EOL
                i = something(findnexttoken(Tokens.COMMA, ts, i), length(ts)+1)
                val = restore(ts[i0:i-1])
            else
                val = ""
            end
            k = uppercase(varname)
            if haskey(vars, k)
                vars[k].decl *= (isempty(vars[k].decl) ? "" : ",") * decltype
                #isempty(vars[k].dim) && (vars[k].dim = dim)
                !isempty(dim) && (vars[k].dim = dim)
                isempty(vars[k].val) && (vars[k].val = val)
            else
                push!(vars, k => Var(decltype,varname,dim,val))
            end
        else
            @error "unexpected kind \"$(tt.kind)\" of token \"$tt\""
            exit(1)
        end
    end

    return vars
end

function collectcommon(lines::AbstractVector)
    rx = r"^\h*common\h*\/\h*(\w+)\h*\/\h*(.+)"mi
    matched = Dict{String,String}()
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            matched[m.captures[1]] = m.captures[2]
        end
    end
    rxbr = r"(\(((?>[^()]++|(?1))*)\))" # https://regex101.com/r/inyyeW/2
    rxsqbr = r"(\[((?>[^\[\]]++|(?1))*)\])"
    for i in keys(matched)
        matched[i] = replace(matched[i], r"[ \t\n\r]" => s"")   # drop spaces
        matched[i] = replace(matched[i], rxsqbr => s"")         # drop braces
        matched[i] = replace(matched[i], rxbr => s"")           # drop braces
    end
    return matched
end

function collectlabels(lines::AbstractVector)

    # capture 'DO's and 'GOTO's labels
    dolabels = Accumulator{String,Int}()
    gotolabels = Set{String}()
    rx = r"^(?:\h*\d+:?\h*|\h*)(?:do|for)\h*(\d+)(?:\h*,|)\h*\w+\h*="i
    #rx = r"^\h*(?:do|for)\h+(\d+)\h+"i
    for i in axes(lines,1)
        (m = match(rx, lines[i])) !== nothing && push!(dolabels, m.captures[1])
    end
    rx = r"(?:^\h*|\W)go\h*?to\h+(\d+)$"i
    for i in axes(lines,1)
        (m = match(rx, lines[i])) !== nothing && push!(gotolabels, m.captures[1])
    end

    # GO TO( 20, 40, 70, 110, 140 )ISAVE( 1 )
    # https://regex101.com/r/ba4wNU/2
    rx = r"(?:^\h*|\W)go\h*?to\h*\((\h*\d+\h*(?:\h*,\h*\d+\h*)*)\)\h*.*$"i
    #rx = r"(?:^\h*|\W)go\h*?to\h*\((\h*\d+\h*(?:\h*,(?:\n[ ]{5}[^ ])?\h*\d+\h*)*)\)\h*.*$"i
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            labels = map(strip, split(m.captures[1], ','))
            foreach(l->push!(gotolabels, l), labels)
        end
    end

    # GOTO EXPR,(10,20,30,40)
    # https://regex101.com/r/iYG7BH/6
    rx = r"^(\h*\d+:?\h*|\h*)(.*?|)(\h*|)go\h*?to\h+(\w+)\h*(?:,\h*|)\((\h*\d+\h*(?:,\h*\d+\h*)*)\)(.*)$"mi
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            labels = map(strip, split(m.captures[5], ','))
            foreach(l->push!(gotolabels, l), labels)
        end
    end

    # 113 IF (IERROR) 119,114,119
    # https://regex101.com/r/5R0qtY/1/
    rx = r"\bif\b\h*(\(((?>[^\(\)]++|(?1))*)\))\h*(\d+\h*,\h*\d+\h*,\h*\d+)"mi
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) !== nothing
            labels = map(strip, split(m.captures[3], ','))
            foreach(l->push!(gotolabels, l), labels)
        end
    end

    @debug "dolabels: \"$dolabels\""
    @debug "gotolabels: \"$gotolabels\""
    return dolabels, gotolabels
end

"""
replace array's braces with square brackets
"""
function replacearraysbrackets(code::AbstractString, arrays)

    # remove whitespace between array name and left brace
    for a in arrays
        rx = Regex("(\\b$(a)\\b)\\h+\\(","i")
        code = replace(code, rx => s"\1(")
    end

    # lookup "ARRAYNAME(" and replace
    braces = r"\(((?>[^()]|(?R))*)\)" # complementary braces
    brackets = s"[\1]"
    for a in arrays
        rx = Regex("\\b$(a)\\b\\(","i")
        for m in collect(eachmatch(rx, code))
            o = m.offset
            code = code[1:prevind(code,o+ncodeunits(a))] *
                   replace(code[thisind(code,o+ncodeunits(a)):end], braces => brackets, count=1)
            #code = code[1:prevind(code,o)] * replace(code[o:end], braces => brackets, count=1)
        end
    end

    # catch second braces, for arrays of strings, https://regex101.com/r/2WIc50/1
    rx = r"(\[((?>[^\[\]]++|(?1))*)\])\h*(\(((?>[^\(\)]++|(?1))*)\))"
    code = replace(code, rx => s"\1[\4]")

    # fix uncompleted ranges in square braces like [1:] and [:1] and so on
    brackets = r"\[((?>[^\[\]]++|(?0))*)\]" # complementary brackets

    # trim line breaks after ':'
    rx = r":\s*(#CMMNT\d{10}|)\t\r\t\h*"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s":") * code[o+length(m.match):end]
        end
    end

    # trim spaces around ':'
    rx = r"\[(.*|)([ ]+|):([ ]+|)(.*|)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s"[\1:\4]") * code[o+length(m.match):end]
        end
    end

    # uppend with "end": `A[2:]` => `A[2:end]`
    rx = r"\[(.*[^,]):(,.*|)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => SS("[\\1:$(mask("end"))\\2]")) *
                   code[o+length(m.match):end]
        end
    end

    # prepend with '1': `A[:2]` => `A[1:2]`
    rx = r"\[:([^,].*)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s"[1:\1]") * code[o+length(m.match):end]
        end
    end

    return code
end

function stripwhitesinbrackets(code::AbstractString)
    brackets = r"\[((?>[^\[\]]++|(?0))*)\]" # complementary brackets
    # trim spaces inside
    for m in reverse(collect(eachmatch(brackets, code)))
        o = m.offset
        code = code[1:prevind(code,o)] *
               replace(m.match, " " => "") *
               code[thisind(code,o+ncodeunits(m.match)):end]
    end
    return code
end

function convertstrings(code::AbstractString)
    # `S = 'A'`       => `S[1] = 'A'; S[2:end] .= 0`
    # `S[1:2] = 'A'`  => `S[1] = 'A'; S[2] = ' '`
    # `S = "ABC"`     => `S .= collect("ABC")[1:min(length(S),length("ABC"))]`
    # ###`S[2] = "ABC"`  => `S[2] = collect("ABC")[1]`
    # `S[:2] = "ABC"` => `S[1:2] .= collect("ABC")[1:2]`
    # `S = "ABC"`     => `S .= collect("ABC")`
    # `S = 'A'`       => `S .= 'A'`
    # `S = 'A'`       => `S .= 'A'`
    #

    # convert `A = "ABC"` into `A .= collect("ABC")`
    # capture string https://regex101.com/r/KTFUCj/3
    rx = r"(\w+\s*)=(\s*)(\"(?:[^\"\\]|\\.)*\")"
    for m in reverse(collect(eachmatch(rx, code)))
        o = m.offset
        code = code[1:prevind(code,o)] *
               replace(m.match, rx => SS("\1.=\2collect(\3)")) *
               code[thisind(code,o+ncodeunits(m.match)):end]
    end

    # convert `A[2] = "ABC"` into `A[2] = collect("ABC")[1]`
    # capture string https://regex101.com/r/KTFUCj/3
    rx = r"(\w+\[\w+\]\s*)=(\s*)(\"(?:[^\"\\]|\\.)*\")"
    for m in reverse(collect(eachmatch(rx, code)))
        o = m.offset
        code = code[1:prevind(code,o)] *
               replace(m.match, rx => SS("\1=\2collect(\3)[1]")) *
               code[thisind(code,o+ncodeunits(m.match)):end]
    end

    return code
end

"""
save and replace all strings with its hash
"""
function savestrings(code::AbstractString, strings = OrderedDict{String, String}(); stringtype = "")

    # replace 'somestring' with STRINGTYPEPREFIX"somestring" in code
    # https://regex101.com/r/q4fzK9/1
    rxstr = r"('((?>[^'\n]*|(?1))*)')"m
    for m in reverse(collect(eachmatch(rxstr, code)))
        key = mask(@sprintf "STR%s" mod(hash("$(m.captures[1])"), 2^20))
        if length(m.captures[1]) > 3 || length(m.captures[1]) == 2
            strings[key] = stringtype * '"' * m.captures[2] * '"'
            code = replace(code, m.match => key)
        else # 'X' -- some single symbol, but in F77 it also string
            #strings[key] = string(m.captures[1])
            strings[key] = stringtype * '"' * m.captures[2] * '"'
            code = replace(code, m.match => key)
        end
    end
    rxstr = r"(\"((?>[^\"\n]*|(?1))*)\")"m
    for m in reverse(collect(eachmatch(rxstr, code)))
        key = mask(@sprintf "STR%s" mod(hash("$(m.captures[1])"), 2^20))
        strings[key] = stringtype * string(m.captures[1])
        code = replace(code, m.match => key)
    end

    # save empty strings ''
    key = mask(@sprintf "STR%s" mod(hash("\"\""), 2^20))
    strings[key] = stringtype * "\"\""
    code = replace(code, mask("''") => key)

    return code, strings
end

function restorestrings(code, strings)
    return foldl(replace, reverse(collect(strings)), init=code)
end

function processiostatements(code::AbstractString, formatstrings)

    # collect LABEL FORMAT(FMT897348)
    rx = r"^\h*(\d+)\h+format\h*\((FMT\d{6})\)"mi
    formats = Dict{String,String}()
    for m in collect(eachmatch(rx, code))
        formats[m.captures[1]] = m.captures[2]
    end

    # replace READ/WRITE with READ/println/@printf
    rxhead = r"(\bread\b|\bwrite\b)\h*\("mi
    rxargs = r"^((?:.*[\n])*.*?)(\h*#[^\n]*)"m  # https://regex101.com/r/PqvHQD/1
    for (i,m) in enumerate(reverse(collect(eachmatch(rxhead, code))))
        o = m.offset
        readwrite = lowercase(replace(m.match, rxhead => s"\1"))
        iolength, io, label, fmt, pmt, paramsstr, args = parsereadwrite(code, o)
        #@show iolength, io, label, fmt, pmt, paramsstr, args
        if length(label) > 0 && haskey(formats, label)
            fmt = formatstrings[formats[label]]
            #fmt = replace(formatstrings[formats[label]], mask("''") => "'")
            #formats[label] = fmt
        elseif length(label) > 0 && haskey(formatstrings, label)
            fmt = replace(formatstrings[label], r"^'\((.*)\)'$" => s"\1")
            #fmt = replace(replace(formatstrings[label], r"^'\((.*)\)'$" => s"\1"), mask("''")=>"'")
        elseif length(label) > 0
            @warn("IO format string in variable $label")
        else
            formats[@sprintf "FMT%06d" mod(hash(fmt), 2^19)] = fmt
        end
        fmt = convertformat(fmt)
        if  readwrite == "read" && (io == "5" || io == "*")
            io = "stdin"
        elseif  readwrite == "write" && (io == "6" || io == "*")
            io = "stdout"
        end
        args = occursin(rxargs, args) ? replace(args, rxargs=>s", \1)\2") : ", $(strip(args)))"
        str  =  readwrite == "read" ? "READ(" :
                fmt       == ""    ? "println(" : "$(mask('@'))printf("
        if str == "$(mask('@'))printf("
            fmt = occursin(r"\$$", fmt) ? replace(fmt, @r_str("\\\$\$") => "") : fmt*"\\n"
        end

        #length(fmt) > 0 && @show fmt
        #fmt = fmt != "" ? ", \"$fmt\"" : ", $label"
        fmt = fmt != "" ? ", \"$fmt\"" :
                          label != "" ? ", $label" : ""
        str *= io * fmt * args
        code = code[1:prevind(code,o)] * str * code[thisind(code,o+iolength):end]
    end

    return code
end

function parsereadwrite(str, pos)
    pos0 = pos
    label = fmt = IU = formatvar = ""
    params = Dict{String,String}()
    cmdtaken = false
    inparams = false
    nextparam = 1
    startsparam = Int[]
    endsparam = Int[]
    startargs = 0
    endargs = -1
    inbraces  = 0

    while true
        ttype = picktoken(str, pos)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            pos = skipspaces(str, pos)
        elseif !cmdtaken && ttype == 'l'
            lex, pos = taketoken(str, pos)
            if lowercase(lex) == "read" || lowercase(lex) == "write"
                cmdtaken = true
                pos = skipspaces(str, pos)
            else
                # it is not an 'READ/WRITE' statement
                return 0, "", "", "", "", "", "", ""
            end
        elseif cmdtaken && inparams && ttype == 'l'
            lex, pos = taketoken(str, pos)
            LEX = uppercase(lex)
            if LEX == "FMT" && picktoken(str, skipspaces(str, pos)) == '='
                label, pos = taketoken(str, skipspaces(str, skiptoken(str, skipspaces(str, pos))))
            elseif (LEX == "ERR" || LEX == "END" || LEX == "REC") &&
                   picktoken(str, skipspaces(str, pos)) == '='
                val, pos = taketokenexpr(str, skipspaces(str, skiptoken(str, skipspaces(str, pos))))
                #val, pos = taketoken(str, skipspaces(str, skiptoken(str, skipspaces(str, pos))))
                params[uppercase(lex)] = val
            elseif nextparam == 1
                IU = lex
            else
                # this is the format string in the variable
                formatvar = lex
                #occursin(r"^FMT\d{6}$", lex) ||
                #@warn("parsereadwrite(): $(@__LINE__): unknown FORMAT-parameter = \"$lex\"" *
                #      " in FORTRAN line:\n\"$(strip(str[thislinerange(str, pos)]))\"\n")
            end
        elseif cmdtaken && inparams && ttype == 'd' && nextparam == 1
            IU, pos = taketoken(str, pos)
        elseif cmdtaken && inparams && ttype == 'd'
            label, pos = taketoken(str, pos)
        elseif cmdtaken && inparams && ttype == '*' && nextparam == 1
            IU, pos = taketoken(str, pos)
        elseif cmdtaken && inparams && ttype == ''' || ttype == '*'
            fmt, pos = taketoken(str, pos)
            fmt = concatcontinuedlines(fmt)
            fmt = replace(fmt, mask("''") => "'")
        elseif ttype == '('
            inbraces += 1
            pos = skiptoken(str, pos)
            if length(startsparam) == 0
                push!(startsparam, pos)
                inparams = true
            end
        elseif ttype == ')'
            inbraces -= 1
            if inbraces == 0 && inparams
                push!(endsparam, prevind(str, pos))
                startargs = nextind(str, pos)
                inparams = false
                endargs = pos = skipupto('\n', str, pos)
                endargs = prevind(str, endargs)
                break
            end
            pos = skiptoken(str, pos)
        elseif ttype == ',' && inparams && inbraces == 1
            push!(endsparam, prevind(str, pos))
            pos = skiptoken(str, pos)
            push!(startsparam, thisind(str, pos))
            nextparam += 1
        else
            pos = skiptoken(str, pos)
        end
    end

    IU = strip(str[startsparam[1]:endsparam[1]])

    label == "" && (label = formatvar)

    if fmt == "*"
        fmt = ""
    else
        fmt = replace(fmt, r"^\((.*)\)$"=>s"\1") # it is unenclose "()"
    end
    fmt = foldl((a,b) -> a*','*b, map(strip, split(fmt, ",")))

    return ncodeunits(str[pos0:prevind(str,pos)]), IU, label, fmt, params,
           str[startsparam[1]:endsparam[end]], str[startargs:endargs]
end

function splitformat(str)
    pos = 1
    lexemstarts = Int[1]
    lexemends = Int[]
    inbraces  = 0

    while true
        ttype = picktoken(str, pos)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            pos = skipwhitespaces(str, pos)
        elseif ttype == ',' && inbraces == 0
            push!(lexemends, prevind(str, pos))
            pos = skiptoken(str, pos)
            push!(lexemstarts, thisind(str, pos))
        elseif ttype == '/' && inbraces == 0
            if pos > lexemstarts[end]
                push!(lexemends, prevind(str, pos))
                push!(lexemstarts, thisind(str, pos))
            end
            push!(lexemends, thisind(str, pos))
            pos = skiptoken(str, pos)
            push!(lexemstarts, thisind(str, pos))
        elseif ttype == '('
            inbraces += 1
            pos = skiptoken(str, pos)
        elseif ttype == ')'
            inbraces -= 1
            pos = skiptoken(str, pos)
        else
            pos = skiptoken(str, pos)
        end
    end
    push!(lexemends, ncodeunits(str))
    @assert length(lexemstarts) == length(lexemends)

    fmt = String[]
    for (i,j) in zip(lexemstarts, lexemends)
        push!(fmt, strip(str[i:j]))
    end

    return fmt
end

function convertformat(formatstring)
    FMT = OrderedDict(
        r"^\/$"                => "\\n",             # / -> START NEW RECORD
        r"^(\d*)P$"i           => "",                # Scale Factor (P): scale from/to mantissa
        r"^A$"i                => s"%s",             # A -> %s
        r"^A(\d*)$"i           => s"%\1s",           # A9 -> %9s
        r"^I(\d*)$"i           => s"%\1i",           # I5 -> %5i
        r"^I(\d+)\.(\d+)$"i    => s"%0\1i",          # I5.5 -> %05i there can be mistake if \2>\1
        r"^[ED](\d*\.\d*)$"i   => s"%\1E",           # E7.2 -> %7.2E
        r"^[ED]S(\d*\.\d*)$"i  => s"%\1E",           # ES7.2 -> %7.2E ???
        r"^P[ED](\d*\.\d*)$"i  => s"%\1E",           # E7.2 -> %7.2E Scientific format with Scale Factor P
        r"^P1[ED](\d*\.\d*)$"i => s"%\1E",           # E7.2 -> %7.2E
        r"^P2[ED](\d*\.\d*)$"i => s"%\1E%\1E",       # E7.2 -> %7.2E ??? is it right?
        r"^P3[ED](\d*\.\d*)$"i => s"%\1E%\1E%\1E",   # E7.2 -> %7.2E
        r"^F(\d*\.\d*)$"i      => s"%\1F",           # F7.2 -> %7.2F
        r"^X$"i                => s" ",              # X -> ' '
        r"^T\d*$"i             => s" ",              # Tx -> ' ' : move to absolute position (column) x
        r"^([-]?\d+)P$"i       => s"",               # 1P -> '' https://docs.oracle.com/cd/E19957-01/805-4939/z4000743a6e2/index.html
        r"^([^']*)'(.*)'([^']*)$" => s"\1\2\3",      # unenclose ''
    )
    # split on tokens and apply repeats
    format = parseformat(formatstring)
    # converting
    for rs in FMT
        format = map( a->replace(a, rs), format)
    end
    return foldl(*, format)
end

"""
split format string on tokens and apply repeats
"""
function parseformat(formatstring)
    #@show formatstring
    rep = r"^(\d+)(.*)$"
    formatparts = splitformat(formatstring)
    #@show formatparts

    format = Vector{String}(undef, 0)
    for (i,el) in enumerate(formatparts)
        if occursin(rep, el)
            num = parse(Int, replace(el, rep=>s"\1"))
            str = replace(el, rep => s"\2")
            if occursin(r"\(.*\)", str)
                str = replace(str,  r"\((.*)\)"=> s"\1")
                fmt = parseformat(str)
                [append!(format, fmt) for j=1:num]
            else
                [push!(format, str) for j=1:num]
            end
        else
            push!(format, el)
        end
    end
    return format
end

function insertabsentreturn(code::AbstractString)
    # find last return if exist. https://regex101.com/r/evx2Lu/13
    # else insert new one
    # Note: make a possessive regex for "ERROR: LoadError: PCRE.exec error: match limit exceeded"
    rx = r"((?<=\n|\r)\h*\d+\h+|(?<=\n|\r)\h*)(return)((?:\h*#?CMMNT\d+\s*|\s*)++)(\h*end\h*(?:function|(?:recursive\h+|)subroutine|program|module|block|)(?:\h*#?CMMNT\d+\s*|\s*))$"mi
    if (m = match(rx, code)) !== nothing
        code = replace(code, rx => SS("\\1$(mask("lastreturn"))\\3\\4"))
    elseif occursin(r"\bend[\h\r\n]*$"mi, code)
        m = collect(eachmatch(r"(?<=\n|\r)(\h+|)end"mi, code))[end]
        o = m.offset; l = ncodeunits(m.captures[1])
        code = code[1:prevind(code,o)] * "$(' '^l)$(mask("lastreturn"))\n" *
                code[thisind(code,o):end]
    end
    # mask all other `return`
    code = replace(code, r"\breturn\b"i => SS("$(mask("return"))"))
    return code
end

function processcommon(code::AbstractString, commons, arrays)
    if length(commons) == 0
        #code = replace(code, Regex("\\b$(mask("lastreturn"))\\b") => SS("$(mask("return"))"))
        return code
    end
    #lines = splitonlines(code)

    # @unpack COMMONs before use
    rx = r"^([ ]*)common\h*\/\h*(\w+)\h*\/\h*"mi
    for m in reverse(collect(eachmatch(rx, code)))
        if haskey(commons, m.captures[2])
            c1 = m.captures[1]; c2 = m.captures[2]
            vars = commons[c2]
            s = SS("$(c1)global $c2\n$c1$(mask('@'))unpack $(vars) = $c2\n$(m.match)")
            code = replace(code, m.match => s)
        end
    end

    # @pack! 'COMMON's back before return
    packstr = "      $(mask('@'))label Lreturn\n"
    for k in keys(commons)
        # arrays don't needed because they should not be reallocated
        v = split(commons[k], ',')
        v = v[findall(a->lowercase(a) ∉ map(lowercase,arrays), v)]
        if length(v) > 0
            v = foldl((a,b)->a*','*' '*b, v)
            packstr *= "      $(mask('@'))pack$(mask('!')) $(k) = $(v)\n"
        end
    end

    # insert `@pack!` line before final return
    rx = Regex("^(\\h*\\d+\\h+|\\h+)$(mask("lastreturn"))", "mi")
    m = match(rx, code)
    if isempty(strip(m.captures[1])) # without leading label
        code = replace(code, m.match => SS("$(m.captures[1])\n$packstr      $(mask("lastreturn"))"))
    else
        label = strip(m.captures[1])
        l = length(m.captures[1])
        code = replace(code, m.match =>
                        SS(' '^l * "$(mask('@'))label L$(label)\n" *
                           ' '^l * "$(lstrip(packstr))      $(mask("lastreturn"))"))
    end

    # replace each RETURN with "@goto Lreturn" and restore last return
    code = replace(code, Regex("\\b$(mask("return"))\\b") => SS("$(mask('@'))goto Lreturn"))
    code = replace(code, Regex("\\b$(mask("lastreturn"))\\b") => SS("$(mask("return"))"))

    return code
end

function commentoutdeclarations(code, commentstrings)

    patterns = [
        r"^\h*common\h*\/\h*\w+\h*\/\h*"mi,
        r"^\h*dimension\h+"mi,
        r"^\h*implicit\h*\w+"mi,
        r"^\h*real(?!.*parameter).*::"mi, # https://regex101.com/r/7YxbnB/1
        r"^\h*real\h*$"mi,
        r"^\h*real\h+(?!.*function)"mi, # https://regex101.com/r/CFiiOx/1
        r"^\h*real\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*double\h*precision(?!.*parameter).*::"mi,
        r"^\h*double\h*precision\h*$"mi,
        r"^\h*double\h*precision\h+(?!.*function)"mi,
        r"^\h*complex(?!.*parameter).*::"mi,
        r"^\h*complex\h*$"mi,
        r"^\h*complex\h+(?!.*function)"mi,
        r"^\h*complex\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer(?!.*parameter).*::"mi,
        r"^\h*integer\h*$"mi,
        r"^\h*integer\h+(?!.*function)"mi,
        r"^\h*integer\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer\h*\(\h*KIND\h*=\h*\w+\h*\)\h*(?!.*function)"mi,
        r"^\h*character(?!.*parameter).*::"mi,
        r"^\h*character\h*$"mi,
        r"^\h*character\h+(?!.*function)"mi,
        r"^\h*character\h*\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*character\h*\(\h*\*\h*\)\h*(?!.*function)"mi,
        r"^\h*character\h*\(\h*\d+\h*\)\h*(?!.*function)"mi,
        r"^\h*character\h*\*?\h*\(\h*\*\h*\)\h*(?!.*function)"mi,
        r"^\h*character\h*\*?\h*\(\h*\d+\h*\)\h*(?!.*function)"mi,
        r"^\h*character\h*\*?\h*\([\h\w_*+-]+\)\h*(?!.*function)"mi,
        r"^\h*external\h+"mi,
        r"^\h*logical(?!.*parameter).*::"mi,
        r"^\h*logical\h*$"mi,
        r"^\h*logical\h+(?!.*function)"mi,
        r"^\h*allocatable\h+"mi,
        r"^\h*entry\h+.*"mi,
        r"^\h*save\h+.*"mi,
        r"^\h*equivalence\h+.*"mi,
        r"^\h+\d+\h+format\h*\("mi,
        # Import statements
        r"^\h*use\h.*"mi,
    ]

    for rx in patterns
        for m in reverse(collect(eachmatch(rx, code)))
            r = continuedlinesrange(code, m.offset)
            str = '#' * replace(code[r], r"(\r\t|\n\t)"m => s"\1#")
            code = code[1:prevind(code,r[1])] * str * code[nextind(code,r[end]):end]
        end
    end

    return savecomments(code, commentstrings, identity)
end

function processdostatements(code::AbstractString)
    rx = r"^(\h*\d*\h*|\h*)do\h*(?:\d+\h*,?|)\h*[\w_]+\h*="mi
    #rx = r"^(\h*\d*\h*|)do\h*(?:\d+\h*,?|)\h*[\w_]+\h*=.*$"i
    i = 1
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        success,s1,lbl,var,ex1,ex2,ex3,com = parsedostatement(code[r])
        #@show success,s1,lbl,var,ex1,ex2,ex3,com
        if success
            length(lbl) > 0 && (lbl *= " ")
            length(ex3) > 0 && (ex3 *= ":")
            str = "$(m.captures[1])for $lbl$var = $ex1:$ex3$ex2$com\n"
            #str = m.captures[1] * "for " * lbl * var * " = " * ex1 * ":" * ex3 * ex2 * com * '\n'
            code = code[1:prevind(code,r[1])] * str * code[nextind(code,r[end]):end]
        end
    end
    return code
end

function parsedostatement(line)
    # "DO10,I=10,A(1,1)*TONUM('10,000'),-2" is the valid fortran 'DO' statement
    # "DO10I=10,A" is also valid
    # https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vn8c/index.html

    pos     = 1
    inbraces = nextarg = 0
    reflabel = label = countervar = ""
    starts   = [0,0,0]

    while true
        ttype = picktoken(line, pos)
        #print("_$ttype")
        if ttype == 'e' || ttype == '#'
            break
        elseif ttype == ' '
            pos = skipspaces(line, pos)
        elseif nextarg == 0 && ttype == 'd'
            reflabel, pos = taketoken(line, pos)
        elseif nextarg == 0 && ttype == 'l'
            lex, pos = taketoken(line, skipspaces(line, pos))
            if lowercase(lex) == "do"
                lex, pos = taketoken(line, skipspaces(line, pos))
                if isdigit(lex[1]) # LABEL
                    label = lex
                    pos = skiptoken(',', line, skipspaces(line, pos))
                    lex, pos = taketoken(line, skipspaces(line, pos))
                end
                countervar = lex
            elseif occursin(r"^do\d+$"i, lex)
                label = replace(lex, r"^do(\d+)$"i=>s"\1")
                pos = skiptoken(',', line, skipspaces(line, pos))
                countervar, pos = taketoken(line, skipspaces(line, pos))
            elseif occursin(r"^do\d+[a-z].*$"i, lex)
                label = replace(lex, r"^do(\d+)([a-zA-Z].*)$"i=>s"\1")
                countervar = replace(lex, r"^do(\d+)([a-zA-Z].*)$"i=>s"\2")
            elseif occursin(r"^do[a-z].*$"i, lex)
                countervar = replace(lex, r"^do([a-z].*)$"i=>s"\1")
            else # something goes wrong
                return false,"","","","","","",""
            end
            pos = skipspaces(line, pos)
            if picktoken(line, pos) != '='
                # it is not the 'DO' statement
                return false,"","","","","","",""
            end
            pos = skiptoken(line, pos)
            starts[nextarg+=1] = pos
        elseif ttype == '(' || ttype == '['
            inbraces += 1
            pos = skiptoken(line, pos)
        elseif ttype == ')' || ttype == ']'
            inbraces -= 1
            pos = skiptoken(line, pos)
        elseif ttype == ','
            pos = skiptoken(line, pos)
            if inbraces == 0
                starts[nextarg+=1] = pos
            end
        elseif ttype == 'l'
                lex, pos1 = taketoken(line, pos)
                occursin(r"#?CMMNT\d{10}", lex) && break
                pos = pos1
        else
            pos = skiptoken(line, pos)
        end
    end

    if starts[1] == 0 || starts[2] == 0
        return false,"","","","","","",""
    elseif starts[3] == 0
        starts[3] = pos + 1
    end

    # cross our fingers to have simple ascii without multibyte characters
    return true, reflabel, label, countervar,
           strip(line[starts[1]:starts[2]-2]),
           strip(line[starts[2]:starts[3]-2]),
           strip(line[starts[3]:pos-1]),
           rstrip(line[pos:end])                   # trailing comment
end

function replacedocontinue(code, dolabels, gotolabels; dontfixcontinue=false)
    # https://regex101.com/r/FG2iyI/4
    #rxlabeled = r"^(?=[ ]{0,4}\d[ ]{0,4})([\d ]{0,5})"
    rxlabeled = r"^(?:\h*(\d+)\h+)\b\w+\b"im
    rxcont = r"^(\h*(\d+)\h+)continue\h*"im
    rxlabel = r"^(\h*(\d+)\h+)(.*?)$"im

    # insert absent 'CONTINUE' in 'LABEL SOMETHING'
    if !dontfixcontinue
        for mx in reverse(collect(eachmatch(rxlabeled, code)))
            r = continuedlinesrange(code, mx.offset)
            isnothing(match(rxcont, code[r])) || continue
            str = stripcomments(rstrip(concatcontinuedlines(code[r])))
            (m = match(rxlabel, str)) !== nothing || continue
            if dolabels[m.captures[2]] > 0
                res = " "^length(m.captures[1]) * replace(str, rxlabel=>s"\3")
                res *= "\n" * m.captures[1] * "continue\n"
                code = code[1:prevind(code,r[1])] *
                       replace(res, r"\t\t"=>"\t\r\t") * code[nextind(code,r[end]):end]
            end
        end
    end

    lines, comments = splitoncomment(code)
    # replace 'LABEL CONTINUE' with 'end do'
    for i in axes(lines,1)
        m = match(rxcont, lines[i])
        if !isnothing(m) && dolabels[m.captures[2]] > 0
            lbl = m.captures[2]
            #if lbl in gotolabels
            #    # 'CYCLE' emulation in the old fortran,
            #    # beware to have 'GOTO LABEL' somewhere
            #    rxgoto = Regex("go\\h*to\\h+"*lbl, "i")
            #    lines = map(a->replace(a, rxgoto => s"cycle"), lines)
            #    # there also can be conditional GOTOs with CYCLE and BREAK
            #    # it should be fixed by hands
            #end
            lines[i] = " "^length(m.captures[1]) * "end do"
            dolabels[lbl] -= 1
            for _ = dolabels[lbl]:-1:1
                lines[i] *= '\n' * ' '^length(m.captures[1]) * "end do"
            end
            dolabels[lbl] = 0
            # save LABEL for another GOTOs
            if lbl in gotolabels
                lines[i] *= '\n' * lbl * " "^(length(m.captures[1])-length(lbl)) * "continue"
            end
        end
    end

    # not replaced LABELs should be fixed by hands
    rx = r"^(\h*(?:do|for)\h+)(\d+)(\h+.*)$"i
    for i in axes(lines,1)
        m = match(rx, lines[i])
        if !isnothing(m) && dolabels[m.captures[2]] > 0
            lines[i] = replace(lines[i], rx => s"\1FIXME: L\2\3")
        end
    end

    lines = map(*, lines, comments)
    return sameasarg(code, lines)
end

"""
replace conditional GOTOs with GOTOs in branch of julia's `if else` statements
"""
function processconditionalgotos(code)
    # GOTO (10,20,30,40) SOMEEXPR
    # https://regex101.com/r/nKZmku/3
    rx = r"(\h*\d+:?\h*|\h*)go\h*?to\h*\(\h*\d+\h*,.*"i
    rxfull = r"^(\h*\d+:?\h*|\h*)(.*?|)(\h*|)go\h*?to\h*\((\h*\d+\h*(?:,\h*\d+\h*)*)\)(?:\h*,|)\h*([^#\n]*)\h*(#.*|)$"mi
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        str = stripcommentsbutlast(rstrip(concatcontinuedlines(code[r])))
        c = match(rxfull, str).captures
        hl = length(c[1]); head = ' '^hl
        if isempty(strip(c[2]))
            ifexpr = c[1]
            tail = isempty(strip(c[6])) ? "" : "\n$(head)$(c[6])"
        else # there can prepend only the one-line-if, thus expand it to multiline
            @assert occursin(r"^if"i, c[2]) == true
            head *= "    "
            ifexpr = "$(c[1])$(c[2]) then\n$head"
            tail = isempty(strip(c[6])) ? "\n$(' '^hl)end" : "\n$(head)$(c[6])\n$(' '^hl)end"
        end
        expr = c[5]; gotos = map(strip, split(c[4], ','))
        if occursin(r"^[A-Za-z][A-Za-z0-9_]*$", expr)
            cond = expr
        else
            cond = @sprintf "GTO%s" mod(hash("$(expr)"), 2^16)
            ifexpr *= head * "$cond = $expr\n"
        end
        for (i,g) = enumerate(gotos)
            ifexpr *= "if ($(cond) == $(i)) then\n$(head)    $(mask('@'))goto L$(g)\n$(head)else"
        end
        ifexpr *= "\n$(head)    $(mask('@'))error(\"non-exist label " *
                  "$(cond)=$(mask('"'))\$($(cond))$(mask('"')) " *
                  "in goto list\")\n$(head)end$(tail)\n"
        code = code[1:prevind(code,r[1])] * ifexpr * code[nextind(code,r[end]):end]
    end

    # GOTO EXPR,(10,20,30,40)
    # https://regex101.com/r/nHwb10/4
    # https://regex101.com/r/iYG7BH/6
    rx = r"(\h*\d+:?\h*|\h*)go\h*?to\h*\w+\h*(?:,\h*|)\(\h*\d+\h*,.*"i
    rxfull = r"^(\h*\d+:?\h*|\h*)(.*?|)(\h*|)go\h*?to\h+(\w+)\h*(?:,\h*|)\((\h*\d+\h*(?:,\h*\d+\h*)*)\)(.*)$"mi
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        str = stripcommentsbutlast(rstrip(concatcontinuedlines(code[r])))
        c = match(rxfull, str).captures
        head = ' '^length(c[1])
        tail = isempty(strip(c[6])) ? "" : "\n$(head)$(c[6])"
        @assert c[2] == ""
        ifexpr = isempty(strip(c[2])) ? c[1] : c[1] * c[2] * '\n' * head
        var = c[4]; gotos = map(strip, split(c[5], ','))
        for g in gotos
            ifexpr *= "if ($(var) == $(g)) then\n$(head)    $(mask('@'))goto L$(g)\n$(head)else"
        end
        ifexpr *= "\n$(head)    $(mask('@'))error(\"non-exist label " *
                  "$(var)=$(mask('"'))\$($(var))$(mask('"')) " *
                  "in goto list\")\n$(head)end$(tail)\n"
        code = code[1:prevind(code,r[1])] * ifexpr * code[nextind(code,r[end]):end]
    end
    return code
end

function processifstatements(code)

    rx = r"(^\h{0,4}\d*\h*|;)(\h*else|)\h*if\h*\(.*"mi
    for m in reverse(collect(eachmatch(rx, code)))
        o = m.offset
        s1 = m.captures[1]; s2 = lowercase(m.captures[2])
        iflength, wothen, ifcond, ifexec, afterthencomment = parseifstatement(code, o)
        length(ifexec) > 0 && isletter(ifexec[1]) && (ifexec = " " * ifexec)
        length(afterthencomment) > 0 && (afterthencomment = " " * afterthencomment)
        #@debug o, iflength, ifcond, ifexec, afterthencomment
        if wothen == true
            str = s1 * s2 * "$(mask("if")) (" * strip(ifcond) * ")"
            if isoneline(ifcond) && isoneline(ifexec)
                # insert "end" before comment or at tail
                str *= replace(ifexec, r"^\h*([^#]*)\h*(#.*|)$" => s" \1 end \2")
            else
                str *= ifexec * "\n" * " "^length(s1) * "end "
            end
        elseif wothen == false && picktoken(ifcond, skipwhitespaces(ifcond, 1)) in ('#','\n',' ')
            str = s1 * s2 * "$(mask("if")) (" * ifcond * ")" * afterthencomment
        elseif wothen == false
            str = s1 * s2 * "$(mask("if")) " * strip(ifcond) * afterthencomment
        else
            str = ""; iflength = 0
        end
        code = code[1:prevind(code,o)] * str * code[thisind(code,o+iflength):end]
    end

    return replace(code, Regex("$(mask("if"))") => "if")
end

function parseifstatement(str, pos)
    pos0 = pos
    reflabel = ""
    afterthencomment = ""
    iftaken = false
    wothen = nothing
    startcond = startexec = 0
    endcond = endexec = -1
    inbraces  = 0

    while true
        ttype = picktoken(str, pos)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            pos = skipspaces(str, pos)
        elseif !iftaken && ttype == 'd'
            reflabel, pos = taketoken(str, pos)
        elseif !iftaken && ttype == 'l'
            lex, pos = taketoken(str, pos)
            if lowercase(lex) == "else"
                lex, pos = taketoken(str, skipspaces(str, pos))
            end
            if lowercase(lex) == "if" || lowercase(lex) == "elseif"
                iftaken = true
            else
                # that is not an "if" statement
                return 0, wothen, "", "", ""
            end
        elseif iftaken && inbraces == 0 && ttype == 'l'
            lex, pos = taketoken(str, pos)
            if lowercase(lex) == "then"
                wothen = false
                if picktoken(str, skipwhitespaces(str, pos)) == '#'
                    pos, pold = skipcomment(str, skipwhitespaces(str,pos)), pos
                    afterthencomment = str[pold:prevind(str, pos)]
                else
                    pos = skipwhitespaces(str, pos)
                end
                break
            else
                wothen = true
                endexec = pos = skipupto('\n', str, pos)
                endexec = prevind(str, endexec)
            end
        elseif ttype == '('
            inbraces += 1
            pos = skiptoken(str, pos)
            startcond == 0 && (startcond = pos)
        elseif ttype == ')'
            inbraces -= 1
            if inbraces == 0 && endcond == -1
                endcond = prevind(str, pos)
                startexec = nextind(str, pos)
            end
            pos = skiptoken(str, pos)
        else
            pos = skiptoken(str, pos)
        end
    end

    return ncodeunits(str[pos0:prevind(str,pos)]), wothen,
           str[startcond:endcond], str[startexec:endexec], afterthencomment
end

function processarithmif(code::AbstractString)
    # 113 IF (IERROR) 119,114,119
    # https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vn9p/index.html
    # https://regex101.com/r/5R0qtY/1/
    rx = r"\bif\b\h*(\(((?>[^\(\)]++|(?1))*)\))\h*(\d+)\h*,\h*(\d+)\h*,\h*(\d+)"mi
    for m in reverse(collect(eachmatch(rx, code)))
        line = m.match
        str = "\n"
        expr = m.captures[2]
        L1, L2, L3 = m.captures[3:5]
        if occursin(r"^[A-Za-z][A-Za-z0-9_]*$", expr)
            cond = expr
        else
            cond = @sprintf "CND%s" mod(hash("$(expr)"), 2^16)
            str *= "      $cond = $expr\n"
        end
        str *= "      if $cond < 0\n" *
               "        goto $L1\n" *
               "      elseif $cond == 0\n" *
               "        goto $L2\n" *
               "      else\n" *
               "        goto $L3\n" *
               "      end\n      "
        o = m.offset; l = ncodeunits(m.match)
        code = code[1:prevind(code,o)] * str * code[thisind(code,o+l):end]
    end
    return code
end

function shapinglabels(code::AbstractString)
    rx = r"^(((\h+|)(\d+)(\h+))(\w+|))"m
    for m in reverse(collect(eachmatch(rx, code)))
        label, tail = m.captures[4], m.captures[6]
        o = m.offset; l = ncodeunits(m.captures[1])
        head = ' '^ncodeunits(m.captures[2])
        str = "$(head)$(mask('@'))label L$label"
        uppercase(tail) != "CONTINUE" && (str *= "\n$head$tail")
        code = code[1:prevind(code,o)] * str * code[thisind(code,o+l):end]
    end
    return code
end

function processimplieddoloops(code)
    W(s) = replace("$s", " " => "")

    # should be:
    # `A = [(real(I), I=3, 5)]`   => `A = [real(I) for I=3, 5]`
    # TODO: not released: `/(real(I), I=3, 5)/` => `[real(I) for I=3, 5]`
    # `(WRITE(*,*)(A(I), I=3, 5)` => `print(view(A, 3:5)...)` # splatting in splatprintviews()
    # `(READ(*,*)(A(I), I=3, 5)`  => `READ(view(A, 3:5))`
    # `DATA (A(I), I=3, 5)`       => `DATA view(A, 3:5)`

    # https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/io.html
    # ( item-1, item-2, ...., item-n, DO-var = initial, final, step )
    # https://regex101.com/r/RW27vC/2
    rxbr = r"(?<=\W)(\(((?>[^()]++|(?1))*)\))"
    rxsqbr = r"(?<=\W)\[\h*(\(((?>[^()]++|(?1))*)\))\h*\]"
    rxslbr = r"(?<=\W)\/\h*(\(((?>[^()]++|(?1))*)\))\h*\/"

    for m in reverse(collect(eachmatch(rxsqbr, code)))
        str = m.captures[2]
        if occursin('=', str)
            str = processimplieddoloops(str)
            ex = nothing
            if length(collect(eachmatch(r"=", str))) == 1
                try ex = Meta.parse('('*str*')') catch end
            end
            #(ex == nothing || ex.head == :incomplete) && continue
            if (p = matchex(:(_E, _I = _1, _2, _INC), ex)) !== nothing
                iex, var, r1, r2, inc = p[:_E], p[:_I], p[:_1], p[:_2], p[:_INC]
                code = replace(code, m.match => "[$iex for $var = $(W(r1)):$(W(inc)):$(W(r2))]")
            elseif (p = matchex(:(_E, _I = _1, _2), ex)) !== nothing
                iex, var, r1, r2 = p[:_E], p[:_I], p[:_1], p[:_2]
                code = replace(code, m.match => "[$iex for $var = $(W(r1)):$(W(r2))]")
            else
                code = replace(code, m.match => "[("*str*")]")
            end
        end
    end

    for m in reverse(collect(eachmatch(rxbr, code)))
        str = m.captures[2]
        if occursin('=', str)
            str = processimplieddoloops(str)
            ex = nothing
            if length(collect(eachmatch(r"=", str))) == 1
                try ex = Meta.parse('('*str*')') catch end
            end
            #@show str, ex
            if (p = matchex(:(_NAME(_E), _I = _1, _2, _INC), ex)) !== nothing ||
               (p = matchex(:(_NAME[_E], _I = _1, _2, _INC), ex)) !== nothing
                nam, iex, var, r1, r2, inc = p[:_NAME], p[:_E], p[:_I], p[:_1], p[:_INC], p[:_2]
                r1, r2 = map(r -> rewrite_all(iex, [var => r]), (r1, r2))
                str = "view($nam, $(W(r1)):$(W(inc)):$(W(r2)))"
            elseif (p = matchex(:(_NAME(_E1,_E2), _I = _1, _2, _INC), ex)) !== nothing ||
                   (p = matchex(:(_NAME[_E1,_E2], _I = _1, _2, _INC), ex)) !== nothing
                nam, iex1, iex2, var, r1, r2, inc =
                    p[:_NAME], p[:_E1], p[:_E2], p[:_I], p[:_1], p[:_INC], p[:_2]
                if matchingex(var, iex1)
                    r1, r2 = map(r -> rewrite_all(iex1, [var => r]), (r1, r2))
                    str = "view($nam, $(W(r1)):$(W(inc)):$(W(r2)), $(W(iex2)))"
                else
                    r1, r2 = map(r -> rewrite_all(iex2, [var => r]), (r1, r2))
                    str = "view($nam, $(W(iex1)), $(W(r1)):$(W(inc)):$(W(r2)))"
                end
            elseif (p = matchex(:(_NAME(_E), _I = _1, _2), ex)) !== nothing ||
                   (p = matchex(:(_NAME[_E], _I = _1, _2), ex)) !== nothing
                nam, iex, var, r1, r2 = p[:_NAME], p[:_E], p[:_I], p[:_1], p[:_2]
                r1, r2 = map(r -> rewrite_all(iex, [var => r]), (r1, r2))
                str = "view($nam, $(W(r1)):$(W(r2)))"
            elseif (p = matchex(:(_NAME(_E1,_E2), _I = _1, _2), ex)) !== nothing ||
                   (p = matchex(:(_NAME[_E1,_E2], _I = _1, _2), ex)) !== nothing
                nam, iex1, iex2, var, r1, r2 = p[:_NAME], p[:_E1], p[:_E2], p[:_I], p[:_1], p[:_2]
                if matchingex(var, iex1)
                    r1, r2 = map(r -> rewrite_all(iex1, [var => r]), (r1, r2))
                    str = "view($nam, $(W(r1)):$(W(r2)), $(W(iex2)))"
                else
                    r1, r2 = map(r -> rewrite_all(iex2, [var => r]), (r1, r2))
                    str = "view($nam, $(W(iex1)), $(W(r1)):$(W(r2)))"
                end
            elseif (p = matchex(:(_E, _I = _1, _2), ex)) !== nothing
                iex, var, r1, r2 = p[:_E], p[:_I], p[:_1], p[:_2]
                str = "($iex for $var = $(W(r1)):$(W(r2)))"
            else
                str = '(' * str * ')'
            end
            code = replace(code, m.match => str)
        end
    end

    return code
end

function splatprintviews(code::AbstractString)
    rx = r"view(\(((?>[^()]++|(?1))*)\))"
    rxprint = Regex("\\bprintln\\b|\\b$(mask('@'))printf\\b")
    for m in reverse(collect(eachmatch(rx, code)))
        line = code[continuedlinesrange(code, m.offset)]
        if occursin(rxprint, line)
        #if occursin(r"\bprintln\b|\bat_iZjAcpPokMprintf\b", line)
            o = m.offset
            l = ncodeunits(m.match)
            code = code[1:prevind(code,o)] * m.match * "..." * code[thisind(code,o+l):end]
        end
    end
    return code
end

function processselectcase(code)
    lines = splitonlines(code)
    var = ""
    case1 = false
    for i in axes(lines,1)
        if occursin(r"select\h+case"i, lines[i])
            expr = replace(lines[i], r"^\h*select\h+case\h*\((.*)\)(\h*#.*|\h*)$"mi => s"\1")
            if occursin(r"^[A-Za-z][A-Za-z0-9_]*$", expr)
                var = expr
                lines[i] = replace(lines[i], r"^(\h*)(select\h+case.*)$"mi => SS("\\1#\\2"))
            else
                var = @sprintf "EX%s" mod(hash("$(expr)"), 2^16)
                lines[i] = replace(lines[i], r"^(\h*)(select\h+case.*)$"mi =>
                                             SS("\\1"*var*" = "*expr*" #\\2") )
            end
            case1 = true
        end
        # TODO: process 'CASE (2:4)' with `if var in (2:4)` or `if 2 <= var <= 4`
        if !isempty(var) && occursin(r"^\h*case\h*\("mi, lines[i])
            exprs = replace(lines[i], r"^\h*case\h*\((.*)\)(\h*#.*|\h*)$"mi => s"\1")
            if occursin(',', exprs)
                exprs = replace(exprs, r"(\h*,\h*)" => SS(" || "*var*" == "))
            end
            exprs = "if $(var) == " * exprs
            if !case1 exprs = "else" * exprs end
            lines[i] = replace(lines[i], r"^(\h*)case\h*\(.*\)(\h*#.*|\h*)$"mi =>
                                         SS("\\1"*exprs*"\\2"))
            case1 = false
        elseif !isempty(var) && occursin(r"^\h*case\h+default"mi, lines[i])
            lines[i] = replace(lines[i], r"^(\h*)case\h+default(\h*#.*|\h*)$"mi => s"\1else\2")
            var = ""
        elseif !isempty(var) && occursin(r"^\h*end\h+select"mi, lines[i])
            var = ""
        end
    end

    return sameasarg(code, lines)
end

function processparameters(code::AbstractString)

    rx = r"^\h*parameter\h*\("mi
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        line = code[r]
        str = lex = ""
        statementtaken = false
        inbraces = 0; p = 1
        while true
            ttype = picktoken(line, p)
            if ttype == 'e'
                str *= '\n'
                break
            elseif inbraces == 0 && ttype == '('
                inbraces += 1
                p = nextind(line, p)
                str *= ' '
            elseif ttype == '('
                p1 = p; p = skipbraces(line, p)
                str *= "Complex" * line[p1:prevind(line, p)]
            elseif inbraces == 1 && ttype == ')'
                inbraces -= 1
                p = nextind(line, p)
                str *= ' ' * line[p:end]
                break
            elseif !statementtaken && ttype == 'l'
                lex, p = taketoken(line, p)
                @assert uppercase(lex) == "PARAMETER"
                str *= ' '^length(lex)
                statementtaken = true
            elseif ttype == ','
                p = skiptoken(line, p)
                str *= ';'
            else
                p1 = p; p = skiptoken(line, p)
                str *= line[p1:prevind(line, p)]
            end
        end
        code = code[1:prevind(code,r[1])] * str * code[nextind(code,r[end]):end]
    end

    rx = r"\bparameter\b[^\n#]*::"mi
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        line = code[r]
        str = ""
        statementtaken = false
        inbraces = 0; p = 1
        while true
            ttype = picktoken(line, p)
            if ttype == 'e'
                str *= '\n'
                break
            elseif !statementtaken && ttype == ':' && picknexttoken(line, p) == ':'
                p = nextind(line, nextind(line, p))
                str *= "      "
                statementtaken = true
            elseif !statementtaken
                p1 = p; p = skiptoken(line, p)
                #str *= ' '^(p-p1)
            elseif ttype == '(' || ttype == '['
                inbraces += 1
                p = nextind(line, p)
                str *= ttype
            elseif ttype == ')' || ttype == ']'
                inbraces -= 1
                p = nextind(line, p)
                str *= ttype
            elseif inbraces == 0 && ttype == ','
                p = nextind(line, p)
                str *= ';'
            else
                p1 = p; p = skiptoken(line, p)
                str *= line[p1:prevind(line, p)]
            end
        end
        code = code[1:prevind(code,r[1])] * str * code[nextind(code,r[end]):end]
    end

    return code
end

function processdatastatement(code::AbstractString, arrays)

    rx = r"^\h*\bdata\b"mi
    for m in reverse(collect(eachmatch(rx, code)))
        r = continuedlinesrange(code, m.offset)
        line = code[r]
        str = lex = ""
        p = 1; statementtaken = itscal = itvector = afterlist = false
        while true
            ttype = picktoken(line, p)
            if ttype == 'e'
                str *= '\n'
                break
            elseif ttype == '['
                p1 = p; p = skipbrackets(line, p)
                str *= line[p1:prevind(line, p)]
                itscal = true
            elseif ttype == '('
                p1 = p; p = skipbraces(line, p)
                str *= line[p1:prevind(line, p)]
            elseif !statementtaken && ttype == 'l'
                lex, p = taketoken(line, p)
                @assert uppercase(lex) == "DATA"
                str *= ' '^length(lex)
                statementtaken = true
            elseif ttype == 'l'
                p1 = p; lex, p = taketoken(line, p)
                str *= line[p1:prevind(line, p)]
            elseif ttype == '/'
                p1 = p; p = nextind(line, p)
                p = skipupto('/', line, p)
                itvector = (lex in arrays && !itscal)
                str[end:end] != ' ' && (str *= ' ')
                str *= itvector ? ".= (" : "= "
                # TODO: here the repeats like `8*0.0` should be caught
                str *= replace(line[nextind(line,p1):prevind(line,p)], r"(?<=\d) ([\dDE])"i=>s"\1")
                p = skiptoken(line, p)
                itvector && (str *= ')')
                afterlist = true
                itscal = itvector = false
            elseif ttype == ',' && afterlist
                p = skiptoken(line, p)
                str *= ';'
                afterlist = false
            else
                p1 = p; p = skiptoken(line, p)
                str *= line[p1:prevind(line, p)]
            end
        end
        code = code[1:prevind(code,r[1])] * str * code[nextind(code,r[end]):end]
    end

    return code
end

function processlinescontinuation(code::AbstractString)
    # fix absent tabs
    code = replace(code, r"\t?(\n|\r)\t|\t(\n|\r)\t?" => SS("\t\r\t"))
    # move arithmetic operator from start of continuator to tail of previous line
    rx = r"(\h*)(#?CMMNT\d{10}|)\t\r\t(\h*|)(\/\/|[+*\/,=(-]|==|<=|>=|!=|<|>|&&|\|\|)"
    code = replace(code, rx => SS("\\4\\1\\2\t\r\t\\3"))
    return code
end

function concatcontinuedlines(code)
    # Note: trailing comments are not allowed here
    return code isa AbstractVector ?
           map(a->replace(a, r"\t(\n|\r)\t?|\t?(\n|\r)\t"=>"\t\t"), code) :
           replace(code, r"\t(\n|\r)\t?|\t?(\n|\r)\t"=>"\t\t")
end

function concatcontinuedlines(lines::AbstractVector, i)
    # Note: trailing comments are not allowed here
    result = lines[i]; n = 1
    k = i
    while k > 1 && iscontinuedlines(lines[k-1], lines[k])
        result = lines[k-1] * result
        n += 1; k -= 1
    end
    k = i
    while k < length(lines) && iscontinuedlines(lines[k], lines[k+1])
        result *= lines[k+1]
        n += 1; k += 1
    end
    return result, n
end

iscontinuedlines(line1, line2) = occursin(r"\t$", line1) || occursin(r"^\t", line2)
function iscontinuedline(str, i)
    i = thisind(str, i)
    while i < ncodeunits(str) && str[i] !='\n'
        i = nextind(str, i)
    end
    if i == ncodeunits(str) || str[prevind(str, i)] != '\t'
        return false
    else
        return true
    end
end


function concatlines(line1, line2)
    # comments will be stripped
    return replace(line1, r"^(.*)(?:#?CMMNT\d{10}|)\t?$"m => s"\1") *
           replace(line2, r"^\t?(.*)$"m => s"\1")
end


function splitoncomment(code)
    rx = r"(?:#?CMMNT\d{10})"
    #rx = r"(?=#?CMMNT\d{10})"
    #rx = r"\h*#.*$"
    lines = splitonlines(code)
    comments  = ["" for i=axes(lines,1)]
    for i in axes(lines,1)
        (m = match(rx, lines[i])) !== nothing && (comments[i] = m.match)
        lines[i] = replace(lines[i], rx => s"")
    end
    # cut out all whitespaces not masked by the comments
    lines = map(a->replace(a, r"[ ]*(\t|)$" => s"\1"), lines)
    return lines, comments
end

stripcomments(code::AbstractString) = replace(code, r"(?<=[ ])[ ]*#?CMMNT\d{10}"m => "")
stripcommentsbutlast(code::AbstractString) = replace(code, r"(?<=[ ])[ ]*#?CMMNT\d{10}(?!\t?$)"m => "")

tostring(code::AbstractVector) =
    replace(foldl((a,b) -> a*'\n'*b, code),
            r"\t?(\n|\r)\t|\t(\n|\r)\t?" => "\t\r\t")
tostring(code::AbstractString) = code

splitonlines(code::AbstractVector) = deepcopy(code)
splitonlines(code) = split(code, r"\n|\r")
splitonfulllines(code) = split(code, r"\n")

sameasarg(like::AbstractString, code::AbstractString)  = code
sameasarg(like::AbstractVector, lines::AbstractVector) = lines
sameasarg(like::AbstractVector, code::AbstractString)  = splitonlines(code)
sameasarg(like::AbstractString, lines::AbstractVector) = tostring(lines)


const repairreplacements = OrderedDict(
    #??? # replace 'PRINT' => 'WRITE'
    #??? r"(PRINT)\h*([*]|'\([^\n\)]+\)')\h*,"mi => s"write(*,\2)",
    # Goto
    r"^(\s*)GO\s*TO\s+(\d+)"mi => s"\1goto \2"  ,
    r"(\h)GO\s*TO\s+(\d+)"mi   => s"\1goto \2"  ,
    r"(\W)GO\s*TO\s+(\d+)"mi   => s"\1 goto \2" ,
    # NOOP
    r"^\h*CONTINUE\h*\n"mi     => "",
    # "**" => '^'
    r"(?<!\*)\*\*(?!\*)"m => s"^",
    # spaces in numbers into underscores "1 000 E3" => "1_000E3"
    r"(\b\d+)[ ](\d+\b)"m                      => s"\1@\2",  # '@' will be restored later
    r"(\b\d+)[ ](\d+(?:[eEdD][-+]?\d+)?\b)"m   => s"\1@\2",
    r"(\b\d+)[ ]([eEdD][-+]?\d+\b)"m           => s"\1\2",
    # fix floating point numbers with exponent
    # https://regex101.com/r/RRHEyN/1
    r"([-+]?[0-9]*\.?[0-9]+)[eE]([-+]?[0-9]+)"m => s"\1f\2",
    r"([-+]?[0-9]*\.?[0-9]+)[d]([-+]?[0-9]+)"m  => s"\1e\2",
    r"([-+]?[0-9]*\.?[0-9]+)[D]([-+]?[0-9]+)"m  => s"\1E\2",
    # fix inclomplete floating point numbers: "1./VAR" => "1.0/VAR"
    r"(\d)\.([/*+-])"               => s"\1.0\2",
    # custom precision floats: "1.0e+3_MY_SUPER_FLOAT" => "MY_SUPER_FLOAT(1.0e+3)"
    r"([-+]?[0-9]*\.?[0-9]+(?:[fFeE][-+]?[0-9]+)?)_(\w+)"m => s"\2(\1)",
    r"(\d)@(\d)"                                           => s"\1_\2",
    # spaces in keywords
    r"\bDOUBLE\b\h+\bPRECISION\b"mi => "DOUBLEPRECISION",
    r"\bDO\b\h+\bWHILE\b"mi         => "DOWHILE",
    r"\bELSE\b\h+\bIF\b"mi          => "ELSEIF",
    r"\bSELECT\b\h+\bCASE\b"mi      => "SELECTCASE",
    r"\bCASE\b\h+\bDEFAULT\b"mi     => "CASEDEFAULT",
    r"\bBLOCK\b\h+\bDATA\b"mi       => "BLOCKDATA", # also need to catch corresponding end
    r"^(\h*|)(\bEND\b)\h+(\b(?:IF|DO|SELECT|FUNCTION|SUBROUTINE|PROGRAM|MODULE|INTERFACE|TYPE)\b)"mi => s"\1\2\3",
    # modern logical symbols
    r"\.true\."i           => "true"    ,
    r"\.false\."i          => "false"   ,
    r"(\h+)\.or\.(\h+)"i   => s"\1||\2" ,
    r"(\h+)\.and\.(\h+)"i  => s"\1&&\2" ,
    r"(\h+)\.not\.(\h+)"i  => s"\1!"    ,
    r"(\h+)\.eq\.(\h+)"i   => s"\1==\2" ,
    r"(\h+)\.ne\.(\h+)"i   => s"\1!=\2" ,
    r"(\h+)\/=(\h+)"i      => s"\1!=\2" ,
    r"(\h+)\.le\.(\h+)"i   => s"\1<=\2" ,
    r"(\h+)\.ge\.(\h+)"i   => s"\1>=\2" ,
    r"(\h+)\.gt\.(\h+)"i   => s"\1>\2"  ,
    r"(\h+)\.lt\.(\h+)"i   => s"\1<\2"  ,
    r"(\h+)\.or\.\h*"i     => s"\1|| "  ,
    r"(\h+)\.and\.\h*"i    => s"\1&& "  ,
    r"(\h+)\.not\.\h*"i    => s"\1!"    ,
    r"(\h+)\.eq\.\h*"i     => s"\1== "  ,
    r"(\h+)\.ne\.\h*"i     => s"\1!= "  ,
    r"(\h+)\/=\h*"i        => s"\1!= "  ,
    r"(\h+)\.le\.\h*"i     => s"\1<= "  ,
    r"(\h+)\.ge\.\h*"i     => s"\1>= "  ,
    r"(\h+)\.gt\.\h*"i     => s"\1> "   ,
    r"(\h+)\.lt\.\h*"i     => s"\1< "   ,
    r"\h*\.or\.\h*"i       => s" || "   ,
    r"\h*\.and\.\h*"i      => s" && "   ,
    r"\h*\.not\.\h*"i      => s" !"     ,
    r"\h*\.eq\.\h*"i       => s" == "   ,
    r"\h*\.ne\.\h*"i       => s" != "   ,
    r"\h*\/=\h*"i          => s" != "   ,
    r"\h*\.le\.\h*"i       => s" <= "   ,
    r"\h*\.ge\.\h*"i       => s" >= "   ,
    r"\h*\.gt\.\h*"i       => s" > "    ,
    r"\h*\.lt\.\h*"i       => s" < "    ,
    # remove whitespace between function name and left brace https://regex101.com/r/CxV232/3
    r"(\b(?!elseif|if|while|data|parameter|selectcase|case)[\w_]+\b)[\h]+\("mi => s"\1(",
    # "a%b" => "a.b"
    r"(\b\w+)\h*%\h*(\w+\b)"m => s"\1.\2",
)


const subscriptreplacements = OrderedDict(
    r"\b(\w+)_0([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₀\2",
    r"\b(\w+)_1([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₁\2",
    r"\b(\w+)_2([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₂\2",
    r"\b(\w+)_3([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₃\2",
    r"\b(\w+)_4([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₄\2",
    r"\b(\w+)_5([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₅\2",
    r"\b(\w+)_6([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₆\2",
    r"\b(\w+)_7([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₇\2",
    r"\b(\w+)_8([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₈\2",
    r"\b(\w+)_9([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1₉\2",
    r"\b(\w+)_a([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₐ\2",
    r"\b(\w+)_e([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₑ\2",
    r"\b(\w+)_h([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₕ\2",
    r"\b(\w+)_i([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ᵢ\2",
    r"\b(\w+)_j([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ⱼ\2",
    r"\b(\w+)_k([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₖ\2",
    r"\b(\w+)_l([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₗ\2",
    r"\b(\w+)_m([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₘ\2",
    r"\b(\w+)_n([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₙ\2",
    r"\b(\w+)_o([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₒ\2",
    r"\b(\w+)_p([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₚ\2",
    r"\b(\w+)_r([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ᵣ\2",
    r"\b(\w+)_s([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₛ\2",
    r"\b(\w+)_t([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₜ\2",
    r"\b(\w+)_u([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ᵤ\2",
    r"\b(\w+)_v([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ᵥ\2",
    r"\b(\w+)_x([₀₁₂₃₄₅₆₇₈₉ᵢᵣᵤᵥᵦᵧᵨᵩᵪₐₑₒₓₕₖₗₘₙₚₛₜⱼ]*)\b"  => s"\1ₓ\2",
)
const greeksubscriptreplacements = OrderedDict(
    r"\b(\w+)_schwa\b" => s"\1ₔ",
    r"\b(\w+)_beta\b"  => s"\1ᵦ",
    r"\b(\w+)_gamma\b" => s"\1ᵧ",
    r"\b(\w+)_rho\b"   => s"\1ᵨ",
    r"\b(\w+)_phi\b"   => s"\1ᵩ",
    r"\b(\w+)_chi\b"   => s"\1ᵪ",
)
const greeksreplacements = OrderedDict(
    # Not all greek letters are distinguished from the latin letters thus skipped.
    # Skipped letters have case insensitive lowercase counterparts.

    #r"\bALPHA([_\d]*)\b"       => s"Α\1" ,
    #r"\bBETA([_\d]*)\b"        => s"Β\1" ,
    r"\bGAMMA([_\d]*)\b"       => s"Γ\1" ,
    r"\bDELTA([_\d]*)\b"       => s"Δ\1" ,
    #r"\bEPSILON([_\d]*)\b"     => s"Ε\1" ,
    #r"\bZETA([_\d]*)\b"        => s"Ζ\1" ,
    #r"\bETA([_\d]*)\b"         => s"Η\1" ,
    r"\bTHETA([_\d]*)\b"       => s"Θ\1" ,
    #r"\bIOTA([_\d]*)\b"        => s"Ι\1" ,
   # r"\bKAPPA([_\d]*)\b"       => s"Κ\1" ,
    r"\bLAMBDA([_\d]*)\b"      => s"Λ\1" ,
    r"\bXI([_\d]*)\b"          => s"Ξ\1" ,
    r"\bPI([_\d]*)\b"          => s"Π\1" ,
    #r"\bRHO([_\d]*)\b"         => s"Ρ\1" ,
    r"\bSIGMA([_\d]*)\b"       => s"Σ\1" ,
    #r"\bTAU([_\d]*)\b"         => s"Τ\1" ,
    #r"\bUPSILON([_\d]*)\b"     => s"Υ\1" ,
    r"\bPHI([_\d]*)\b"         => s"Φ\1" ,
    #r"\bCHI([_\d]*)\b"         => s"Χ\1" ,
    r"\bPSI([_\d]*)\b"         => s"Ψ\1" ,
    r"\bOMEGA([_\d]*)\b"       => s"Ω\1" ,
    r"\balpha([_\d]*)\b"i       => s"α\1" ,
    r"\bbeta([_\d]*)\b"i        => s"β\1" ,
    r"\bgamma([_\d]*)\b"       => s"γ\1" ,
    r"\bdelta([_\d]*)\b"       => s"δ\1" ,
    r"\bzeta([_\d]*)\b"i        => s"ζ\1" ,
    r"\beta([_\d]*)\b"i         => s"η\1" ,
    r"\btheta([_\d]*)\b"       => s"θ\1" ,
    r"\biota([_\d]*)\b"i        => s"ι\1" ,
    r"\bkappa([_\d]*)\b"i       => s"κ\1" ,
    r"\blambda([_\d]*)\b"      => s"λ\1" ,
    r"\bmu([_\d]*)\b"          => s"μ\1" ,
    r"\bnu([_\d]*)\b"          => s"ν\1" ,
    r"\bxi([_\d]*)\b"          => s"ξ\1" ,
    r"\bpi([_\d]*)\b"          => s"π\1" ,
    r"\brho([_\d]*)\b"i         => s"ρ\1" ,
    r"\bvarsigma([_\d]*)\b"i    => s"ς\1" ,
    r"\bsigma([_\d]*)\b"       => s"σ\1" ,
    r"\btau([_\d]*)\b"i         => s"τ\1" ,
    r"\bupsilon([_\d]*)\b"i     => s"υ\1" ,
    r"\bvarphi([_\d]*)\b"i      => s"φ\1" ,
    r"\bchi([_\d]*)\b"i         => s"χ\1" ,
    r"\bpsi([_\d]*)\b"         => s"ψ\1" ,
    r"\bomega([_\d]*)\b"       => s"ω\1" ,
    r"\bvartheta([_\d]*)\b"    => s"ϑ\1" ,
    r"\bphi([_\d]*)\b"         => s"ϕ\1" ,
    r"\bvarpi([_\d]*)\b"i       => s"ϖ\1" ,
    r"\bStigma([_\d]*)\b"i     => s"Ϛ\1" ,
    r"\bDigamma([_\d]*)\b"i    => s"Ϝ\1" ,
    r"\bKoppa([_\d]*)\b"       => s"Ϟ\1" ,
    #r"\bSampi([_\d]*)\b"i      => s"Ϡ\1" ,
    r"\bvarkappa([_\d]*)\b"i    => s"ϰ\1" ,
    r"\bvarrho([_\d]*)\b"i      => s"ϱ\1" ,
    r"\bVARTHETA([_\d]*)\b"    => s"ϴ\1" ,
    r"\bepsilon([_\d]*)\b"i     => s"ϵ\1" ,
    r"\bbackepsilon([_\d]*)\b"i => s"϶\1" ,

    #r"\bforall([_\d]*)\b"i     => s"∀\1" ,
    r"\bpartial([_\d]*)\b"i    => s"∂\1" ,
    r"\bexists([_\d]*)\b"i     => s"∃\1" ,
    r"\bnexists([_\d]*)\b"i    => s"∄\1" ,
    r"\bnabla([_\d]*)\b"i      => s"∇\1" ,
)

const trivialreplacements = OrderedDict(
    #r"(\bif\b|\belseif\b|\bend\b|\bdo\b|\bwhile\b)"i => s"\L\1",
    r"(\bif\b)"i     => s"if"     ,
    r"(\belseif\b)"i => s"elseif" ,
    r"(\bend\b)"i    => s"end"    ,
    r"(\bdo\b)"i     => s"do"     ,
    r"(\bwhile\b)"i  => s"while"  ,
    # Some specific functions
    # https://gcc.gnu.org/onlinedocs/gcc-9.2.0/gfortran/Intrinsic-Procedures.html
    # https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc8/index.html
    r"\bmod\b\h*\("i       => s"mod("       ,
    r"\bd?sqrt\b\h*\("i    => s"sqrt("      ,
    r"\bd?exp\b\h*\("i     => s"exp("       ,
    r"\bd?log\b\h*\("i     => s"log("       ,
    r"\bd?log10\b\h*\("i   => s"log10("     ,
    r"\bd?conjg\b\h*\("i   => s"conj("      ,
    r"\b[da]?imag\b\h*\("i => s"imag("      ,
    r"\bd?sign\b\("i       => "copysign("   ,
    r"\bd?max1?\b\h*\("i   => "max("        ,
    r"\bd?min1?\b\h*\("i   => "min("        ,
    r"\bd?abs\b\h*\("i     => "abs("        ,
    r"\bd?sin\b\h*\("i     => "sin("        ,
    r"\bd?asin\b\h*\("i    => "asin("       ,
    r"\bd?cos\b\h*\("i     => "cos("        ,
    r"\bd?acos\b\h*\("i    => "acos("       ,
    r"\bd?tan\b\h*\("i     => "tan("        ,
    r"\bd?atan2?\b\h*\("i  => "atan("       ,
    r"\b(i[dq])?nint\b\h*\("i => "round(Int, ",
    r"\bint\b\h*\("i       => "trunc(Int, " ,
    r"\breal\b\h*\("i      => "Float32("    , # ?
    r"\bfloat\b\h*\("i     => "Float32("    , # ?
    r"\bsngl\b\h*\("i      => "Float32("    , # ?
    r"\bdble\b\h*\("i      => "float("      ,
    r"\bdfloat\b\h*\("i    => "float("      ,
    r"\bcmplx\b\h*\("i     => "ComplexF32(" , # ?
    r"\bdcmplx\b\h*\("i    => "complex("    ,
    r"\bceiling\b\h*\("i   => "ceil(Int, "  ,
    r"\bkind\b\h*\("i      => "sizeof("     ,
    r"\brand\b\h*\(\)"i    => "rand()"      ,
    # https://regex101.com/r/whrGry/2
    #r"(?:\brand\b\h*)(\(((?>[^()\n]+|(?1))+)\))"i => s"rand(MersenneTwister(round(Int,\2)))",
    r"\blen\b\h*\("i       => "length(",
    # https://regex101.com/r/whrGry/1
    r"(?:\blen_trim\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"length(rstrip(\2))",
    r"(?:\blnblnk\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"length(rstrip(\2))",
    # TODO: index, scan, verify, fdate, ctime
    r"\bi?char\b\h*\("i  => "Char("    ,
    r"\brepeat\b\h*\("i  => "repeat("  ,
    r"\badjustr\b\h*\("i => "ADJUSTR(" ,
    r"\badjustl\b\h*\("i => "ADJUSTL(" ,
    # TODO: reshape, shape, lbound, ubound
    r"\bbackspace\b\h+(\w+)"i  => s"seek(\1, position(\1)-1) # BACKSPACE \1",
    r"\brewind\b\h+(\w+)"i     => s"seek(\1, 0) # REWIND \1",
    r"\bend\h*file\b\h+(\w+)"i => s"seekend(\1) # ENDFILE \1",

    # TODO: to complete
    r"\bMPI_SEND\b\h*\("i  => "MPI.Send("   ,
    r"\bMPI_ISEND\b\h*\("i => "MPI.Isend("  ,
    r"\bMPI_RECV\b\h*\("i  => "MPI.Recv!("  ,
    r"\bMPI_IRECV\b\h*\("i => "MPI.Irecv!(" ,
    r"\bMPI_WTIME\b\h*\("i => "MPI.Wtime("  ,

    # LSAME
    r"\bLSAME\b\h*\("i     => "LSAME("      ,
    #r"(?:\blsame\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"((a,b)->==(uppercase(Char(first(a))),uppercase(Char(first(b))))\1",

    # Custom functions
    #r"(\b[\w_]+\b)\(1:J_LEN\(\1\)\)"i => s"rstrip(\1)",
    #r"(?:\bJ_LEN\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"length(rstrip\1)",
)


# Regex/substitution pairs for replace(). Order matters here.
const replacements = OrderedDict(
    # Powers use ^ not **
    r"([^*\n])([*]{2})([^*])" => s"\1^\3",
    # constant arrays
    r"\(/" => s"[",
    r"/\)" => s"]",
    # Space after commas in expressions
    r"(,)([^\s,\[\]\(\)]+\])" => s"@\2",
    r"(,)(\S)" => s"\1 \2",
    r"@(?!LABEL)" => ",",
    # Replace "return\nend" => "return nothing\nend"
    r"(\h*)\bRETURN\b(\h+end|)$"i => s"\1return nothing\2",
    # include ESCAPEDSTR
    r"^(\s*)INCLUDE\h+([\w]+)$"mi => s"\1include(\2)",
    # Replace do LABEL ... -> do ...
    r"for\h+(?:\d+)(\h+.*)$" => s"for\1",
    # Replace do while -> while
    r"\bDOWHILE\b\h*"i => s"while ",
    r"^(\h*)\bDO\b\h*$"i => s"\1while true",
    # Replace ELSE with else
    r"^(\s*)\bELSE\b"m => s"\1else",
    # Relace END XXXX with end
    r"^(\h*)END\h*(?:FUNCTION|(?:RECURSIVE\h+|)SUBROUTINE|PROGRAM|MODULE|BLOCK|DO|IF|SELECT)\h*\w*$"mi => s"\1end",
    r"^(\h*)END$"mi => s"\1end",
    # Don't need CALL
    r"([^ ])\bCALL\b(\h+)"i => s"\1 ",
    r"\bCALL\b(\h+)"i => s"",
    # Fix assignments
    r"(?<=[^\s=<>!/\\])=(?=[^\s=])" => " = ",
    # Remove expression's brackets after if/elseif/while https://regex101.com/r/eF9bK5/17
    r"^(\h*)(IF|ELSEIF|WHILE)(\h*)(\(((?>[^()]++|(?4))*)\))\h*$"m => s"\1\2\3\5",
    # breaks
    r"\bEXIT\b"mi => s"break",
    #r"\breturn\b"mi => s"return",
    r"\bCYCLE\b"mi => s"continue",
    r"\bSTOP\b\h+(\d+)"mi => s"exit(\1)",
    r"\bSTOP\b"mi => s"exit(1)",
    # Strings concatenation operator
    r"//" => " * ",
    # Format floats as "5.0" not "5."
    r"(\W\d+)\.(?=\D)" => s"\1.0",
    r"(\W\d+)\.$"m => s"\1.0",
    r"(?<=\W)\.(\d+)"m => s"0.\1",
    # Goto
    r"^(\s*)GO\s*TO\s+(\d+)"mi => s"\1@goto L\2",
    r"(\h)GO\s*TO\s+(\d+)"mi => s"\1@goto L\2",
    r"(\W)GO\s*TO\s+(\d+)"mi => s"\1 @goto L\2",
    # ASSIGN statement
    r"^(\h*\d+:?\h*|\h*)ASSIGN\h*(\d+)\h*TO\h*(\w+)$"i => s"\1\3 = \2",
    # fix end
    r"^(.+[^ ])[ ]([ ]+)END$" => s"\1 end\2",
    # remove whitespace between function name and left brace https://regex101.com/r/CxV232/3
    r"(\b(?!ELSEIF|IF)[\w_]+\b)[\h]+\("i => s"\1(",
    # Implicit declaration
    r"^\h*IMPLICIT(.*)"mi => s"",
    # returns
    Regex("\\b$(mask("lastreturn"))\\b") => "return",
    Regex("\\b$(mask("return"))\\b") => "return",
    r"\breturn\b$" => "return nothing",
    # main program
    r"^(\h*)PROGRAM\h*$"mi => s"\1PROGRAM NAMELESSPROGRAM",
    # rstrip()
    r"\h*$" => s"",
)


const headersprocessing = OrderedDict(
    # Reorganise functions and doc strings. This may be very project specific.
    # https://regex101.com/r/DAIHhl/1
    r"^(\h+)subroutine(\h+)(\w+)(\(([^)]*)\))(\h*#.*?|\h*)\n(#\h*\n)?#\h*function:\h*(.*?)#\h*\n"is =>
    SubstitutionString("\"\"\"\n    \\3(\\5)\n\n\\8\"\"\"\nfunction \\3(\\5)\n#\n"),
    # Simple subroutine
    r"^\h*(?:recursive\h+|)subroutine\h+(\w+)\h*\("mi => s"function \1(",
    r"^\h*(?:recursive\h+|)subroutine\h+(\w+)\h*"mi => s"function \1()",
    #r"^\h*subroutine\h+"mi => s"function ",
    r"^\h*program\h+(\w+)"mi => s"function \1()",
    r"^\h*block\h*data\h+(\w+)"mi => s"function \1()",
    r"^\h*module\h+"mi => s"module ",
    r"^\h*real(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
    r"^\h*complex(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
    r"^\h*double\h*precision\h*function\h+"mi => s"function ",
    r"^\h*integer(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
    r"^\h*character(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
    r"^\h*logical\h*function\h+"mi => s"function ",
    r"^\h*function\h+"mi => s"function ",
)


# Patterns to remove
const removal = [
    # Declarations
    r"\n\s*implicit none"i,
    r"\n\s*real,\s*external\s.*"i,
    r"\n\s*external\s.*"i,
    r"\n\s*intrinsic\s.*"i,
    r"\n\s*contains\s.*"i,
]


# Some routines for Tokenize.jl

restore(x) = join(untokenize.(x))

function matchedbrace(k::Tokens.Kind)
    if k == Tokens.LPAREN
        return Tokens.RPAREN
    elseif k == Tokens.RPAREN
        return Tokens.LPAREN
    elseif k == Tokens.LSQUARE
        return Tokens.RSQUARE
    elseif k == Tokens.RSQUARE
        return Tokens.LSQUARE
    elseif k == Tokens.LBRACE
        return Tokens.RBRACE
    elseif k == Tokens.RBRACE
        return Tokens.LBRACE
    end
end
matchedbrace(t::Token) = matchedbrace(t.kind)

function isleftbrace(k::Tokens.Kind)
    if k == Tokens.LPAREN || k == Tokens.LSQUARE || k == Tokens.LBRACE
        return true
    else
        return false
    end
end
isleftbrace(t::Token) = isleftbrace(t.kind)

function isrightbrace(k::Tokens.Kind)
    if k == Tokens.RPAREN || k == Tokens.RSQUARE || k == Tokens.RBRACE
        return true
    else
        return false
    end
end
isrightbrace(t::Token) = isrightbrace(t.kind)


skipbraces(ts::AbstractArray{<:Token}, i) =
    isleftbrace(ts[i]) ? nextind(ts, findmatchedbrace(ts, i)) :
        isrightbrace(ts[i]) ? prevind(ts, findmatchedbrace(ts, i)) : nothing

Base.nextind(ts::AbstractArray{<:Token}, ::Nothing) = nothing
Base.prevind(ts::AbstractArray{<:Token}, ::Nothing) = nothing
Base.thisind(ts::AbstractArray{<:Token}, ::Nothing) = nothing
Base.nextind(ts::AbstractArray{<:Token}, i::Integer) = i < length(ts) ? i + 1 : nothing
Base.prevind(ts::AbstractArray{<:Token}, i::Integer) = i > 1 ? i - 1 : nothing
Base.thisind(ts::AbstractArray{<:Token}, i::Integer) = 1 <= i <= length(ts) ? i : nothing

function findmatchedbrace(ts::AbstractArray{<:Token}, i)
    checkbounds(ts, i)
    if ts[i].kind == Tokens.LPAREN || ts[i].kind == Tokens.LSQUARE || ts[i].kind == Tokens.LBRACE
        left = ts[i].kind; right = matchedbrace(ts[i])
    else
        left = matchedbrace(ts[i]); right = ts[i].kind
    end
    inbraces = 1
    if ts[i].kind == left # '('
        for i = i+1:length(ts)
            if ts[i].kind == left
                inbraces += 1
            elseif inbraces == 1 && ts[i].kind == right
                return i
            elseif ts[i].kind == right
                inbraces -= 1
            end
        end
    elseif ts[i].kind == right # ')'
        for i = i-1:-1:1
            if ts[i].kind == right
                inbraces += 1
            elseif inbraces == 1 && ts[i].kind == left
                return i
            elseif ts[i].kind == left
                inbraces -= 1
            end
        end
    else
        @error "Unknown brace kind of Token: \"$(ts[i].kind)\" in \"$(ts[i])\""
    end
    return nothing
end


"""
$(SIGNATURES)
Return the position of requested token
There is loop up to C or EOL or # or )|]|}
All types of brackets aware.
"""
function findnexttoken(c, ts, i)
    checkbounds(ts, i)
    while i <= length(ts)
        t = ts[i]
        if t.kind == c
            return i
        elseif t.kind == Tokens.ENDMARKER || t.kind == Tokens.COMMENT || isrightbrace(t)
            break
        elseif isleftbrace(t)
            i = findmatchedbrace(ts, i) + 1
        else
            i += 1
        end
    end
    return nothing
end

# Some parsing routines


function picktoken(str, i)
    # may be there should be used `Base.Unicode.category_abbrev()` or `Base.Unicode.category_code()`
    len = ncodeunits(str)
    i = thisind(str, i)
    i > len && return 'e'
    c = str[i]
    if isspace(c) && c != '\n'
        return ' '
    elseif isdigit(c)
        return 'd'
    elseif isletter(c)
        return 'l'
    elseif c == '\t' || c == '\r'
        return ' '
    elseif c == '\n' && iscontinuedline(str, i)
        return ' '
    elseif c == '\n'
        return 'e'
    else
        return c
    end
end

function catchtokenvar(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    start = i
    while true
        i = nextind(str, i)
        if eos() || !(picktoken(str, i) in ('l', 'd', '_'))
            return start, prevind(str,i)
        end
    end
end
function catchtokenstring(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    @assert str[i] == '''
    start = nextind(str, i)
    while true
        i = nextind(str, i)
        eos() && return start, prevind(str, i)
        if picktoken(str, i) == ''' && picktoken(str, nextind(str, i)) == ''' # double ' -- escaped
            i = nextind(str, i)
        elseif picktoken(str, i) == '''
            return start, prevind(str, i)
        end
    end
end
" Return the position next to '(' and position of previous symbol to companion ')' "
function catchbrackets(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    @assert str[i] == '['
    inbraces = 1
    start = nextind(str, i)
    while true
        i = nextind(str, i)
        eos() && return start, prevind(str, i)
        if picktoken(str, i) == '['
            inbraces += 1
        elseif inbraces == 1 && picktoken(str, i) == ']'
            return start, prevind(str, i)
        elseif picktoken(str, i) == ']'
            inbraces -= 1
        end
    end
end
function catchbraces(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    @assert str[i] == '('
    inbraces = 1
    start = nextind(str, i)
    while true
        i = nextind(str, i)
        eos() && return start, prevind(str, i)
        if picktoken(str, i) == '('
            inbraces += 1
        elseif inbraces == 1 && picktoken(str, i) == ')'
            return start, prevind(str, i)
        elseif picktoken(str, i) == ')'
            inbraces -= 1
        end
    end
end
function catchtoken(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    start = i
    ttype = picktoken(str, i)
    while true
        i = nextind(str, i)
        if eos() || picktoken(str, i) != ttype
            return start, prevind(str,i)
        end
    end
end
function catchspaces(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return len, prevind(str,len)
    start = i
    while !eos() && (isspace(str[i]) && (str[i] != '\n' || iscontinuedline(str, i)))
        i = nextind(str, i)
    end
    return start, prevind(str,i)
end
function catchcomment(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return len, prevind(str,len)
    start = i
    while !eos() && str[i] ∉ ('\r', '\n')
    #while !eos() && str[i] != '\n'
        i = nextind(str, i)
    end
    return start, prevind(str, i)
end

function marktoken(str, i)
    len = ncodeunits(str); i = thisind(str, i)
    i > len && return len, prevind(str, len), len+1
    ttype = picktoken(str, i)
    if ttype == 'l'      # catch keywords and identifiers
        startpos, endpos = catchtokenvar(str, i)
        return startpos, endpos, nextind(str, endpos)
    elseif ttype == '''  # strings
        startpos, endpos = catchtokenstring(str, i)
        return startpos, endpos, nextind(str, nextind(str, endpos))
    elseif ttype == ' '  # catch spaces
        startpos, endpos = catchspaces(str, i)
        return startpos, endpos, nextind(str, endpos)
    elseif ttype == '#'  # catch comment
        startpos, endpos = catchcomment(str, i)
        return startpos, endpos, nextind(str, endpos)
    #elseif ttype == 'd'  # catch integer|float numbers
    elseif ttype == 'd'  # catch integer number
        startpos, endpos = catchtoken(str, i)
        return startpos, endpos, nextind(str, endpos)
    elseif ttype == '*'  # catch '**' operator
        startpos, endpos = catchtoken(str, i)
        return startpos, endpos, nextind(str, endpos)
    else                 # all other tokens should be singles
        return i, i, nextind(str, i)
    end
end

"""
take symbols till end of expression: ',' or ')'
"""
function taketokenexpr(str, i)
    eos() = i>len; len = ncodeunits(str); i1 = i = thisind(str, i)
    eos() && return str[len:len-1], len+1
    while true
        picktoken(str, skipspaces(str,i)) in (',', ')') && break
        i = skiptoken(str, skipspaces(str,i))
        eos() && break
    end
    return str[i1:prevind(str,i)], i
end

taketoken(str, i)     =  ( (p1,p2,p3) -> (str[p1:p2], p3) )(marktoken(str, i)...)
picknexttoken(str, i) =  picktoken(str, marktoken(str, i)[3])
skiptoken(str, i)     =  marktoken(str, i)[3]
skiptoken(c::AbstractChar, str, i) = picktoken(str,i) == c ? skiptoken(str,i) : i
skipspaces(str, i)    = catchspaces(str, i)[2] + 1
skipcomment(str, i)   = catchcomment(str, i)[2] + 1
skipbraces(str, i)    = nextind(str, nextind(str, catchbraces(str, i)[2]))
skipbrackets(str, i)  = nextind(str, nextind(str, catchbrackets(str, i)[2]))
isoneline(str)        = !('\n' in str || '\r' in str)
existind(str, i)      = thisind(str, min(max(1,i), ncodeunits(str)))

function skipwhitespaces(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return len+1
    while !eos() && str[i] in (' ', '\t', '\r')
        i = nextind(str, i)
    end
    return i
end

"""
Return the last position of token of requested char
"""
function skipupto(c, str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    while true
        eos() && return len+1
        if str[i] == ''' != c # skip string
            i = skiptoken(str, i)
        elseif str[i] == c == '\n'
            if iscontinuedline(str, i)
                i = nextind(str, i)
            else
                break
            end
        elseif str[i] == '#' != c
            i = skipcomment(str, i)
        elseif str[i] == c
            break
        else
            i = nextind(str, i)
        end
    end
    return i
end

function startoflineind(str, i)
    i = thisind(str, i)
    i < 1 && return 0
    while i>1 && str[prevind(str, i)] !='\n'
        i = prevind(str, i)
    end
    return i
end
function endoflineind(str, i)
    i = thisind(str, i)
    i > ncodeunits(str) && return ncodeunits(str)+1
    while i < ncodeunits(str) && str[i] !='\n'
        i = nextind(str, i)
    end
    return i
end
thislinerange(str, i) = UnitRange(startoflineind(str, i):endoflineind(str, i))
function nextlinerange(str, i)
    i = thisind(str, i)
    firstip1 = nextind(str, endoflineind(str, i))
    if firstip1 >= ncodeunits(str)
        return UnitRange(firstip1:ncodeunits(str))
    else
        return UnitRange(firstip1:endoflineind(str, firstip1))
    end
end
function prevlinerange(str, i)
    i = thisind(str, i)
    lastim1 = prevind(str, startoflineind(str, i))
    if lastim1 < 1
        return UnitRange(1:lastim1)
    else
        return UnitRange(startoflineind(str, lastim1):lastim1)
    end
end

function continuedlinesrange(str, i)
    i = thisind(str, i)
    firsti = startoflineind(str, i)
    while iscontinuedlines(str[prevlinerange(str, firsti)], str[thislinerange(str, firsti)])
        firsti = startoflineind(str, prevind(str, firsti))
    end
    lasti = endoflineind(str, i)
    while iscontinuedlines(str[thislinerange(str, lasti)], str[nextlinerange(str, lasti)])
        lasti = endoflineind(str, nextind(str, lasti))
    end
    return UnitRange(firsti, lasti)
end

function printlines(lines)
    println()
    #for l in lines println(l) end
    for (i,l) in enumerate(lines) print("[$i]: \"$l\"\n") end
end


end # of module FortranTransplier

isinteractive() || FortranTransplier.CLI.main(ARGS)
