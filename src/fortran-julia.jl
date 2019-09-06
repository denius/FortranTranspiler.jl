#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

using DataStructures, Printf, JuliaFormatter

#=

# KNOWN ISSUES
    * julia's `for` loop counter is local scope variable, unlike fortran.
      Thus the codes that use the value of counter after the loop will be broken
      and should be fixed via `while` loop.
    * due to julia's GC, some variables initialized inside the loop/if blocks may disappear
      upon they exit these blocks. Consider initilizing them with the appropriate value in
      the initial part of the function code.
    * all strings in fortran are `Vector{Char}` and should stay same type to preserve
      assignment statements. (use `join` and `collect` to convert to strings and back)
    * whitespaces matter despite of the old fortran, which can ignore them all:
    ** in julia whitespaces between name of functions/arrays and braces are not allowed.
    * julia can't propagate back changed values of scalars via functions args, unlike fortran.
      Thus such changed scalars should be returned via functions `return` statement.
    * some unserviceable comments are cutted off (eg. in expanded lines continuations).
    * this script is mostly for fixed-form fortran and free-form is not tested yet.
    * in the `DATA` statement can occur uncatched repetitions like `DATA SOMEARRAY/8*0,1,2/`
    * `READ`/`WRITE` is ugly
    * can't substitude dynamic - created at runtime FORMAT string
    * `FORMAT` is incorrect sometimes
    * fortran functions returns result via assignment the returned to function name,
      in julia it is the reassignment, thus should be fixed

# TODO
    return intent(out):
    LSAME
    XERBLA
    P_VALUE
    RWK/IWK
    M_GETMEM
=#

function printusage()
    println("fortran-julia.jl

    This julia script converts fortran 77, 90 code into julia.
    It uses naive regex replacements to do as much as possible,
    but the output may need further cleanup.
    Based on https://gist.github.com/rafaqz/fede683a3e853f36c9b367471fde2f56

    Usage:
        fortran-julia.jl -h | --help
        fortran-julia.jl [--lowercase | --uppercase] ... [--] <filename1.f> <filename2.f> ...

    Samples:
        fortran-julia.jl *.f*
        fortran-julia.jl **/*.[fF]*

    Options:
        -h, --help       Show this screen.
        --version        Show version.
        -q, --quiet      Suppress all console output.
        -v, --verbose    Be verbose.
        -v, --verbose    Be more verbose.
        --uppercase      Convert all identifiers to upper case.
        --lowercase      Convert all identifiers to lower case.
        --greeks         Replace the greeks letters names thats starts the vars names
                         with corresponding unicode symbols, like DELTA1 -> δ1.
        --subscripts     Replace tail suffixes in vars names like SOMEVAR_? with unicode subscripts, if exist.
        --               The rest of the args is only the filenames.
        --formatting     Try to format with JuliaFormatter package.
        --returnnothing
        --packarrays
        --strings
")
end


function main(args)
    args, fnames = parse_args(args)
    casetransform = args["--uppercase"] ? uppercase :
                    args["--lowercase"] ? lowercase : identity

    for fname in fnames
        isfile(fname) || continue
        args["--verbose"] && print("$(fname):\n")
        code = open(fname) |> read |> String
        result = convertfromfortran(code, casetransform=casetransform, quiet=args["--quiet"],
                                    verbose=args["--verbose"], veryverbose=args["--veryverbose"])
        write(splitext(fname)[1] * ".jl", result)
        if args["--formatting"]
            try
                result = JuliaFormatter.format_text(result, margin = 192)
                write(splitext(fname)[1] * ".jl", result)
            catch
                @info("$(splitext(fname)[1] * ".jl")\nFORMATTING IS NOT SUCCESSFULL\n")
            end
        end
    end
    return nothing
end

isfixedformfortran = true

function convertfromfortran(code; casetransform=identity, quiet=true,
                            verbose=false, veryverbose=false)

    # should be '\n' newlines only
    for rx in newlinereplacements
        code = replace(code, rx)
    end

    # expand tabs
    code = replace(code, r"^[ ]{0,5}\t"m => "      ")
    code = replace(code, r"\t"m => "  ")

    # replace 'PRINT' => 'WRITE'
    code = replace(code, r"(print)\h*([*]|'\([^\n\)]+\)')\h*,"mi => s"write(*,\2)")

    global isfixedformfortran = !occursin(r"^[^cC*!\n][^\n]{0,3}[^\d\h\n]"m, code)

    code = savespecialsymbols(code)

    # convert comments
    code = convertfortrancomments(code)

    # mark lines of code occupied by each subroutine
    subs, subnames = markbysubroutine(code)

    alllines = splitonlines(code)

    # converted code will be stored in string `result`
    result = "using Printf, FortranFiles, Parameters, OffsetArrays"

    for i = 1:length(subs)-1

        code = foldl((a,b) -> a*'\n'*b, view(alllines, subs[i]:subs[i+1]-1))
        verbose && print("[$(subs[i]):$(subs[i+1]-1)] $(strip(subnames[i]))\n")

        scalars, arrays = collectvars(code)
        veryverbose && length(scalars) > 0 && println("    scalars : $scalars")
        veryverbose && length(arrays) > 0  && println("    arrays  : $arrays")

        # replace array's brackets with square braces
        code = replacearraysbrackets(code, arrays)

        commons = collectcommon(code)
        veryverbose && length(commons) > 0 && println("    COMMONs : $commons")

        code, formats = collectformat(code)
        veryverbose && length(formats) > 0 && println("    FORMATs : $formats")
        #formats = map(convertformat, formats)

        # comment out all vars declarations
        code = commentoutdeclarations(code)

        code = processcommon(code, commons, arrays)

        code = processdostatements(code)

        dolabels, gotolabels = collectlabels(code)
        code = replacedocontinue(code, dolabels, gotolabels)
        code = processconditionalgotos(code)

        code = processiostatements(code, formats)

        code = processifstatements(code)

        #write("test.jl", code)
        # process multiline replacements
        for rx in multilinereplacements
            code = replace(code, rx)
        end

        #write("test.jl", code)
        code = processselectcase(code)

        code = processparameters(code)

        code = processdata(code)

        lines = splitonlines(code)

        lines, comments = splitoncomment(lines)

        lines = map(casetransform, lines)

        for rx in singlewordreplacements
            lines = map(a->replace(a, rx), lines)
        end

        lines = processlinescontinuation(lines)

        # straight syntax conversions
        for rx in replacements
            lines = map(a->replace(a, rx), lines)
        end
        #write("test.jl", foldl((a,b) -> a*'\n'*b, lines))

        # append comments back
        lines = map(*, lines, comments)

        #lines = includeformat(lines, formats)

        # concat string lines back together and
        # create functions declaration with custom header
        code = foldl((a,b) -> a*'\n'*b, lines)
        for rx in headersprocessing
            code = replace(code, rx)
        end

        # removing unnecessaries
        for rx in removal
            code = replace(code, rx => "")
        end

        code = restorespecialsymbols(code)

        try
            Meta.parse(code, 1)
        catch e
            quiet || @error("$(strip(subnames[i]))\n$e\n")
            #quiet || @error("$(strip(subnames[i]))\n$e\n$(stacktrace(catch_backtrace()))\n")
        end

        # concat all subroutines together back
        result = result * '\n' * code

    end # of loop over subroutines

    #result = restorespecialsymbols(result)

    return result
end

function parse_args(lines)
    function therecanbeonlyone!(args, v1, v2)
        args[v1] = args[v1] || args[v2]
        delete!(args, v2)
    end

    fnames = Vector{String}()
    args = OrderedDict{String, Any}()
    args["--quiet"]         = args["-q"]  = false
    args["--verbose"]       = args["-v"]  = false
    args["--veryverbose"]   = args["-vv"] = false
    args["--uppercase"]     = false
    args["--lowercase"]     = false
    args["--greeks"]        = false
    args["--subscripts"]    = false
    args["--formatting"]    = false
    args["--returnnothing"] = false
    args["--packarrays"]    = false
    args["--strings"]       = false
    #args["--"]  = false

    while length(lines) > 0
        if first(args) == "--help" || first(args) == "-h"
            printusage()
            exit(0)
        elseif haskey(args, first(lines))
            args[first(lines)] += 1
            popfirst!(lines)
        elseif first(args) == "--" # end of options
            popfirst!(lines)
            break
        elseif occursin(r"^--\w\w+$", first(lines)) || occursin(r"^-\w$", first(lines))
            @error("unrecognized option in ARGS: '$(first(lines))'")
            exit(0)
        else
            push!(fnames, first(lines))
            popfirst!(lines)
        end
    end
    while length(lines) > 0
        push!(fnames, first(lines))
        popfirst!(lines)
    end

    if args["--verbose"] == 2 || args["-v"] == 2 || args["-vv"] == 1
        args["--veryverbose"] = args["--verbose"] = args["-vv"] = args["-v"] = true
    end
    for (k,v) in args
        !isa(args[k], Bool) && args[k] == 1 && (args[k] = true)
    end

    therecanbeonlyone!(args, "--veryverbose", "-vv")
    therecanbeonlyone!(args, "--verbose", "-v")
    therecanbeonlyone!(args, "--quiet", "-q")

    if args["--veryverbose"]
        for (key,val) in args @printf("  %-13s  =>  %-5s\n", key, repr(val)); end
        println("  fnames: $fnames")
    end

    return args, fnames
end

function markbysubroutine(code)
    lines, comments = splitoncomment(code)

    rx = r"^(\h*)((?:[\w*\h]+\h+)function|(?:recursive\h+|)subroutine|program|block\h*data|module)\h*.*$"mi
    marked = Vector{Int}(undef, 0)
    names  = Vector{String}(undef, 0)
    for i in axes(lines,1)
        if occursin(rx, lines[i])
            push!(marked, i)
            push!(names, match(rx, lines[i]).match)
        end
    end

    rx = r"^\h*end\h*(?!do|if)"i
    #rx = r"^\h*end\h*(?:(function|subroutine|program|module|)\h*\w*)"i
    marked = sort(collect(marked))
    for i = length(marked):-1:2
        for j = marked[i]:-1:marked[i-1]
            if occursin(rx, lines[j])
                marked[i] = j+1
                break
            end
        end
    end
    marked[1] = 1
    push!(marked, length(lines)+1)

    return marked, names
end


function convertfortrancomments(code)
    rx = r"(^[cC*]|!)([^\n]*)"m
    code1 = code isa AbstractVector ? foldl((a,b) -> a*'\n'*b, code) : code
    code1 = replace(code1, rx => s"#\2")
    code1 = replace(code1, r"#=" => "# ") # accident multiline comment
    return code isa AbstractVector ? splitonlines(code1) : code1
end

function stripcomments(code)
    # https://regexr.com/4j4vi via
    # https://stackoverflow.com/questions/9203774/regular-expression-for-comments-but-not-within-a-string-not-in-another-conta
    juliacomments = r"(#=([^*]|\*(?!#))*?=#|#[^\"\n\r]*(?:\"[^\"\n\r]*\"[^\"\n\r]*)*[\r\n])(?=[^\"]*(?:\"[^\"]*\"[^\"]*)*$)"m
    dlms = [r"#"m]
    lines = splitonlines(code)
    for i in axes(lines,1)
        for dlm in dlms
            lines[i] = split(lines[i], dlm, limit=2, keepempty=true)[1]
        end
        lines[i] = rstrip(lines[i])
    end
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function splitoncomment(code)
    rx = r"\h*#.*$"
    lines = splitonlines(code)
    comments  = ["" for i=axes(lines,1)]
    for i in axes(lines,1)
        (m = match(rx, lines[i])) != nothing && (comments[i] = m.match)
        lines[i] = replace(lines[i], rx => s"")
    end
    # cut out all whitespaces not masked by the comments
    rx = r"\h*$"
    lines = map(a->replace(a, rx => s""), lines)
    return lines, comments
end

function processlinescontinuation(lines::AbstractVector)
    # `lines` should be lines without comments
    rx1 = r"^(.*[^,&])(&|)$"
    rx2 = r"^([ ]{5}[^ ]\h*)([+*\/,-]|==|<=|>=|!=|<|>|&&|\|\||=)(.*)$"
    # move arithmetic operator from start of continuator to tail of previous line
    for i = 2:length(lines)
        if (m2 = match(rx2, lines[i])) != nothing
            if occursin(rx1, lines[i-1])
                lines[i-1] = replace(lines[i-1], rx1 => SubstitutionString("\\1"* m2.captures[2] *"\\2"))
                lines[i]   = replace(lines[i],   rx2 => s"\1\3")
            end
        end
    end
    rx3 = r"^([ ]{5}[^ ]|[ ]{6})(.*)$"
    rx4 = r"^(.*)(&|)$"
    lines = map(a->replace(a, rx3 => s"      \2"), lines)
    lines = map(a->replace(a, rx4 => s"\1"), lines)
    return lines
end

function concatlinecontinuation(code)
    # Note: trailing comments are not allowed here
    replacements = OrderedDict(
        r"\n     [.+&$\w]\h*" => "",  # fixed-form startline continuations
        r"&\n"                => ""   # free-form continuations
    )
    code1 = deepcopy(code)
    for rx in replacements
        code1 = replace(code1, rx)
    end
    return code1
end

function splitonlines(code)
    if code isa AbstractVector
        lines = deepcopy(code)
    else
        lines = split(code, '\n')
    end
    return lines
end

function printlines(lines)
    print("\n")
    #for (i,l) in enumerate(lines) print("$l\n") end
    for (i,l) in enumerate(lines) print("line[$i]: \"$l\"\n") end
end

function commentoutdeclarations(code)

    patterns = [
        # Match declarations
        r"^\h*common\h*\/\h*\w+\h*\/\h*"mi,
        r"^\h*dimension\h+"mi,
        r"^\h*implicit\h*\w+"mi,
        r"^\h*real.*::"mi,
        r"^\h*real\h+(?!.*function)"mi, # https://regex101.com/r/CFiiOx/1
        r"^\h*real\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*double\h*precision.*::"mi,
        r"^\h*double\h*precision\h+(?!.*function)"mi,
        r"^\h*complex.*::"mi,
        r"^\h*complex\h+(?!.*function)"mi,
        r"^\h*complex\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer.*::"mi,
        r"^\h*integer\h+(?!.*function)"mi,
        r"^\h*integer\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer\h*\(\h*KIND\h*=\h*\w+\h*\)\h*(?!.*function)"mi,
        r"^\h*character.*::"mi,
        r"^\h*character\h+(?!.*function)"mi,
        r"^\h*character\h*\(\h*\*\h*\)\h*(?!.*function)"mi,
        r"^\h*character\*\h*\(\h*\*\h*\)\h*(?!.*function)"mi,
        r"^\h*character\h*\(\h*\d+\h*\)\h*(?!.*function)"mi,
        r"^\h*character\*\h*\(\h*\d+\h*\)\h*(?!.*function)"mi,
        r"^\h*character\*\h*\([\h\w_*+-]+\)\h*(?!.*function)"mi,
        r"^\h*character\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*external\h+"mi,
        r"^\h*logical.*::"mi,
        r"^\h*logical\h+(?!.*function)"mi,
        r"^\h*entry\h+.*"mi,
        r"^\h+\d+\h+format\h*\("mi
    ]

    lines = splitonlines(code)
    for i in markbypattern(patterns, code)
        lines[i] = "#" * lines[i]
    end

    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function collectvars(code)
    code1 = stripcomments(code)
    code1 = concatlinecontinuation(code1)
    lines = splitonlines(code1)

    patterns = [
        # Match declarations
        r"^\h*common\h*/\h*\w+\h*/\h*"mi,
        r"^\h*dimension\h+"mi,
        r"^\h*real\h+(?!.*function)"mi,
        r"^\h*real\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*real\h*\(\h*KIND\h*=\h*\w+\h*\)\h*(?!.*function)"mi,
        r"^\h*double\h+precision\h+(?!.*function)"mi,
        r"^\h*complex\h+(?!.*function)"mi,
        r"^\h*complex\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer\h+(?!.*function)"mi,
        r"^\h*integer\*[0-9]+\h+(?!.*function)"mi,
        r"^\h*integer\h*\(\h*KIND\h*=\h*\w+\h*\)\h*(?!.*function)"mi,
        r"^\h*logical\h+(?!.*function)"mi,
        r"^\h*character\h+"mi,
        r"^\h*character\*\h+"mi
    ]

    matched = Vector{String}()
    for rx in patterns
        for l in lines
            occursin(rx, l) && push!(matched, replace(replace(l, rx => s""), r" " => s""))
        end
    end

    matchedvars = map(a->replace(a, r"\*\(" => s"("), deepcopy(matched))  # drop '*' in 'character*()'
    matchedvars = map(a->replace(a, r"\*[\d]+" => s""), matchedvars)      # clean '*7' in 'A*7'
    matchedvars = map(a->replace(a, r"(\(((?>[^\(\)]++|(?1))*)\))" => s""), matchedvars)    # clean braces
    matchedvars = map(a->replace(a, r"(\/((?>[^\/]*|(?1))*)\/)" => s""), matchedvars)  # clean /.../ -- statics

    # catch only arrays names
    for i in 1:length(matched)
        matched[i] = replace(matched[i], r"\*([\d]+)" => s"(\1)") # character A*7 -> A(7)
        matched[i] = replace(matched[i], r"(\(((?>[^\(\)]++|(?1))*)\))" => s"()")   # clean inside braces
        # https://regex101.com/r/weQOSe/4
        matched[i] = replace(matched[i], r"(\/((?>[^\/]*|(?1))*)\/)" => s"")   # clean /.../ -- statics
        matched[i] = replace(matched[i], r"[^(),]+," => s"")      # drop non-arrays (without braces)
        matched[i] = replace(matched[i], r",[^(),]+$"m => s"")    # drop last non-array (without braces)
        matched[i] = replace(matched[i], r"^[^(),]+$"m => s"")    # drop single non-array (without braces)
        matched[i] = replace(matched[i], r"\*\(" => s"(")         # drop '*' in character*()
    end

    # add undoubted arrays
    # Note: to any fortran 'CHARACTER' can be applied get index operator: 'C(:1)'
    patterns = [r"^\h*.*dimension\h*\(.+\)\h*::\h*"mi,
               r"^\h*character\h+"mi,
               r"^\h*character\*\h+"mi,
               r"^\h*character\*\d+\h+"mi,
               r"^\h*character\(\d+\)\h*"mi,
               r"^\h*character\(\*\)\h*"mi,
               r"^\h*character\*\(\*\)\h*"mi,
               r"^\h*character\*\([\h\w_*+-]+\)\h*"mi ]
    for rx in patterns
        for l in lines
            if occursin(rx, l)
                line = replace(l, rx => "")
                line = replace(line, r" " => s"")
                line = replace(line, r"(\/((?>[^\/]*|(?1))*)\/)" => s"")   # clean /.../ -- statics
                line = replace(line, r"(\(((?>[^()]++|(?1))*)\))" => s"")  # braces and its containment
                line = replace(line, r"\*" => s"")
                push!(matched, line)
            end
        end
    end

    arrays = Vector{String}()
    for l in matched
        length(l) > 0 && append!(arrays, split(l, ","))
    end
    for i in 1:length(arrays)
        arrays[i] = replace(arrays[i], r"^(\w+)\(.*\)$"m => s"\1") # drop braces
    end
    arrays = unique(arrays)

    # all other will be scalars
    scalars = Vector{String}()
    for l in matchedvars
        length(l) > 0 && append!(scalars, split(l, ","))
    end
    scalars = unique(vcat(scalars, arrays))
    scalars = scalars[findall(a->a ∉ arrays, scalars)]

    return scalars, arrays
end

function collectcommon(code)
    code1 = stripcomments(code)
    code1 = concatlinecontinuation(code1)
    lines = splitonlines(code1)

    rx = r"^\h*common\h*\/\h*(\w+)\h*\/\h*(.+)"mi
    matched = Dict{String,String}()
    for i in axes(lines,1)
        m = match(rx, lines[i])
        if !isnothing(m) matched[m.captures[1]] = m.captures[2] end
    end
    for i in keys(matched)
        matched[i] = replace(matched[i], r" " => s"")               # drop spaces
        matched[i] = replace(matched[i], r"\[[^\]]+\]" => s"")      # drop braces
    end
    return matched
end

function iscontinued(line1, line2)
    if occursin(r"^[^#&]*&\h*(#.*|)$"m, line1) ||  # free-form continuation
       occursin(r"^     [.+&$\w].*$"m, line2)  #=||  # fixed-form startline continuation
       occursin(r"^[#]", line2)                  =#  # commented out line
        return true
    else
        return false
    end
end

function concatlines(line1, line2)
    # comments will be skipped
    return replace(line1, r"^([^#&]*)(&|)(\h*#.*|\h*)$"m => s"\1") *  # free-form continuation
           replace(line2, r"^(     [.+&$\w]|)(.*)$"m => s"\2")        # fixed-form startline continuation
end

function processdostatements(code)
    lines = splitonlines(code)
    rx = r"^(\h*\d*\h*|)do\h*(?:\d+\h*,?|)\h*[\w_]+\h*=.*$"i
    i = 1
    while true
        if occursin(rx, lines[i])
            # cut out this line with all continuations
            head = replace(lines[i], rx=>s"\1")
            str = "" * lines[i]
            while iscontinued(lines[i], i<length(lines) ? lines[i+1] : "")
                splice!(lines, i)
                str = concatlines(str,lines[i])
            end
            splice!(lines, i)
            # replace with new julia's `for` statement
            success,s1,s2,s3,s4,s5,s6,sc = parsedostatement(str)
            length(s2) > 0 && (s2 *= " ")
            length(s6) > 0 && (s6 *= ":")
            success && (str = head * "for " * s2 * s3 * " = " * s4 * ":" * s6 * s5 * sc)
            insert!(lines, i, str)
        end
        i += 1
        i > length(lines) && break
    end
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function parsedostatement(line)
    # "DO10,I=10,A(1,1)*TONUM('10,000'),-2" is the valid fortran 'DO' statement
    # https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vn8c/index.html

    reflabel   = ""
    label      = ""
    countervar = ""
    position   = 1
    inbraces   = 0
    nextarg    = 0
    starts     = [0,0,0]

    while true
        ttype = picktoken(line, position)
        #print("_$ttype")
        if ttype == 'e' || ttype == '#'
            break
        elseif ttype == ' '
            position = skipspaces(line, position)
        #elseif ttype == '#'
        #    break #position = skipspaces(line, position)
        elseif nextarg == 0 && ttype == 'd'
            reflabel, position = taketoken(line, position)
        elseif nextarg == 0 && ttype == 'l'
            lex, position = taketoken(line, skipspaces(line, position))
            if lowercase(lex) == "do"
                lex, position = taketoken(line, skipspaces(line, position))
                if isdigit(lex[1]) # LABEL
                    label = lex
                    position = skiptoken(',', line, skipspaces(line, position))
                    lex, position = taketoken(line, skipspaces(line, position))
                end
                countervar = lex
            elseif occursin(r"^do\d+$"i, lex)
                label = replace(lex, r"^do(\d+)$"i=>s"\1")
                position = skiptoken(',', line, skipspaces(line, position))
                countervar, position = taketoken(line, skipspaces(line, position))
            elseif occursin(r"^do\d+[a-z].*$"i, lex)
                label = replace(lex, r"^do(\d+)([a-zA-Z].*)$"i=>s"\1")
                countervar = replace(lex, r"^do(\d+)([a-zA-Z].*)$"i=>s"\2")
            elseif occursin(r"^do[a-z].*$"i, lex)
                countervar = replace(lex, r"^do([a-z].*)$"i=>s"\1")
            else # something goes wrong
                return false,"","","","","","",""
            end
            position = skipspaces(line, position)
            if picktoken(line, position) != '='
                # it is not the 'DO' statement
                return false,"","","","","","",""
            end
            position = skiptoken(line, position)
            starts[nextarg+=1] = position
        elseif ttype == '('
            inbraces += 1
            position = skiptoken(line, position)
        elseif ttype == ')'
            inbraces -= 1
            position = skiptoken(line, position)
        elseif ttype == ','
            position = skiptoken(line, position)
            if inbraces == 0
                starts[nextarg+=1] = position
            end
        else
            position = skiptoken(line, position)
        end
    end

    if starts[1] == 0 || starts[2] == 0
        return false,"","","","","","",""
    elseif starts[3] == 0
        starts[3] = position + 1
    end

    # cross our fingers to have simple ascii without multibyte characters
    return true, reflabel, label, countervar,
           strip(line[starts[1]:starts[2]-2]),
           strip(line[starts[2]:starts[3]-2]),
           strip(line[starts[3]:position-1]),
           line[position:end]                   # trailing comment
end

function processifstatements(code)
    code1 = code isa AbstractVector ? foldl((a,b) -> a*'\n'*b, code) : deepcopy(code)

    rx = r"^(\h{0,4}\d*\h*|)(else|)\h*if\h*\(.*"mi
    for (i,m) in enumerate(reverse(collect(eachmatch(rx, code1))))
        o = m.offset
        s1 = replace(m.match, rx => s"\1")
        s2 = lowercase(replace(m.match, rx => s"\2"))
        iflength, wothen, ifcond, ifexec, afterthencomment = parseifstatement(code1, o)
        length(ifexec) > 0 && isletter(ifexec[1]) && (ifexec = " " * ifexec)
        length(afterthencomment) > 0 && (afterthencomment = " " * afterthencomment)
        #@show o, iflength, ifcond, ifexec, afterthencomment
        if wothen == true
            str = s1 * s2 * "if_iZjAcpPokM (" * strip(ifcond) * ")"
            if isoneline(ifcond) && isoneline(ifexec)
                # insert "end" before comment or at tail
                str *= replace(ifexec, r"^\h*([^#]*)\h*(#.*|)$" => s" \1 end \2")
            else
                str *= ifexec * "\n" * " "^length(s1) * "end "
            end
        elseif wothen == false && picktoken(ifcond, skipwhitespaces(ifcond, 1)) in ('#','\n',' ')
            str = s1 * s2 * "if_iZjAcpPokM (" * ifcond * ")" * afterthencomment
        elseif wothen == false
            str = s1 * s2 * "if_iZjAcpPokM " * strip(ifcond) * afterthencomment
        else
            str = ""
        end
        code1 = code1[1:prevind(code1,o)] * str * code1[thisind(code1,o+iflength):end]
    end

    code1 = replace(code1, r"if_iZjAcpPokM" => "if")
    return code isa AbstractVector ? split(code1, '\n') : code1
end

function parseifstatement(str, position)
    position0 = position
    reflabel = ""
    afterthencomment = ""
    iftaken = false
    wothen = nothing
    startcond = startexec = 0
    endcond = endexec = -1
    inbraces  = 0

    while true
        ttype = picktoken(str, position)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            position = skipspaces(str, position)
        elseif !iftaken && ttype == 'd'
            reflabel, position = taketoken(str, position)
        elseif !iftaken && ttype == 'l'
            lex, position = taketoken(str, position)
            if lowercase(lex) == "else"
                lex, position = taketoken(str, skipspaces(str, position))
            end
            if lowercase(lex) == "if" || lowercase(lex) == "elseif"
                iftaken = true
            else
                # that is not an "if" statement
                return 0, wothen, "", "", ""
            end
        elseif iftaken && inbraces == 0 && ttype == 'l'
            lex, position = taketoken(str, position)
            if lowercase(lex) == "then"
                wothen = false
                if picktoken(str, skipwhitespaces(str, position)) == '#'
                    afterthencomment = str[position:prevind(str, (position = skipcomment(str, skipwhitespaces(str,position)))...)]
                else
                    position = skipwhitespaces(str, position)
                    #position = skipspaces(str, position)
                end
                break
            else
                wothen = true
                endexec = position = skipupto('\n', str, position)
                endexec = prevind(str, endexec)
            end
        elseif ttype == '('
            inbraces += 1
            position = skiptoken(str, position)
            startcond == 0 && (startcond = position)
        elseif ttype == ')'
            inbraces -= 1
            if inbraces == 0 && endcond == -1
                endcond = prevind(str, position)
                startexec = nextind(str, position)
            end
            position = skiptoken(str, position)
        else
            position = skiptoken(str, position)
        end
    end

    return ncodeunits(str[position0:prevind(str,position)]), wothen,
           str[startcond:endcond], str[startexec:endexec], afterthencomment
end

function processiostatements(code, formats)
    code1 = code isa AbstractVector ? foldl((a,b) -> a*'\n'*b, code) : deepcopy(code)
    rxhead = r"(\bread\b|\bwrite\b)\h*\("mi
    rxargs = r"^((?:.*[\n])*.*?)(\h*#[^\n]*)"m  # https://regex101.com/r/PqvHQD/1
    for (i,m) in enumerate(reverse(collect(eachmatch(rxhead, code1))))
        o = m.offset
        readwrite = lowercase(replace(m.match, rxhead => s"\1"))
        iolength, io, label, fmt, pmt, paramsstr, args = parsereadwrite(code1, o)
        #@show iolength, io, label, fmt, pmt, paramsstr, args
        length(label) > 0 && haskey(formats, label) && (fmt = formats[label])
        #length(fmt) > 0 && @show fmt
        fmt = convertformat(fmt)
        if  readwrite == "read" && (io == "5" || io == "*")
            io = "stdin"
        elseif  readwrite == "write" && (io == "6" || io == "*")
            io = "stdout"
        end
        args = occursin(rxargs, args) ? replace(args, rxargs=>s"\1)\2") : rstrip(args)*')'
        str  =  readwrite == "read" ? "READ(" :
                fmt       == ""    ? "println(" : "at_iZjAcpPokMprintf("
        if str == "at_iZjAcpPokMprintf("
            fmt = occursin(r"\$$", fmt) ? replace(fmt, @r_str("\\\$\$") => "") : fmt*"\\n"
        end

        #length(fmt) > 0 && @show fmt
        fmt != "" && (fmt = '"' * fmt * "\", " )
        str *= io * ", " * fmt * args
        code1 = code1[1:prevind(code1,o)] * str * code1[thisind(code1,o+iolength):end]
    end
    return code isa AbstractVector ? split(code1, '\n') : code1
end

function parsereadwrite(str, position)
    position0 = position
    label = ""
    fmt = ""
    IU = ""
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
        ttype = picktoken(str, position)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            position = skipspaces(str, position)
        elseif !cmdtaken && ttype == 'l'
            lex, position = taketoken(str, position)
            if lowercase(lex) == "read" || lowercase(lex) == "write"
                cmdtaken = true
                position = skipspaces(str, position)
            else
                # that is not an 'READ/WRITE' statement
                return 0, "", "", "", "", "", "", ""
            end
        elseif cmdtaken && inparams && ttype == 'l'
            lex, position = taketoken(str, position)
            LEX = uppercase(lex)
            if LEX == "FMT" && picktoken(str, skipspaces(str, position)) == '='
                label, position = taketoken(str, skipspaces(str, skiptoken(str, skipspaces(str, position))))
            elseif (LEX == "ERR" || LEX == "END") && picktoken(str, skipspaces(str, position)) == '='
                val, position = taketoken(str, skipspaces(str, skiptoken(str, skipspaces(str, position))))
                params[uppercase(lex)] = val
            elseif nextparam == 1
                IU = lex
            else
                @warn("parsereadwrite(): $(@__LINE__): unknown FORMAT-parameter = \"$lex\"" *
                      " in FORTRAN line:\n\"$(strip(thisline(str, position)))\"\n")
            end
        elseif cmdtaken && inparams && ttype == 'd' && nextparam == 1
            IU, position = taketoken(str, position)
        elseif cmdtaken && inparams && ttype == 'd'
            label, position = taketoken(str, position)
        elseif cmdtaken && inparams && ttype == '*' && nextparam == 1
            IU, position = taketoken(str, position)
        elseif cmdtaken && inparams && ttype == ''' || ttype == '*'
            fmt, position = taketoken(str, position)
        elseif ttype == '('
            inbraces += 1
            position = skiptoken(str, position)
            if length(startsparam) == 0
                push!(startsparam, position)
                inparams = true
            end
        elseif ttype == ')'
            inbraces -= 1
            if inbraces == 0 && inparams
                push!(endsparam, prevind(str, position))
                startargs = nextind(str, position)
                inparams = false
                endargs = position = skipupto('\n', str, position)
                endargs = prevind(str, endargs)
                break
            end
            position = skiptoken(str, position)
        elseif ttype == ',' && inparams && inbraces == 1
            push!(endsparam, prevind(str, position))
            position = skiptoken(str, position)
            push!(startsparam, thisind(str, position))
            nextparam += 1
        else
            position = skiptoken(str, position)
        end
    end

    IU = strip(str[startsparam[1]:endsparam[1]])

    if fmt == "*"
        fmt = ""
    else
        fmt = replace(fmt, r"^\((.*)\)$"=>s"\1") # it is unenclose "()"
    end
    fmt = foldl((a,b) -> a*','*b, map(strip, split(fmt, ",")))

    return ncodeunits(str[position0:prevind(str,position)]), IU, label, fmt, params,
           str[startsparam[1]:endsparam[end]], str[startargs:endargs]
end

function splitformat(str)
    position = 1
    lexemstarts = Int[1]
    lexemends = Int[]
    inbraces  = 0

    while true
        ttype = picktoken(str, position)
        #print("_$ttype")
        if ttype == 'e' || ttype == '\n'
            break
        elseif ttype == ' '
            position = skipwhitespaces(str, position)
        elseif ttype == ',' && inbraces == 0
            push!(lexemends, prevind(str, position))
            position = skiptoken(str, position)
            push!(lexemstarts, thisind(str, position))
        elseif ttype == '/' && inbraces == 0
            if position > lexemstarts[end]
                push!(lexemends, prevind(str, position))
                push!(lexemstarts, thisind(str, position))
            end
            push!(lexemends, thisind(str, position))
            position = skiptoken(str, position)
            push!(lexemstarts, thisind(str, position))
        elseif ttype == '('
            inbraces += 1
            position = skiptoken(str, position)
        elseif ttype == ')'
            inbraces -= 1
            position = skiptoken(str, position)
        else
            position = skiptoken(str, position)
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

function picktoken(str, i)
    islbrace(c::AbstractChar) = c=='(' || c=='['
    isrbrace(c::AbstractChar) = c==')' || c==']'
    len = ncodeunits(str)
    i = thisind(str, i)
    i > len && return 'e'
    c = str[i]
    if isspace(c) && c != '\n'
        return ' '
    elseif islbrace(c)
        return '('
    elseif isrbrace(c)
        return ')'
    elseif isdigit(c)
        return 'd'
    elseif isletter(c)
        return 'l'
    elseif c == '\n' && catchcontinuedlines(str, i)[1]
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
    while !eos() && isspace(str[i])
        if str[i] != '\n'
            i = nextind(str, i)
        else
            catched, _, iend = catchcontinuedlines(str, i)
            if catched
                i = nextind(str, iend)
            else
                break
            end
        end
    end
    return start, prevind(str,i)
end
function catchcomment(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return len, prevind(str,len)
    start = i
    while !eos() && str[i] != '\n'
        i = nextind(str, i)
    end
    return start, prevind(str, i)
end
function catchcontinuedlines(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return false, len, prevind(str,len)
    str[i] != '\n' && return false, i, prevind(str,i)
    if isfixedformfortran && len>nextind(str,i,6) &&
       str[nextind(str,i):nextind(str,i,5)] == "     " && str[nextind(str,i,6)] != ' '
        return true, i, nextind(str,i,6)
    elseif !isfixedformfortran && iscontinuedlinefreeform(str, i)
        # look behind for '&' at end of line
        return true, i, i
    else
        return false, i, prevind(str,i)
    end
end

function iscontinuedlinefreeform(str, i)
    i = thisind(str, i)
    len = i
    while len < ncodeunits(str) && str[len] !='\n'
        len = nextind(str, len)
    end
    # look from start of line
    while i>1 && str[prevind(str, i)] !='\n'
        i = prevind(str, i)
    end
    while i <= len
        if str[i] == ' '
            i = skipwhitespaces(str, i)
        elseif str[i] == '''
            i = skiptoken(str, i)
        elseif str[i] == '#'
            return false
        elseif str[i] == '&'
            t, i = '&', nextind(str, i)
            i = skipwhitespaces(str, i)
            if i > ncodeunits(str) || str[i] == '#' || str[i] == '\n'
                return true
            elseif str[i] == '&' # catch "&&" operator
                i = nextind(str, i)
            else
                @warn("unexpected '$(str[i])' after '$(t)'")
                i = skiptoken(str, i)
            end
        else
            i = skiptoken(str, i)
        end
    end
    return false
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

taketoken(str, i)     =  ( (p1,p2,p3) -> (str[p1:p2], p3) )(marktoken(str, i)...)
picknexttoken(str, i) =  picktoken(str, marktoken(str, i)[3])
skiptoken(str, i)     =  marktoken(str, i)[3]
skiptoken(c::AbstractChar, str, i) = picktoken(str,i) == c ? skiptoken(str,i) : i
skipspaces(str, i)    = catchspaces(str, i)[2] + 1
skipcomment(str, i)   = catchcomment(str, i)[2] + 1
isoneline(str)        = !('\n' in str)
existind(str, i)      = thisind(str, min(max(1,i), ncodeunits(str)))

function skipwhitespaces(str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    eos() && return len+1
    while !eos() && str[i] == ' '
        i = nextind(str, i)
    end
    return i
end

function skipupto(c, str, i)
    eos() = i>len; len = ncodeunits(str); i = thisind(str, i)
    while true
        eos() && return len+1
        if str[i] == ''' != c # skip string
            i = skiptoken(str, i)
        elseif str[i] == c == '\n'
            catched, _, iend = catchcontinuedlines(str, i)
            if catched
                i = nextind(str, iend)
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

function thisline(str, i)
    i = thisind(str, i)
    len = i
    while len < ncodeunits(str) && str[len] !='\n'
        len = nextind(str, len)
    end
    # look from start of line
    while i>1 && str[prevind(str, i)] !='\n'
        i = prevind(str, i)
    end
    return str[i:len]
end

function collectlabels(code)
    code1 = stripcomments(code)
    code1 = concatlinecontinuation(code1)
    lines = splitonlines(code1)

    # capture 'DO's and 'GOTO's labels
    dolabels = Set{String}()
    gotolabels = Set{String}()
    rx = r"^\h*(?:do|for)\h+(\d+)\h+"i
    for i in axes(lines,1)
        (m = match(rx, lines[i])) != nothing && push!(dolabels, m.captures[1])
    end
    rx = r"(?:^\h*|\W)go\h*?to\h+(\d+)$"i
    for i in axes(lines,1)
        (m = match(rx, lines[i])) != nothing && push!(gotolabels, m.captures[1])
    end
    # https://regex101.com/r/ba4wNU/2
    rx = r"(?:^\h*|\W)go\h*?to\h*\((\h*\d+\h*(?:\h*,\h*\d+\h*)*)\)\h*.*$"i
    #rx = r"(?:^\h*|\W)go\h*?to\h*\((\h*\d+\h*(?:\h*,(?:\n[ ]{5}[^ ])?\h*\d+\h*)*)\)\h*.*$"i
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) != nothing
            labels = map(strip, split(m.captures[1], ','))
            foreach(l->push!(gotolabels, l), labels)
        end
    end
    return dolabels, gotolabels
end

function processconditionalgotos(code)
    lines = splitonlines(code)
    # https://regex101.com/r/nKZmku/3
    rx = r"^(\h*|)go\h*?to\h*\(\h*\d+\h*,.*$"i
    rxfull = r"^(\h*|)go\h*?to\h*\((\h*\d+\h*(?:\h*,\h*\d+\h*)*)\)\h*([^#]*)\h*(#.*|)$"i
    i = 1
    while true
        if occursin(rx, lines[i])
            # cut out this line with all continuations
            head = replace(lines[i], rx=>s"\1")
            str = "" * lines[i]
            while iscontinued(lines[i], i<length(lines) ? lines[i+1] : "")
                splice!(lines, i)
                str = concatlines(str,lines[i])
            end
            splice!(lines, i)
            # replace with branch of julia's `if else` statements
            gotos = map(strip, split(replace(str, rxfull => s"\2"), ','))
            var = replace(str, rxfull => s"\3")
            ifexpr = "" * head
            for (i,g) = enumerate(gotos)
                ifexpr *= "if ($(var) .eq. $(i)) then\n$(head)    at_iZjAcpPokMgoto L$(g)\n$(head)else"
                #ifexpr *= "if ($(var) == $(i))\n$(head)    at_iZjAcpPokMgoto L$(g)\n$(head)else"
            end
            ifexpr *= "\n$(head)    at_iZjAcpPokMerror(\"non-exist label $(var)=\\'\$($(var))\\' in goto list\")\n$(head)end\n"
            #println("ifexpr: $ifexpr")
            insert!(lines, i, ifexpr)
        end
        i += 1
        i > length(lines) && break
    end
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function replacedocontinue(code, dolabels, gotolabels)
    lines, comments = splitoncomment(code)

    # insert absent CONTINUE in LABEL SOMETHING
    # https://regex101.com/r/FG2iyI/4
    rxlabeled = r"^(?=[ ]{0,4}\d[ ]{0,4})([\d ]{0,5})"
    rxcont = r"^(\h*(\d+)\h+)continue\h*"i
    rxlabel = r"^(\h*(\d+)\h+)(.*?)$"i
    for i in axes(lines,1)
        if !isnothing(match(rxlabeled, lines[i])) && isnothing(match(rxcont, lines[i]))
            m = match(rxlabel, lines[i])
            if m.captures[2] in dolabels
                lines[i] = " "^length(m.captures[1]) * replace(lines[i], rxlabel=>s"\3")
                lines[i] = lines[i] * "\n" * m.captures[1] * "continue\n"
            end
        end
    end

    code1 = foldl((a,b) -> a*'\n'*b, map(*, lines, comments))
    lines, comments = splitoncomment(code1)

    # replace 'LABEL CONTINUE' with 'end do'
    for i in axes(lines,1)
        m = match(rxcont, lines[i])
        if !isnothing(m) && m.captures[2] in dolabels
            if m.captures[2] in gotolabels
                # 'CYCLE' emulation in the old fortran
                # here some risk to have 'GOTO LABEL' outside of the loop
                rxgoto = Regex("go\\h*to\\h+"*m.captures[2], "i")
                lines = map(a->replace(a, rxgoto => s"cycle"), lines)
                pop!(gotolabels, m.captures[2])
            end
            lines[i] = " "^length(m.captures[1]) * "end do"
            pop!(dolabels, m.captures[2])
        end
    end

    # non replaced LABELs should be fixed by hands
    rx = r"^(\h*(?:do|for)\h+)(\d+)(\h+.*)$"i
    for i in axes(lines,1)
        m = match(rx, lines[i])
        if !isnothing(m) && m.captures[2] in dolabels
            lines[i] = replace(lines[i], rx => s"\1FIXME: L\2\3")
        end
    end

    lines = map(*, lines, comments)
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

"""
replace array's braces with square brackets
"""
function replacearraysbrackets(code, arrays)
    braces = Regex("\\(((?>[^()]|(?R))*)\\)") # complementary braces
    brackets = SubstitutionString("[\\1]")
    for a in arrays
        w = Regex("\\b$(a)\\b","i") # whole word
        for m in collect(eachmatch(w, code))
            o = m.offset
            if code[min(o+length(a),end)] == '('
                code = code[1:o-1] * replace(code[o:end], braces => brackets, count=1)
            end
        end
    end

    # catch second braces, for arrays of strings, https://regex101.com/r/2WIc50/1
    rx = r"(\[((?>[^\[\]]++|(?1))*)\])\h*(\(((?>[^\(\)]++|(?1))*)\))"
    code = replace(code, rx => s"\1[\4]")

    # fix uncompleted ranges in square braces like [1:] and [:1]
    brackets = Regex("\\[((?>[^\\[\\]]++|(?0))*)\\]") # complementary brackets
    # trim spaces around ':'
    rx = r"\[(.*|)([ ]+|):([ ]+|)(.*|)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s"[\1:\4]") * code[o+length(m.match):end]
        end
    end
    # uppend with "end"
    rx = r"\[(.*[^,]):(,.*|)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s"[\1:lastpos_iZjAcpPokM\2]") *
                   code[o+length(m.match):end]
        end
    end
    # prepend with '1'
    rx = r"\[:([^,].*)\]"
    for m in reverse(collect(eachmatch(brackets, code)))
        if occursin(rx, m.match)
            o = m.offset
            code = code[1:o-1] * replace(m.match, rx => s"[1:\1]") * code[o+length(m.match):end]
        end
    end
    return code
end

function collectformat(code)
    code1 = stripcomments(code)
    code1 = concatlinecontinuation(code1)
    lines = splitonlines(code1)

    lines = splitonlines(code)
    rx = r"^\h*(\d+)\h+format\h*\("mi
    rxextr = r"^\h*(\d+)\h+format\h*\((.+)\)"mi
    formats = Dict{String,String}()
    i = 1
    while true
        if occursin(rx, lines[i])
            # collect the lines that make up the format statement
            str = "" * lines[i]
            while iscontinued(lines[i], i<length(lines) ? lines[i+1] : "")
                str = concatlines(str, lines[i+=1])
            end
            if (m = match(rxextr, str)) != nothing
                formats[m.captures[1]] = m.captures[2]
            end
        end
        (i += 1) > length(lines) && break
    end

    return code, formats
end

function convertformat(formatstring)
    rep = r"^(\d+)(.*)$"
    FMT = OrderedDict(
        r"^\/$"           => "\\n",               # / -> START NEW RECORD
        r"^(\d*)P$"i      => "",                  # Scale Factor (P): scale from/to mantissa
        r"^A$"i           => s"%s",               # A -> %s
        r"^A(\d*)$"i      => s"%\1s",             # A9 -> %9s
        r"^I(\d*)$"i      => s"%\1i",             # I5 -> %5i
        r"^E(\d*\.\d*)$"i => s"%\1E",             # E7.2 -> %7.2E
        r"^PE(\d*\.\d*)$"i => s"%\1E",            # E7.2 -> %7.2E Scientific format with Scale Factor P
        r"^F(\d*\.\d*)$"i => s"%\1F",             # F7.2 -> %7.2F
        r"^X$"i           => s" ",                # X -> ' '
        r"^([-]?\d+)P$"i  => s"",                 # 1P -> '' https://docs.oracle.com/cd/E19957-01/805-4939/z4000743a6e2/index.html
        r"^([^']*)'(.*)'([^']*)$" => s"\1\2\3",   # unenclose ''
    )
    format = splitformat(formatstring)

    # eval repeats in the format
    format2 = Vector{String}(undef, 0)
    for (i,el) in enumerate(format)
        if occursin(rep, el)
            num = parse(Int, replace(el, rep=>s"\1"))
            str = replace(el, rep => s"\2")
            occursin(r"\(.*\)", str) && (str = replace(str,  r"\((.*)\)"=> s"\1"))
            [push!(format2, str) for j=1:num]
        else
            push!(format2, el)
        end
    end
    format = format2
    # make conversion
    for rs in FMT
        format = map( a->replace(a, rs), format)
    end
    return foldl(*, format)
end

function processcommon(code, commons, arrays)
    length(commons) == 0 && return code
    lines = splitonlines(code)

    # @unpack COMMONs before use
    rx = r"^(#(\h*)common\h*\/\h*(\w+)\h*\/\h*.+)"i
    for i in axes(lines,1)
        if (m = match(rx, lines[i])) != nothing && haskey(commons, m.captures[3])
            ss =  SubstitutionString("\\2global \\3\n\\2at_iZjAcpPokMunpack " *
                                     commons[m.captures[3]] * " = \\3\n\\1")
            lines[i] = replace(lines[i], rx=>ss)
        end
    end

    # @pack! 'COMMON's back before return
    packstr = "      at_iZjAcpPokMlabel Lreturn\n"
    for k in keys(commons)
        # arrays are not needed because they should not be reallocated
        v = split(commons[k], ',')
        v = v[findall(a->lowercase(a) ∉ map(lowercase,arrays), v)]
        if length(v) > 0
            v = foldl((a,b)->a*','*b, v)
            packstr *= "      at_iZjAcpPokMpackbang_iZjAcpPokM $(k) = $(v)\n"
        end
    end
    lines = foldl((a,b) -> a*'\n'*b, lines)
    # Regex for 'END' of subroutine https://regex101.com/r/evx2Lu/4
    rx = r"(\W)(\w+|)((?:\h*#[^\n]*|\h*)\n+)(\h*end\h*(?:function|(?:recursive\h+|)subroutine|program|module|block|)(?:\h*#[^\n]*|\h*))$"mi
    vm = collect(eachmatch(rx, lines))
    if length(vm) == 0
        @info("something wrong")
        return code
    end

    #length(vm) > 1 && @error("too many end of subroutine")
    m = vm[end]
    if lowercase(m.captures[2]) != "return"
        lines = replace(lines, rx => SubstitutionString("\\1\\2\\3" * packstr * "      ret_iZjAcpPokM\n" * "\\4"))
    else
        lines = replace(lines, rx => SubstitutionString("\\1\n" * packstr * "      ret_iZjAcpPokM\\3\\4"))
    end
    lines, comments = splitoncomment(lines)

    # replace each RETURN with "@goto Lreturn"
    lines = map(a->replace(a, r"\b(return)\b(\h+end|)$"mi => s"at_iZjAcpPokMgoto Lreturn\2"), lines)
    lines = foldl((a,b) -> a*'\n'*b, map(*, lines, comments))
    lines = replace(lines, r"ret_iZjAcpPokM"m => "return")

    return code isa AbstractVector ? splitonlines(lines) : lines
end

function markbypattern(patterns, code)
    lines = stripcomments(code)
    lines = splitonlines(lines)
    marked = Set{Int}()
    for rx in patterns
        for i in axes(lines,1)
            occursin(rx, lines[i]) && push!(marked, i)
        end
    end
    # mark lines with continuations
    marked2 = copy(marked)
    for i in marked
        while occursin(r",$", lines[i]) || occursin(r"&$", lines[i]) ||
              occursin(r"^$", lines[i]) || # full comment-line can occur inside the declaration lines
              occursin(r"^[ ]{5}[^ ]", lines[min(i+1,end)])
            push!(marked2, i+=1)
        end
    end
    return sort(collect(marked2))
end

function markcontinued(code, i)
    lines = splitonlines(code)
    marked = [i]
    # mark lines with continuations
    while i < length(lines) && iscontinued(lines[i], lines[i+1])
        push!(marked, i+=1)
    end
    return marked
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
                lines[i] = replace(lines[i], r"^(\h*)(select\h+case.*)$"mi => SubstitutionString("\\1" * "#\\2") )
            else
                var = @sprintf "EX%s" mod(hash("$(expr)"), 2^16)
                lines[i] = replace(lines[i], r"^(\h*)(select\h+case.*)$"mi => SubstitutionString("\\1"*var*" = "*expr*" #\\2") )
            end
            case1 = true
        end
        # TODO: I also want to process 'CASE (2:4)' with `if var in (2:4)` or `if 2 <= var <= 4`
        if !isempty(var) && occursin(r"^\h*case\h*\("mi, lines[i])
            exprs = replace(lines[i], r"^\h*case\h*\((.*)\)(\h*#.*|\h*)$"mi => s"\1")
            if occursin(',', exprs)
                exprs = replace(exprs, r"(\h*,\h*)" => SubstitutionString(" || "*var*" == "))
            end
            exprs = "if $(var) == " * exprs
            if !case1 exprs = "else" * exprs end
            lines[i] = replace(lines[i], r"^(\h*)case\h*\(.*\)(\h*#.*|\h*)$"mi => SubstitutionString("\\1"*exprs*"\\2"))
            case1 = false
        elseif !isempty(var) && occursin(r"^\h*case\h+default"mi, lines[i])
            lines[i] = replace(lines[i], r"^(\h*)case\h+default(\h*#.*|\h*)$"mi => s"\1else\2")
            var = ""
        elseif !isempty(var) && occursin(r"^\h*end\h+select"mi, lines[i])
            var = ""
        end
    end

    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function processparameters(code)
    patterns = [
        r"^\h*parameter\h*\("mi
    ]
    marked = markbypattern(patterns, code)
    lines, comments = splitoncomment(code)
    for i in marked
        lines[i] = replace(lines[i], patterns[1] => s"      ")
        lines[i] = replace(lines[i], r"," => s";")
        lines[i] = replace(lines[i], r"\)$" => s"")
        lines[i] = replace(lines[i], r"^     [.+&$]" => "      ")
    end

    # TODO: F90 with :: declarations
    patterns = [
        r"\h*[^,#]+\h*,\h*parameter\h*::\h*"mi
    ]

    lines = map(*, lines, comments)
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function processdata(code)
    patterns = [
        r"^\h*\bdata\b\h*[.\w]+(?:\h*\[\h*\d+\h*\]|)\h*\/"mi
    ]
    marked = markbypattern(patterns, code)

    # 'DATA A/1/' -- scalars initializations
    rx = r"(\b\w+\b(?:\h*\[\h*\d+\h*\]|)\h*)(\/((?>[^,\/\n*#]*|(?1))*)\/)"mi
    ss = s"\1 = \3"
    lines, comments = splitoncomment(code)
    for i in marked
        lines[i] = replace(lines[i], rx => ss)
    end
    lines = map(*, lines, comments)
    code1 = foldl((a,b) -> a*'\n'*b, lines)

    replacements = OrderedDict(
    # DATA statements https://regex101.com/r/Z0neJ8/1 https://regex101.com/r/XAol49/2
    r"^(\h*)(data\h*)([.\w]+)\h*(\/((?>[^\/,]++|(?3))*)\/)"mi => @s_str("\\1\\3 = \\5"),
    r"^(\h*)(data\h*)([.\w]+)\h*(\/((?>[^\/]++|(?3))*)\/)"mi => @s_str("\\1\\3 .= (\\5)"),
    r"^(\h*)(data\h*)([,.\w]+)\h*(\/((?>[^\/]++|(?3))*)\/)"mi => @s_str("\\1\\3 = (\\5)"),
    r"^(\h*)(data\h*)(view\(((?>[^()]++|(?3))*)\))\h*(\/((?>[^\/]++|(?3))*)\/)"mi =>
        @s_str("\\1\\3 .= (\\6)"),
    )
    for rx in replacements
        code1 = replace(code1, rx)
    end

    # gather all remaining 'DATA'
    rx = r"(\b\w+\b\h*)(\/((?>[^\/\n#]*|(?1))*)\/)"mi
    ss = s"\1 .= (\3)"
    lines, comments = splitoncomment(code1)
    for i in marked
        lines[i] = replace(lines[i], rx => ss)
        lines[i] = replace(lines[i], r"^(\h*)\bdata\b"mi => s"\1")
        # eat spaces in the numbers 0.111 222 333
        lines[i] = replace(lines[i], r"([^ \[=]) ([^ ])"mi => s"\1\2")
    end
    lines = map(*, lines, comments)

    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function savespecialsymbols(code)
    # save '#', '@', '!', "''"
    # https://regex101.com/r/YrkC9C/1
    rx = r"('((?>[^'\n]*|(?1))*)'){2,}"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = m.match
        len = ncodeunits(code)
        while true
            str = replace(str, r"('((?>[^'\n]*|(?1))*)')('((?>[^'\n]*|(?1))*)')" => s"'\2GhUtwoap_iZjAcpPokM\4'")
            len == (len = ncodeunits(str)) && break
        end
        code = replace(code, m.match => str)
    end
    # save '!' inside strings
    rx = r"('((?>[^'\n]*|(?1))*)')"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"!" => "bang_iZjAcpPokM")
        code = replace(code, m.match => str)
    end
    # save '"' inside strings
    rx = r"('((?>[^'\n]*|(?1))*)')"m
    for m in reverse(collect(eachmatch(rx, code)))
        str = replace(m.match, r"\"" => "quote_iZjAcpPokM")
        code = replace(code, m.match => str)
    end
    code  = replace(code, r"#" => "sha_iZjAcpPokM")
    code  = replace(code, r"@" =>  "at_iZjAcpPokM")
    return code
end
function restorespecialsymbols(code)
    # restore saved symbols
    code = replace(code,      r"sha_iZjAcpPokM" => "#")
    code = replace(code,       r"at_iZjAcpPokM" => "@")
    code = replace(code,     r"bang_iZjAcpPokM" => "!")
    code = replace(code,    r"quote_iZjAcpPokM" => "\\\"")
    code = replace(code, r"GhUtwoap_iZjAcpPokM" => "'")
    code = replace(code,  r"lastpos_iZjAcpPokM" => "end")
    return code
end

const newlinereplacements = OrderedDict(
    # \r\n -> \n, \r -> \n
    r"\r\n" => @s_str("\n"),
    r"\r" => @s_str("\n")
)

const multilinereplacements = OrderedDict(
    # Labels
    r"^(\h*)(\d)(\h+)continue(.*)$"mi     => @s_str(" \\1\\3@label L\\2 \\4"),
    r"^(\h*)(\d\d)(\h+)continue(.*)$"mi   => @s_str("  \\1\\3@label L\\2 \\4"),
    r"^(\h*)(\d\d\d)(\h+)continue(.*)$"mi => @s_str("   \\1\\3@label L\\2 \\4"),
    r"^(\h*)(\d+)(\h+)continue(.*)$"mi    => @s_str("    \\1\\3@label L\\2 \\4"),
    # https://regex101.com/r/J5ViSG/1
    r"^([\h]{0,4})([\d]{1,5})(\h*)(.*)"m             => @s_str("      \\1@label L\\2\n\n     \\3\\4"),
    # array repeating statement https://regex101.com/r/R9g9aU/3 for READ/WRITE and DATA
    # complex expressions like "(A(I)=1,SIZE(A,1))" are ommited
    r"\(([\w]+)\(\h*([^()]+)\h*\),\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\4:\\5)"),
    r"\(([\w]+)\[\h*([^\[\]]+)\h*\],\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\4:\\5)"),
    r"\(([\w]+)\(\h*([^()]+)\h*,\h*([^()]+)\h*\),\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\5:\\6, \\3)"),
    r"\(([\w]+)\[\h*([^\[\]]+)\h*,\h*([^\[\]]+)\h*\],\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\5:\\6, \\3)"),
    r"\(([\w]+)\(\h*([^()]+)\h*,\h*([^()]+)\h*\),\h*(\3)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\2, \\5:\\6)"),
    r"\(([\w]+)\[\h*([^\[\]]+)\h*,\h*([^\[\]]+)\h*\],\h*(\3)\h*\=\h*([^()]+)\h*,\h*(.*?)\)"mi => @s_str("view(\\1, \\2, \\5:\\6)"),
)

const singlewordreplacements = OrderedDict(
    # Use greek symbols, is it need?
    r"\balpha([_\d]*)\b"i   => s"α\1",
    r"\bbeta([_\d]*)\b"i    => s"β\1",
    r"\bgamma([_\d]*)\b"i   => s"γ\1",
    r"\bdelta([_\d]*)\b"i   => s"δ\1",
    r"\bepsilon([_\d]*)\b"i => s"ϵ\1",
    r"\blambda([_\d]*)\b"i  => s"λ\1",
    r"\btheta([_\d]*)\b"i   => s"θ\1",
    # Swap logical symbols
    r"\.true\."i => "true",
    r"\.false\."i => "false",
    r"(\s+)\.or\.(\s+)"i => s"\1||\2",
    r"(\s+)\.and\.(\s+)"i => s"\1&&\2",
    r"(\s+)\.not\.(\s+)"i => s"\1!",
    r"(\s+)\.eq\.(\s+)"i => s"\1==\2",
    r"(\s+)\.ne\.(\s+)"i => s"\1!=\2",
    r"(\s+)\/=(\s+)"i    => s"\1!=\2",
    r"(\s+)\.le\.(\s+)"i => s"\1<=\2",
    r"(\s+)\.ge\.(\s+)"i => s"\1>=\2",
    r"(\s+)\.gt\.(\s+)"i => s"\1>\2",
    r"(\s+)\.lt\.(\s+)"i => s"\1<\2",
    r"(\s+)\.or\.\s*"i => s"\1|| ",
    r"(\s+)\.and\.\s*"i => s"\1&& ",
    r"(\s+)\.not\.\s*"i => s"\1!",
    r"(\s+)\.eq\.\s*"i => s"\1== ",
    r"(\s+)\.ne\.\s*"i => s"\1!= ",
    r"(\s+)\/=\s*"i    => s"\1!= ",
    r"(\s+)\.le\.\s*"i => s"\1<= ",
    r"(\s+)\.ge\.\s*"i => s"\1>= ",
    r"(\s+)\.gt\.\s*"i => s"\1> ",
    r"(\s+)\.lt\.\s*"i => s"\1< ",
    r"\s*\.or\.\s*"i => s" || ",
    r"\s*\.and\.\s*"i => s" && ",
    r"\s*\.not\.\s*"i => s" !",
    r"\s*\.eq\.\s*"i => s" == ",
    r"\s*\.ne\.\s*"i => s" != ",
    r"\s*\/=\s*"i    => s" != ",
    r"\s*\.le\.\s*"i => s" <= ",
    r"\s*\.ge\.\s*"i => s" >= ",
    r"\s*\.gt\.\s*"i => s" > ",
    r"\s*\.lt\.\s*"i => s" < ",
    # Some specific functions
    # https://gcc.gnu.org/onlinedocs/gcc-9.2.0/gfortran/Intrinsic-Procedures.html
    # https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc8/index.html
    r"\bmod\b\h*\("i => s"mod(",
    r"\bd?sqrt\b\h*\("i => s"sqrt(",
    r"\bd?exp\b\h*\("i => s"exp(",
    r"\bd?log\b\h*\("i => s"log(",
    r"\bd?log10\b\h*\("i => s"log10(",
    r"\bd?conjg\b\h*\("i => s"conj(",
    r"\b[da]?imag\b\h*\("i => s"imag(",
    r"\bd?sign\b\("i => "copysign(",
    r"\bd?max1?\b\h*\("i => "max(",
    r"\bd?min1?\b\h*\("i => "min(",
    r"\bd?abs\b\h*\("i => "abs(",
    r"\bd?sin\b\h*\("i => "sin(",
    r"\bd?asin\b\h*\("i => "asin(",
    r"\bd?cos\b\h*\("i => "cos(",
    r"\bd?acos\b\h*\("i => "acos(",
    r"\bd?tan\b\h*\("i => "tan(",
    r"\bd?atan2?\b\h*\("i => "atan(",
    r"\b(i[dq])?nint\b\h*\("i => "round(Int,",
    r"\bint\b\h*\("i => "trunc(Int,",
    r"\breal\b\h*\("i => "Float32(", # ?
    r"\bfloat\b\h*\("i => "Float32(", # ?
    r"\bsngl\b\h*\("i => "Float32(", # ?
    r"\bdble\b\h*\("i => "float(",
    r"\bdfloat\b\h*\("i => "float(",
    r"\bcmplx\b\h*\("i => "ComplexF32(", # ?
    r"\bdcmplx\b\h*\("i => "complex(",
    r"\brand\b\h*\(\)"i => "rand()",
    # https://regex101.com/r/whrGry/2
    r"(?:\brand\b\h*)(\(((?>[^()\n]+|(?1))+)\))"i => s"rand(MersenneTwister(round(Int,\2)))",
    r"\blen\b\h*\("i => "length(",
    # https://regex101.com/r/whrGry/1
    r"(?:\blen_trim\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"length(rstrip(\2))",

    r"\bMPI_SEND\b"i  => "MPI.Send",
    r"\bMPI_ISEND\b"i => "MPI.Isend",
    r"\bMPI_RECV\b"i  => "MPI.Recv!",
    r"\bMPI_IRECV\b"i => "MPI.Irecv!",

    # LSAME
    r"\bLSAME\b\h*\("i => "LSAME(",
    #r"(?:\blsame\b\h*)(\(((?>[^()]++|(?1))*)\))"i => s"((a,b)->==(uppercase(first(a)),uppercase(b)))\1",

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
    # Escape quote "
    r"\"'" => "\\\"'",
    r"'\"" => "'\\\"",
    # Only double quotes allowed for strings
    r"'([^',])'" => s"@\1@",
    "'" => "\"",
    r"@([^\\])@" => s"'\1'",
    r"@([\\])@" => "'\\\\'",
    # unescape "\/" in strings
    r"[\\][\/]" => s"/",
    # Spaces around math operators
    #r"([\*\+\/=])(?=\S)" => s"\1 ",
    #r"(?<=\S)([\*\+\/=])" => s" \1",
    # Spaces around - operators, except after e
    # r"([^e][\-])(\S)" => s"\1 \2",
    #r"(?<!\W\de)(\h*\-\h*)" => s" - ",
    # Space after commas in expressions
    r"(,)([^\s,\[\]\(\)]+\])" => s"@\2",
    r"(,)(\S)" => s"\1 \2",
    #r"@" => ",",
    r"@(?!label)" => ",",
    # Replace XXXXXX with xxxxxx
    r"(\h*)\breturn\b(\h+end|)$"i => s"\1return nothing\2",
    # include
    r"^(\s*)include\h+\"(.*)([.][^.]+)\"(\h*#.*|\h*)$"mi => s"\1include(\"\2.jl\")\4",
    # Replace do LABEL ... -> do ...
    r"for\h+(?:\d+)(\h+.*)$" => s"for\1",
    # Replace do while -> while
    r"do\h+while"i => s"while",
    r"^(\h*)\bdo\b\h*$"i => s"\1while true",
    # Replace ELSE with else
    #r"^(\s*)ELSE"m => s"\1else",
    r"^(\s*)\bELSE\b"m => s"\1else",
    # Relace END XXXX with end
    r"^(\h*)end\h*(?:function|(?:recursive\h+|)subroutine|program|module|block|do|if|select)\h*\w*$"mi => s"\1end",
    r"^(\h*)END$"mi => s"\1end",
    # Don't need CALL
    r"([^ ])\bcall\b(\h+)"i => s"\1 ",
    r"\bcall\b(\h+)"i => s"",
    # Fix assignments
    r"(?<=[^\s=<>!/\\])=(?=[^\s=])" => " = ",
    #r"(?<=[^\s=])=(?=[^\s=])" => " = ",
    # splatting view() in print: https://regex101.com/r/SPfeco/3
    r"((?:\bat_iZjAcpPokMprintf\b|\bprintln\b).*?(?:\bview\b\h*))(\(((?>[^()]++|(?2))*)\))" => s"\1\2...",
    # Remove expression's brackets after if/elseif/while https://regex101.com/r/eF9bK5/17
    r"^(\h*)(if|elseif|while)(\h*)(\(((?>[^()]++|(?4))*)\))\h*$"m => s"\1\2\3\5",
    # Process PARAMETER
    r"^(\h*)parameter(\h+)\((.*)\)(\h*?#.*|\h*?)$"mi => s"\1\2\3\4",
    r"^(\h*).*;\h*parameter\h*::\h*(.*)(\h*?#.*|\h*?)$"mi => s"\1\2\3",
    # breaks
    r"\bexit\b"mi => s"break",
    #r"\breturn\b"mi => s"return",
    r"\bcycle\b"mi => s"continue",
    r"\bstop\b\h+(\d+)"mi => s"exit(\1)",
    r"\bstop\b"mi => s"exit(1)",
    r"//" => " * ",
    # Format floats as "5.0" not "5."
    r"(\W\d+)\.(?=\D)" => s"\1.0",
    r"(\W\d+)\.$"m => s"\1.0",
    r"(?<=\W)\.(\d+)"m => s"0.\1",
    # Floats: 0E0 -> 0.0f+0, 0D0 -> 0.0e+0
    #r"(^|(?<=\W))([-+]?[0-9]+\.?[0-9]*)(([eE])([-+]?[0-9]+))" => s"\2f\5", # is it need?
    r"(^|(?<=\W))([-+]?[0-9]+\.?[0-9]*)(([dD])([-+]?[0-9]+))" => s"\2e\5",
    # Goto
    r"^(\s*)go\h*to\s+(\d+)"mi => s"\1@goto L\2",
    r"(\h)go\h*to\s+(\d+)"mi => s"\1@goto L\2",
    r"(\W)go\h*to\s+(\d+)"mi => s"\1 @goto L\2",
    # fix end
    r"^(.+[^ ])[ ]([ ]+)end$" => s"\1 end\2",
    # remove whitespace between function name and left brace https://regex101.com/r/CxV232/3
    r"(\b(?!elseif|if)[\w_]+\b)[\h]+\("i => s"\1(",
    # Replace suberror with error and mark for fixup
    r"(\W)suberror\((.*?),.*?\)" => s"\1 error(\2)",
    # Implicit declaration
    r"^\h*implicit(.*)"mi => s"",
)

const headersprocessing = OrderedDict(
    # Reorganise functions and doc strings. This may be very project specific.
    # https://regex101.com/r/DAIHhl/1
    r"(\h+)subroutine(\h+)(\w+)(\(([^)]*)\))(\h*#.*?|\h*)\n(#\h*\n)?#\h*function:\h*(.*?)#\h*\n"is =>
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
    ## Trailing whitespace
    #r"\h*$"m,
    # Variable declarations
    #r"\n\s*real\s.*"i,
    #r"\n\s*integer\*?\d*\s.*"i,
    #r"\n\s*logical\s.*"i,
    #r"\n\s*character\s.*"i,
    #r"\n\s*character[*].*"i,
    #r"\n\s*parameter\s.*"i,
    r"\n\s*implicit none"i,
    r"\n\s*real,\s*external\s.*"i,
    r"\n\s*external\s.*"i,
    r"\n\s*intrinsic\s.*"i,
    r"\n\s*contains\s.*"i,
    # Import statements
    r"\n\s*use\s.*"i,
]

main(ARGS)
