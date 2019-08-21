#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

using DataStructures, Printf, DocOpt, JuliaFormatter

const doc =

"""fortran-julia.jl

    This julia script converts fortran 77, 90 code into julia.
    It uses naive regex replacements to do as much as possible,
    but the output WILL need further cleanup.
    From https://gist.github.com/rafaqz/fede683a3e853f36c9b367471fde2f56

    Usage:
        fortran-julia.jl -h | --help
        fortran-julia.jl [--lowercase | --uppercase] [--formatting] <filename.f>

    Options:
        -h --help     Show this screen.
        --version     Show version.
        --uppercase   Convert all identifiers to upper case.
        --lowercase   Convert all identifiers to lower case.
        --formatting  You have luck enough.

"""
#=
KNOWN ISSUES:
    * COMMON, READ/WRITE is very ugly
    * FORMAT is incorrect
    * this script is mostly for fixed-form fortran

TODO:
    return intent(out):
    LSAME
    XERBLA
    P_VALUE
    RWK/IWK
    M_GETMEM
=#

function main()


    args = docopt(doc, version=v"0.1")

    println("Used: $args")

    if args["--uppercase"]
        casetransform = uppercase
    elseif args["--lowercase"]
        casetransform = lowercase
    else
        casetransform = identity
    end

    isnothing(args["<filename.f>"]) && exit

    f = open(args["<filename.f>"])
    code = read(f, String)

    result = convertfromfortran(code, casetransform)

    write(splitext(args["<filename.f>"])[1] * ".jl", result)

    if args["--formatting"]
        result = JuliaFormatter.format_text(result)
    end

    write(splitext(args["<filename.f>"])[1] * ".jl", result)

    return nothing
end

function convertfromfortran(code, casetransform=identity)

    # process newlines
    for rx in newlinereplacements
        code = replace(code, rx)
    end

    # expand tabs
    code = replace(code, r"^[ ]{0,5}\t"m => "      ")
    code = replace(code, r"\t"m => "  ")

    # replace PRINT => WRITE
    code = replace(code, r"(print)\h*([*]|'\([^\n\)]+\)')\h*,"mi => s"write(*,\2)")

    # replace non - continue/format labels
    # CAN'T because this goto's spaghetti logic can broke all!
    # Solution: it is necessary to replace "DO LABEL I=" loops without "CONTINUE"
    # with julia's "while" loops
    #code = replace(code, r"(\n)([ ]{0,4})(\d{1,5})(\h+)(?!.*continue|.*format)([^\n]*)"mi =>
    #                     @s_str("\\1     \\2\\4\\5\n\\2\\3\\4continue"))

    isfixedform = !occursin(r"^[^cC*!\n][^\n]{0,3}[^\d\h\n]"m, code)

    # save '#' and '!'
    code  = replace(code, r"#" => "iZjAcpPokMoshUYgRvpYxh")
    # https://regex101.com/r/eF9bK5/13
    code  = replace(code, r"('[^'!\n]*)(!)([^'\n]*')" => s"\1iZjAcpPGjkERcWRlKIpYxh\3")

    # convert comments
    code = convertfortrancomments(code)
    code = replace(code, r"#=" => "# ") # accident multiline comment

    # split text onto subroutines
    subs, subnames = markbysubroutine(code)
    alllines = splitonlines(code)
    result = "using Printf, FortranFiles, OffsetArrsys"

    for i = 1:length(subs)-1

        code = foldl((a,b) -> a*'\n'*b, view(alllines, subs[i]:subs[i+1]-1))
        print("\nsub name: $(strip(subnames[i]))\n")
        print("$(subs[i]) : $(subs[i+1]-1)\n")

        # replace array's brackets with square braces
        code = replacearraysbrackets(code)

        commons = collectcommon(code)
        code, formats = collectformat(code)
        formats = convertformat(formats)

        # comment out all vars declarations
        code = commentoutdeclarations(code)

        code = processcommon(code, commons)

        code = replacedocontinue(code)

        #printlines(splitonlines(code))
        # process multiline replacements
        for rx in multilinereplacements
            code = replace(code, rx)
        end
        #printlines(splitonlines(code))

        # select case
        code = processselectcase(code)

        # parameters
        code = processparameters(code)

        lines = splitonlines(code)

        lines = enclosereadwrite(lines)

        lines, comments = splitoncomment(lines)
        lines = map(casetransform, lines)

        # straight syntax conversions
        for rx in replacements
            lines = map(a->replace(a, rx), lines)
        end
        lines = map(*, lines, comments)

        lines = includeformat(lines, formats)

        # concate string lines back together and
        # create functions declaration with custom header
        code = foldl((a,b) -> a*'\n'*b, lines)
        for rx in headersprocessing
            len = length(code)
            while true
                code = replace(code, rx)
                len == length(code) && break
                len = length(code)
            end
        end

        # some removing
        for rx in removal
            code = replace(code, rx => "")
        end

        result = result * '\n' * code

    end # of loop over subroutines

    # restore saved '#', '!'
    result = replace(result, r"iZjAcpPokMoshUYgRvpYxh" => "#")
    result = replace(result, r"iZjAcpPGjkERcWRlKIpYxh" => "!")

    return result
end

function markbysubroutine(code)
    lines, comments = splitoncomment(code)

    rx = r"^(\h*)((?:[\w*\h]+\h+)function|subroutine|program|block\h*data|module)\h*.*$"mi
    marked = Vector{Int}(undef, 0)
    names  = Vector{String}(undef, 0)
    for i = 1:length(lines)
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

    #println("subs ranges: $marked")
    return marked, names
end


function convertfortrancomments(code)
    # https://regex101.com/r/AwMX32/2
    #rx = r"^[cC*][^\n]*[\n]|(![^'\n]*(?:'[^'\n]*'[^'\n]*)*[\n])(?=[^']*(?:'[^']*'[^']*)*$)"m
    #rx = r"^[cC*][^\n]*|![^\n]*"m
    rx = r"(^[cC*]|!)[^\n]*"m
    code1 = code isa AbstractVector ? foldl((a,b) -> a*'\n'*b, code) : code
    chars = collect(code1)
    for m in collect(eachmatch(rx, code1))
        chars[m.offset] = '#'
    end
    return code isa AbstractVector ? splitonlines(join(chars)) : join(chars)
end

function stripcomments(code)
    # https://regexr.com/4j4vi via
    # https://stackoverflow.com/questions/9203774/regular-expression-for-comments-but-not-within-a-string-not-in-another-conta
    juliacomments = r"(#=([^*]|\*(?!#))*?=#|#[^\"\n\r]*(?:\"[^\"\n\r]*\"[^\"\n\r]*)*[\r\n])(?=[^\"]*(?:\"[^\"]*\"[^\"]*)*$)"m
    dlms = [r"#"m]
    lines = splitonlines(code)
    for i = 1:length(lines)
        for dlm in dlms
            lines[i] = split(lines[i], dlm, limit=2, keepempty=true)[1]
        end
        lines[i] = rstrip(lines[i])
    end
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function splitoncomment(code)
    # https://regex101.com/r/NohlNS/1
    #rx = r"[\h]*(#[^'\n]*(?:'[^'\n]*'[^'\n]*)*)(?=[^']*(?:'[^']*'[^']*)*$)"m
    rx = r"\h*#.*$"
    lines = splitonlines(code)
    comments  = ["" for i=1:length(lines)]
    for i = 1:length(lines)
        (m = match(rx, lines[i])) != nothing && (comments[i] = m.match)
        lines[i] = replace(lines[i], rx => s"")
    end
    return lines, comments
end

function expandlinecontinuation(code)
    # Note: trailing comment is not allowed here
    replacements = OrderedDict(
        r"\n     [.+&$\w]\h*" => "",  # fixed-form startline continuations
        r"&\n"              => ""   # free-form continuations
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

function printlines(list)
    print("\n")
    #for (i,l) in enumerate(list) print("$l\n") end
    for (i,l) in enumerate(list) print("list[$i]:$l\n") end
end

function commentoutdeclarations(code)

    matches = [
        # Match declarations
        r"^\h*common\h*\/\h*\w+\h*\/\h*"mi,
        r"^\h*dimension\h+"mi,
        r"^\h*implicit\h*\w+"mi,
        r"^\h*real.*::"mi,
        r"^\h*real\h+(?!function)"mi,
        r"^\h*real\*[0-9]+\h+(?!function)"mi,
        r"^\h*double\h*precision\h+(?!function)"mi,
        r"^\h*integer.*::"mi,
        r"^\h*integer\h+(?!function)"mi,
        r"^\h*integer\*[0-9]+\h+(?!function)"mi,
        r"^\h*character.*::"mi,
        r"^\h*character\h+"mi,
        r"^\h*character\*"mi,
        r"^\h*external\h+"mi,
        r"^\h*logical.*::"mi,
        r"^\h*logical\h+(?!function)"mi,
        r"^\h+\d+\h+format\h*\("mi
    ]

    list = splitonlines(code)
    for i in markbypattern(matches, code)
        list[i] = "#" * list[i]
    end

    return code isa AbstractVector ? list : foldl((a,b) -> a*'\n'*b, list)
end

function collectarrays(code)
    lines = stripcomments(code)
    lines = expandlinecontinuation(lines)
    list = splitonlines(lines)

    matches = [
        # Match declarations
        r"^\h*common\h*/\h*\w+\h*/\h*"mi,
        r"^\h*dimension\h+"mi,
        r"^\h*real\h+(?!function)"mi,
        r"^\h*real\*[0-9]+\h+(?!function)"mi,
        r"^\h*double\h+precision\h+(?!function)"mi,
        r"^\h*integer\h+(?!function)"mi,
        r"^\h*integer\*[0-9]+\h+(?!function)"mi,
        r"^\h*character\h+"mi,
        r"^\h*character\*\h+"mi
    ]

    matched = Vector{String}()
    for rx in matches
        for i in list
            occursin(rx, i) && push!(matched, strip(i))
        end
    end
    # match only arrays names
    for i in 1:length(matched)
        for rx in matches
            matched[i] = replace(matched[i], rx => s"")               # drop declaration
        end
        matched[i] = replace(matched[i], r" " => s"")               # drop spaces
        matched[i] = replace(matched[i], r"\([^)]+\)" => s"()")     # clean inside braces
        matched[i] = replace(matched[i], r"[^(),]+," => s"")        # drop non-arrays (without braces)
        matched[i] = replace(matched[i], r",[^(),]+$"m => s"")      # drop last non-array (without braces)
        matched[i] = replace(matched[i], r"^[^(),]+$"m => s"")      # drop single non-array (without braces)
        matched[i] = replace(matched[i], r"\*\(" => s"(")           # drop '*' in character*()
    end
    matches = [r"^\h*.*dimension\h*\(.+\)\h*::\h*"mi,
               r"^\h*character\*\(\*\)\h*"mi,
               r"^\h*character\*\(\d+\)\h*"mi ]
    for rx in matches
        for i in list
            occursin(rx, i) && push!(matched, strip(replace(i, rx=>"")))
        end
    end
    vars = Vector{String}()
    for i in 1:length(matched)
        length(matched[i]) > 0 && append!(vars, split(matched[i], ","))
    end
    vars = unique(vars)
    for i in 1:length(vars)
        vars[i] = replace(vars[i], r" " => s"")
        vars[i] = replace(vars[i], r"^(\w+)\(.*\)$"m => s"\1")
    end

    return vars
end

function collectcommon(code)
    lines = stripcomments(code)
    lines = expandlinecontinuation(lines)
    list = splitonlines(lines)
    rx = r"^\h*common\h*\/\h*(\w+)\h*\/\h*(.+)"mi
    matched = Dict{String,String}()
    for i in 1:length(list)
        m = match(rx, list[i])
        if !isnothing(m) matched[m.captures[1]] = m.captures[2] end
    end
    # match all var names
    for i in keys(matched)
        matched[i] = replace(matched[i], r" " => s"")               # drop spaces
        matched[i] = replace(matched[i], r"\[[^\]]+\]" => s"")       # drop braces
    end
    #println("common: $matched")
    return matched
end

function replacedocontinue(code)
    # FIXME: there can be dropped 'continue'
    lines, comments = splitoncomment(code)
    dolabels = Set{String}()
    rx = r"^\h*do\h+(\d+)\h+"i
    for i in 1:length(lines)
        m = collect(eachmatch(rx, lines[i]))
        length(m) > 0 && push!(dolabels, m[1].captures[1])
    end
    rx = r"^(\h*(\d+)\h+)continue\h*"i
    for i in 1:length(lines)
        m = collect(eachmatch(rx, lines[i]))
        if length(m) > 0
            if m[1].captures[2] in dolabels
                #println("capture: $(m[1])")
                lines[i] = " "^length(m[1].captures[1]) * "end do"
                pop!(dolabels, m[1].captures[2])
            end
        end
    end
    lines = map(*, lines, comments)
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

# replace array's brackets with square braces
function replacearraysbrackets(code)
    braces = Regex("\\(((?>[^()]|(?R))*)\\)") # complementary braces
    sqbraces = SubstitutionString("[\\1]")
    arrays = collectarrays(code)
    println("arrays: $arrays")
    for a in arrays
        w = Regex("\\b$(a)\\b","i") # whole word
        for m in collect(eachmatch(w, code))
            o = m.offset
            if code[min(length(code), o+length(a))] == '('
                code = code[1:o-1] * replace(code[o:end], braces => sqbraces, count=1)
            end
        end
    end
    return code
end

function collectformat(code)
    lines = stripcomments(code)
    lines = expandlinecontinuation(lines)
    list = splitonlines(lines)
    rx = r"^\h*(\d+)\h+format\h*\((.+)\)"mi
    format = Dict{String,String}()
    for i in 1:length(list)
        m = match(rx, list[i])
        if !isnothing(m) format[m.captures[1]] = m.captures[2] end
    end
    FMT = OrderedDict(
        r"^([^']*)'(.*)'([^']*)$"=>s"\1\2\3" # unenclose ''
       )
    for i in keys(format)
        format[i] = replace(format[i], r"[/]" => ",\\n,")
        v = split(format[i], ",")
        #v = map(strip, v)
        for rs in FMT
            v = map( a->replace(a, rs), v)
        end
        format[i] = foldl((a,b) -> a*','*b, v)
    end
    # collect all other format strings
    lines = splitonlines(code)
    idx = 100001
    rx = r"^(\h*\d*\h*(?:write|read)\h*\(\h*[\*\w]+\h*,\h*)('\(.+\)')(\h*(?:,\h*[ =\w]+\h*)?\)\h*.*)"mi
    for i in 1:length(lines)
        m = match(rx, lines[i])
        if !isnothing(m)
            format[string(idx)] = replace(m.captures[2], r"'\((.*)\)'" => s"\1")
            lines[i] = m.captures[1] * string(idx) * m.captures[3]
            idx += 1
        end
    end
    #println("FORMATs: $format")
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines), format
end

function convertformat(formatstrings)
    rep = r"^(\d+)(.*)$"
    FMT = OrderedDict(
        r"^A$"i           => s"%s",               # A -> %s
        r"^A(\d*)$"i      => s"%\1s",             # A9 -> %9s
        r"^I(\d*)$"i      => s"%\1i",             # I5 -> %5i
        r"^E(\d*\.\d*)$"i => s"%\1E",             # E7.2 -> %7.2E
        r"^F(\d*\.\d*)$"i => s"%\1F",             # F7.2 -> %7.2F
        r"^X$"i           => s" ",                # X -> ' '
    )
    for k in keys(formatstrings)
        format = split(formatstrings[k], ",") # FIXME: format is recursive
        format = map(strip, format)
        format2 = Vector{String}(undef, 0)
        for (i,el) in enumerate(format) # eval repeats in format
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
        for rs in FMT
            format = map( a->replace(a, rs), format)
        end
        formatstrings[k] = foldl(*, format)
        if occursin(r"\$$", formatstrings[k])
            formatstrings[k] = replace(formatstrings[k], @r_str("\\\$\$") => "")
        else
            formatstrings[k] *= "\\n"
        end
    end
    return formatstrings
end

function includeformat(code, format)
    # https://regex101.com/r/auYywX/1
    rx = r"^(\h*\d*\h*)(write|read)(\h*\(\h*)([\*\w]+)(\h*,\h*)([*]|\d+)(\h*(?:,\h*[ =\w]+\h*)?\))(\h*.*)"mi
    lines = splitonlines(code)
    for i in 1:length(lines)
        m = match(rx, lines[i])
        if !isnothing(m)
            io = m.captures[4]
            if  lowercase(m.captures[2]) == "read" && (m.captures[4] == "5" || m.captures[4] == "*")
                io = "stdin"
            elseif  lowercase(m.captures[2]) == "write" && (m.captures[4] == "6" || m.captures[4] == "*")
                io = "stdout"
            end
            if  lowercase(m.captures[2]) == "read"
                cmd = "READ"
            elseif m.captures[6] == "*"
                cmd = "println"
            else
                cmd = "@printf"
            end
            if m.captures[6] == "*"
                ss = m.captures[1] * cmd * m.captures[3] * io * m.captures[5] * m.captures[8]
            else
                ss = m.captures[1] * cmd * m.captures[3] * io * m.captures[5] *
                          '"' * format[m.captures[6]] * "\"," * m.captures[8]
            end
            lines[i] = replace(lines[i], rx=>ss)
        end
    end
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function enclosereadwrite(code)
    rx = r"^\h*\d*\h*(write|read)\h*\(\h*[\*\w]+\h*,[^)]+\).*"mi
    lines, comments = splitoncomment(code)
    for i in 1:length(lines)
        if occursin(rx, lines[i])
            j = i
            while occursin(r",$", lines[j]) || occursin(r"&$", lines[j]) ||
                  occursin(r"^[ ]{5}[^ ]", lines[min(j+1,end)])
                j += 1
            end
            lines[j] = lines[j] * ")"
        end
    end
    lines = map(*, lines, comments)
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end


function processcommon(code, commons)
    lines, comments = splitoncomment(code)
    for k in keys(commons)
        for var in split(commons[k], ",")
            w = Regex("\\b$(var)\\b","i") # whole word
            ss = SubstitutionString("$(k).$(var)")
            lines = map( a->replace(a, w=>ss), lines)
        end
    end
    lines = map(*, lines, comments)
    rx = r"^(#(\h*)common\h*\/\h*(\w+)\h*\/\h*.+)"i
    ss =  SubstitutionString("\\2global \\3\n\\1")
    lines = map( a->replace(a, rx=>ss), lines)
    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
end

function markbypattern(patterns, code)
    lines = stripcomments(code)
    list = splitonlines(lines)
    marked = Set{Int}()
    for rx in patterns
        for (i,l) in enumerate(list)
            occursin(rx, l) && push!(marked, i)
        end
    end
    # mark lines with continuations
    marked2 = copy(marked)
    for i in marked
        while occursin(r",$", list[i]) || occursin(r"&$", list[i]) ||
              occursin(r"^[ ]{5}[^ ]", list[min(i+1,end)])
            push!(marked2, i+=1)
        end
    end
    return sort(collect(marked2))
end

function processselectcase(code)
    lines = splitonlines(code)
    var = ""
    case1 = false
    for i = 1:length(lines)
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
    matches = [
        r"^\h*parameter\h*\("mi
    ]
    marked = markbypattern(matches, code)
    lines, comments = splitoncomment(code)
    for i in marked
        lines[i] = replace(lines[i], matches[1] => s"      ")
        lines[i] = replace(lines[i], r"," => s";")
        lines[i] = replace(lines[i], r"\)$" => s"")
        lines[i] = replace(lines[i], r"^     [.+&$]" => "      ")
    end

    # TODO: F90 with :: declarations
    matches = [
        r"\h*[^,#]+\h*,\h*parameter\h*::\h*"mi
    ]

    lines = map(*, lines, comments)

    return code isa AbstractVector ? lines : foldl((a,b) -> a*'\n'*b, lines)
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
    r"^([\h]{0,4})([\d]{1,5})(\h*)(.*)"m             => @s_str("     \\1@label L\\2\n\n   \\3\\4"),
    #r"(?=^(\h*)(\d+)(\h+)(.*))(?:!^.{5}\w)"m             => @s_str("    \\1@label L\\2\n\n    \\3\\4"),
    #r"^(\h*)(\d+)(\h+)(.*)$"m             => @s_str("    \\1@label L\\2\n\n    \\3\\4")
    # IF one-liners: https://regex101.com/r/eF9bK5/16 via
    # https://stackoverflow.com/questions/36357344/how-to-match-paired-closing-bracket-with-regex
    r"^(\h*)(if\h*)(\(((?>[^()]++|(?3))*)\))(\h*(?!.*then)[^#\n]+(?:\n[ ]{5}[.$&\w][^\n]+)*)"mi =>
    SubstitutionString("\\1\\2\\3\\5 end "),
    # https://regex101.com/r/eF9bK5/15
    r"^(\h*)(if\h*)(\(((?>[^()]++|(?3))*)\))(\h*)((#[^\n]*|)*\n[ ]{5}[.$&\w])+(\h*(?!.*then)[^\n]+)"mi =>
    SubstitutionString("\\1\\2\\3\\5\\7\n      \\8\n\\1end"),
    # array repeating statement https://regex101.com/r/R9g9aU/2
    # skip complex expressions like "(A(I)=1,SIZE(A,1))"
    r"\(([^)]+)\(([^()]+)\),\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\4:\\5)"),
    r"\(([^)]+)\[([^\[\]]+)\],\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\4:\\5)"),
    r"\(([^)]+)\(([^()]+),([^()]+)\),\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\5:\\6, \\3)"),
    r"\(([^)]+)\[([^\[\]]+),([^\[\]]+)\],\h*(\2)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\5:\\6, \\3)"),
    r"\(([^)]+)\(([^()]+),([^()]+)\),\h*(\3)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\5:\\6, \\2)"),
    r"\(([^)]+)\[([^\[\]]+),([^\[\]]+)\],\h*(\3)\h*\=\h*([^()]+)\h*,\h*(.*)\)"mi => @s_str("view(\\1, \\5:\\6, \\2)"),
    # DATA statements https://regex101.com/r/XAol49/1 https://regex101.com/r/XAol49/2
    r"^(\h*)(data\h*)([.\w]+)\h*(\/((?>[^\/,]++|(?3))*)\/)"mi => @s_str("\\1\\3 = \\5"),
    r"^(\h*)(data\h*)([.\w]+)\h*(\/((?>[^\/]++|(?3))*)\/)"mi => @s_str("\\1\\3 .= (\\5)"),
    r"^(\h*)(data\h*)([,.\w]+)\h*(\/((?>[^\/]++|(?3))*)\/)"mi => @s_str("\\1\\3 = (\\5)"),
    r"^(\h*)(data\h*)(view\(((?>[^()]++|(?3))*)\))\h*(\/((?>[^\/]++|(?3))*)\/)"mi => @s_str("\\1\\3 .= (\\6)"),
)

# Regex/substitution pairs for replace(). Order matters here.
const replacements = OrderedDict(
    # Remove '&' / '.' multiline continuations
    r"^     [.+&$\w](.*)"m => s"      \1",
    # Powers use ^ not **
    r"\*\*" => "^",
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
    # DO loop to for loop, https://regex101.com/r/e2dw0k/13
    r"do(\h+\d+)?\h+(.+)(?<!\(|\[),(?![\w\s]*[\)\]])(.+)(?<!\(|\[),(?![\w\s]*[\)\]])(.+)"i => s"for \2:\4:\3",
    #r"do(\h+\d+)?\h+(.+)(?<!\(),(?![\w\s]*[\)])(.+)(?<!\(),(?![\w\s]*[\)])(.+)"i => s"for \2:\4:\3",
    r"do(\h+\d+)?\h+(.+)(?<!\(|\[),(?![\w\s]*[\)\]])(.+)"i => s"for \2:\3",
    #r"do(\h+\d+)?\h+(.+)(?<!\(),(?![\w\s]*[\)])(.+)"i => s"for \2:\3",
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
    r"(\s*)RETURN" => s"\1return",
    # include
    r"^(\s*)include\h+\"(.*)([.][^.]+)\"(\h*#.*|\h*)$"mi => s"\1include(\"\2.jl\")\4",
    # Replace do while -> while
    r"do\h+while"i => s"while",
    # Replace ELSEIF/ELSE IF with elseif
    r"^(\s*)else\s*if"mi => s"\1elseif",
    # Replace "IF(" to "if ("
    r"^(\s*)(elseif|if)\("mi => s"\1\2 (",
    # Remove THEN
    r"([)\s])then(\s*)"i => s"\1\2",
    # Replace IF with if
    r"^(\s*)IF(?=\W)"m => s"\1if",
    # Replace ELSE with else
    #r"^(\s*)ELSE(?=\W)"m => s"\1else",
    r"^(\s*)ELSE"m => s"\1else",
    # Relace END XXXX with end
    r"^(\h*)end\h*(?:function|subroutine|program|module|block|do|if)\h*\w*$"mi => s"\1end",
    r"^(\h*)END$"mi => s"\1end",
    # Replace expnent function
    r"(\W)exp\("i => s"\1exp(",
    # Don't need CALL
    r"(\s*)call(\h+)"i => s"\1",
    # Use real math symbols, is it need?
    r"\bgamma\b"i => "γ",
    r"\btheta\b"i => "θ",
    r"\bepsilon\b"i => "ϵ",
    r"\blambda\b"i => "λ",
    r"\balpha\b"i => "α",
    r"\bbeta\b"i => "β",
    # Swap logical symbols
    r"\.true\."i => "true",
    r"\.false\."i => "false",
    r"(\s+)\.or\.(\s+)"i => s"\1||\2",
    r"(\s+)\.and\.(\s+)"i => s"\1&&\2",
    r"(\s+)\.not\.(\s+)"i => s"\1!",
    r"(\s+)\.eq\.(\s+)"i => s"\1==\2",
    r"(\s+)\.ne\.(\s+)"i => s"\1!=\2",
    r"(\s+)\.le\.(\s+)"i => s"\1<=\2",
    r"(\s+)\.ge\.(\s+)"i => s"\1>=\2",
    r"(\s+)\.gt\.(\s+)"i => s"\1>\2",
    r"(\s+)\.lt\.(\s+)"i => s"\1<\2",
    r"(\s+)\.or\.\s*"i => s"\1|| ",
    r"(\s+)\.and\.\s*"i => s"\1&& ",
    r"(\s+)\.not\.\s*"i => s"\1!",
    r"(\s+)\.eq\.\s*"i => s"\1== ",
    r"(\s+)\.ne\.\s*"i => s"\1!= ",
    r"(\s+)\.le\.\s*"i => s"\1<= ",
    r"(\s+)\.ge\.\s*"i => s"\1>= ",
    r"(\s+)\.gt\.\s*"i => s"\1> ",
    r"(\s+)\.lt\.\s*"i => s"\1< ",
    r"\s*\.or\.\s*"i => s" || ",
    r"\s*\.and\.\s*"i => s" && ",
    r"\s*\.not\.\s*"i => s" !",
    r"\s*\.eq\.\s*"i => s" == ",
    r"\s*\.ne\.\s*"i => s" != ",
    r"\s*\.le\.\s*"i => s" <= ",
    r"\s*\.ge\.\s*"i => s" >= ",
    r"\s*\.gt\.\s*"i => s" > ",
    r"\s*\.lt\.\s*"i => s" < ",
    # Fix assignments
    r"(?<=[^\s=<>!/\\])=(?=[^\s=])" => " = ",
    #r"(?<=[^\s=])=(?=[^\s=])" => " = ",
    # Add end after single line if: https://regex101.com/r/eF9bK5/6
    #r"^(\h*)if(\h*)(\(((?>[^()]++|(?3))*)\))\h*((?!then)[^\h].+)"i => s"\1if\2\3 \5 end",
    #### Single-line IF statement with various terminations
    ###r"^(\h*)if\s*\((.*?)\)\s*exit(\s*#.*|\s*)$"mi => s"\1\2 && break\3",
    ###r"^(\h*)if\s*\((.*?)\)\s*return(\s*#.*|\s*)$"mi => s"\1\2 && return\3",
    ###r"^(\h*)if\s*\((.*?)\)\s*cycle(\s*#.*|\s*)$"mi => s"\1\2 && continue\3",
    ###r"^(\h*)if\s*\((.*?)\)\s*go\h*to\s+(\d+\h*)(?!end)$"i => s"\1\2 && @goto L\3",
    ####r"(^\h*)if\s*\((.*?)\)\s*go\h*to\s+(\d+\h*)(?!end)(.*)"i => s"\1\2 && @goto L\3\4",
    # Remove expression's brackets after if/elseif/while https://regex101.com/r/eF9bK5/17
    r"^(\h*)(if|elseif|while)(\h*)(\(((?>[^()]++|(?4))*)\))\h*$"m => s"\1\2\3\5",
    # Process PARAMETER
    r"^(\h*)parameter(\h+)\((.*)\)(\h*?#.*|\h*?)$"mi => s"\1\2\3\4",
    r"^(\h*).*;\h*parameter\h*::\h*(.*)(\h*?#.*|\h*?)$"mi => s"\1\2\3",
    # breaks
    r"\bexit\b"mi => s"break",
    r"\breturn\b"mi => s"return",
    r"\bcycle\b"mi => s"continue",
    r"\bstop\b\h+(\d+)"mi => s"exit(code=\1)",
    r"\bstop\b"mi => s"exit(code=1)",
    # Some specific functions
    r"\bsign\b\("i => "copysign(",
    r"(?<=\W)MAX\(" => "max(",
    r"(?<=\W)MIN\(" => "min(",
    r"(?<=\W)ABS\(" => "abs(",
    r"//" => " * ",
    # Format floats as "5.0" not "5."
    r"(\W\d+)\.(\D)" => s"\1.0\2",
    r"(\W\d+)\.$"m => s"\1.0",
    r"(?<=\W)\.(\d+)"m => s"0.\1",
    # Floats: 0E0 -> 0.0f+0, 0D0 -> 0.0e+0
    #r"(^|(?<=\W))([-+]?[0-9]+\.?[0-9]*)(([eE])([-+]?[0-9]+))" => s"\2f\5", # is it need?
    r"(^|(?<=\W))([-+]?[0-9]+\.?[0-9]*)(([dD])([-+]?[0-9]+))" => s"\2e\5",
    # Goto
    r"^(\s*)go\h*to\s+(\d+)"mi => s"\1@goto L\2",
    r"(\h)go\h*to\s+(\d+)"mi => s"\1@goto L\2",
    r"(\W)go\h*to\s+(\d+)"mi => s"\1 @goto L\2",
    # Tab to 4 spaces
    r"\t" => "    ",
    # fix end
    r"^(.+[^ ])[ ]([ ]+)end$" => s"\1 end\2",
    # Relace suberror with error and mark for fixup
    r"(\W)suberror\((.*?),.*?\)" => s"\1 error(\2)",
    # Mark #FIXME the various things this script can't handle
    #r"(write|while\s)" => s"#FIXME \1",
    # Implicit declaration
    r"^\h*implicit(.*)"mi => s"",
    # Custom functions
    r"\bMPI_SEND\b"i => "MPI.Send",
    r"\bMPI_ISEND\b"i => "MPI.Isend",
    r"\bMPI_RECV\b"i => "MPI.Recv!",
    r"\bMPI_IRECV\b"i => "MPI.Irecv!",
    # J_LEN(a) -> length(strip(a))
    r"(\w+)\(1:J_LEN\(\1\)\)"i => s"strip(\1)",
    #r"J_LEN\(([^()]|(?R))*\)"i => s"length(strip(\1))",
    r"J_LEN\("i => s"length(",
)

const headersprocessing = OrderedDict(
    # Reorganise functions and doc strings. This may be very project specific.
    # https://regex101.com/r/DAIHhl/1
    r"(\h+)subroutine(\h+)(\w+)(\(([^)]*)\))(\h*#.*?|\h*)\n(#\h*\n)?#\h*function:\h*(.*?)#\h*\n"is =>
    SubstitutionString("\"\"\"\n    \\3(\\5)\n\n\\8\"\"\"\nfunction \\3(\\5)\n#\n"),
    # Simple subroutine
    r"^\h*subroutine\h+"mi => s"function ",
    r"^\h*program\h+(\w+)"mi => s"function \1()",
    r"^\h*block\h*data\h+(\w+)"mi => s"function \1()",
    r"^\h*module\h+"mi => s"module ",
    r"^\h*real(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
    r"^\h*double\h*precision\h*function\h+"mi => s"function ",
    r"^\h*integer(?:\*\d{1,2})?\h*function\h+"mi => s"function ",
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
    r"\n\s*real, external\s.*"i,
    r"\n\s*external\s.*"i,
    r"\n\s*intrinsic\s.*"i,
    r"\n\s*contains\s.*"i,
    # Import statements
    r"\n\s*use\s.*"i,
]

main()
