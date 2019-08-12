#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' \
    "${BASH_SOURCE[0]}" "$@"
=#

#=
This julia script converts fortran 77, 90 code into julia.
It uses naive regex replacements to do as much as possible,
but the output WILL need further cleanup.

To run from a shell: 

julia fortran-julia.jl filename.f90

Output is written to filename.jl.
=#

using DataStructures, Printf

function convertfortrancomments(code)
    dlms = [
        r"^\*"m,
        r"^c"mi,
        "!" # there is will be wrong case on escaped "!"
    ]
    lines = splitonlines(code)
    for i = 1:length(lines)
        for r in dlms
            lines[i] = replace(lines[i], r => s"#", count=1)
        end
    end
    if !(code isa AbstractVector)
        lines = foldl((a,b) -> a*'\n'*b, lines)
    end
    return lines
end

function stripcomments(code)
    dlms = [
        r"#"m
    ]
    lines = splitonlines(code)
    for i = 1:length(lines)
        for dlm in dlms
            lines[i] = split(lines[i], dlm, limit=2, keepempty=true)[1]
        end
        lines[i] = rstrip(lines[i])
    end
    if !(code isa AbstractVector)
        lines = foldl((a,b) -> a*'\n'*b, lines)
    end
    return lines
end

function splitonfortrancomment(lines)
    dlms = [
        r"^\*"m,
        r"^c"mi,
        "!" # there is will be wrong case on escaped "!"
    ]
    code = splitonlines(lines)
    comments  = ["" for i=1:length(code)]
    for i = 1:length(code)
        for dlm in dlms
            if length(code[i]) > 0
                r = split(code[i], dlm, limit=2, keepempty=true)
                code[i] = r[1]
                length(r) > 1 && (comments[i] = r[2] * comments[i])
            end
        end
    end
    for i = 1:length(code)
        if length(comments[i]) > 0 comments[i] = "#" * comments[i] end
    end
    # absorb tail spaces to comment
    for i = 1:length(code)
        r = collect(eachmatch(r"\h+$", code[i]))
        if length(r) > 0
            code[i] = replace(code[i], r"\h+$" => "")
            comments[i] = r[1].match * comments[i]
        end
    end
    return code, comments
end

function splitoncomment(lines)
    dlm = r"#"m
    code = splitonlines(lines)
    comments  = ["" for i=1:length(code)]
    for i = 1:length(code)
        if length(code[i]) > 0
            r = split(code[i], dlm, limit=2, keepempty=true)
            code[i] = r[1]
            length(r) > 1 && (comments[i] = "#" * r[2])
        end
    end
    # absorb tail spaces to comment
    for i = 1:length(code)
        r = collect(eachmatch(r"\h+$", code[i]))
        if length(r) > 0
            code[i] = replace(code[i], r"\h+$" => "")
            comments[i] = r[1].match * comments[i]
        end
    end
    return code, comments
end

function expandlinecontinuation(code)
    # Note: trailing comment is not allowed here
    replacements = OrderedDict(
        r"\n     [.+&$]\h*" => " ",  # fixed-form startline continuations
        r"&\n"              => " "   # free-form continuations
    )
    code1 = deepcopy(code)
    for r in replacements
        code1 = replace(code1, r)
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
    for (i,l) in enumerate(list) print("list[$i]:$l\n") end
end

function collectarrays(code)
    lines = stripcomments(code)
    lines = expandlinecontinuation(lines)
    list = splitonlines(lines)

    matches = [
        # Match declarations
        r"^\h*real\h+"mi,
        r"^\h*real\*[0-9]+\h+"mi,
        r"^\h*double\h+precision\h+"mi,
        r"^\h*integer\h+"mi,
        r"^\h*integer\*[0-9]+\h+"mi,
        r"^\h*character\h+"mi,
        r"^\h*character\*"mi
    ]

    matched = Vector{String}()
    for r in matches
        for i in list
            occursin(r, i) && push!(matched, strip(i))
        end
    end
    # match only arrays names
    for i in 1:length(matched)
        for r in matches
            matched[i] = replace(matched[i], r => s"")               # drop declaration
        end
        matched[i] = replace(matched[i], r" " => s"")               # drop spaces
        matched[i] = replace(matched[i], r"\([^)]+\)" => s"()")     # clean inside braces
        matched[i] = replace(matched[i], r"[^(),]+," => s"")        # drop non-arrays (without braces)
        matched[i] = replace(matched[i], r",[^(),]+$"m => s"")      # drop last non-array (without braces)
        matched[i] = replace(matched[i], r"^[^(),]+$"m => s"")      # drop single non-array (without braces)
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

function markbypattern(patterns, code)
    lines = stripcomments(code)
    list = splitonlines(lines)
    marked = Set{Int}()
    for r in patterns
        for (i,l) in enumerate(list)
            occursin(r, l) && push!(marked, i)
        end
    end
    # mark lines with continuations
    marked2 = copy(marked)
    for i in marked
        while occursin(r",$", list[i]) || occursin(r"&$", list[i])
            push!(marked2, i+=1)
        end
    end
    return sort(collect(marked2))
end

function commentoutdeclarations(code)

    matches = [
        # Match declarations
        r"^\h*real\h+"mi,
        r"^\h*real\*[0-9]+\h+"mi,
        r"^\h*double\h+precision\h+"mi,
        r"^\h*integer\h+"mi,
        r"^\h*integer\*[0-9]+\h+"mi,
        r"^\h*character\h+"mi,
        r"^\h*character\*"mi,
        r"^\h*external\h+"mi,
        r"^\h*logical\h+"mi
    ]

    list = splitonlines(code)
    for i in markbypattern(matches, code)
        list[i] = "#" * list[i]
    end

    if !(code isa AbstractVector)
        list = foldl((a,b) -> a*'\n'*b, list)
    end

    return list
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
        # FIXME: I also want to process 'CASE (2:4)'
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
    if !(code isa AbstractVector)
        lines = foldl((a,b) -> a*'\n'*b, lines)
    end
    return lines
end

function processparameters(code)
    matches = [
        r"^\h*parameter\h+\("mi
    ]
    marked = markbypattern(matches, code)
    lines, comments = splitoncomment(code)
    for i in marked
        lines[i] = replace(lines[i], matches[1] => s"      ")
        lines[i] = replace(lines[i], r"," => s";")
        lines[i] = replace(lines[i], r"\)$" => s"")
        lines[i] = replace(lines[i], r"^     [.+&$]" => "      ")
    end

    # FIXME: F90
    matches = [
        r"\h*[^,#]+\h*,\h*parameter\h*::\h*"mi
    ]

    lines = map(*, lines, comments)
    if !(code isa AbstractVector)
        lines = foldl((a,b) -> a*'\n'*b, lines)
    end
    return lines
end

function processonelineif(code)
    # https://regex101.com/r/eF9bK5/3 via
    # https://stackoverflow.com/questions/36357344/how-to-match-paired-closing-bracket-with-regex
    r = r"(\h*)(if\h*|IF\h*)(\(((?>[^()]++|(?3))*)\))\h*\n     [.+&$]([^\n]+)"m
    s = SubstitutionString("\\1\\2\\3\n      \\5\n\\1end")
    code = replace(code, r => s)
    return code
end


const newlinereplacements = OrderedDict(
    # \r\n -> \n, \r -> \n
    r"\r\n" => @s_str("\n"),
    r"\r" => @s_str("\n")
)

const multilinereplacements = OrderedDict(
    #### Remove '&' / '.' multiline continuations in fixed-form
    ###r"^     [.+&$]" => "      ",
    # Labels
    r"^(\h*)(\d+)\h+(.*)$"m => @s_str("    @label L\\2\n\n    \\3"),
)

# Regex/substitution pairs for replace(). Order matters here.
const replacements = OrderedDict(
    # Remove '&' / '.' multiline continuations
    r"^     [.+&$]"m => "      ",
    # Powers use ^ not **
    #r"\*\*([^*\n]+)" => s"^\1",
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
    # DO loop to for loop
    r"do (.*),(.*)"i => s"for \1:\2",
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
    # Replace IF followed by ( to if (
    r"^(\s*)(elseif|if)\("mi => s"\1\2 (",
    # Remove THEN
    r"([)\s])then(\s*)"i => s"\1\2",
    # Replace IF with if
    r"^(\s*)IF(?=\W)"m => s"\1if",
    # Replace ELSE with else
    #r"^(\s*)ELSE(?=\W)"m => s"\1else",
    r"^(\s*)ELSE"m => s"\1else",
    # Relace END XXXX with end
    r"^(\h*)end\h*?\w*?(\h*#.*|\h*)$"mi => s"\1end\2",
    # Replace expnent function
    r"(\W)exp\("i => s"\1exp(",
    # Don't need CALL
    r"(\s*)call(\h+)" => s"\1",
    r"(\s*)CALL(\h+)" => s"\1",
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
    r"(?<=[^\s=])=(?=[^\s=])" => " = ",
    # Add end after single line if with an = assignment
    r"(?<=\s)if\s*([\(].*?) = (.*?)(\s*#.*|)$"mi => s"if \1 = \2 end\3",
    # Single-line IF statement with various terminations
    r"(^\h*)if\s*(.*?)\s*return(\s*#.*|\s*)$"mi => s"\1\2 && return\3",
    r"(^\h*)if\s*(.*?)\s*cycle(\s*#.*|\s*)$"mi => s"\1\2 && continue\3",
    r"(^\h*)if\s*(.*?)\s*goto\s+(\d+\h*)(?!end)(.*)"i => s"\1\2 && @goto L\3\4",
    # Remove expression's brackets after if/elseif/while
    r"^(\h*)(if|elseif|while)(\h*)\(([^#]+)\)\h*$"m => s"\1\2\3\4",
    r"^(\h*)(if|elseif|while)(\h*)\(([^#]+)\)(\h*#.*)$"m => s"\1\2\3\4\5",
    # Process PARAMETER
    r"^(\h*)parameter(\h+)\((.*)\)(\h*?#.*|\h*?)$"mi => s"\1\2\3\4",
    r"^(\h*).*;\h*parameter\h*::\h*(.*)(\h*?#.*|\h*?)$"mi => s"\1\2\3",
    # Some specific functions
    r"\bsign\b\("i => "copysign(",
    r"(?<=\W)MAX\(" => "max(",
    r"(?<=\W)MIN\(" => "min(",
    r"(?<=\W)ABS\(" => "abs(",
    r"//" => " * ",
    # Format floats as "5.0" not "5."
    r"(\W\d+)\.(\D)" => s"\1.0\2",
    # Goto
    r"^(\s*)goto\s+(\d+)"mi => s"\1@goto L\2",
    r"(\h)goto\s+(\d+)"mi => s"\1@goto L\2",
    r"(\W)goto\s+(\d+)"mi => s"\1 @goto L\2",
    # Tab to 4 spaces
    r"\t" => "    ",
    # Relace suberror with error and mark for fixup
    r"(\W)suberror\((.*?),.*?\)" => s"\1 error(\2)",
    # Mark #FIXME the various things this script can't handle
    #r"(write|while\s)" => s"#FIXME \1",
    # Implicit declaration
    r"$\h*implicit none"mi => s"",
    # Custom functions
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
)

# Patterns to remove
const removal = [
    ## Trailing whitespace
    #r"\h*$"m,
    # Variable declarations
    r"\n\s*real\s.*"i,
    r"\n\s*real, external\s.*"i,
    r"\n\s*integer\*?\d*\s.*"i,
    r"\n\s*implicit none"i,
    r"\n\s*logical\s.*"i,
    r"\n\s*character\s.*"i,
    r"\n\s*character[*].*"i,
    r"\n\s*external\s.*"i,
    r"\n\s*intrinsic\s.*"i,
    #r"\n\s*parameter\s.*"i,
    # Import statements
    r"\n\s*use\s.*"i,
]

function main()

    # Load the file from the first command line argument
    filename = string(ARGS[1])
    f = open(filename)
    code = read(f, String)

    # process newlines
    for r in newlinereplacements
        code = replace(code, r)
    end

    # convert comments
    code = convertfortrancomments(code)

    # replace array's brackets with square braces
    code = replacearraysbrackets(code)

    # comment out all vars declarations
    code = commentoutdeclarations(code)

    # process multiline replacements
    for r in multilinereplacements
        code = replace(code, r)
    end

    # IF one-liners
    code = processonelineif(code)

    # select case
    code = processselectcase(code)

    # parameters
    code = processparameters(code)

    lines = splitonlines(code)

    # most syntax conversions
    for r in replacements
        for i = 1:length(lines)
            lines[i] = replace(lines[i], r)
        end
    end

    # concate string lines back together and
    # create functions declaration with custom header
    code = foldl((a,b) -> a*'\n'*b, lines)
    for r in headersprocessing
        len = length(code)
        while true
            code = replace(code, r)
            len == length(code) && break
            len = length(code)
        end
    end

    # some removing
    for r in removal
        code = replace(code, r => "")
    end

    # TODO
    # :: declarations
    # CHARACTER *NN
    # PARAMETER: Done
    # SWITCH CASE: Done
    # include: Done
    # write(*,*)
    # write(*,*) && STOP
    # INPUT:
    # OUTPUT:
    # P_VALUE
    # RWK/IWK
    # M_GETMEM
    # MPI_
    # LSAME
    # XERBLA
    # DO NN I=

    # Write the output to a .jl file with the same filename stem.
    stem = split(filename, ".")[1]
    outfile = stem * ".jl"
    write(outfile, code)

end

main()

