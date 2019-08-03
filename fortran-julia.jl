#!/usr/bin/env julia

#=
This julia script converts fortran 90 code into julia.
It uses naive regex replacements to do as much as possible,
but the output WILL need further cleanup.

Known conversion problems such as GOTO are commented and marked with FIXME

Most variable declaration lines are entirely deleted, which may or
may not be useful. 


To run from a shell: 

julia fortran-julia.jl filename.f90

Output is written to filename.jl.
=#

using DataStructures, Printf

function collect_arrays(code)
    replacements = OrderedDict(
        # Skip comments
        r"\h*!.*" => s"",
        r"^\*.*$"m => s"",
        r"^c.*$"m => s"",
        r"^C.*$"m => s"",
        # Remove '&' / '.' multiline continuations
        #r"^\h*&\h*"m => " ",
        #r"^     \S\h*"m => " ",
        ##r"^\h*\S\h*"m => " ",
        r"\n\h*[.+&$]\h*" => " ",
    )
    matches = [
        # Match declarations
        r"^\h*real\h+[^(]+\("mi,
        r"^\h*integer\h+[^(]+\("mi,
        r"^\h*integer\*.\h+[^(]+\("mi,
        r"^\h*character\h+"mi,
        r"^\h*character[*]"mi,
    ]
    for r in replacements
        len = length(code)
        while true
            code = replace(code, r)
            len == length(code) && break
            len = length(code)
        end
    end
    list = split(code, "\n")
    matched = Vector{String}()
    for r in matches
        #append!(matched, list[findall(i->occursin(r, i), list)])
        for i in list
            occursin(r, i) && push!(matched, strip(i))
        end
    end
    # match only vars names
    for i in 1:length(matched)
        matched[i] = replace(matched[i], r"^\S+\h+(.*)$"m => s"\1") # drop declaration
        matched[i] = replace(matched[i], r" " => s"")               # drop spaces
        matched[i] = replace(matched[i], r"\([^)]+\)" => s"()")     # clean inside braces
        matched[i] = replace(matched[i], r"[^(),]+," => s"")        # drop non-arrays (without braces)
        matched[i] = replace(matched[i], r",[^(),]+$"m => s"")      # drop last non-array (without braces)
    end
    vars = Vector{String}()
    for i in 1:length(matched)
        append!(vars, split(matched[i], ","))
    end
    for i in 1:length(vars)
        vars[i] = replace(vars[i], r" " => s"")
        vars[i] = replace(vars[i], r"^(\w+)\(.*\)$"m => s"\1")
    end

    @show vars

    return vars

end

function processselectcase(lines)
    var = ""
    case1 = false
    for i = 1:length(lines)
        if occursin(r"select\h+case"i, lines[i])
            expr = replace(lines[i], r"^\h*select\h+case\h*\((.*)\)(\h*#.*|\h*)$"mi => s"\1")
            var = @sprintf "EX%s" mod(hash("$(expr)"), 2^16)
            lines[i] = replace(lines[i], r"^(\h*)(select\h+case.*)$"mi => SubstitutionString("\\1"*var*" = "*expr*" #\\2") )
            case1 = true
        end
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
    return lines
end

function processparameters(lines)
    for i = 1:length(lines)
        if occursin(r"^\h*parameter\h+\("mi, lines[i]) ||
           occursin(r"\h*[^,#]+\h*,\h*parameter\h*::\h*"mi, lines[i])
            lines[i] = replace(lines[i], r"," => s";")
        end
    end
    return lines
end

# Regex/substitution pairs for replace(). Order matters here.
const multilinereplacements = OrderedDict(
    # \r\n -> \n, \r -> \n
    r"\r\n" => @s_str("\n"),
    r"\r" => @s_str("\n"),
    # Lowercase everything not commented
    #r"^(?!.*!).*"m => lowercase,
    # Lowercase start of lines with comments
    #r"^.*!"m => lowercase,
    # Remove '&' / '.' multiline continuations
    r"\n\h*[.+&$]\h*" => " ",
    #r"^\h*&\h*"m => " ",
    #r"^     \S\h*"m => " ",
    ##r"(?<=#)\s*&\s*" => "",
    ##r"\n\s*\.\s*" => " ",
    # Labels
    r"^(\h*)(\d+)\h+(.*)$"m => @s_str("    @label L\\2\n\n    \\3"),
)

const replacements = OrderedDict(
    # Comments use # not ! or *
    "!" => "#",
    r"^\*"m => "#",
    r"^c"m => "#",
    r"^C"m => "#",
    #r"^#(\*+)"m => s"#",
    r"====" => s"----",
    r"-===$"m => s"----",
    r"--==$"m => s"----",
    r"---=$"m => s"----",
    r"^#(=+)"m => s"#",
    # Powers use ^ not **
    r"\*\*([^*\n]+)" => s"^\1",
    #r"\*\*" => "^",
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
    r"^(\s*)ELSE(?=\W)"m => s"\1else",
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
    r"\s*\.or\.\s*"i => " || ",
    r"\s*\.and\.\s*"i => " && ",
    r"\s*\.not\.\s*"i => " !",
    r"\s*\.eq\.\s*"i => " == ",
    r"\s*\.ne\.\s*"i => " != ",
    r"\s*\.le\.\s*"i => " <= ",
    r"\s*\.ge\.\s*"i => " >= ",
    r"\s*\.gt\.\s*"i => " > ",
    r"\s*\.lt\.\s*"i => " < ",
    # Fix assignments
    r"(?<=[^\s=])=(?=[^\s=])" => " = ",
    # Add end after single line if with an = assignment
    r"(?<=\s)if\s*([\(].*?) = (.*?)(\s*#.*|)$"mi => s"if \1 = \2 end\3",
    # Single-line IF statement with various terminations
    #r"(?<=\s)if\s*(.*?)\s*return\s*"i => s"\1 && return\2",
    r"(^\h*)if\s*(.*?)\s*return(\s*#.*|\s*)$"mi => s"\1\2 && return\3",
    r"(^\h*)if\s*(.*?)\s*cycle(\s*#.*|\s*)$"mi => s"\1\2 && continue\3",
    r"(^\h*)if\s*(.*?)\s*goto\s*(.*)"i => s"\1\2 && @goto L\3",
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
    # Trailing whitespace
    r"\h*$"m,
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

    # process includes
    #includes = r"^\h*include\h+'([^']+)'(\h*!.*|\h*)$"mi
    #for m in collect(eachmatch(includes, code))
    #    file2include = dirname(filename) * m.captures[1]
    #    code2include = read(open(file2include), String)
    #    code = replace(code,  includes => SubstitutionString(code2include), count=1)
    #end

    # replace with square braces for arrays
    braces = Regex("\\(((?>[^()]|(?R))*)\\)") # complementary braces
    sqbraces = SubstitutionString("[\\1]")
    arrays = collect_arrays(deepcopy(code))
    for a in arrays
        w = Regex("\\b$(a)\\b") # whole word
        for m in collect(eachmatch(w, code))
            o = m.offset
            if code[min(length(code), o+length(a))] == '('
                code = code[1:o-1] * replace(code[o:end], braces => sqbraces, count=1)
            end
        end
    end

    # Process multiline replacements.
    for r in multilinereplacements
        code = replace(code, r)
    end

    # split string on vector of string lines and process each separately
    lines = split(code, '\n')
    lines = processselectcase(lines)
    lines = processparameters(lines)
    for r in replacements
        for i = 1:length(lines)
            #@show lines[i], r
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
    #println(code)

    for r in removal
        code = replace(code, r => "")
    end

    # TODO
    # CHARACTER *NN
    # PARAMETER
    # SWITCH CASE: Done
    # include
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
