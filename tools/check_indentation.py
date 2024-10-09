import re
import sys
import glob


def check_if_match(actual_indent, expected_indent, continued_indent, continuation_line, line_num, file_path):
    if continuation_line:
        if actual_indent != continued_indent:
            print(f"Indentation error in {file_path}, line {line_num}: "
                f"Expected {continued_indent} spaces, found {actual_indent}")
            return False
    else:
        if actual_indent != expected_indent:
            print(f"Indentation error in {file_path}, line {line_num}: "
                f"Expected {expected_indent} spaces, found {actual_indent}")
            return False
    return True

def check_indentation(file_path):
    procedure_indent = 2
    module_program_indent = 2
    loop_conditional_indent = 3
    continuation_indent = 5

    inside_module_program = False
    inside_procedure = False
    inside_loop_conditional = False
    inside_select = False
    specifier_line = False
    inside_derived_type = False
    on_continued_if_line = False
    inside_procedure_arguments = False
    inside_associate_arguments = False

    expected_indent = 0  # Default expected indentation
    continuation_line = False  # Flag to indicate if the previous line was a continuation
    open_bracket_count = 0  # Count of unbalanced open brackets
    close_bracket_count = 0  # Count of unbalanced close brackets
    unbalanced_brackets = 0
    continued_indent = 0

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            stripped_line = line.rstrip()

            # Skip empty lines
            if not stripped_line:
                continuation_line = False
                continue

            # If comment line !###, skip
            if re.match(r'^\s*!###', stripped_line):
                continue

            # Check if preprocessing directive
            if re.match(r'^\s*#', stripped_line):
                continue

            # Replace all numbers at the start of the line with the same number of spaces
            stripped_line = re.sub(r'^\d+', lambda x: ' ' * len(x.group()), stripped_line)

            # Check if line starts with comment
            if re.match(r'^\s*!', stripped_line):
                actual_indent = len(stripped_line) - len(stripped_line.lstrip())
                if not check_if_match(actual_indent, expected_indent, continued_indent, continuation_line, line_num, file_path):
                    return False
                continue

            # Check if line starts with close bracket, if so, update the indentation
            if re.match(r'^\s*[)\]]', stripped_line):
                continued_indent = expected_indent + ( unbalanced_brackets - 1 ) * continuation_indent

            # Count open and close brackets
            open_bracket_count += stripped_line.count('(')
            open_bracket_count += stripped_line.count('[')
            close_bracket_count += stripped_line.count(')')
            close_bracket_count += stripped_line.count(']')
            if stripped_line.count('(') + stripped_line.count('[') < stripped_line.count(')') + stripped_line.count(']'):
                unbalanced_brackets -= 1
            elif stripped_line.count('(') + stripped_line.count('[') > stripped_line.count(')') + stripped_line.count(']'):
                unbalanced_brackets += 1


            # Detect end of do loop, if statement, or where statement
            if re.match(r'^\s*end\s*(do|if|where|select)\b', stripped_line, re.IGNORECASE):
                expected_indent -= loop_conditional_indent

            # Detect else statements in if and where blocks, can be "PATTERN", "PATTERN\s*if", or "PATTERN\s*where"
            if inside_loop_conditional and re.match(r'^\s*else\s*(if|where)?\b', stripped_line, re.IGNORECASE):
                prior_indent = expected_indent
                expected_indent -= loop_conditional_indent
                specifier_line = True

            # Detect case, type, and rank statements within select, can be "PATTERN(", "PATTERN (" or "PATTERN default"
            if ( inside_select and re.match(r'^\s*(case|type is|rank)\s*\(', stripped_line, re.IGNORECASE) ) or \
               ( inside_select and re.match(r'^\s*(case|type|rank)\s+default\b', stripped_line, re.IGNORECASE) ):
                prior_indent = expected_indent
                expected_indent -= loop_conditional_indent
                specifier_line = True

            # Detect if contains line
            if re.match(r'^\s*contains\b', stripped_line, re.IGNORECASE):
                prior_indent = expected_indent
                specifier_line = True
                if inside_derived_type:
                    expected_indent -= loop_conditional_indent - 1
                else:
                    expected_indent -= module_program_indent

            # Detect end of associate block
            if re.match(r'^\s*end\s*associate\b', stripped_line, re.IGNORECASE):
                expected_indent -= loop_conditional_indent

            # Detect end of interface block
            if re.match(r'^\s*end\s*interface\b', stripped_line, re.IGNORECASE):
                expected_indent -= loop_conditional_indent

            # Detect end of derived type block
            if inside_derived_type and re.match(r'^\s*end\s*type\b', stripped_line, re.IGNORECASE):
                expected_indent -= loop_conditional_indent
                inside_derived_type = False

            # Detect end of procedure block
            if inside_procedure and re.match(r'^\s*end\s*(function|subroutine|procedure)\b', stripped_line, re.IGNORECASE):
                expected_indent -= procedure_indent
                inside_procedure = False

            # Detect end of module or program
            if re.match(r'^\s*end\s*(module|program)\b', stripped_line, re.IGNORECASE):
                expected_indent -= module_program_indent



            # Check actual indentation
            actual_indent = len(stripped_line) - len(stripped_line.lstrip())
            if not check_if_match(actual_indent, expected_indent, continued_indent, continuation_line, line_num, file_path):
                return False
            

            # strip comments from end of line
            stripped_line = re.sub(r'!.*', '', stripped_line).strip()

            if continuation_line:
                continued_indent = expected_indent + unbalanced_brackets * continuation_indent

            # Check for line continuation character
            if stripped_line.endswith('&'):
                if not continuation_line:
                    continuation_line = True
                    # Set expected indentation for next line
                    continued_indent = expected_indent + continuation_indent
                    if unbalanced_brackets == 0:
                        unbalanced_brackets = 1
            else:
                # If it was a continuation line, reset to normal expected indentation
                open_bracket_count = 0
                close_bracket_count = 0
                unbalanced_brackets = 0
                if continuation_line:
                    continuation_line = False
                if inside_procedure_arguments:
                    inside_procedure_arguments = False
                    expected_indent += procedure_indent
                if inside_associate_arguments:
                    inside_associate_arguments = False
                    expected_indent += loop_conditional_indent
                    

            # Reset from contains line
            if specifier_line:
                specifier_line = False
                expected_indent = prior_indent
            

            # Detect module or program blocks (specifically avoid module procedure/function)
            if re.match(r'^\s*module\b(?!\s+(procedure|function))', stripped_line, re.IGNORECASE) or \
                re.match(r'^\s*program\b', stripped_line, re.IGNORECASE):
                expected_indent += module_program_indent

            # Detect procedure blocks, can be "module (function|subroutine|procedure)" or "(function|subroutine|procedure)" but not "procedure," or "procedure ::"
            if re.match(r'^\s*(module\s+)?(function|subroutine|procedure)\b', stripped_line, re.IGNORECASE) and \
                not re.match(r'^\s*(function|subroutine|procedure)\s*(,|::)', stripped_line, re.IGNORECASE):
                inside_procedure = True
                if stripped_line.lower().endswith("&"):
                    inside_procedure_arguments = True
                else:
                    expected_indent += procedure_indent


            # Detect derived type block
            if re.match(r'^\s*type\s*(::|,)', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent
                inside_derived_type = True

            # Detect interface block, can be "abstract interface" or "interface"
            if re.match(r'^\s*(abstract\s+)?interface\b', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent

            # Detect associate block
            if re.match(r'^\s*associate\b', stripped_line, re.IGNORECASE):
                if stripped_line.lower().endswith("&"):
                    inside_associate_arguments = True
                else:
                    expected_indent += loop_conditional_indent

            # Detect do loop, and where statement with optional "NAME:"
            if re.match(r'^\s*\w+\s*:\s*(do|where)\b', stripped_line, re.IGNORECASE) or \
               re.match(r'^\s*(do|where)\b', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent
                inside_loop_conditional = True

            # Detect "if RANDOM then" statement with optional "NAME:"
            if re.match(r'^\s*\w*\s*:\s*if\s*\(.*\)\s*then\b', stripped_line, re.IGNORECASE) or \
               re.match(r'^\s*if\s*\(.*\)\s*then\b', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent
                inside_loop_conditional = True

            # Detect line ends with "then" from unfinished if statement 
            if on_continued_if_line and not continuation_line:
                if stripped_line.lower().endswith("then"):
                    expected_indent += loop_conditional_indent
                on_continued_if_line = False

            # Detect if linebreak statement with optional "NAME:"
            if continuation_line and \
               ( re.match(r'^\s*\w+\s*:\s*(if)\b', stripped_line, re.IGNORECASE) or \
                 re.match(r'^\s*(if)\b', stripped_line, re.IGNORECASE) ):
                on_continued_if_line = True

            # Detect select type, select case, and select rank
            if re.match(r'^\s*select\s+(type|case|rank)\b', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent
                inside_select = True



    return True


def main():
    # Check all Fortran files in the directory
    for file_path in glob.glob("**/*.f90", recursive=True):
        if not check_indentation(file_path):
            sys.exit(1)
        else:
            print(f"{file_path} passed indentation check.")
    print("All files passed indentation check.")

if __name__ == "__main__":
    main()
