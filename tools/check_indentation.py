import os
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
    inside_do_concurrent_limits = False
    interface_block = False

    expected_indent = 0  # Default expected indentation
    continuation_line = False  # Flag to indicate if the previous line was a continuation
    open_bracket_count = 0  # Count of unbalanced open brackets
    close_bracket_count = 0  # Count of unbalanced close brackets
    unbalanced_brackets = 0
    continued_indent = 0
    num_single_quotes = 0
    num_double_quotes = 0
    unbalanced_quotes = False

    success = True

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            
            # check length of line does not exceed
            if len(line) > 81:
                print(f"Line {line_num} in {file_path} exceeds 80 characters: {len(line)}")
                success = False

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

            # Check for lines with ! in first column and skip
            if re.match(r'^!', stripped_line):
                continue

            # If unbalanced quotes, check line starts with ampersand
            if unbalanced_quotes:
                if not re.match(r'^\s*&', stripped_line):
                    print(f"Unbalanced quotes in {file_path}, line {line_num}")
                    print(stripped_line)
                    print("Number of single quotes: ", num_single_quotes)
                    print("Number of double quotes: ", num_double_quotes)
                    return False

            # Replace all numbers at the start of the line with the same number of spaces
            stripped_line = re.sub(r'^\d+', lambda x: ' ' * len(x.group()), stripped_line)

            # Check if line starts with comment
            if re.match(r'^\s*!', stripped_line):
                actual_indent = len(stripped_line) - len(stripped_line.lstrip())
                if not check_if_match(actual_indent, expected_indent, continued_indent, continuation_line, line_num, file_path):
                    return False
                continue

            # Handle previous line with unbalanced quotes
            first_single_quote_pos = 0
            first_double_quote_pos = 0
            if num_single_quotes % 2 != 0:
                first_single_quote_pos = stripped_line.find("'") + 1
            if num_double_quotes % 2 != 0:
                first_double_quote_pos = stripped_line.find('"') + 1
            stripped_line_excld_quote = stripped_line[max(first_single_quote_pos, first_double_quote_pos):]

            # Check if line contains unbalanced quotes by counting the number of quotes
            if re.search(r"(?<!\\)'", stripped_line):
                num_single_quotes += stripped_line.count("'") - stripped_line.count("\\'")
            if re.search(r'(?<!\\)"', stripped_line):
                num_double_quotes += stripped_line.count('"') - stripped_line.count('\\"')

            # Check for quoted quotes and remove from the count, i.e. "'" or '"'
            quoted_quotes = re.findall(r'"[^"]*\'[^"]*"', stripped_line)
            if len(quoted_quotes) > 0:
                count = 0
                for q in quoted_quotes:
                    count += q.count(r'"')
                num_single_quotes -= count
            quoted_quotes = re.findall(r"'[^']*\"[^']*'", stripped_line)
            if len(quoted_quotes) > 0:
                count = 0
                for q in quoted_quotes:
                    count += q.count(r'"')
                num_double_quotes -= count

            # Handle unbalanced quotes
            last_single_quote_pos = len(stripped_line_excld_quote)
            last_double_quote_pos = len(stripped_line_excld_quote)
            if num_single_quotes % 2 != 0:
                last_single_quote_pos = stripped_line_excld_quote.rfind("'")
            if num_double_quotes % 2 != 0:
                last_double_quote_pos = stripped_line_excld_quote.rfind('"')
            stripped_line_excld_quote = stripped_line_excld_quote[:min(last_single_quote_pos, last_double_quote_pos)+1]
            
            # Check if line contains unbalanced quotes by counting the number of quotes
            if num_single_quotes % 2 != 0 or num_double_quotes % 2 != 0:
                unbalanced_quotes = True
            else:
                unbalanced_quotes = False
                num_single_quotes = 0
                num_double_quotes = 0

            # Check if line starts with close bracket, if so, update the indentation
            if re.match(r'^\s*[)\]]', stripped_line_excld_quote):
                continued_indent = expected_indent + ( unbalanced_brackets - 1 ) * continuation_indent

            # Count open and close brackets
            open_bracket_count += stripped_line_excld_quote.count('(')
            open_bracket_count += stripped_line_excld_quote.count('[')
            close_bracket_count += stripped_line_excld_quote.count(')')
            close_bracket_count += stripped_line_excld_quote.count(']')
            if stripped_line_excld_quote.count('(') + stripped_line_excld_quote.count('[') < \
                    stripped_line_excld_quote.count(')') + stripped_line_excld_quote.count(']'):
                unbalanced_brackets -= 1
            elif stripped_line_excld_quote.count('(') + stripped_line_excld_quote.count('[') > \
                    stripped_line_excld_quote.count(')') + stripped_line_excld_quote.count(']'):
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
               ( inside_select and re.match(r'^\s*(case|class|rank)\s+default\b', stripped_line, re.IGNORECASE) ):
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

            # Detect end of block block
            if re.match(r'^\s*end\s*block\b', stripped_line, re.IGNORECASE):
                expected_indent -= procedure_indent

            # Detect end of associate block
            if re.match(r'^\s*end\s*associate\b', stripped_line, re.IGNORECASE):
                expected_indent -= loop_conditional_indent

            # Detect end of interface block
            if re.match(r'^\s*end\s*interface\b', stripped_line, re.IGNORECASE):
                interface_block = False
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
                if unbalanced_quotes:
                    print(f"Unbalanced quotes in {file_path}, line {line_num}")
                    return False
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
            if not continuation_line and specifier_line:
                specifier_line = False
                expected_indent = prior_indent
            

            # Detect module or program blocks (specifically avoid module procedure/function)
            if re.match(r'^\s*module\b(?!\s+(procedure|function))', stripped_line, re.IGNORECASE) or \
                re.match(r'^\s*program\b', stripped_line, re.IGNORECASE):
                expected_indent += module_program_indent

            # Detect procedure blocks, can be "module (function|subroutine|procedure)" or "(function|subroutine|procedure)" but not "procedure," or "procedure ::"
            if re.match(r'^\s*(integer\s+|logical\s+)?(pure\s+)?(module\s+)?(recursive\s+)?(function|subroutine|procedure)\b', stripped_line, re.IGNORECASE) and \
                not re.match(r'^\s*(function|subroutine|procedure)\s*(,|::)', stripped_line, re.IGNORECASE) and \
                not ( interface_block and re.match(r'^\s*procedure\s+\w+,', stripped_line, re.IGNORECASE) ):
                # print(stripped_line)
                # print(re.match(r'^\s*(integer\s+)', stripped_line, re.IGNORECASE), stripped_line)
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
                interface_block = True
                expected_indent += loop_conditional_indent

            # Detect associate block
            if re.match(r'^\s*associate\b', stripped_line, re.IGNORECASE):
                if stripped_line.lower().endswith("&"):
                    inside_associate_arguments = True
                else:
                    expected_indent += loop_conditional_indent
            
            # Detect block block
            if re.match(r'^\s*block\b', stripped_line, re.IGNORECASE):
                expected_indent += procedure_indent

            # Detect do loop, and where statement with optional "NAME:"
            if re.match(r'^\s*\w+\s*:\s*(do|where)\b', stripped_line, re.IGNORECASE) or \
               re.match(r'^\s*(do|where)\b', stripped_line, re.IGNORECASE):
                expected_indent += loop_conditional_indent
                inside_loop_conditional = True

            # Detect do concurrent linebreak statement with optional "NAME:"
            if re.match(r'^\s*\w+\s*:\s*do\s+concurrent\b', stripped_line, re.IGNORECASE) or \
               re.match(r'^\s*do\s+concurrent\b', stripped_line, re.IGNORECASE):
                inside_do_concurrent_limits = True
                expected_indent -= loop_conditional_indent
                inside_loop_conditional = False
            
            if inside_do_concurrent_limits and not continuation_line:
                if stripped_line.lower().endswith(")"):
                    expected_indent += loop_conditional_indent
                inside_do_concurrent_limits = False
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



    return success


def main():
    # Check directory of script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir+"/..")
    
    # Get list of .f90 and .F90 files from ../src/fortran and ../app
    fortran_files = glob.glob("src/fortran/**/*.f90", recursive=True) + glob.glob("src/fortran/**/*.F90", recursive=True)
    app_files = glob.glob("app/**/*.f90", recursive=True) + glob.glob("app/**/*.F90", recursive=True)
    all_files = fortran_files + app_files

    # Check indentation of all files
    success = True
    for file_path in all_files:
        if not check_indentation(file_path):
            success = False
        else:
            print(f"{file_path} passed indentation check.")

    # Exit with error if any file failed indentation check
    if not success:
        sys.exit("Indentation check failed.")
    else:
        print("All files passed indentation check.")

if __name__ == "__main__":
    main()
