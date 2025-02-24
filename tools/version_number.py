import re

def update_version(new_version):
    # Update fpm.toml
    with open('fpm.toml', 'r') as file:
        content = file.read()
    content = re.sub(r'version = "\d+\.\d+\.\d+.*"', f'version = "{new_version}"', content)
    with open('fpm.toml', 'w') as file:
        file.write(content)

    # Update Fortran module
    with open('src/fortran/lib/mod_io_utils.F90', 'r') as file:
        content = file.read()
    content = re.sub(r'character\(len=\*\), parameter :: raffle__version__ = "\d+\.\d+\.\d+.*"', f'character(len=*), parameter :: raffle__version__ = "{new_version}"', content)
    with open('src/fortran/lib/mod_io_utils.F90', 'w') as file:
        file.write(content)

def get_version():
    # get the version number from fpm.toml
    with open('fpm.toml', 'r') as file:
        content = file.read()
        match = re.search(r'version = "(\d+\.\d+\.\d+.*)"', content)
        print(match.group(1))
        if match:
            return match.group(1)

if __name__ == '__main__':
    update_version(get_version())
