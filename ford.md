project:
summary: A Fortran library and executable for structure prediction at material interfaces
src_dir: ./src/fortran
    ./app
exclude: **/f90wrap_*.f90
output_dir: doc/html
preprocess: false
predocmark: !!
fpp_extensions: f90
display: public
         protected
         private
source: true
graph: true
md_extensions: markdown.extensions.toc
coloured_edges: true
sort: permission-alpha
author: Ned Thaddeus Taylor
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/nedtaylor/raffle
project_download: https://github.com/nedtaylor/raffle/releases
github: https://github.com/nedtaylor

{!README.md!}