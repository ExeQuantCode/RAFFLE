from ase.io import write
from ase.build import bulk

hcp_carbon = bulk("C", "hcp", a=3.567)
write("hcp_carbon.xyz", hcp_carbon)
