from pyscf import gto, scf, ao2mo, lo
from pyscf.lo.boys import Boys
import sys
import numpy as np
import tempfile

'''
A simple example to call integral transformation for given orbitals
'''

mol = gto.Mole()

mol.basis = {'Ar': gto.basis.parse('''
Ar    S
      5.450000E+05           4.558280E-05          -1.295510E-05           0.000000E+00           4.049900E-06           0.000000E+00
      8.164000E+04           3.541080E-04          -1.004280E-04           0.000000E+00           3.136910E-05           0.000000E+00
      1.858000E+04           1.857970E-03          -5.295830E-04           0.000000E+00           1.656460E-04           0.000000E+00
      5.261000E+03           7.768510E-03          -2.213960E-03           0.000000E+00           6.916620E-04           0.000000E+00
      1.717000E+03           2.742320E-02          -7.968450E-03           0.000000E+00           2.497900E-03           0.000000E+00
      6.199000E+02           8.238360E-02          -2.458030E-02           0.000000E+00           7.710740E-03           0.000000E+00
      2.416000E+02           2.012300E-01          -6.577980E-02           0.000000E+00           2.087140E-02           0.000000E+00
      9.979000E+01           3.567810E-01          -1.379420E-01           0.000000E+00           4.439650E-02           0.000000E+00
      4.315000E+01           3.495630E-01          -2.016300E-01           0.000000E+00           6.802240E-02           0.000000E+00
      1.914000E+01           1.182660E-01          -4.128340E-02           0.000000E+00           1.413500E-02           0.000000E+00
      7.488000E+00           5.601900E-03           4.846800E-01           0.000000E+00          -2.074890E-01           0.000000E+00
      3.205000E+00           4.834730E-04           5.792240E-01           0.000000E+00          -4.250450E-01           0.000000E+00
      1.196000E+00          -1.480850E-04           8.790830E-02           1.000000E+00           7.048550E-02           0.000000E+00
      5.204000E-01           2.920250E-05          -7.275530E-03           0.000000E+00           7.336270E-01           0.000000E+00
      1.954000E-01          -2.316040E-05           2.328840E-03           0.000000E+00           3.960050E-01           1.000000E+00
Ar    S
      0.0685000              1.0000000
Ar    P
      7.618000E+02           2.369760E-03           0.000000E+00          -6.672110E-04           0.000000E+00
      1.802000E+02           1.901990E-02           0.000000E+00          -5.327170E-03           0.000000E+00
      5.750000E+01           8.808070E-02           0.000000E+00          -2.554940E-02           0.000000E+00
      2.124000E+01           2.563770E-01           0.000000E+00          -7.571970E-02           0.000000E+00
      8.388000E+00           4.387110E-01           0.000000E+00          -1.411330E-01           0.000000E+00
      3.416000E+00           3.475690E-01           0.000000E+00          -9.327680E-02           0.000000E+00
      1.206000E+00           5.667400E-02           1.000000E+00           2.828720E-01           0.000000E+00
      4.523000E-01          -5.238820E-03           0.000000E+00           5.624500E-01           0.000000E+00
      1.545000E-01           1.643760E-03           0.000000E+00           3.250590E-01           1.000000E+00
Ar    P
      0.0487000              1.0000000
Ar    D
      1.254000E+00           1.000000E+00           0.000000E+00
      4.100000E-01           0.000000E+00           1.000000E+00
Ar    D
      0.1690000              1.0000000
Ar    F
      8.900000E-01           1.0000000
Ar    F
      0.4060000              1.0000000
Ar    S
      0.245645              1.0000000
Ar    P
      0.430082              1.0000000
Ar    D
      0.622557              1.0000000
Ar    S
      0.098496              1.0000000
Ar    P
      0.169341              1.0000000
Ar    D
      0.242160              1.0000000
Ar    S
      0.052725              1.0000000
Ar    P
      0.089894              1.0000000
Ar    D
      0.127840              1.0000000
'''),
 'K': gto.basis.parse('''
Ar    S
      5.450000E+05           4.558280E-05          -1.295510E-05           0.000000E+00           4.049900E-06           0.000000E+00
      8.164000E+04           3.541080E-04          -1.004280E-04           0.000000E+00           3.136910E-05           0.000000E+00
      1.858000E+04           1.857970E-03          -5.295830E-04           0.000000E+00           1.656460E-04           0.000000E+00
      5.261000E+03           7.768510E-03          -2.213960E-03           0.000000E+00           6.916620E-04           0.000000E+00
      1.717000E+03           2.742320E-02          -7.968450E-03           0.000000E+00           2.497900E-03           0.000000E+00
      6.199000E+02           8.238360E-02          -2.458030E-02           0.000000E+00           7.710740E-03           0.000000E+00
      2.416000E+02           2.012300E-01          -6.577980E-02           0.000000E+00           2.087140E-02           0.000000E+00
      9.979000E+01           3.567810E-01          -1.379420E-01           0.000000E+00           4.439650E-02           0.000000E+00
      4.315000E+01           3.495630E-01          -2.016300E-01           0.000000E+00           6.802240E-02           0.000000E+00
      1.914000E+01           1.182660E-01          -4.128340E-02           0.000000E+00           1.413500E-02           0.000000E+00
      7.488000E+00           5.601900E-03           4.846800E-01           0.000000E+00          -2.074890E-01           0.000000E+00
      3.205000E+00           4.834730E-04           5.792240E-01           0.000000E+00          -4.250450E-01           0.000000E+00
      1.196000E+00          -1.480850E-04           8.790830E-02           1.000000E+00           7.048550E-02           0.000000E+00
      5.204000E-01           2.920250E-05          -7.275530E-03           0.000000E+00           7.336270E-01           0.000000E+00
      1.954000E-01          -2.316040E-05           2.328840E-03           0.000000E+00           3.960050E-01           1.000000E+00
Ar    S
      0.0685000              1.0000000
Ar    P
      7.618000E+02           2.369760E-03           0.000000E+00          -6.672110E-04           0.000000E+00
      1.802000E+02           1.901990E-02           0.000000E+00          -5.327170E-03           0.000000E+00
      5.750000E+01           8.808070E-02           0.000000E+00          -2.554940E-02           0.000000E+00
      2.124000E+01           2.563770E-01           0.000000E+00          -7.571970E-02           0.000000E+00
      8.388000E+00           4.387110E-01           0.000000E+00          -1.411330E-01           0.000000E+00
      3.416000E+00           3.475690E-01           0.000000E+00          -9.327680E-02           0.000000E+00
      1.206000E+00           5.667400E-02           1.000000E+00           2.828720E-01           0.000000E+00
      4.523000E-01          -5.238820E-03           0.000000E+00           5.624500E-01           0.000000E+00
      1.545000E-01           1.643760E-03           0.000000E+00           3.250590E-01           1.000000E+00
Ar    P
      0.0487000              1.0000000
Ar    D
      1.254000E+00           1.000000E+00           0.000000E+00
      4.100000E-01           0.000000E+00           1.000000E+00
Ar    D
      0.1690000              1.0000000
Ar    F
      8.900000E-01           1.0000000
Ar    F
      0.4060000              1.0000000
Ar    S
      0.245645              1.0000000
Ar    P
      0.430082              1.0000000
Ar    D
      0.622557              1.0000000
Ar    S
      0.098496              1.0000000
Ar    P
      0.169341              1.0000000
Ar    D
      0.242160              1.0000000
Ar    S
      0.052725              1.0000000
Ar    P
      0.089894              1.0000000
Ar    D
      0.127840              1.0000000
'''),
'ghost': gto.basis.load('aug-cc-pvdz', 'Ar'),
'Ar@2': gto.basis.load('aug-cc-pvdz', 'Ar'),
#'K': gto.basis.load('aug-cc-pvdz', 'Ar'),
             }

mol.build(
     #atom = 'K 0 0 0.0; ghost 0 0 5.0; Ar 0 0 10.0',  # in Angstrom
     atom = 'K 0 0 0.0; Ar@2 0 0.0 4.5; Ar 0 0 9.0',  # in Angstrom
     symmetry = True,
     charge = 1,
     spin = 0,
 )

