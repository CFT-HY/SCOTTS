#!/usr/bin/env python3

import sys, string, numpy

if len(sys.argv) in [1,4,5] or len(sys.argv) > 6:
    sys.stderr.write('''
Usage:
  python demo.py <silo file 1> <silo file 2> <mesh> <var> <output silo file>

Reads silo files 1 and 2, and diffs var in both files (laid out on  mesh)
then writes the diff (and mesh) to the output file.

If you specify just one file, lists meshes and vars.
If you specify two files, lists meshes and vars common to both.

Requires pyvisfile: https://mathema.tician.de/software/pyvisfile/
''')
    sys.exit(1)

from pyvisfile.silo import SiloFile, DB_READ, DB_CLOBBER, DB_NODECENT

db1 = SiloFile(sys.argv[1], create=False, mode=DB_READ)

qmesh1 = db1.get_toc().qmesh_names
qvar1 = db1.get_toc().qvar_names

sys.stderr.write('Quadmeshes in db1: '
                 + ' '.join(qmesh1) + '\n')
sys.stderr.write('Vars in db1: '
                 + ' '.join(qvar1) + '\n')

if len(sys.argv) == 2:
    sys.exit(0)

sys.stderr.write('\n')
    
db2 = SiloFile(sys.argv[2], create=False, mode=DB_READ)

qmesh2 = db2.get_toc().qmesh_names
qvar2 = db2.get_toc().qvar_names

sys.stderr.write('Quadmeshes in db2: '
                 + ' '.join(qmesh2) + '\n')
sys.stderr.write('Vars in db2: '
                 + ' '.join(qvar2) + '\n')

sys.stderr.write('\n')

sys.stderr.write('Common to both:\n')

qmesh_common = list(set(qmesh1) & set(qmesh2))
qvar_common = list(set(qvar1) & set(qvar2))

sys.stderr.write('Quadmeshes: ' + ' '.join(qmesh_common) + '\n')
sys.stderr.write('Quadvars: ' + ' '.join(qvar_common) + '\n')

if len(sys.argv) == 3:
    sys.exit(0)

mdiff = sys.argv[3]
vdiff = sys.argv[4]
outfile = sys.argv[5]

sys.stderr.write('\n')

sys.stderr.write('Going to diff mesh "%s" with var "%s" into file "%s" \n'
                 % (mdiff, vdiff, outfile))

if not mdiff in qmesh_common:
    sys.stderr.write("Chosen mesh not in both files!\n")
    sys.exit(1)

if not vdiff in qvar_common:
    sys.stderr.write("Chosen var not in both files!\n")
    sys.exit(1)

f = SiloFile(outfile, mode=DB_CLOBBER)

meshin = db1.get_quadmesh(mdiff)
f.put_quadmesh(mdiff, meshin.coords)

qvar1 = db1.get_quadvar(vdiff)
qvar2 = db2.get_quadvar(vdiff)


qvar_diff = qvar1.vals[0] - qvar2.vals[0]
f.put_quadvar1(vdiff + "diff", mdiff,
              numpy.asarray(qvar_diff, order="F"),
              qvar_diff.shape,
              centering=DB_NODECENT)

f.close()
sys.stderr.write('Written diff to file "%s" as variable "%sdiff"\n'
                 % (outfile, vdiff))
