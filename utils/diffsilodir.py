import sys, string, numpy, os

from pyvisfile.silo import SiloFile, DB_READ, DB_CLOBBER, DB_NODECENT

if len(sys.argv) != 4:
    sys.stderr.write('''
Usage:
  python demo.py <dir1> <dir 2> <output_dir>

For all silo files with common names in dir1 and dir2 finds the diff of all common var 
(laid out on common mesh) then writes the diff (and mesh) to the output_dir.

Requires pyvisfile: https://mathema.tician.de/software/pyvisfile/
''')
    sys.exit(1)


# Get file list from directory paths

dir1=sys.argv[1]
dir2=sys.argv[2]

dir1_files = os.listdir(dir1)

dir2_files = os.listdir(dir2)

if not os.path.exists(sys.argv[3]):
    sys.stderr.write("Error, {:} does not exist, attempting".format(sys.argv[3])
                      + "to make directory.\n")
    os.mkdir(sys.argv[3])

# Find common silo files

common_files=list(set(dir1_files).intersection(dir2_files))

common_silo=[fname for fname in common_files if ".silo" in fname]

prog = 0



for i in range(len(common_silo)):
    # Find common meshes and variables.    

    db1=SiloFile(dir1+common_silo[i], create=False, mode=DB_READ)

    db2=SiloFile(dir2+common_silo[i], create=False, mode=DB_READ)

    qmesh1 = db1.get_toc().qmesh_names
    qvar1= db1.get_toc().qvar_names

    qmesh2 = db1.get_toc().qmesh_names
    qvar2 = db1.get_toc().qvar_names

    cm_mesh = list(set(qmesh1) & set(qmesh2))
    cm_var = list(set(qvar1) & set(qvar2))
    
    # Find diffs for each file, then save to output directory

    
    
    silo_out = SiloFile(sys.argv[3]+common_silo[i], mode=DB_CLOBBER)
    meshnames=[]
    for var in cm_var:
        var1 = db1.get_quadvar(var)
        var2 = db2.get_quadvar(var)

        if var1.meshname != var2.meshname:
            sys.stderr.write("Error, common var not on same mesh.\n")
            break
        meshin = db1.get_quadmesh(var1.meshname)
        if var1.meshname not in meshnames:
            silo_out.put_quadmesh(var1.meshname,meshin.coords)
            meshnames.append(var1.meshname)
            
        var_diff = var1.vals[0] - var2.vals[0]

        silo_out.put_quadvar1(var + "diff", var1.meshname,
                                  numpy.asarray(var_diff, order = "F"),
                                  var_diff.shape, centering = DB_NODECENT)

    silo_out.close()

    if (i*100/(len(common_silo))%10==0) and (i*100/(len(common_silo))>prog):
        sys.stderr.write("{:}% done\n".format(i*100/(len(common_silo))))
        prog = i*100/(len(common_silo))

                  
