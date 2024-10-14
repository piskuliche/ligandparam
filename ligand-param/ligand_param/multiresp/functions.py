def WriteFitSh(base):
    import subprocess
    sh = open("%s.resp.sh"%(base),"w")
    sh.write("""#!/bin/bash
out=%s.resp
cat <<EOF > ${out}.qwt
1
1.

EOF

resp -O -i ${out}.inp -o ${out}.out -p ${out}.punch -t ${out}.qout -w ${out}.qwt -e ${out}.esp

rm -f ${out}.punch ${out}.qout ${out}.qwt ${out}.esp

"""%(base))
    sh.close()
    subprocess.call(["bash","%s.resp.sh"%(base)])