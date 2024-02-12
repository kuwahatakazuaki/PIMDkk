import subprocess
import gzip
import os
import parse_molcas_out
import sys
import shutil
import re

# configuration
PYMOLCAS        = "/home/osamu/.opt/openmolcas/omp_gen1int/pymolcas"
MOLCAS_MEM      = '4096'
OMP_NUM_THREADS = '2'

def gzipfile(infile, outfile):
    with open(infile, r"rb") as f:
        content = f.read()
    with gzip.open(outfile, r"w") as f:
        f.write(content)

#print(f"PWD: {os.environ['PWD']}")
#print(f"sys.argv={sys.argv}")
# scr_parent is not used, for compatibility.
(input_file, output_file, scr_parent, istepsv) = sys.argv[1:]

#print("<< os.environ >>")
#for k in sorted(os.environ.keys()):
#    print(f"{k}={os.environ[k]}")

input_file  = os.path.abspath(input_file)
output_file = os.path.abspath(output_file)
istepsv     = int(istepsv)

dirname     = os.path.dirname(input_file)
if dirname == '':
    dirname = '.'
basename    = os.path.basename(input_file)
basebody    = os.path.splitext(basename)[0]
scr_child   = dirname + "/molcas/"
jobiph      = scr_child + basebody + ".JobIph"
molcas_err  = dirname + "/molcas.err"

#os.environ['OMPI_MCA_initial_workdir'] = dirname
#os.environ['PBS_O_WORKDIR'] = dirname

#if "/00001/" in input_file:
#    print("$ ls -la " + dirname+"/../..")
#    subprocess.run(['ls', '-l', dirname+"/../.."])

#print(f"input_file={input_file}")
#print(f"output_file={output_file}")
#print(f"istepsv={istepsv}")
#print(f"dirname={dirname}")
#print(f"basename={basename}")
#print(f"basebody={basebody}")
#print(f"scr_child={scr_child}")
#print(f"jobiph={jobiph}")

os.environ['WorkDir']         = scr_child
os.environ['MOLCAS_MEM']      = MOLCAS_MEM
os.environ['OMP_NUM_THREADS'] = OMP_NUM_THREADS

os.makedirs(scr_child, exist_ok=True)

# if "/00001/" in input_file:
#     print("ls -la " + dirname)
#     subprocess.run(['ls', '-l', dirname])
#     print("ls -la " + scr_child)
#     subprocess.run(['ls', '-l', scr_child])

# molcas.xyzを使用するかどうかに関わらず設置する。
# symlinkはtargetが存在しなくても作成できる。
symlink2xyz = scr_child + '/molcas.xyz'
if os.path.exists(symlink2xyz):
    os.remove(symlink2xyz)
os.symlink('../molcas.xyz', symlink2xyz)
#if "/00001/" in input_file:
#    print("ls -la " + dirname)
#    subprocess.run(['ls', '-l', dirname])
#    print("ls -la " + scr_child)
#    subprocess.run(['ls', '-l', scr_child])

symlink2jobiph = scr_child + '/molcas.JobIph'
if os.path.exists(symlink2jobiph):
    os.remove(symlink2jobiph)
os.symlink('../molcas.JobIph', symlink2jobiph)

#print("symlink2xyz:    ", symlink2xyz)
#print("symlink2jobiph: ", symlink2jobiph)

#p = subprocess.Popen(['cat', 'grad.asd.out'], stdout=subprocess.PIPE)
# 少なくともversion 3.6.8でos.chdir()はos.environ['PWD']を
# 変更しないので手動で修正する必要がある。
os.environ['OLDPWD'] = os.environ['PWD']
os.environ['PWD'] = scr_child
os.chdir(scr_child)
#print(f'input_file: {input_file}')
with open(molcas_err, r"w") as fe:
    #print(f"pwd: {os.environ['PWD']}")
    #print(f"pwd: {os.getcwd()}")
    p = subprocess.Popen([PYMOLCAS, input_file], stdout=subprocess.PIPE, stderr=fe)
    # ファイルサイズは1MBもないようなので一括でreadする。
    molcas_stdout = p.stdout.read()
with gzip.open(output_file, r"w") as fo:
    fo.write(molcas_stdout)
parse_molcas_out.readlog(molcas_stdout.decode('utf-8').split("\n"))
re_molden = re.compile("\.rasscf\.molden(?:\.\d+)?$")

#if "/00001/" in input_file:
#    print("$ ls -la " + dirname+"/../..")
#    subprocess.run(['ls', '-l', dirname+"/../.."])

for f in os.listdir():
    if re_molden.search(f):
        gzipfile(f,f"{dirname}/{f}.{istepsv:08d}.gz")
        os.remove(f)

