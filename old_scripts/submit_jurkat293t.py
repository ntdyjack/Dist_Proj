#script to submit splatter tests looping over a shell script
import os
import sys
import time
#what log-fold changes do I want to use
#length = int(len(range(min,max)) * (1/by))
#lS_vec = [min + float(x) * by for x in range(length)]
info_vec = range(0,10)

for i in info_vec:
    i = str(i)
    out = '-o /users/ndyjack/Dist_Proj/outputs/jurkat293t.nC.' + i + '.out '
    err = '-e /users/ndyjack/Dist_Proj/outputs/jurkat293t.nC.' + i + '.err '
    script ='/users/ndyjack/Dist_Proj/scripts/submit_jurkat293t.sh ' + i 
    cmd = 'qsub ' + err + out + script
#    print(cmd)
    os.system(cmd)
    time.sleep(2)
