import subprocess
import os


args = ['git', 'clone', '--depth=1', 'https://github.com/acycliq/SfN2018.git', 'D:\crap\SfN2018']
res = subprocess.Popen(args, stdout=subprocess.PIPE)
output, _error = res.communicate()

if not _error:
    print(output)
else:
    print(_error)