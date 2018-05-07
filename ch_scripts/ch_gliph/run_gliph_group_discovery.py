import os


for inpfile in ['a', 'b']:
    input_file = 'cmvtet_gliphed_{}.txt'.format(inpfile)
    cmd = 'perl bin/gliph-group-discovery.pl --tcr={}'.format(input_file)
    os.system(cmd)