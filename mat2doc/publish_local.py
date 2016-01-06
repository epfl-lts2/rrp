#!/usr/bin/python

import platform

if platform.system() == 'Darwin':
    homefolder = '/Users/nati'
else:
    homefolder = '/home/nati'


project=homefolder+'/work/git/toolbox/r2p/'

# Configure HTML placement at remote server
user = 'nati'

# -------- Configuration of mat2doc ------------
mat2docpath=homefolder+'/work/git/mat2doc'

# Configure matlab paths
unlocxpath = homefolder+'/work/git/toolbox/unlocbox';
ltfatpath = homefolder+'/work/git/toolbox/ltfat';
gsppath = homefolder+'/work/git/toolbox/gspbox';
rrppath = homefolder+'/work/git/toolbox/r2p';



