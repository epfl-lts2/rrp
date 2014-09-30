import sys,os


def dos2unix(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name=os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('dos2unix '+name)

def unix2dos(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name=os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('unix2dos '+name)    

def safe_mkdir(path):
    if  1 != os.path.isdir(path):
        os.system('mkdir '+path)


def do_rebuild_file(source,dest,mode):
    if mode=='rebuild':
        return True

    if not os.path.exists(dest):
        if mode=='cached':
            print 'Error: Cached version of '+dest+ ' does not exist'
            sys.exit()

        print dest +' missing'

        return True

    is_newer = os.path.getmtime(source)>os.path.getmtime(dest)

    if mode=='auto':
        return is_newer

    if mode=='cached':
        return False