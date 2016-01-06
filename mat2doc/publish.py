#!/usr/bin/python

import sys,os,platform
import publishfunctions

curdir = os.path.dirname(os.path.realpath(__file__))

# ------- Configuration parameters -------------
projectname='r2p'

if platform.system()=='Darwin':
    homefolder='/Users/user'
else:
    homefolder='/home/user'

project=homefolder+'/'+projectname+'/'

# Configure HTML placement at remote server
host='nperraud,r2p@web.sourceforge.net'
www='/home/project-web/r2p/htdocs/'
outputdirweb= '~/work/git/website/rrp-website/'

# -------- Configuration of mat2doc ------------
mat2docpath='~/mat2doc'

try :
    from publish_local import *
except:
    pass



# -------- Automatique configuration ------------
import conf
outputdir=conf.outputdir
outputdirphp=outputdir+'/'+projectname+'-php/'
outputdirmat=outputdir+'/'+projectname+'-mat/'
outputdirtex=outputdir+'/'+projectname+'-tex/'
outputdirpackage=outputdir+'/'+projectname+'-package/'


f=file(project+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

# ------- do not edit below this line ----------
f = open(project + 'mat2doc/startup.m', 'w')
f.write('addpath ' + unlocxpath + '\n')
f.write('init_unlocbox;\n')
f.write('addpath ' + ltfatpath + '\n')
f.write('ltfatstart;\n');
f.write('addpath ' + gsppath + '\n')
f.write('gsp_start;\n');
f.write('addpath ' + rrppath + '\n')
f.write('init_unlocboxrr;\n');
f.close()



todo=sys.argv



#  Optional parameters

if 'fast' in todo:
    plot='--no-execplot '    
else:
    plot='--execplot '


if 'rebuild' in todo:
    build='--rebuild '
    mode = 'rebuild'
elif 'cached' in todo:
    build='--cached '
    mode = 'cached'
else:
    build='--auto '
    mode = 'auto'


if 'phpall' in todo:
    todo.append('mat')
    todo.append('php')
    todo.append('tex')
    todo.append('compiletex')
    todo.append('package')
    todo.append('sendphp')



# #  Publish
# if 'mat' in todo:
#     s=mat2docpath+'/mat2doc.py '+project+' mat'
#     os.system(s)

    
# if 'php' in todo:
#     s=mat2docpath+'/mat2doc.py '+plot+build+' '+project+' php'
#     os.system(s)


# if 'html' in todo:
#     s=mat2docpath+'/mat2doc.py '+plot+build+' '+project+' html'
#     os.system(s)

# if 'tex' in todo:
#     s=mat2docpath+'/mat2doc.py '+plot+build+' '+project+' tex'
#     os.system(s)
#  Publish
for mode in  ['mat', 'php', 'html', 'tex']:
    if mode in todo:
        s = '%s %s/mat2doc.py %s%s %s %s' % ('PYTHONPATH="%s:$PYTHONPATH"' % (curdir,), mat2docpath, plot if mode != 'mat' else '', build if mode != 'mat' else '', project, mode,)
        os.system(s)

# if 'compiletex' in todo: 
        


        # files =  os.listdir('./'+projectname+'-tex')

        # for el in files:
        #     if os.path.isdir('./'+projectname+'-tex/'+el):
        #         print el
        #         #  Copy necessary files
        #         os.system('cp '+project+'/mat2doc/experiments.tex ./'+projectname+'-tmp/')
        #         os.system('cp '+project+'/mat2doc/project.bib ./'+projectname+'-tmp/')
        #         os.system('cp '+project+'/mat2doc/Makefile ./'+projectname+'-tmp/')

        #         #  Copy all the files in the tmp directory
        #         os.system('cp -r ./'+projectname+'-tex/'+el+ '/ ./'+projectname+'-tmp/')

        #         files2 = os.listdir('./'+projectname+'-tex/'+el)
        #         for el2 in files2:
        #             if ((len( el2 )>4) and (el2[-4:]=='.tex')):
        #                 if el2 != 'Contents.tex' and el2 != 'experiments.tex':
        #                     print el2[:-4]
        #                     os.system('mv ./'+projectname+'-tmp/'+el2+' ./'+projectname+'-tmp/temp.tex')
        #                     os.system('cd '+projectname+'-tmp/ && make clean')
        #                     os.system('cd '+projectname+'-tmp/ && make')
        #                     os.system('mv ./'+projectname+'-tmp/experiments.pdf ./'+projectname+'-tex/'+el+'/'+el2[:-4]+'.pdf' )


        #         #  Clear all the files in the tmp directories
        #         os.system('rm ./'+projectname+'-tmp/*')
                
if 'compiletex' in todo: 
        
    files =  os.listdir('./'+projectname+'-tex')

    path1 = 'unlocbox-rr-compile-tex/'
    publishfunctions.safe_mkdir(path1)

    for el in files:
        if os.path.isdir('./'+projectname+'-tex/'+el):
            print 'Cheking ' +el +' directory'

            #  Make the directory
            path2 = path1+el+'/'
            publishfunctions.safe_mkdir(path2)


            files2 = os.listdir('./'+projectname+'-tex/'+el)
            for el2 in files2:
                if ((len( el2 )>4) and (el2[-4:]=='.tex')):
                    if el2 != 'Contents.tex' and el2 != 'experiments.tex':
                        
                        #  Make the directory
                        path3 = path2+el2[:-4]+'/'
                        publishfunctions.safe_mkdir(path3)
                        source = './'+projectname+'-tex/'+el+'/'+el2
                        dest = path3+el2
                        if publishfunctions.do_rebuild_file(source,dest,mode):

                            print '   Rebuild: '+el2[:-4]
                            #  Copy necessary files
                            os.system('cp '+project+'/mat2doc/experiments.tex '+path3)
                            os.system('cp '+project+'/mat2doc/project.bib  '+path3)
                            os.system('cp '+project+'/mat2doc/Makefile  '+path3)

                            #  Copy all the necessary files in the tmp directory
                            os.system('cp -r ./'+projectname+'-tex/'+el+'/'+el2[:-4]+'* '+path3)

                            #  Compile
                            os.system('cp '+path3+el2+' '+path3+'temp.tex')

                            os.system('cd '+path3+' && make clean')
                            os.system('cd '+path3+' && make')

                            os.system('mv '+path3+'experiments.pdf ./'+projectname+'-tex/'+el+'/'+el2[:-4]+'.pdf' )
                        else:
                            print '   Nothing to do for: '+el2[:-4]




#  Packaging
if 'package' in todo:

    files =  os.listdir('./'+projectname+'-mat')

    for el in files:
        if os.path.isdir('./'+projectname+'-mat/'+el):
            fname=outputdirpackage+el
            #  Copy pdf files if they exists
            if os.path.isdir('./'+projectname+'-tex/'+el):
                os.system('ls '+'./'+projectname+'-tex/'+el+'/*.pdf')
                os.system('cp ./'+projectname+'-tex/'+el+'/*.pdf'+' ./'+projectname+'-mat/'+el+'/.')
            # Create the Unix src package
            os.system('tar -zcvf '+fname+'.tgz -C'+projectname+'-mat/ '+el+'/')
            print projectname+'-mat/'+el+'/'
            # Create the Windows src package
            os.system('rm '+fname+'.zip')
            publishfunctions.unix2dos(outputdir+projectname+'-mat/'+el)
            os.system('cd '+projectname+'-mat/ '+'&& zip -r '+fname+'.zip '+el+'/')


 

#  Send to the server


if 'sendphp' in todo:
    os.system('rm '+outputdirphp+'index.php')
    os.system('rm '+outputdirphp+'contentsmenu.php')
    os.system('rm '+outputdirphp+'lookup.php')
    s='rsync -av '+outputdirphp+' '+host+':'+www
    os.system(s)  
    s='rsync -av '+outputdirpackage+' '+host+':'+www+'archive/'
    os.system(s)  


if 'sendweb' in todo:
    s="rsync --verbose --archive --exclude '.git' "+outputdirweb+' '+host+':'+www
    os.system(s)  






