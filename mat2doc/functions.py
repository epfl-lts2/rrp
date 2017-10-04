#!/usr/bin/python


import sys
import os


def process_content_file(path):
    with open(path) as f:
        content = f.readlines()
        text = '<div class="section" id="figures-of-the-paper">'
        text += '\n<h2>Figures of the paper</h2>'
        text += '\n<ul class="nav nav-list doc-sidenav affix-top simple">'
        for line in content:
            if line[0:4] == '<li>':
                text += line
        text += '</ul>\n'
    with open(path, mode='w') as f:
        f.write(text)


if __name__ == '__main__':
    global_path = sys.argv[1]
    files = []
    print(global_path)
    for (dirpath, dirnames, filenames) in os.walk(global_path):
        for filename in filenames:
            if filename == 'contentsmenu.html':
                files.append(dirpath+'/'+filename)
    print(files)
    for file in files:
        process_content_file(file)