"""
This is a setup.py script manually created from the python bool

Usage:
    python setup.py 
"""
from distutils.core import setup
"""
#root = '/dist/RnaEditor.app/Contents'
#pythonhome = os.path.join(root, 'Frameworks/Python.framework/Versions/Current')

APP = ['Main.py']
DATA_FILES = ['ui/icons/rnaEditor_icon.pdf','ui/icons/rnaEditor_icon.png','ui/icons/rnaEditor_icon.svg','configuration.txt']
OPTIONS = {'argv_emulation': True,
'iconfile': '/Users/david/git/rnaEditor/ui/icons/rnaEditor.icns',
'plist': {'CFBundleShortVersionString':'0.1.0'}
}
"""
setup(
    
    name='RnaEditor',
    version = "1.0",
    author = "David John",
    author_email = "john@med.uni-frankfurt.de",
    py_modules = ["CallEditingSites","createDiagrams","Gene","Genome","gtfHandler","Helper","Main","MapFastq","recountReads","RnaEdit","Transcript","VariantSet",   "ui",]
)
