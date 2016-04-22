"""
This is a setup.py script manually created from the python bool

Usage:
    python setup.py 
"""
from distutils.core import setup
"""
#root = '/dist/RnaEditor.app/Contents'
#pythonhome = os.path.join(root, 'Frameworks/Python.framework/Versions/Current')

APP = ['RnaEditor.py']
DATA_FILES = ['ui/icons/rnaEditor_icon.pdf','ui/icons/rnaEditor_icon.png','ui/icons/rnaEditor_icon.svg','configuration.txt']
OPTIONS = {'argv_emulation': True,
'iconfile': '/Users/david/git/rnaEditor/ui/icons/rnaEditor.icns',
'plist': {'CFBundleShortVersionString':'0.1.0'}
}
"""
setup(
    
    name='RnaEditor',
    version = "0.1",
    author = "David John",
    author_email = "john@med.uni-frankfurt.de",
    url = "http://rnaeditor.uni-frankfurt.de",
    download_url = "http://rnaeditor.uni-frankfurt.de/install.php",
    py_modules = ["CallEditingSites","createDiagrams","Gene","Genome","gtfHandler","Helper","MapFastq","recountReads","RnaEditor","Transcript","VariantSet"],
    scripts = [],
    package_dir = {"pysam" : "/usr/local/lib/python2.7/dist-packages/pysam/",
                   "PyQt4" : "/usr/lib/python2.7/dist-packages/PyQt4/",
                   "numpy" : "/usr/lib/python2.7/dist-packages/numpy/",
                   "matplotlib" : "/usr/lib/pymodules/python2.7/matplotlib/"},
      packages= ["ui","pysam","PyQt4","numpy","matplotlib"]
)
