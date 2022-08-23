"""
py2app/py2exe build script for RNAEditor.

Usage (Mac OS X):
    python setup.py py2app

Usage (Windows):
    python setup.py py2exe

Usage (Linux):
    python sytup.py install
"""
from __future__ import print_function
from setuptools import setup, find_packages
import io, sys, codecs, os

mainscript = 'RNAEditor.py'
name="RNAEditor"
OPTIONS = {'iconfile': 'ui/icons/rnaEditor.icns','plist': {'CFBundleShortVersionString':'0.1.0'}
}
DATA_FILES = ['ui/icons/rnaEditor_icon.pdf','ui/icons/RNAeditor.png','ui/icons/rnaEditor_icon.svg','configuration.txt','ui/icons/inputTab_icon.png']

if sys.platform == 'darwin':
    extra_options = dict(
        setup_requires=['py2app'],
        app=[mainscript],
        # Cross-platform applications generally expect sys.argv to
        # be used for opening files.
        options=dict(py2app=OPTIONS),
    )
elif sys.platform == 'win32':
    extra_options = dict(
        setup_requires=['py2exe'],
        app=[mainscript],
    )
    sys.exit("RNAEditor does not support Windows yet")
else:
    extra_options = dict(
        # Normally unix-like platforms will use "setup.py install"
        # and install the main script as such
        scripts=[mainscript],
    )
    sys.exit("Install all the requirements and run RNAEditor.py")

setup(

    name=name,
    version="0.1",
    author="David John",
    packages=find_packages(),
    description="RNAEditor is tool to analyze RNA editing events from RNA-seq data.",
    author_email="john@med.uni-frankfurt.de",
    url="http://rnaeditor.uni-frankfurt.de",
    download_url="http://rnaeditor.uni-frankfurt.de/install.php",
    data_files=DATA_FILES,
    install_requires=['numpy>=1.9.0', 'pysam>=0.9.0', 'matplotlib>=1.4.3'],
    **extra_options
)

