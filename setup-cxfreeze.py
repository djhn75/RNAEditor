#import sys
from cx_Freeze import setup, Executable


# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = [], excludes = [])


#base = 'Win32GUI' if sys.platform=='win32' else None
base=None

executables = [
    Executable('Main.py', base=base, targetName = 'RnaEditor')
]

setup(name='RnaEditor',
      version = '1.0',
      description = 'A tool to detect editing events from rna-seq data.',
      options = dict(build_exe = buildOptions),
      executables = executables)
