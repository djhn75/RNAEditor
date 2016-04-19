#import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = {
"includes": [],
"packages": ['os','sys','re','copy','gc','pysam','__builtin__','array','collection','gzip','itertools','operator','datetime','subprocess','errno','collections','traceback','PyQt4','matplotlib','numpy','multiprocessing','pyface'],
'excludes' : ['boto.compat.sys',
                 'boto.compat._sre',
                 'boto.compat._json',
                 'boto.compat._locale',
                 'boto.compat._struct',
                 'boto.compat.array',
                 'collections.abc',
              'collections.sys'],
"include_files": []}


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
