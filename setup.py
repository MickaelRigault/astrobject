#! /usr/bin/env python
#

DESCRIPTION = "astrobject: Basic tools for astrophysics and cosmology"
LONG_DESCRIPTION = """\
This module gathers the basic tools one usually needs for astrophysics and cosmology usage.
"""

DISTNAME = 'astrobject'
AUTHOR = 'Astrobject Developers'
MAINTAINER = 'Mickael Rigault' 
MAINTAINER_EMAIL = 'mrigault@physik.hu-berlin.de'
URL = 'https://github.com/MickaelRigault/astrobject/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/MickaelRigault/astrobject/tarball/0.6'
VERSION = '0.6.3'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import propobject
    except ImportError:
        install_requires.append('propobject')
    try:
        import astropy
    except ImportError:
        install_requires.append('astropy')

    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = [
            'astrobject', 
            'astrobject.utils', 
            'astrobject.collections', 
            'astrobject.instruments', 
            'astrobject.utils.plot']

    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          install_requires=install_requires,
          packages=packages,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3.5',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
      )
