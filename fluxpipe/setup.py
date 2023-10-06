""" Forward-model coronal magnetic fields and the solar wind using magnetograms.
FLUXpipe
This is a long description
"""

import sys
from setuptools import setup, find_packages
# from src.movie.dep.sunback import versioneer

# short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = False #{'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

# try:
#     with open("README.md", "r") as handle:
#         long_description = handle.read()
# except:
#     long_description = "\n".join(short_description[2:])

packs = find_packages()

setup(
    # Self-descriptive entries which should always be present
    name='FLUXpipe',
    author='Gilly',
    author_email='gilly@swri.org',
    description="Forward-model coronal magnetic fields and the solar wind using magnetograms.",
    long_description="Forward-model coronal magnetic fields and the solar wind using magnetograms.",
    long_description_content_type="text/markdown",
    version=100,  # versioneer.get_version(),
    # cmdclass=100,  # versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    url='https://github.com/GillySpace27/sunback',  # Website

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=packs,

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Required packages, pulls from pip if needed; do not use for Conda deployment
    install_requires=packs,
    # install_requires=["boto3", "matplotlib", "twine", "pillow", "appscript;platform_system=='Darwin'",
    #                   "moviepy", 'pippi', 'parfive', 'playsound', 'opencv-python', 'numba', "bs4",
    #                   "sunpy", "scipy", "astropy", "html5lib", "numpy"],

    platforms=['Windows', 'Linux', 'Mac OS-X'],            # Valid platforms your code works on, adjust to your flavor
    # 'Linux','Mac OS-X','Unix',

    python_requires=">=3.0",          # Python version restrictions

    classifiers=[
        # How mature is this project? Common values are
        #   1 - Planning
        #   2 - Pre-Alpha
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # # Indicate who your project is intended for
        # 'Intended Audience :: End Users/Desktop',
        # 'Intended Audience :: Science/Research',
        # 'Topic :: Desktop Environment',
        # 'Topic :: Desktop Environment :: Screen Savers',
        # 'Topic :: Multimedia :: Graphics :: Viewers',

        # # Pick your license as you wish (should match "license" above)
        # 'License :: OSI Approved :: BSD License',

        # # Specify the Python versions you support here. In particular, ensure
        # # that you indicate whether you support Python 2, Python 3 or both.
        # # 'Programming Language :: Python :: 3',
        # # 'Programming Language :: Python :: 3.2',
        # # 'Programming Language :: Python :: 3.3',
        # 'Programming Language :: Python :: 3 :: Only',

        # # Platforms
        # 'Operating System :: Microsoft :: Windows',
        # 'Operating System :: POSIX',
        # 'Operating System :: MacOS',
    ]

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
