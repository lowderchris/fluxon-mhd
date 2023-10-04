from setuptools import setup, find_packages
import sys

# Set conditional requirement for pytest-runner
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

setup(
    # Self-descriptive entries which should always be present
    name='fluxpipe',
    author='Gilly',
    author_email='gilly@swri.org',
    description="Forward-model coronal magnetic fields and the solar wind using magnetograms.",
    long_description="Forward-model coronal magnetic fields and the solar wind using magnetograms.",
    long_description_content_type="text/markdown",
    version='105',
    license='BSD-3-Clause',
    url='https://github.com/GillySpace27/sunback',

    # Packages to include
    packages=find_packages(),
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # External dependencies
    install_requires=[],

    platforms=['Windows', 'Linux', 'Mac OS-X'],
    python_requires=">=3.0",

    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
    ]
)
