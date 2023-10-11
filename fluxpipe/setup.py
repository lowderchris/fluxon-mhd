from setuptools import setup, find_packages
import versioneer
import sys

with open('requirements.txt') as f:
    required = f.read().splitlines()

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
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    url='https://github.com/lowderchris/fluxon-mhd',

    # Packages to include
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[],
    install_requires=[],
    python_requires=">=3.0",
    classifiers=[
        'Development Status :: 4 - Beta',
    ],
)
