from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='rosewater',
    version='0.1.0',
    description='Expression-aware gene assignment of ROSE output.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/j-andrews7/rosewater",
    author='Jared Andrews',
    author_email='jared.andrews07@gmail.com',
    python_requires='>=3',
    py_modules=['rosewater'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta'
    ],
    install_requires=[
        'Click',
        'pandas',
        'pybedtools'
    ],
    entry_points={
        'console_scripts': [
            'rosewater=rosewater:rosewater'
        ],
    },
)