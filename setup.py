from setuptools import setup, find_packages

setup(
    name='ligandparam',
    version='0.2.0',
    author='Zeke Piskulich',
    author_email='piskulichz@gmail.com',
    description='A ligand parmameterization package for Amber',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/my-python-package',
    packages=find_packages(),
    install_requires=open('requirements.txt').read().splitlines(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)