
import glob, setuptools

setuptools.setup(
    name = 'af2genomics',
    version = '24.dev0',
    author = 'jurgjn',
    author_email = 'jurgjn@users.noreply.github.com',
    url = 'https://github.com/jurgjn/af2genomics',
    scripts = glob.glob('workflow/scripts/*'),
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    packages = ['af2genomics'],
    python_requires='>=3.9', # @functools.cache
)
