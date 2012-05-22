from distutils.core import setup

setup(
	name = 'genompy',
	description='Utilities for DNA copy-number analysis',
	version = '0.1.0',
	author = 'David J H Shih',
	author_email = 'djh.shih@gmail.com',
	packages = ['genompy', 'genompy/tests', 'genompy/plot'],
	scripts = ['bin/gistic2mat.py', 'bin/cnplot.py'],
	url = 'TBD',
	license = 'LICENSE.txt',
	long_description = open('README.txt').read(),
	install_requires = [],
)

