from setuptools import setup

setup(name='pfam',
      version='0.1',
      description='Programmaing interface for Pfam database.',
      url='https://github.com/AlbertoBoldrini/python-pfam',
      author='Alberto Boldrini',
      author_email='alberto.boldrini@studenti.unitn.it',
      license='MIT',
      packages=['pfam'],
      install_requires=[
          'requests',
      ],
      zip_safe=False)