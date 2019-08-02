from setuptools import setup

setup(name='StoatyDive',
      author='Florian Heyl',
      description='StoatyDive is a tool to evaluate and classify predicted peak profiles to assess the binding specificity of a protein to its targets.',
      author_email='heylf@informatik.uni-freiburg.de',
      url='https://github.com/heylf/StoatyDive',
      license="GNU General Public License v3.0  (GPL)",
      version='1.0.4',
      packages=['lib'],
      package_data={'': ['*.r', '*.R']},
      include_package_data=True,
      zip_safe = False,
      scripts=['bin/StoatyDive.py'])