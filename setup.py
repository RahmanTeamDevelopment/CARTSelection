from setuptools import setup

setup(
    name = 'CARTSelection',
    version = '1.7.0',
    description = 'CARTSelection tool',
    url = 'https://github.com/RahmanTeamDevelopment/CARTSelection',
    author = 'Marton Munz',
    author_email = 'munzmarci@gmail.com',
    license = 'MIT',
    packages=['cartselection'],
    scripts=['bin/CARTSelection.py', 'bin/cartselection'],
    zip_safe=False
)
