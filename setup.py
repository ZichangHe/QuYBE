from setuptools import setup, find_namespace_packages

setup(
        name='QuYBE',
        version='0.1',
        description='Compressed circuit for simulating 1d Hesinburg dynamics',
        author='Zichang He',
        author_email='zichanghe@ucsb.edu',
        url='https://github.com/ZichangHe/QuYBE.git',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        install_requires = ['numpy', 'qiskit'],
        #license = '',
        #package_data ={'spode': ['core/model.json']},
        python_requires='>=3',
        packages=find_namespace_packages()
)