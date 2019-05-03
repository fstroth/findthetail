import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="findthetail",
    version="1.1.0",
    author="Frederik Strothmann",
    author_email="frstrothmann@gmail.com",
    description="Package implementing the tail detection method detaild in the paper https://arxiv.org/abs/1805.10040.",
    license="MIT",
    keywords="extreme value statistics tail detection",
    url="https://github.com/fstroth/findthetail",
    packages=['findthetail'],
    package_data={'findthetail': ['templates/*']},
    long_description=read('README.md'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
)
