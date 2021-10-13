
from setuptools import setup, find_packages, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pysmartengine",
    version="0.0.3",
    author="WenpingXie",
    author_email="1303061669@qq.com",
    maintainer="WenpingXie",
    description="a python library for internal combustion engine",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/melan-thompson/pysmartengine",
    package_dir={"":"src"},
    packages=find_packages(
        where='src',
        ),
    # package_data={
    #     "mypkg":["./data/*"]
    # },
    install_requires=[
        'Pillow>=5.1.0',
        'numpy>=1.14.4',
        "scipy>=1.6.2",
        "scikit-learn>=0.24.1",
        "scikit-opt>=0.6.1",
        "matplotlib>=3.3.4",
        ""
        ],
    # py_modules=['Engine'],

    include_package_data=True,
    # ext_modules=[Extension("", [])],
    # data_files=["data"],

    # 添加这个选项，在windows下Python目录的scripts下生成exe文件 # 注意：模块与函数之间是冒号
    # entry_points={
    #     'console_scripts': [
    #         'douyin_image=douyin_image:main'
    #     ],
    # },

    # 程序的所属分类列表
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],

    # 此项需要，否则卸载时报windows error
    zip_safe=False
)