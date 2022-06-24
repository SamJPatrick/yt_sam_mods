from setuptools import setup

# Change the following fields to the name of the
# directory where your code is.

setup(name="yt_thingie",
      version="0.0.0",
      description="A yt thingie.",
      author="It is you",
      author_email="me@me.me",
      license="BSD",
      keywords=["yt", "astronomy"],
      url="https://github.com/brittonsmith/yt_thingie",
      packages=["yt_thingie"],
      include_package_data=True,
      classifiers=[
          "Development Status :: 1 - Planning",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          "License :: OSI Approved :: BSD License",
          "Operating System :: MacOS :: MacOS X",
          "Operating System :: POSIX :: Linux",
          "Operating System :: Unix",
          "Natural Language :: English",
          "Programming Language :: Python :: 3",
      ],
      install_requires=[
          'numpy',
          'yt>=4.0',
      ],
)
