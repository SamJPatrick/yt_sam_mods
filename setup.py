from setuptools import setup

# Change the following fields to the name of the
# directory where your code is.

setup(name="yt_sam_mods",
      version="0.1.0",
      description="My modifications to the yt p2p extention package",
      author="Sam Patrick",
      author_email="sam.patrick@ed.ac.uk",
      license="BSD",
      keywords=["yt", "astronomy"],
      url="https://github.com/SamJPatrick/yt_sam_mods",
      packages=["yt_sam_mods"],
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
