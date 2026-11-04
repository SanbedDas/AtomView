# Builds the PTM extension (atomview._ptmlib); all other metadata lives in pyproject.toml.
import glob, sys
from setuptools import Extension, setup

win = sys.platform == "win32"
args = ["/O2", "/std:c++14", "/EHsc"] if win else ["-O2", "-std=c++11", "-fvisibility=hidden"]
link = ["-pthread"] if sys.platform.startswith("linux") else []
setup(ext_modules=[Extension("atomview._ptmlib", sorted(glob.glob("atomview/_ptm/*.cpp")), include_dirs=["atomview/_ptm"],
                             language="c++", extra_compile_args=args, extra_link_args=link)])
