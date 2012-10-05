try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

setup(name="MTpy",
      scripts=["MTpy/utils/CombineEDIs.py",
               "MTpy/utils/runParalanaMT.py",
               "MTpy/utils/wsmt-pv.py",
               "MTpy/utils/occam2d_gui/run1.py",
               "MTpy/core/RunBIRRPSingleStation.py"]
      )
