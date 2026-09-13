import importlib

import mingl


def test_package_has_version():
    assert mingl.__version__ is not None


def test_legacy_mingle_alias_imports():
    assert importlib.import_module("MINGLE").__version__ == mingl.__version__
