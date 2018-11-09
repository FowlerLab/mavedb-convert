import unittest

from mavedbconvert import disable_logging


if __name__ == "__main__":
    disable_logging()
    loader = unittest.TestLoader()
    tests = loader.discover(start_dir='./', pattern="test_*.py")
    unittest.TextTestRunner().run(tests)
