import os

__author__ = 'fzhang'

TEST_MTPY_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)  # assume tests is on the root level of mtpy

TEST_DIR = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
TEST_TEMP_DIR = os.path.normpath(os.path.join(TEST_DIR, "temp"))

if not os.path.isdir(TEST_TEMP_DIR):
    os.mkdir(TEST_TEMP_DIR)
