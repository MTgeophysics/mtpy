import os

__author__ = 'fzhang'

TEST_DIR = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
TEST_TEMP_DIR = os.path.normpath(TEST_DIR, "temp")

if not os.path.isdir(TEST_TEMP_DIR):
    os.mkdir(TEST_TEMP_DIR)
