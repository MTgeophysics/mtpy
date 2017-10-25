from unittest import TestCase

from qtpy.QtTest import QTest

from mtpy.gui.SmartMT.start import StartGUI


class TestStartGUI(TestCase):
    """
    only testing the loading of the gui at the moment. to make sure all the necessory packages and dependencies are
    correctly loaded
    """
    def setUp(self):
        self.smartMT = StartGUI()
        self.smartMT.show()
        QTest.qWaitForWindowActive(self.smartMT)

    def test_main(self):
        self.assertTrue(self.smartMT.isVisible())
        self.assertTrue(self.smartMT.isFullScreen())

    def tearDown(self):
        self.smartMT.close()

