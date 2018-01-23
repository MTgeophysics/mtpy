from tests.SmartMT import SmartMTGUITestCase


class TestStartGUI(SmartMTGUITestCase):
    """
    only testing the loading of the gui at the moment. to make sure all the necessory packages and dependencies are
    correctly loaded
    """

    def test_main(self):
        self.assertTrue(self.smartMT.isVisible())
        self.assertTrue(self.smartMT.isMaximized())
