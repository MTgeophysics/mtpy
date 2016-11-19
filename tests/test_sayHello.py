from unittest import TestCase
from hello4tests import say_hello

__author__ = 'fzhang'


class TestSayHello(TestCase):

    def setUp(self):
        print ("Calling setUp")
        self.obj= say_hello.SayHello()

    def tearDown(self):
        print ("Calling tearDown")
        del self.obj

    def test_say(self, himsg="Friend"):
        msg=self.obj.say(himsg)
        self.assertEquals(msg, himsg)
        #self.assertEquals(msg, "fail_at_wrong_msg")
