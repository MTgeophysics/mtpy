import unittest
from mtpy.utils import *
import tempfile
import numpy as np

class TestFilehandling(unittest.TestCase):

    def setUp(self):
        #building files in temporary directory

        ts_length = 10
        sampling  = 0.1

        timeaxis = np.arange(ts_length)*sampling

        filelist = []
        for i in range(3):
            FH, F = tempfile.mkstemp()
            print FH,F
            filelist.append(F)
            data = np.zeros((ts_length,2))
            data[:, 0] = timeaxis
            data[:, 1] = np.random.normal(0, 1., ts_length)
            np.savetxt(F,np.array(data))

        self.filelist = filelist


    def test_combine_ts_files(self):
        #check the concatenation of data files - write 3 files to one 
        #and compare first and last rows

        new_FH, new_F = tempfile.mkstemp()

        fulldata = None
        idx =0 
        for i in self.filelist:
            print idx
            idx +=1
            with open(i,'r') as f:
                data = np.loadtxt(f)
                print data 
                if fulldata == None:
                    fulldata = data
                    continue
                np.concatenate((fulldata, data))
                print fulldata
            print

        np.savetxt(new_F, fulldata)

        re_read_data = np.loadtxt(new_F)
        first_line = re_read_data[0,:]
        last_line = re_read_data[-1,:]

        data1 = np.loadtxt(self.filelist[0])[0,:]
        data2 = np.loadtxt(self.filelist[-1])[-1,:]
        
        print last_line, data2
        self.assertTrue(first_line[0] == data1[0])
        self.assertTrue(last_line[1] == data2[1])



    # def test_choice(self):
    #     element = random.choice(self.seq)
    #     self.assertTrue(element in self.seq)

    # def test_sample(self):
    #     with self.assertRaises(ValueError):
    #         random.sample(self.seq, 20)
    #     for element in random.sample(self.seq, 5):
    #         self.assertTrue(element in self.seq)



if __name__ == '__main__':
    unittest.main()
