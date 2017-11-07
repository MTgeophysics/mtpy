from mtpy.core.edi import Edi

def test_read_write():

    #path2edi = 'data/AMT/15125A_spe.edi'
    path2edi='data/AMT/15125A_imp.edi'

    edi_obj = Edi(edi_fn=path2edi)
    # change the latitude
    edi_obj.lat = 45.7869


    new_edi_fn = "%s_2"%(path2edi)

    ret_edi= edi_obj.write_edi_file(new_edi_fn=new_edi_fn)

    print(ret_edi)

if __name__ =="__main__":
    test_read_write()
