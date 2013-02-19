#!/usr/bin/env python

"""
This module contains exceptions for MTpy. 



@UofA, 2013
(LK)

"""

#=================================================================

class MTpyError_float( Exception ): 
    pass


class MTpyError_inputarguments( Exception ):
    pass


class MTpyError_ts_data( Exception ): 
    pass

class MTpyError_config_file( Exception ): 
    pass

class MTpyError_file_handling( Exception ): 
    pass

class MTpyError_edi_file( Exception ): 
    pass

class MTpyError_EDI( Exception ): 
    pass
 
class MTpyError_Z( Exception ):
    pass
