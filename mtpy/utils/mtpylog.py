"""
A more Pythonic way of logging:
Define a class MtPyLog to wrap the python logging module;
Use a (optional) configuration file (yaml, ini, json) to configure the logging,
then return a logger object with the intended config setting.
see also: http://www.cdotson.com/2015/11/python-logging-best-practices/
"""

import logging
import logging.config
import inspect
import os
#import json
import yaml



class MtPyLog():
    def __init__(self, path2configfile=None):
        """
        configure/setup the logging according to configfile
        :param configfile: .ini, .conf, .json, .yaml, dict
        """
        self.configfile = path2configfile

        if self.configfile is None or self.configfile=='':
            logging.basicConfig()
        elif self.configfile.endswith('yaml') or self.configfile.endswith('yml'):

            this_module_file_path = os.path.abspath(__file__)
            print this_module_file_path

            logging.debug("module file: %s", this_module_file_path)

            yaml_path = os.path.join(os.path.dirname(this_module_file_path), path2configfile)
            print yaml_path

            if os.path.exists(yaml_path):
                with open(yaml_path, 'rt') as f:
                    config = yaml.safe_load(f.read())
                logging.config.dictConfig(config)
            else:
                logging.exception("the config yaml file %s does not exist?", yaml_path)

        elif self.configfile.endswith('.conf') or self.configfile.endswith('.ini'):
            logging.config.fileConfig(self.configfile, disable_existing_loggers=False)
            # must change the default disable_existing_loggers=True to False to make this behave 100% OK
        elif self.configfile.endswith('.json'):
            pass
        else:
            raise Exception("logging configuration file %s is not supported" % self.configfile)

    def get_mtpy_logger(self, loggername=''):
        """
        create a named logger (try different)
        :param loggername: the name (key) of the logger object in this Python interpreter.
        :return:
        """

        logger = logging.getLogger(loggername)  # configured explicitely specifically in logging.conf

        return logger


def test_none_configfile():

    this_fun_name=inspect.getframeinfo(inspect.currentframe())[2]

    print ("starting", this_fun_name)
    # 1 user provides config file to use from envar or other methods
    UsersOwnConfigFile = ''  # ''logging.yaml'
    # 2 construct a MtPyLog object
    myobj = MtPyLog(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    #logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use the root's
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # = root config
    # logger = myobj.get_mtpy_logger()  # root


    print(logger, id(logger), logger.level,logger.handlers)


    # 4 use the named-logger
    logger.debug(this_fun_name+' debug message')
    logger.info(this_fun_name+' info message')
    logger.warn(this_fun_name+ ' warn message')
    logger.error(this_fun_name+' error message')
    logger.critical(this_fun_name+' critical message')

    print ("End of: ", this_fun_name)


def test_yaml_configfile():
    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]
    print ("starting", this_fun_name)

    # 1 user provides config file to use from envar or other methods
    UsersOwnConfigFile='logging.yml'
    # 2 construct a MtPyLog object
    myobj = MtPyLog(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use the default or root's??
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # named logger
    # logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # not good, considered as rootLogger; use the above

    # logger.setLevel(logging.DEBUG)
    print(logger, id(logger), logger.name, logger.level, logger.handlers)

    # create console handler and set level to debug
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.INFO)
    # # add ch to logger
    # logger.addHandler(ch)
    # print(logger, id(logger), logger.name,logger.level, logger.handlers)

    # 4 use the named-logger
    logger.debug(this_fun_name + ' debug message')
    logger.info(this_fun_name + ' info message')
    logger.warn(this_fun_name + ' warn message')
    logger.error(this_fun_name + ' error message')
    logger.critical(this_fun_name + ' critical message')

    print ("End of: ", this_fun_name)

def test_ini_configfile(UsersOwnConfigFile='logging.conf'):
    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]
    print ("starting", this_fun_name)

    # 1 user provides config file to use from envar or other methods

    # 2 construct a MtPyLog object
    myobj = MtPyLog(UsersOwnConfigFile)
    # 3 create a named-logger object
    #logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use the default or root's??
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # named logger
    # logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # not good, considered as rootLogger; use the above

    #logger.setLevel(logging.DEBUG)
    print(logger, id(logger), logger.name, logger.level, logger.handlers)

    # create console handler and set level to debug
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.INFO)
    # # add ch to logger
    # logger.addHandler(ch)
    # print(logger, id(logger), logger.name,logger.level, logger.handlers)

    # 4 use the named-logger
    logger.debug(this_fun_name + ' debug message')
    logger.info(this_fun_name + ' info message')
    logger.warn(this_fun_name + ' warn message')
    logger.error(this_fun_name + ' error message')
    logger.critical(this_fun_name + ' critical message')

    print ("End of: ", this_fun_name)

def test_json_configfile():
    pass

####################################################
# Example application code
# quick test of this class

if __name__ == "__main__":


    # before any configuration of logging, behavior different than the end.
    logging.debug("Start: how about the old logging format?")
    logging.warn("Start: how about the old logging format?")  # warn is default logging level

    test_none_configfile()

    test_ini_configfile()

    test_yaml_configfile()

    # 1 user decide what config file to use from envar or other methods
    UsersOwnConfigFile = '' #''logging.yaml'
    # 2 construct a MtPyLog object
    myobj = MtPyLog(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use the root's
    logger = myobj.get_mtpy_logger(__name__)  #__main__  # = root config
    #logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # root

    # 4 use the named-logger
    print(logger, id(logger), logger.name, logger.level, logger.handlers)
    logger.debug('debug message')
    logger.info('info message')
    logger.warn('warn message')
    logger.error('error message')
    logger.critical('critical message')

    ### now what happen to logging default?  logging has been configured
    logging.debug("End: how about the old logging format?")
    logging.warn("End: how about the old logging format?")
    logging.info("End: how about the old logging format?")
