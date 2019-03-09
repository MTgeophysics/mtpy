"""
A more Pythonic way of logging:
Define a class MtPyLog to wrap the python logging module;
Use a (optional) configuration file (yaml, ini, json) to configure the logging,
It will return a logger object with the user-provided config setting.
see also: http://www.cdotson.com/2015/11/python-logging-best-practices/
"""

import os
# import json
import yaml
import logging
import logging.config
import inspect


# DEBUG is good for debug and development
# logging.getLogger().setLevel(logging.DEBUG)


class MtPyLog(object):
    # def __init__(self, path2configfile=None):
    @staticmethod
    def load_configure(path2configfile=None):
    # def load_configure(path2configfile='logging.yml'):
        """
        configure/setup the logging according to the input configfile

        :param configfile: .yml, .ini, .conf, .json, .yaml.
        Its default is the logging.yml located in the same dir as this module.
        It can be modofied to use env variables to search for a log config file.
        """
        configfile = path2configfile

        if configfile is None or configfile == '':
            logging.basicConfig()
        elif configfile.endswith('yaml') or configfile.endswith('yml'):

            this_module_file_path = os.path.abspath(__file__)

            logging.info("module file: %s", this_module_file_path)

            yaml_path = os.path.join(
                os.path.dirname(this_module_file_path),
                path2configfile)
            logging.info('Effective yaml configuration file %s', yaml_path)

            if os.path.exists(yaml_path):
                with open(yaml_path, 'rt') as f:
                    config = yaml.safe_load(f.read())
                logging.config.dictConfig(config)
            else:
                logging.exception(
                    "the config yaml file %s does not exist?", yaml_path)

        elif configfile.endswith('.conf') or configfile.endswith('.ini'):
            logging.config.fileConfig(
                configfile, disable_existing_loggers=False)
            # must change the default disable_existing_loggers=True to False to
            # make this behave 100% OK
        elif configfile.endswith('.json'):
            pass
        else:
            raise Exception(
                "logging configuration file %s is not supported" %
                configfile)

    @staticmethod
    def get_mtpy_logger(loggername=''):
        """
        create a named logger (try different)
        :param loggername: the name (key) of the logger object in this Python interpreter.
        :return:
        """

        # configured explicitely specifically in logging.conf
        logger = logging.getLogger(loggername)

        return logger


def test_none_configfile():
    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]

    print(("starting", this_fun_name))
    # 1 user provides config file to use from envar or other methods
    UsersOwnConfigFile = ''  # ''logging.yaml'
    # 2 construct a MtPyLog object
    myobj = MtPyLog(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use
    # the root's
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # = root config
    # logger = myobj.get_mtpy_logger()  # root

    print((logger, id(logger), logger.level, logger.handlers))

    # 4 use the named-logger
    logger.debug(this_fun_name + ' debug message')
    logger.info(this_fun_name + ' info message')
    logger.warn(this_fun_name + ' warn message')
    logger.error(this_fun_name + ' error message')
    logger.critical(this_fun_name + ' critical message')

    print(("End of: ", this_fun_name))


def test_yaml_configfile(yamlfile='logging.yml'):
    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]
    print(("starting", this_fun_name))

    # 1 user provides config file to use from envar or other methods
    UsersOwnConfigFile = yamlfile
    # 2 construct a MtPyLog object
    myobj = MtPyLog.load_configure(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use
    # the default or root's??
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # named logger
    # logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # not good, considered as
    # rootLogger; use the above

    # logger.setLevel(logging.DEBUG)
    print((logger, id(logger), logger.name, logger.level, logger.handlers))

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

    print(("End of: ", this_fun_name))


def test_ini_configfile(UsersOwnConfigFile='logging.conf'):
    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]
    print(("starting", this_fun_name))

    # 1 user provides config file to use from envar or other methods

    # 2 construct a MtPyLog object
    myobj = MtPyLog.load_configure(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use
    # the default or root's??
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # named logger
    # logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # not good, considered as
    # rootLogger; use the above

    # logger.setLevel(logging.DEBUG)
    print((logger, id(logger), logger.name, logger.level, logger.handlers))

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

    print(("End of: ", this_fun_name))


def test_json_configfile():
    pass


####################################################
# Example application code
# quick test of this class

if __name__ == "__main__":
    # before any configuration of logging, behavior different than the end.
    logging.debug("Start: how about the old logging format?")
    # warn is default logging level
    logging.warn("Start: how about the old logging format?")

    # test_none_configfile()
    #
    # test_ini_configfile()

    test_yaml_configfile(yamlfile='logging.yml')

    # 1 user decide what config file to use from envar or other methods
    UsersOwnConfigFile = ''  # ''logging.yaml'
    # 2 construct a MtPyLog object
    myobj = MtPyLog.load_configure(UsersOwnConfigFile)
    # 3 create a named-logger object
    # logger = myobj.get_mtpy_logger('simpleExample')
    # logger = myobj.get_mtpy_logger('simpleExample2') # not configured, use
    # the root's
    logger = myobj.get_mtpy_logger(__name__)  # __main__  # = root config
    # logger = myobj.get_mtpy_logger()  # root
    # logger = myobj.get_mtpy_logger('')  # root

    # 4 use the named-logger
    print((logger, id(logger), logger.name, logger.level, logger.handlers))
    logger.debug('debug message')
    logger.info('info message')
    logger.warn('warn message')
    logger.error('error message')
    logger.critical('critical message')

    # now what happen to logging default?  logging has been configured
    logging.debug("End: how about the old logging format?")
    logging.warn("End: how about the old logging format?")
    logging.info("End: how about the old logging format?")
