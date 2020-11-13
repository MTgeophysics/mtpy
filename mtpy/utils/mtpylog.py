"""
A more Pythonic way of logging:
Define a class MtPyLog to wrap the python logging module;
Use a (optional) configuration file (yaml, ini, json) to configure the logging,
It will return a logger object with the user-provided config setting.
see also: http://www.cdotson.com/2015/11/python-logging-best-practices/
"""

from pathlib import Path
import yaml
import logging
import logging.config


# DEBUG is good for debug and development
# logging.getLogger().setLevel(logging.DEBUG)


class MtPyLog(object):
    # def __init__(self, path2configfile=None):
    @staticmethod
    def load_configure(config_fn=None):
    # def load_configure(path2configfile='logging.yml'):
        """
        configure/setup the logging according to the input configfile

        :param configfile: .yml, .ini, .conf, .json, .yaml.
        Its default is the logging.yml located in the same dir as this module.
        It can be modofied to use env variables to search for a log config file.
        """ 
        if not isinstance(config_fn, Path):
            config_fn = Path(config_fn)
        
        if not config_fn.exists():
            logging.basicConfig()
            
        elif config_fn.suffix in ['.yaml', '.yml']:

            with open(config_fn, "r") as fid:
                config_dict = yaml.safe_load(fid)
            logging.config.dictConfig(config_dict)
            
            # open root logger
            logger = logging.getLogger(__name__)
        
            # make sure everything is working
            logger.info("Started MTpy")
            logger.debug("Beginning debug mode for MTpy")
            debug_fn = logger.root.handlers[1].baseFilename
            error_fn = logger.root.handlers[2].baseFilename
        
            logger.info("Debug Log file can be found at {0}".format(debug_fn))
            logger.info("Error Log file can be found at {0}".format(error_fn))

        elif config_fn.suffix in ['conf', 'ini']:
            logging.config.fileConfig(
                config_fn, disable_existing_loggers=False)
            # must change the default disable_existing_loggers=True to False to
            # make this behave 100% OK
        else:
            raise Exception(
                "logging configuration file %s is not supported" %
                config_fn)

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

