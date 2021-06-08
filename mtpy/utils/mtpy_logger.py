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
import queue
from concurrent_log_handler import ConcurrentRotatingFileHandler

# =============================================================================
# Global Variables
# =============================================================================
LEVEL_DICT = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL,
}

LOG_FORMAT = logging.Formatter(
    "%(asctime)s [line %(lineno)d] %(name)s.%(funcName)s - %(levelname)s: %(message)s"
)
# Get the configuration file path, should be in same directory as this file
CONF_PATH = Path(__file__).parent
CONF_FILE = Path.joinpath(CONF_PATH, "logging_config.yaml")

# make a folder for the logs to go into.
LOG_PATH = CONF_PATH.parent.parent.joinpath("logs")

if not LOG_PATH.exists():
    LOG_PATH.mkdir()

if not CONF_FILE.exists():
    CONF_FILE = None
    print("No Logging configuration file found, using defaults.")
    
class EvictQueue(queue.Queue):
    def __init__(self, maxsize):
        self.discarded = 0
        super().__init__(maxsize)

    def put(self, item, block=False, timeout=None):
        while True:
            try:
                super().put(item, block=False)
            except queue.Full:
                try:
                    self.get_nowait()
                    self.discarded += 1
                except queue.Empty:
                    pass
                
def speed_up_logs(): 
    rootLogger = logging.getLogger()     
    log_que = EvictQueue(1000)
    queue_handler = logging.handlers.QueueHandler(log_que)
    queue_listener = logging.handlers.QueueListener(log_que, *rootLogger.handlers)
    queue_listener.start()
    rootLogger.handlers = [queue_handler]


def load_configure(config_fn=CONF_FILE):
    # def load_configure(path2configfile='logging.yml'):
    """
    configure/setup the logging according to the input configfile

    :param configfile: .yml, .ini, .conf, .json, .yaml.
    Its default is the logging.yml located in the same dir as this module.
    It can be modofied to use env variables to search for a log config file.
    """
    config_file = Path(config_fn)
    with open(config_file, "r") as fid:
        config_dict = yaml.safe_load(fid)
    logging.config.dictConfig(config_dict)


def get_mtpy_logger(logger_name, fn=None, level="debug"):
    """
    Create a logger, can write to a separate file.  This will write to
    the logs folder in the mt_metadata directory.

    :param logger_name: name of the logger, typically __name__
    :type logger_name: string
    :param fn: file name to write to, defaults to None
    :type fn: TYPE, optional
    :param level: DESCRIPTION, defaults to "debug"
    :type level: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

    logger = logging.getLogger(logger_name)

    # if there is a file name create file in logs directory
    if fn is not None:
        # need to clear the handlers to make sure there is only
        # one call per logger plus stdout
        if logger.hasHandlers():
            logger.handlers.clear()

        logger.propagate = False
        # want to add a stream handler for any Info print statements as stdOut
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(LOG_FORMAT)
        stream_handler.setLevel(LEVEL_DICT["info"])
        logger.addHandler(stream_handler)
        fn = LOG_PATH.joinpath(fn)
        exists = False
        if fn.exists():
            exists = True

        if fn.suffix not in [".log"]:
            fn = Path(fn.parent, f"{fn.stem}.log")

        exists = False
        if fn.exists():
            exists = True

        fn_handler = ConcurrentRotatingFileHandler(fn, maxBytes=2**21)
        # fn_handler = logging.handlers.RotatingFileHandler(fn, maxBytes=2*21)
        fn_handler.setFormatter(LOG_FORMAT)
        fn_handler.setLevel(LEVEL_DICT[level.lower()])
        logger.addHandler(fn_handler)
        if not exists:
            logger.info(f"Logging file can be found {logger.handlers[-1].baseFilename}")
    # commented this out because it leads to a drastic reduction in processing
    # else, give it a null handler, which will go to default logger.
    # else:
    #     null_handler = logging.NullHandler()
    #     null_handler.setFormatter(LOG_FORMAT)
    #     null_handler.setLevel(LEVEL_DICT[level.lower()])
    #     logger.addHandler(null_handler)
    #speed_up_logs()
    return logger
