"""
see also: http://www.cdotson.com/2015/11/python-logging-best-practices/
"""

import logging
import logging.config

logging.config.fileConfig('logging.conf')

# create a named logger (try different)
logger = logging.getLogger('simpleExample')  #configured explicitely specifically in logging.conf
logger = logging.getLogger('simpleExample2') # not configured, use the root's
#logger = logging.getLogger(__name__)  #__main__  # = root
#logger = logging.getLogger()   # root
#logger = logging.getLogger('')  # root

# 'application' code
logger.debug('debug message')
logger.info('info message')
logger.warn('warn message')
logger.error('error message')
logger.critical('critical message')