import os
import sys
import tempfile
import logging.config

__all__ = [
    'tests',
    'base',
    'constants',
    'empiric',
    'enrich',
    'enrich2',
    'exceptions',
    'fasta',
    'utilities',
    'validators',
    'disable_logging',
    'LOGGER',
]

HOMEDIR = os.path.normpath(os.path.expanduser('~/.mave_convert/'))
tempfile.tempdir = HOMEDIR
if not os.path.isdir(HOMEDIR):
    os.mkdir(HOMEDIR)

LOGGER = 'mavedbconvert'

# Initialize the logging via dictionary configuration
logging.config.dictConfig(
    {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'verbose': {
                'format': '[%(levelname)s] %(asctime)s %(module)s %(message)s'
            },
            'simple': {
                'format': '%(levelname)s %(message)s'
            },
        },
        'handlers': {
            'file': {
                'level': 'WARNING',
                'class': 'logging.FileHandler',
                'filename': os.path.join(HOMEDIR, 'info.log'),
                'formatter': 'verbose'
            },
            'console': {
                'level': 'INFO',
                'class': 'logging.StreamHandler',
                'stream': sys.stdout,
                'formatter': 'verbose'
            },
        },
        'loggers': {
            LOGGER: {
                'handlers': ['file', 'console'],
                'level': 'INFO',
                'propagate': True
            },
        },
    }
)


def disable_logging():
    logging.disable(logging.INFO)
    logging.disable(logging.WARN)
    logging.disable(logging.WARNING)
    logging.disable(logging.ERROR)
    logging.disable(logging.CRITICAL)
    logging.disable(logging.DEBUG)
    logging.disable(logging.FATAL)
