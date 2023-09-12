import logging

class CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    blue = "\x1b[34;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    magenta = "\x1b[35;20m"
    green = "\x1b[32;20m"
    reset = "\x1b[0m"
    format = "%(name)s:  %(levelname)s ----- %(message)s)"

    FORMATS = {
        logging.DEBUG: magenta + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


logger = logging.getLogger("API")
logger.setLevel(logging.DEBUG)


logger.handlers = []
logging.getLogger().handlers = []

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.propagate = False
logger.addHandler(ch)