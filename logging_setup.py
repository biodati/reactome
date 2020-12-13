#!/usr/bin/env python
# -*-coding: utf-8 -*-

# Standard Library
import logging
import logging.config
import sys

# Third Party Imports
import structlog

# Local Imports
import settings

urllib3_logger = logging.getLogger("urllib3")
urllib3_logger.setLevel("CRITICAL")


def add_app_context(logger, method_name, event_dict):
    f, name = structlog._frames._find_first_app_frame_and_name(["logging", __name__])
    event_dict["file"] = f.f_code.co_filename.replace(
        "/Users/william/biodati/reactome_to_bel/", "./"
    )
    event_dict["line"] = f.f_lineno
    event_dict["function"] = f.f_code.co_name

    # event_dict["message"] = event_dict.pop("event", "")

    return event_dict


# timestamper = structlog.processors.TimeStamper(fmt="%Y-%m-%d %H:%M:%S")

pre_chain = [
    # Add the log level and a timestamp to the event_dict if the log entry
    # is not from structlog.
    structlog.stdlib.add_log_level,
    # timestamper,
]

log = structlog.get_logger()

logging.config.dictConfig(
    {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "plain": {
                "()": structlog.stdlib.ProcessorFormatter,
                "processor": structlog.dev.ConsoleRenderer(colors=False),
                "foreign_pre_chain": pre_chain,
            },
            "colored": {
                "()": structlog.stdlib.ProcessorFormatter,
                "processor": structlog.dev.ConsoleRenderer(colors=True),
                "foreign_pre_chain": pre_chain,
            },
        },
        "handlers": {
            "default": {
                "level": settings.LOGLEVEL,
                "class": "logging.StreamHandler",
                "formatter": "colored",
            },
            "file": {
                "level": settings.LOGLEVEL,
                "class": "logging.FileHandler",
                "mode": "a",
                "filename": "reactome_conversion.log",
                "formatter": "plain",
            },
        },
        "loggers": {
            "": {"handlers": ["default", "file"], "level": settings.LOGLEVEL, "propagate": True,},
        },
        # "loggers": {"": {"handlers": ["default"], "level": settings.LOGLEVEL, "propagate": True,},},
    }
)

structlog.configure(
    processors=[
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        # timestamper,
        add_app_context,
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
        structlog.processors.ExceptionPrettyPrinter(),
        structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
    ],
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    wrapper_class=structlog.stdlib.BoundLogger,
    cache_logger_on_first_use=True,
)

# log.info("Here 2", test=True, test2=1)
