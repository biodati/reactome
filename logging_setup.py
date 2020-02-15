# Standard Library
import logging
import sys

# Third Party Imports
# from pythonjsonlogger import jsonlogger
import structlog
import structlog._frames

# Studio Imports
import settings

# Setup jsonlogger to print JSON
#   see example of integrating Sentry https://codywilbourn.com/2018/08/23/playing-with-python-structured-logs/
json_handler = logging.StreamHandler(sys.stdout)
# json_handler.setFormatter(jsonlogger.JsonFormatter())  # Breaks structlog.processors.JSONRenderer

log_level = "INFO"

# Add sys.stdout Stream handlers to logging
logging.basicConfig(
    format="%(message)s",
    handlers=[json_handler],
    level=log_level,  # Or whatever the general level should be
)

root_logger = logging.getLogger()
root_logger.setLevel(log_level)


def add_log_level(logger, method_name, event_dict):
    """
    Add the log level to the event dict.
    """
    if method_name == "warn":
        # The stdlib has an alias
        method_name = "warning"

    event_dict["log.level"] = method_name
    return event_dict


def add_app_context(logger, method_name, event_dict):
    f, name = structlog._frames._find_first_app_frame_and_name(["logging", __name__])
    event_dict["file"] = f.f_code.co_filename
    event_dict["line"] = f.f_lineno
    event_dict["function"] = f.f_code.co_name
    return event_dict


def elastic_common_schema_keys(logger, method_name, event_dict):
    """Revise event_dict to use Elastic Common Schema"""

    # mydict[new_key] = mydict.pop(old_key)  replace a key
    return event_dict


processors = [
    structlog.stdlib.filter_by_level,
    structlog.stdlib.add_logger_name,
    add_log_level,
    add_app_context,
    structlog.processors.TimeStamper(fmt="iso", key="@timestamp"),
    structlog.processors.StackInfoRenderer(),  # adds stack if stack_info=True
    structlog.processors.format_exc_info,  # Formats exc_info
    structlog.processors.UnicodeDecoder(),  # Decodes all bytes in dict to unicode
    # elastic_common_schema_keys,
]

processors.append(structlog.processors.ExceptionPrettyPrinter())

processors.append(
    structlog.processors.JSONRenderer(
        indent=2, sort_keys=True
    )  # renders the log JSON obj as a string
)

structlog.configure(
    processors=processors,
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    wrapper_class=structlog.stdlib.BoundLogger,
    cache_logger_on_first_use=True,
)
