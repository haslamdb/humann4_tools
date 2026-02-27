# humann4_tools/logger.py
"""
Logging functionality for HUMAnN3 Tools.
"""

import os
import sys
import logging
import datetime

def setup_logger(log_file=None, log_level=logging.INFO):
    """
    Setup logger for HUMAnN3 Tools.
    
    Args:
        log_file: Path to log file (optional)
        log_level: Logging level (default: INFO)
        
    Returns:
        Logger instance
    """
    # Create logger
    logger = logging.getLogger('humann4_analysis')
    logger.setLevel(log_level)
    
    # Remove any existing handlers to avoid duplication
    if logger.handlers:
        logger.handlers = []
    
    # Create console handler with the specified log level
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Add formatter to console handler
    ch.setFormatter(formatter)
    
    # Add console handler to logger
    logger.addHandler(ch)
    
    # Add file handler if log file is specified
    if log_file:
        # Create directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            
        # Create file handler with the specified log level
        fh = logging.FileHandler(log_file)
        fh.setLevel(log_level)
        
        # Add formatter to file handler
        fh.setFormatter(formatter)
        
        # Add file handler to logger
        logger.addHandler(fh)
    
    return logger

def log_print(message, level='info'):
    """
    Print message to console and log with the specified level.
    
    Args:
        message: Message to print and log
        level: Logging level (info, debug, warning, error, critical)
    """
    # Print to console
    print(message)
    
    # Get logger
    logger = logging.getLogger('humann4_analysis')
    
    # Log with the specified level
    if level.lower() == 'debug':
        logger.debug(message)
    elif level.lower() == 'warning':
        logger.warning(message)
    elif level.lower() == 'error':
        logger.error(message)
    elif level.lower() == 'critical':
        logger.critical(message)
    else:
        logger.info(message)
