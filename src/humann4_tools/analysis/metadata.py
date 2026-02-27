# humann4_tools/analysis/metadata.py
import pandas as pd
import traceback
import logging

def read_and_process_metadata(sample_key, logger):
    """
    Read and process sample metadata file.
    
    Args:
        sample_key: Path to the sample key CSV file
        logger: Logger instance for logging
        
    Returns:
        DataFrame containing the sample metadata
    """
    try:
        df = pd.read_csv(sample_key)
        logger.info(f"Loaded sample key {len(df)} rows, columns: {list(df.columns)}")
        return df
    except Exception as e:
        logger.error(f"Error reading sample key: {str(e)}")
        logger.error(traceback.format_exc())
        raise RuntimeError(e)