import os
import pandas as pd
import logging

def generate_summary_and_log(output_folder):
    """
    Generate a summary DataFrame and log messages to a file and terminal.

    Parameters:
    - output_folder (str): The directory where output files will be saved.

    Returns:
    - df (pandas.DataFrame): The DataFrame containing the summary data.
    """

    # Configure logging
    log_file = os.path.join(output_folder, 'summary.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Create logger
    logger = logging.getLogger()

    # Add a StreamHandler to display log messages in the terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)  # Set the desired verbosity level
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # Log start of function execution
    logger.info('Starting generate_summary_and_log function execution...')

    # Generate summary DataFrame
    data = {'A': [1, 2, 3], 'B': [4, 5, 6]}
    df = pd.DataFrame(data)

    # Log DataFrame summary
    logger.info('Generated DataFrame summary:\n%s', df.to_string())

    # Save DataFrame to a CSV file
    csv_file = os.path.join(output_folder, 'summary.csv')
    df.to_csv(csv_file, index=False)

    # Log end of function execution
    logger.info('generate_summary_and_log function execution completed.')

    return df

# Example usage
summary_df = generate_summary_and_log('./')
print(summary_df)
