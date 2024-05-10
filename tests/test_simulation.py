# TODO: batch_simulation()
import numpy as np
import pandas as pd
import pytest

import src.uspekpy.simulation as sim


class TestReadFileToDataFrame:
    @staticmethod
    def test_read_excel_file(tmp_path):
        # Create a temporary Excel file for testing
        excel_data = {
            'A': [1, 2, 3],
            'B': [4, 5, 6]
        }
        excel_df = pd.DataFrame(excel_data)
        excel_file_path = tmp_path / 'test_excel_file.xlsx'
        excel_df.to_excel(excel_file_path, index=False)

        # Call the function with the Excel file
        result_df = sim.read_file_to_dataframe(str(excel_file_path), sheet_name='Sheet1')

        # Assert that the result is a DataFrame and contains the same data as the original DataFrame
        assert isinstance(result_df, pd.DataFrame)
        assert result_df.equals(excel_df)

    @staticmethod
    def test_read_csv_file(tmp_path):
        # Create a temporary CSV file for testing
        csv_data = {
            'A': [1.0, 2.0, 3.0],
            'B': [4.0, 5.0, 6.0]
        }
        csv_df = pd.DataFrame(csv_data)
        csv_file_path = tmp_path / 'test_csv_file.csv'
        csv_df.to_csv(csv_file_path, index=False)

        # Call the function with the CSV file
        result_df = sim.read_file_to_dataframe(str(csv_file_path))

        # Assert that the result is a DataFrame and contains the same data as the original DataFrame
        assert isinstance(result_df, pd.DataFrame)
        assert result_df.equals(csv_df)

    @staticmethod
    def test_invalid_file_format():
        # Call the function with an invalid file format and expect it to raise a ValueError
        with pytest.raises(ValueError):
            sim.read_file_to_dataframe('invalid_file.txt')


class TestParseBeamParameters:
    @staticmethod
    def test_valid_dataframe():
        # Create a valid DataFrame for testing
        data = {
            'Al filter width (mm)': [1.0],
            'Al filter width rel uncertainty': [0.1],
            'Cu filter width (mm)': [2.0],
            'Cu filter width rel uncertainty': [0.2],
            'Sn filter width (mm)': [1.0],
            'Sn filter width rel uncertainty': [0.1],
            'Pb filter width (mm)': [2.0],
            'Pb filter width rel uncertainty': [0.2],
            'Be filter width (mm)': [1.0],
            'Be filter width rel uncertainty': [0.1],
            'Air gap (mm)': [2.0],
            'Air gap rel uncertainty': [0.2],
            'Peak kilovoltage (kV)': [100],
            'Peak kilovoltage rel uncertainty': [10],
            'Anode angle (deg)': [45],
            'Anode angle rel uncertainty': [5]
        }
        df = pd.DataFrame(data).T
        df.columns = ['label']

        # Call the function with the DataFrame and assert the result
        result = sim.parse_beam_parameters(df, column='label')
        expected_result = {
            'Al': (1.0, 0.1),
            'Cu': (2.0, 0.2),
            'Sn': (1.0, 0.1),
            'Pb': (2.0, 0.2),
            'Be': (1.0, 0.1),
            'Air': (2.0, 0.2),
            'kVp': (100, 10),
            'th': (45, 5)
        }
        assert result == expected_result

    @staticmethod
    def test_invalid_dataframe():
        # Creates an invalid DataFrame with missing columns for testing
        data = {
            'Al filter width (mm)': [1.0],
            'Al filter width rel uncertainty': [0.1]
        }
        df = pd.DataFrame(data).T
        df.columns = ['label']

        # Calls the function with the DataFrame and expect it to raise a KeyError
        with pytest.raises(KeyError):
            sim.parse_beam_parameters(df, column='label')


class TestOutputDigest:
    @staticmethod
    def test_output_digest():
        # Creates input DataFrame containing simulation parameters
        input_data = {
            'Name': ['Label1', 'Label2', 'Label3'],
            'Case1': [1, 2, 3],
            'Case2': [4, 5, 6]
        }
        input_df = pd.DataFrame(input_data)
        input_df.set_index('Name', inplace=True)

        # Creates output DataFrames containing simulation results
        output_data1 = {
            '#': ['Mean', 'Standard deviation', 'Relative uncertainty'],
            'HVL1 Al (mm)': [10.0, 1.0, 0.1],
            'HVL2 Al (mm)': [20.0, 2.0, 0.2],
            'HVL1 Cu (mm)': [30.7, 3.0, 0.3],
            'HVL2 Cu (mm)': [40.0, 4.0, 0.4],
            'Mean energy (keV)': [50.0, 5.0, 0.5],
            'Mean air kerma (keV/g)': [60.0, 6.0, 0.6],
            'Mean conv. coefficient (Sv/Gy).': [70.0, 7.0, 0.7]
        }
        output_df1 = pd.DataFrame(output_data1)

        output_data2 = {
            '#': ['Mean', 'Standard deviation', 'Relative uncertainty'],
            'HVL1 Al (mm)': [70.0, 7.0, 0.7],
            'HVL2 Al (mm)': [70.0, 7.0, 0.7],
            'HVL1 Cu (mm)': [10.0, 1.0, 0.1],
            'HVL2 Cu (mm)': [20.0, 2.0, 0.2],
            'Mean energy (keV)': [30.7, 3.0, 0.3],
            'Mean air kerma (keV/g)': [40.0, 4.0, 0.4],
            'Mean conv. coefficient (Sv/Gy)': [50.0, 5.0, 0.5]
        }
        output_df2 = pd.DataFrame(output_data2)

        # Call the function with the input DataFrame and list of output DataFrames
        result_df = sim.output_digest(input_df, [output_df1, output_df2])

        # Define the expected result DataFrame based on the provided function logic
        expected_result_data = {
            'Label1': [1.0, 4.0],
            'Label2': [2.0, 5.0],
            'Label3': [3.0, 6.0],
            'Results': [np.NaN, np.NaN],
            'HVL1 Al Mean (mm)': [10.0, 70.0],
            'HVL1 Al Standard deviation (mm)': [1.0, 7.0],
            'HVL1 Al Relative uncertainty': [0.1, 0.7],
            'HVL2 Al Mean (mm)': [20.0, 70.0],
            'HVL2 Al Standard deviation (mm)': [2.0, 7.0],
            'HVL2 Al Relative uncertainty': [0.2, 0.7],
            'HVL1 Cu Mean (mm)': [30.7, 10.0],
            'HVL1 Cu Standard deviation (mm)': [3.0, 1.0],
            'HVL1 Cu Relative uncertainty': [0.3, 0.1],
            'HVL2 Cu Mean (mm)': [40.0, 20.0],
            'HVL2 Cu Standard deviation (mm)': [4.0, 2.0],
            'HVL2 Cu Relative uncertainty': [0.4, 0.2],
            'Mean energy Mean (keV)': [50.0, 30.7],
            'Mean energy Standard deviation (keV)': [5.0, 3.0],
            'Mean energy Relative uncertainty': [0.5, 0.3],
            'Mean air kerma (keV/g)': [60.0, 40.0],
            'Mean air kerma Standard deviation (keV/g)': [6.0, 4.0],
            'Mean air kerma Relative uncertainty': [0.6, 0.4],
            'Mean conv. coefficient (Sv/Gy)': [70.0, 50.0],
            'Mean conv. coefficient Standard deviation (Sv/Gy)': [7.0, 5.0],
            'Mean conv. coefficient. Relative uncertainty': [0.7, 0.5],
        }
        expected_result_df = pd.DataFrame(expected_result_data).T
        expected_result_df.columns = ['Case1', 'Case2']

        # Assert that the result is a DataFrame and contains the expected data
        assert isinstance(result_df, pd.DataFrame)
        assert result_df.equals(expected_result_df)
