import os

import pandas as pd


# Function to search a file and delete it if already existing
def buscar_eliminar(ruta):
    """
    Attempt to delete a file specified by the given path.

    Parameters
    ----------
    ruta : str
        The path to the file to be deleted.

    Notes
    -----
    If the file specified by `ruta` exists, it will be removed from the system.
    If the file does not exist, this function does nothing.

    Examples
    --------
    >>> buscar_eliminar('path/to/file.txt')
    # File 'file.txt' is deleted if it exists, otherwise, no action is taken.
    """
    # Step 1: Check if the file specified by `ruta` exists
    if os.path.exists(ruta):
        # Step 2: If the file exists, remove it from the system
        os.remove(ruta)
    else:
        # Step 3: If the file does not exist, do nothing
        # This could be enhanced by providing feedback to the user or logging a message
        print("File does not exist.")


# Function to go from dataframe to txt and save the file
# PAL: mean, standard deviation and var coefficient % which is assumed as relative percentage uncertainty
def guardar_txt(ruta, nombre, calidad, media, angulo, desviacion, v_hpk, energiaMedia, desviacionEnergia,
                coef_variacionEnergia, kermaMedio, desviacionKerma, coef_variacionKerma):
    """
    Append or create a CSV file with the provided data.

    Parameters
    ----------
    ruta : str
        The path to the CSV file.
    nombre : str
        Name of the data.
    calidad : float
        Quality of the data.
    media : float
        Mean of the data.
    angulo : str
        Angle of the data.
    desviacion : float
        Standard deviation of the data.
    v_hpk : float
        Percentile of the data.
    energiaMedia : float
        Mean energy of the data.
    desviacionEnergia : float
        Standard deviation of the mean energy.
    coef_variacionEnergia : float
        Coefficient of variation of the mean energy.
    kermaMedio : float
        Mean kerma of the data.
    desviacionKerma : float
        Standard deviation of the mean kerma.
    coef_variacionKerma : float
        Coefficient of variation of the mean kerma.

    Notes
    -----
    This function appends or creates a CSV file at the specified `ruta` with the provided data. If the file
    does not exist or is empty, it writes the data with headers. If the file already exists and is not empty,
    it appends the data to the existing file without headers.

    Examples
    --------
    >>> guardar_txt('data.csv', 'Sample', 0.9, 25.5, '45°', 3.2, 90, 150.5, 2.0, 1.5, 100.3, 5.0, 2.3)
    # Appends or creates 'data.csv' file with the provided data.
    """
    # Step 1: Create a DataFrame with the provided data and column names
    x = pd.DataFrame({"_Name": [nombre], "_quality": [calidad], "_mean_hK_" + angulo: [media],
                      "______std dev hK______": [desviacion], "_______ur%(hK)_______": [v_hpk],
                      "_______Mean energy_______": [energiaMedia], "_______std dev Emean_______": [desviacionEnergia],
                      "_______ur%(Emean)______": [coef_variacionEnergia],
                      "_______Mean Kair_______": [kermaMedio], "_______std dev Kair_______": [desviacionKerma],
                      "_______ur%(Kair)_______": [coef_variacionKerma]})

    # Step 2: Write the DataFrame to the CSV file
    # If the file does not exist or is empty, write with headers (mode="w")
    # If the file already exists and is not empty, append without headers (mode="a")
    x.to_csv(ruta, header=not (os.path.isfile(ruta) and os.stat(ruta).st_size != 0), index=False, mode="a", sep=",")


def calcular_energia_media(fluencias, energias):
    """
    Calculate the average energy given a list of fluences and their corresponding energies.

    Parameters
    ----------
    fluencias : list of float
        List of fluences.
    energias : list of float
        List of energies corresponding to the fluences.

    Returns
    -------
    float
        The calculated average energy.

    Notes
    -----
    The function calculates the average energy by dividing the sum of the product of each fluence and its corresponding
    energy by the total sum of fluences. This calculation gives a weighted average, where each energy value is weighted
    by its corresponding fluence.

    If the total sum of fluences is zero, indicating no data, the function returns 0 to avoid division by zero errors.

    Examples
    --------
    >>> fluences = [100, 200, 300]
    >>> energies = [0.5, 1.0, 1.5]
    >>> calcular_energia_media(fluences, energies)
    # Returns the average energy calculated from the provided fluences and energies.
    """
    # Step 1: Calculate the sum of the product of each fluence and its corresponding energy
    sumatorio_producto = sum(fluencia * energia for fluencia, energia in zip(fluencias, energias))

    # Step 2: Calculate the total sum of fluences
    sumatorio_fluencias = sum(fluencias)

    # Step 3: Check if the total sum of fluences is zero to avoid division by zero errors
    if sumatorio_fluencias == 0:
        return 0

    # Step 4: Calculate the average energy by dividing the sum of the product by the total sum of fluences
    energia_media = sumatorio_producto / sumatorio_fluencias

    return energia_media


def calcular_kerma_medio(fluencias, energias, coeficientes_mutr_rho):
    """
    Calculate the mean kerma given lists of fluences, energies, and coefficients.

    Parameters
    ----------
    fluencias : list of float
        List of fluences.
    energias : list of float
        List of energies corresponding to the fluences.
    coeficientes_mutr_rho : list of float
        List of coefficients corresponding to the fluences and energies.

    Returns
    -------
    float
        The calculated mean kerma.

    Notes
    -----
    The function calculates the mean kerma by summing the product of each fluence, energy, and coefficient,
    and then dividing this sum by the total sum of fluences.

    If the total sum of fluences is zero, indicating no data, the function returns 0 to avoid division by zero errors.

    Examples
    --------
    >>> fluences = [100, 200, 300]
    >>> energies = [0.5, 1.0, 1.5]
    >>> coefficients = [0.1, 0.2, 0.3]
    >>> calcular_kerma_medio(fluences, energies, coefficients)
    # Returns the mean kerma calculated from the provided fluences, energies, and coefficients.
    """
    # Step 1: Calculate the sum of the product of fluences, energies, and coefficients
    sumatorio_producto_kerma = sum(fluencia * energia * coeficiente_mutr_rho
                                   for fluencia, energia, coeficiente_mutr_rho
                                   in zip(fluencias, energias, coeficientes_mutr_rho))

    # Step 2: Calculate the total sum of fluences
    sumatorio_fluencias = sum(fluencias)

    # Step 3: Check if the total sum of fluences is zero to avoid division by zero errors
    if sumatorio_fluencias == 0:
        return 0

    # Step 4: Calculate the mean kerma by dividing the sum of the product by the total sum of fluences
    kerma_medio = sumatorio_producto_kerma / sumatorio_fluencias

    return kerma_medio


def calcular_hpk(fluencias, energias, coeficientes_mutr_rho, hk):
    """
    Calculate the high-kinetic energy per unit mass (hpk) given lists of fluences, energies, coefficients, and hk values.

    Parameters
    ----------
    fluencias : list of float
        List of fluences.
    energias : list of float
        List of energies corresponding to the fluences.
    coeficientes_mutr_rho : list of float
        List of coefficients corresponding to the fluences and energies.
    hk : list of float
        List of hk values.

    Returns
    -------
    float
        The calculated high-kinetic energy per unit mass (hpk).

    Notes
    -----
    The function calculates the hpk by summing the product of each fluence, energy, coefficient, and hk value,
    and then dividing this sum by the total sum of the product of fluences, energies, and coefficients.

    If the total sum of the product of fluences, energies, and coefficients is zero, indicating no data,
    the function returns 0 to avoid division by zero errors.

    Examples
    --------
    >>> fluences = [100, 200, 300]
    >>> energies = [0.5, 1.0, 1.5]
    >>> coefficients = [0.1, 0.2, 0.3]
    >>> hk_values = [0.01, 0.02, 0.03]
    >>> calcular_hpk(fluences, energies, coefficients, hk_values)
    # Returns the high-kinetic energy per unit mass (hpk) calculated from the provided fluences, energies,
    # coefficients, and hk values.
    """
    # Step 1: Calculate the sum of the product of fluences, energies, coefficients, and hk values
    sumatorio_numerador = sum(fluencia * energia * coeficiente_mutr_rho * hk
                              for fluencia, energia, coeficiente_mutr_rho, hk
                              in zip(fluencias, energias, coeficientes_mutr_rho, hk))

    # Step 2: Calculate the total sum of the product of fluences, energies, and coefficients
    sumatorio_denominador = sum(fluencia * energia * coeficiente_mutr_rho
                                for fluencia, energia, coeficiente_mutr_rho
                                in zip(fluencias, energias, coeficientes_mutr_rho))

    # Step 3: Check if the total sum of the product of fluences, energies, and coefficients is zero
    # to avoid division by zero errors
    if sumatorio_denominador == 0:
        return 0

    # Step 4: Calculate the hpk by dividing the sum of the numerator by the sum of the denominator
    hpk = sumatorio_numerador / sumatorio_denominador

    return hpk
