# List of magnitudes independent from the angle
magnitudes = [(media_energia, desviacion_energia, coef_variacion_energia),
              (media_kerma, desviacion_kerma, coef_variacion_kerma),
              (mean_hvl, sd_hvl, v_hvl),
              (mean_hvl2, sd_hvl2, v_hvl2),
              (mean_hvlCu, sd_hvlCu, v_hvlCu),
              (mean_hvl2Cu, sd_hvl2Cu, v_hvl2Cu)]

# List of angles and their respective values and uncertainties of hpk
angles = [(0, hpk, sd_hpk, v_hpk),
           (15, hpk15, sd_hpk_15, v_hpk_15),
           (30, hpk30, sd_hpk_30, v_hpk_30),
           (45, hpk45, sd_hpk_45, v_hpk_45),
           (60, hpk60, sd_hpk_60, v_hpk_60),
           (75, hpk75, sd_hpk_75, v_hpk_75),
           (90, hpk90, sd_hpk_90, v_hpk_90),
           (180, hpk180, sd_hpk_180, v_hpk_180)]

# Build path of output file
file_path = os.path.join(ruta, f"{f_m}_{calidad}_results.txt")

# Open output file in append mode
with open(file_path, 'a') as file:
    for angle, hpk, sd_hpk, v_hpk in angles:
        file.write(f"Quality: {calidad}; Angle: {angle} deg\n")
        archivo.write(f"Mean value of conversion coeff: {hpk_actual}\n")
        archivo.write(f"Std deviation: {sd_hpk_actual}\n")
        archivo.write(f"ur%(hK): {v_hpk_actual}\n")
        archivo.write(
            f"Mean E keV: {media_energia}, std dev keV: {desviacion_energia}, ur%: {coef_variacion_energia}\n")
        archivo.write(
            f"Kerma Medio - Valor medio: {media_kerma}, Desviación estándar: {desviacion_kerma}, ur%: {coef_variacion_kerma}\n\n")
        archivo.write(
            f"1st HVL mm Al- mean value: {mean_hvl}, std deviation mm Al: {sd_hvl}, ur%: {v_hvl}\n\n")
        archivo.write(
            f"2nd HVL mm Al- mean value: {mean_hvl2}, std deviation mm Al: {sd_hvl2}, ur%: {v_hvl2}\n\n")
        archivo.write(
            f"1st HVL mm Cu- mean value: {mean_hvlCu}, std deviation mm Cu: {sd_hvlCu}, ur%: {v_hvlCu}\n\n")
        archivo.write(
            f"2nd HVL mm Cu- mean value: {mean_hvl2Cu}, std deviation mm Cu: {sd_hvl2Cu}, ur%: {v_hvl2Cu}\n\n")

tiempo_final_columns_9 = time()

tiempo_ejecucion_columns_9 = tiempo_final_columns_9 - tiempo_inicial_columns_9

print('El tiempo de ejecucion para hktable columns 9 en (s) fue:', tiempo_ejecucion_columns_9)
