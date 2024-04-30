import pandas as pd

if __name__ == '__main__':
    angle = 45
    txt_file = "raw_code/hp_0.07_slab.csv_H-60_resultados_fix_seed_10sim.txt"
    csv_file = 'uspekpy/output_fix_seed_10sim.csv'

    with open(txt_file, 'r') as file:
        lines = file.readlines()

    raw_values = {}
    for idx, line in enumerate(lines):
        if f"angle: {angle} deg" in line:
            next_line = lines[idx + 1]
            hk_mean = next_line.split(":")[-1].strip()
            hk_mean = hk_mean.strip("[]")
            raw_values['Mean conv. coefficient. Mean'] = hk_mean

            next_line = lines[idx + 2]
            hk_std = next_line.split(":")[-1].strip()
            raw_values['Mean conv. coefficient. Standard deviation'] = hk_std

            next_line = lines[idx + 3]
            hk_ru = next_line.split(":")[-1].strip()
            raw_values['Mean conv. coefficient. Relative uncertainty'] = float(hk_ru) / 100

            next_line = lines[idx + 4]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            ru = float(next_line.split(':')[3].strip()) / 100
            keys = ['Mean energy Mean', 'Mean energy Standard deviation', 'Mean energy Relative uncertainty']
            raw_values.update(dict(zip(keys, [mean, std, ru])))

            next_line = lines[idx + 5]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            ru = float(next_line.split(':')[3].strip()) / 100
            keys = ['Mean kerma Mean', 'Mean kerma Standard deviation', 'Mean kerma Relative uncertainty']
            raw_values.update(dict(zip(keys, [mean, std, ru])))

            next_line = lines[idx + 7]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            keys = ['HVL1 Al Mean', 'HVL1 Al Standard deviation']
            raw_values.update(dict(zip(keys, [mean, std])))

            next_line = lines[idx + 9]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            keys = ['HVL2 Al Mean', 'HVL2 Al Standard deviation']
            raw_values.update(dict(zip(keys, [mean, std])))

            next_line = lines[idx + 11]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            keys = ['HVL1 Cu Mean', 'HVL1 Cu Standard deviation']
            raw_values.update(dict(zip(keys, [mean, std])))

            next_line = lines[idx + 13]
            mean = next_line.split(':')[1].split(',')[0].strip()
            std = next_line.split(':')[2].split(',')[0].strip()
            keys = ['HVL2 Cu Mean', 'HVL2 Cu Standard deviation']
            raw_values.update(dict(zip(keys, [mean, std])))

    df_r = pd.DataFrame([raw_values.keys(), raw_values.values()]).T
    df_r.columns = ['Name', 'Raw']
    df_r['Raw'] = df_r['Raw'].astype(float)

    df_u = pd.read_csv(csv_file)
    df_u = df_u.iloc[27:, :]
    df_u.columns = ['Name', 'Uspek']
    df_u['Uspek'] = df_u['Uspek'].astype(float)

    df = pd.merge(df_u, df_r, on='Name')
    df['(Uspek-Raw)/Raw (%)'] = (df['Uspek'] - df['Raw']) / df['Raw'] * 100

# Difference when taking the random values:
# Raw code:  first interpolate, the take random values
# XCB: first take random values, then interpolate

# Difference in spectrum bins:
# Raw code: some values are masked
# XCB: all values considered
