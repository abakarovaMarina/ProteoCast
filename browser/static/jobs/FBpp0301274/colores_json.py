import pandas as pd
import matplotlib.pyplot as plt
import json

# Ruta del archivo CSV
file_path = 'res_FBpp0301274.csv'

# Leer el archivo CSV
df = pd.read_csv(file_path, header=None, names=["number", "value"])

# Normalización de valores para determinar el rango
norm = plt.Normalize(vmin=df['value'].min(), vmax=df['value'].max())
cmap = plt.cm.viridis  # Colormap deseado

# Función para mapear los valores de la columna 'value' a colores RGB
def map_color(value):
    rgba = cmap(norm(value))  # Devuelve valores RGBA en el rango 0-1
    rgb = tuple(int(c * 255) for c in rgba[:3])  # Convertir a RGB (0-255)
    return {"r": rgb[0], "g": rgb[1], "b": rgb[2]}  # Formato deseado

# Aplicar la función map_color a la columna 'value' y crear la columna de colores
df['color'] = df['value'].apply(map_color)

print(df.head())
# Crear el JSON de colores de residuos
residue_colors = df[['number', 'color']].to_dict(orient='records')

# Guardar el JSON
output_file = open('residue_colors.json','w')
output_file.write('[')
#with open(output_file, 'w') as f:
#    json.dump(residue_colors, f, indent=4)  # Indent para mejor legibilidad en el archivo JSON
for i in range(0,len(df)):
   res = df['number'][i]
   color = df['color'][i]
   if i == 0:
      output_file.write("{'number': '"+str(res)+"', 'color': "+str(color)+"}")
   else:
      output_file.write(",{'number': '"+str(res)+"', 'color': "+str(color)+"}")

output_file.write(']')
