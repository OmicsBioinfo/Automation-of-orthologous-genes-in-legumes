#!/bin/bash

# Verifica que se haya proporcionado un archivo como argumento
if [[ $# -ne 1 ]]; then
  echo "Uso: $0 <archivo_de_coordenadas>"
  exit 1
fi

# Archivo de coordenadas pasado como argumento
COORDENADAS_FILE="$1"

# Verifica que el archivo de coordenadas existe
if [[ ! -f "$COORDENADAS_FILE" ]]; then
  echo "Error: El archivo '$COORDENADAS_FILE' no existe."
  exit 1
fi

# Lee cada línea del archivo de coordenadas
while IFS= read -r line; do
  # Extrae los valores de ID del gen, referencia, inicio y fin usando sed
  gene_id=$(echo "$line" | sed -n 's/.*Gene ID: \([^,]*\).*/\1/p')
  ref_id=$(echo "$line" | sed -n 's/.*Reference: \([^,]*\).*/\1/p')
  start=$(echo "$line" | sed -n 's/.*Start: \([^,]*\).*/\1/p')
  end=$(echo "$line" | sed -n 's/.*End: \([^,]*\).*/\1/p')

  # Verifica que los datos necesarios estén presentes
  if [[ -n "$gene_id" && -n "$ref_id" && -n "$start" && -n "$end" ]]; then
    # Construye el nombre del archivo de salida con el ID del gen
    output_file="${gene_id}_gene.fasta"

    # Usa curl para descargar la secuencia con las coordenadas especificadas
    curl -o "$output_file" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ref_id&rettype=fasta&retmode=text&seq_start=$start&seq_stop=$end"

    # Verifica si la descarga fue exitosa
    if [[ $? -eq 0 ]]; then
      echo "Secuencia para $gene_id descargada en $output_file"
    else
      echo "Error al descargar la secuencia para $gene_id"
    fi
  else
    echo "Advertencia: Datos incompletos en la línea: $line"
  fi
done < "$COORDENADAS_FILE"

echo "Descarga completada para todos los genes en el archivo."
