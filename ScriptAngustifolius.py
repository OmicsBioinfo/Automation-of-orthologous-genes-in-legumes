import os
import tkinter as tk
from tkinter import messagebox
from Bio import Entrez
import subprocess

# Configura tu correo electrónico para acceder a la API de NCBI
Entrez.email = "frankeusebiofernandez@gmail.com"

# Función para obtener coordenadas de un ID de gen usando NCBI
def buscar_coordenadas_ncbi(gene_id):
    try:
        # Fetch information from the gene database
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        records = Entrez.read(handle)
        handle.close()

        # Explorar cada parte de la respuesta hasta encontrar las coordenadas
        for gene in records:
            for locus in gene.get("Entrezgene_locus", []):
                ref_id = locus.get("Gene-commentary_accession")
                # Verifica que contenga las coordenadas de inicio y fin
                if "Gene-commentary_seqs" in locus:
                    seq_info = locus["Gene-commentary_seqs"][0]
                    if "Seq-loc_int" in seq_info:
                        start = seq_info["Seq-loc_int"]["Seq-interval"]["Seq-interval_from"]
                        end = seq_info["Seq-loc_int"]["Seq-interval"]["Seq-interval_to"]
                        return ref_id, int(start), int(end)

        # Si no se encontraron coordenadas, retornar None
        print(f"Advertencia: No se encontraron coordenadas para el gen {gene_id}.")
        return None

    except Exception as e:
        print(f"Error al obtener coordenadas de {gene_id}: {e}")
        return None

def obtener_secuencias():
    ids_genes = entry_ids.get().strip().split(",")  # Obtener IDs de genes de la entrada
    if not ids_genes:
        messagebox.showwarning("Advertencia", "Por favor ingresa los IDs de los genes.")
        return

    # Crear el nombre del archivo usando los IDs ingresados
    ids_concatenados = "_".join(id.strip() for id in ids_genes)
    coordenadas_file_path = f"CoordenadasAngostifolius{ids_concatenados}.txt"

    with open(coordenadas_file_path, 'w') as coord_file:
        for gene_id in ids_genes:
            gene_id = gene_id.strip()
            # Intentar obtener coordenadas desde NCBI
            coordenadas = buscar_coordenadas_ncbi(gene_id)
            if coordenadas:
                ref_id, start, end = coordenadas
                coord_file.write(f"Gene ID: {gene_id}, Reference: {ref_id}, Start: {start}, End: {end}\n")
            else:
                messagebox.showwarning("Advertencia", f"ID de gen {gene_id} no encontrado o no disponible.")

    messagebox.showinfo("Éxito", f"Coordenadas guardadas en {coordenadas_file_path}")

    # Ejecuta el script de Bash para descargar las secuencias
    try:
        subprocess.run(["./descargar_secuencias1.sh", coordenadas_file_path], check=True)
        messagebox.showinfo("Descarga Completa", "Las secuencias han sido descargadas exitosamente.")

        # Ejecuta el script de Python Soya_Phaseolus.py al finalizar descargar_secuencias1.sh
        subprocess.run(["python3", "BlastAngostifolius.py"], check=True)
        messagebox.showinfo("Proceso Completo", "Lupinus.py se ha ejecutado exitosamente.")

    except subprocess.CalledProcessError:
        messagebox.showerror("Error", "Hubo un error al ejecutar descargar_secuencias1.sh o Lupinus.py.")
# Función para obtener las coordenadas y generar el archivo de coordenadas



# Configuración de la ventana principal
root = tk.Tk()
root.title("Obtener Secuencias Genómicas")

# Etiquetas y entradas
tk.Label(root, text="IDs de Genes (separados por comas):").grid(row=0, column=0, padx=10, pady=10)
entry_ids = tk.Entry(root, width=100)
entry_ids.grid(row=0, column=1, padx=10, pady=10)

tk.Button(root, text="Obtener Secuencias", command=obtener_secuencias).grid(row=1, column=0, columnspan=2, padx=10, pady=10)

# Iniciar la aplicación
root.mainloop()
