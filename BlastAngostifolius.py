import os
from Bio import SeqIO

# Archivos FASTA originales de donde se crearían las bases de datos
proteinas_fasta = "ProteinasAngostifolius.fasta"
genes_fasta = "GenesAngostifolius.fasta"

# Directorios de salida
output_dir_proteinas = "./Proteinas_Seleccionadas_Angostifolius"
output_dir_genes = "./Genes_Ortologos_Resultados_Angostifolius"

# Crear directorios de salida si no existen
os.makedirs(output_dir_proteinas, exist_ok=True)
os.makedirs(output_dir_genes, exist_ok=True)

# Función para guardar una secuencia en formato FASTA
def guardar_secuencia_fasta(fasta_path, seq_id, output_path):
    found = False
    with open(output_path, "w") as out_fasta:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id == seq_id:
                SeqIO.write(record, out_fasta, "fasta")
                print(f"Secuencia {seq_id} guardada en {output_path}.")
                found = True
                break
    if not found:
        print(f"Secuencia {seq_id} no encontrada en {fasta_path}.")

# Procesar archivos gene.fasta en el directorio actual
for archivo in os.listdir("."):
    if archivo.endswith("_gene.fasta"):
        nombre_gene = archivo.split("_")[0]
        print(f"Procesando {archivo}...")

        # Ejecutar BLASTx contra ProteinasPhaseolusDB
        blastx_out = f"{output_dir_proteinas}/{nombre_gene}_blastx.out"
        comando_blastx = f"blastx -query {archivo} -db ProteinasAngostifoliusDB -outfmt '6 sseqid pident length evalue' -max_target_seqs 1 -out {blastx_out}"
        os.system(comando_blastx)

        # Leer el mejor hit de BLASTx
        try:
            with open(blastx_out) as f:
                best_hit = f.readline().strip().split()
                if not best_hit:
                    print(f"No se encontraron hits significativos para {archivo}.")
                    continue
                id_proteina, pident, length, evalue = best_hit

            # Guardar la secuencia proteica del mejor hit en el directorio de salida
            secuencia_proteina_path = f"{output_dir_proteinas}/{id_proteina}.fasta"
            guardar_secuencia_fasta(proteinas_fasta, id_proteina, secuencia_proteina_path)

            # Ejecutar TBLASTN con la secuencia proteica contra GenesPhaseolusDB
            tblastn_out = f"{output_dir_genes}/{id_proteina}_tblastn.out"
            comando_tblastn = f"tblastn -query {secuencia_proteina_path} -db GenesAngostifoliusDB -outfmt '6 sseqid pident length evalue' -max_target_seqs 1 -out {tblastn_out}"
            os.system(comando_tblastn)

            # Leer el mejor hit de TBLASTN
            with open(tblastn_out) as f:
                best_hit_gene = f.readline().strip().split()
                if not best_hit_gene:
                    print(f"No se encontraron hits significativos para {id_proteina} en TBLASTN.")
                    continue
                id_gene, pident_gene, length_gene, evalue_gene = best_hit_gene

            # Guardar la secuencia del gen homólogo en el directorio de salida
            gen_homologo_path = f"{output_dir_genes}/Ortólogoa{seq_id}_{id_gene}.fasta"
            guardar_secuencia_fasta(genes_fasta, id_gene, gen_homologo_path)

        except FileNotFoundError:
            print(f"No se encontró el archivo de salida {blastx_out} o {tblastn_out}.")
