import csv
from datetime import datetime, timezone
from pathlib import Path

from evaluate_clade import analyse_clade

clades = {
    "Insecta": "50557",
    "Eukaryota": "2759",
    "Vertebrata": "7742",
    "Fungi": "4751",
}

summary = []

annotation_time_humans = 3080000000 / 92.
annotation_time_drosophila = 144000000 / 5.8
annotation_time_arabidopsis = 136000000 / 4.2
annotation_time_avg = (annotation_time_humans + annotation_time_drosophila + annotation_time_arabidopsis) / 3

finetuning_time = (158.6 + 158.4 + 158.4) / 3


for clade_name, taxid in clades.items():
    result = analyse_clade(
        taxid,
        taxon_level="species",
        min_assembly_level="chromosome",
    )
    rows = result["rows"]
    lengths = [int(row["total_sequence_length"]) for row in rows]
    n_species = len(rows)
    n_genera = len({row["genus"] for row in rows if row["genus"]})
    n_families = len({row["family"] for row in rows if row["family"]})
    finetuning_species_h = finetuning_time * n_species / 60
    finetuning_genera_h = finetuning_time * n_genera / 60
    finetuning_families_h = finetuning_time * n_families / 60
    annotation_time_h = sum(lengths) / annotation_time_avg / 60

    summary.append({
        "clade": clade_name,
        "n_species": n_species,
        "avg_genome_length_bp": sum(lengths) / len(lengths) if lengths else "",
        "total_genome_length_bp": sum(lengths),
        "n_genera": n_genera,
        "n_families": n_families,
        "finetuning_time_h (species)": finetuning_species_h,
        "finetuning_time_h (genera)": finetuning_genera_h,
        "finetuning_time_h (families)": finetuning_families_h,
        "annotation_time_h": annotation_time_h,
        "total_runtime_h (species)": finetuning_species_h + annotation_time_h,
        "total_runtime_h (genera)": finetuning_genera_h + annotation_time_h,
        "total_runtime_h (families)": finetuning_families_h + annotation_time_h,
    })

timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
output_path = Path("results") / f"clade_summary_{timestamp}.tsv"
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as output:
    columns = [
        "clade", "n_species", "avg_genome_length_bp",
        "total_genome_length_bp", "n_genera", "n_families",
        "finetuning_time_h (species)", "finetuning_time_h (genera)",
        "finetuning_time_h (families)", "annotation_time_h",
        "total_runtime_h (species)",
        "total_runtime_h (genera)", "total_runtime_h (families)",
    ]
    writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"Saved summary table: {output_path}")