from datetime import timedelta

from evaluate_clade import analyse_clade
import sys



if __name__ == "__main__":
    taxid = sys.argv[1]    # E.g. 9443
    tax_level = sys.argv[2] 
    

    finetuning_time = (158.6 + 158.4 + 158.4) / 3
    
    annotation_time_humans = 3080000000 / 92.1
    annotation_time_drosophila = 144000000 / 5.8
    annotation_time_arabidopsis = 136000000 / 4.2
    annotation_time_avg = (annotation_time_humans + annotation_time_drosophila + annotation_time_arabidopsis) / 3
    #print("Annotation runtime (homo sapiens): " + str(annotation_time_humans) + " bp/min")
    #print("Annotation runtime (drosophila): " + str(annotation_time_drosophila) + " bp/min")
    #print("Annotation runtime (arabidopsis): " + str(annotation_time_arabidopsis) + " bp/min")

    res = analyse_clade(taxid, taxon_level=tax_level, min_assembly_level="chromosome")

    total_runtime = finetuning_time * res["n_taxa"] + sum(res["lengths"]) / annotation_time_avg


    print("Total runtime: " + str(divmod(total_runtime, 60)[0])  + " (hours, mins)")

    # aggregated table with date
    # insecta, eukaryota, vertabrata (clades), fungi
    # # species, avg genome length, whole genome length, # (genera,ivmod(total_runtime, 60) families), time in h (and individually finetuning/inference)

    # table per clade
    # species, genus, family, accession number, assembly size
