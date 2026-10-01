from evaluate_clade import analyse_clade
import sys



if __name__ == "__main__":
    taxid = sys.argv[1]    # E.g. 9443
    tax_level = sys.argv[2] 
    

    res = analyse_clade(taxid, taxon_level=tax_level, min_assembly_level="chromosome")
    print(res["n_taxa"])    # number of unique families
    print(res["lengths"])   # list of genome lengths in bp
    res["rows"]      # the table, as a list of dicts