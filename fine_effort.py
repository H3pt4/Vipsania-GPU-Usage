# Input: clade (e. g. NCBI Taxonomy ID)
# Input: taxnomic level (species, genus, family)
# Input: min assembly state (default: chromosome)
# Output: Computing time for finetuning / gene annotation
# Output: Average genome size

# For each species S (with seleced assembly) and taxonomic level L, finetune Vipsania once for each L, using all genomes of level L


# 1. For selected clade, get all clades of level L
# 2. For each level L, get all species S and their genomes with selected assembly state
# 3. For each genome G, get the genome size and calculate the average size
# 4. Compute the time used for finetuning Vipsania


from ete3 import NCBITaxa

# Get subclades of level L

def get_subclades(ncbi, taxid, level):
    """Return list of level-rank taxids that are descendants of taxid (or taxid itself)."""
    all_desc = ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
    all_desc.append(taxid)
    # Check for applicable level !!!


    ranks = ncbi.get_rank(all_desc)
    
    return [t for t, r in ranks.items() if r == level]


def get_species_in_clade(ncbi, taxid):
    """Return list of species-rank taxids that are descendants of taxid (or taxid itself)."""
    all_desc = ncbi.get_descendant_taxa(taxid)
    if not all_desc:
        ranks = ncbi.get_rank([taxid])
        if ranks.get(taxid) == 'species':
            return [taxid]
        return []
    ranks = ncbi.get_rank(all_desc)
    return [t for t, r in ranks.items() if r == 'species']


if __name__ == "__main__":
    ncbi = NCBITaxa()
    taxid = 9606  # Example: Homo sapiens
    level = "genus"
    subclades = get_subclades(ncbi, taxid, level)
    print(f"Subclades of taxid {taxid} at level '{level}': {subclades}")

    species = get_species_in_clade(ncbi, taxid)
    print(f"Species in clade of taxid {taxid}: {species}")