from collections import defaultdict

# Using defaultdict to optimize dictionary creation
def populate_taxa_dict(taxa):
    """
    Populates a dictionary mapping order names to sets of taxa.

    Parameters:
        taxa (list of tuples): A list where each tuple contains a species and its category.

    Returns:
        dict: A dictionary where keys are categories and values are sets of species.
    """
    taxa_dict = defaultdict(set)
    for species, category in taxa:
        taxa_dict[category].add(species)  # Automatically initializes the set if the key doesn't exist
    
    return taxa_dict

# List comprehension version
def populate_taxa_dict_comprehension(taxa):
    """
    Populates a dictionary mapping order names to sets of taxa using list comprehension.

    Parameters:
        taxa (list of tuples): A list where each tuple contains a species and its category.

    Returns:
        dict: A dictionary where keys are categories and values are sets of species.
    """
    return {category: {species for species, cat in taxa if cat == category} for _, category in taxa}

# Main execution
if __name__ == "__main__":
    taxa = [
        ('Myotis lucifugus', 'Chiroptera'),
        ('Gerbillus henleyi', 'Rodentia'),
        ('Peromyscus crinitus', 'Rodentia'),
        ('Mus domesticus', 'Rodentia'),
        ('Cleithrionomys rutilus', 'Rodentia'),
        ('Microgale dobsoni', 'Afrosoricida'),
        ('Microgale talazaci', 'Afrosoricida'),
        ('Lyacon pictus', 'Carnivora'),
        ('Arctocephalus gazella', 'Carnivora'),
        ('Canis lupus', 'Carnivora'),
    ]

    # Using defaultdict
    taxa_dict_defaultdict = populate_taxa_dict(taxa)
    print("Using defaultdict:")
    print(dict(taxa_dict_defaultdict))  # Convert to regular dict for display

    # Using list comprehension
    taxa_dict_comprehension = populate_taxa_dict_comprehension(taxa)
    print("\nUsing list comprehension:")
    print(taxa_dict_comprehension)



