taxa = [ ('Myotis lucifugus','Chiroptera'),
         ('Gerbillus henleyi','Rodentia',),
         ('Peromyscus crinitus', 'Rodentia'),
         ('Mus domesticus', 'Rodentia'),
         ('Cleithrionomys rutilus', 'Rodentia'),
         ('Microgale dobsoni', 'Afrosoricida'),
         ('Microgale talazaci', 'Afrosoricida'),
         ('Lyacon pictus', 'Carnivora'),
         ('Arctocephalus gazella', 'Carnivora'),
         ('Canis lupus', 'Carnivora'),
        ]

# Write a python script to populate a dictionary called taxa_dic derived from
# taxa so that it maps order names to sets of taxa and prints it to screen.
taxa_dict = {}
for species, category in taxa:#check if the catagory already exist 
    if taxa_dict.get(category) is None: 
        taxa_dict[category] = set() # set a new set to reviece new category 
        taxa_dict[category].add(species) 
print(taxa_dict)
# An example output is:
#  
# 'Chiroptera' : set(['Myotis lucifugus']) ... etc. 
# OR, 
# 'Chiroptera': {'Myotis  lucifugus'} ... etc

#### Your solution here #### 

# Now write a list comprehension that does the same (including the printing after the dictionary has been created)  
 
#### Your solution here #### 
taxa_dict = {category: {species for species, cat in taxa if cat == category} for _, category in taxa}
#Outer dictionary comprehension:
#{category: ... for _, category in taxa}: This is a dictionary comprehension 
# that iterates through all the category values in taxa and creates a key-value pair for each category.

#Inner set comprehension:
#{species for species, cat in taxa if cat == category}: 
# This is a set comprehension that iterates through all the entries (species and cat) in taxa and checks if cat is equal to the outer category. 
# If they match, the species is added to the set corresponding to that category



