import csv
import sys
from fuzzywuzzy import fuzz 

#Define function
def is_an_oak(name):
    """ Returns True if name is starts with 'quercus' """
    """
    Returns True if name starts with 'quercus'
    
    >>> is_an_oak('Quercus robur')
    True
    >>> is_an_oak('Fagus sylvatica')
    False
    >>> is_an_oak('Quercuss')
    True
    """
    return fuzz.partial_ratio(name.lower(), 'quercus') > 85  
#This line of code checks whether the given string name,
#  when compared to 'quercus' (the genus name for oaks), 
# has a similarity score greater than 85%. If the similarity score exceeds 85%, 
# it considers the input name to be a close enough match to 'quercus' 
# (indicating that it's likely an oak tree), 
# and the function returns True. Otherwise, it returns False, 
# meaning the string is not considered a match for the oak genus.
#of course you can change the ratio 85 with your own choice 

def main(argv): 
    f = open('../data/TestOaksData.csv','r')
    g = open('../results/oak_result.csv', 'w', newline='')
    taxa = csv.reader(f)
    csvwrite = csv.writer(g)
    oaks = set()
    for row in taxa:
        print(row)
        print ("The genus is: ") 
        print(row[0] + '\n')
        if is_an_oak(row[0]):
            print('FOUND AN OAK!\n')
            csvwrite.writerow([row[0], row[1]])    

    return 0
    
if (__name__ == "__main__"):
    status = main(sys.argv)