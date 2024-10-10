# Average UK Rainfall (mm) for 1910 by month
# http://www.metoffice.gov.uk/climate/uk/datasets
rainfall = (('JAN',111.4),
            ('FEB',126.1),
            ('MAR', 49.9),
            ('APR', 95.3),
            ('MAY', 71.8),
            ('JUN', 70.2),
            ('JUL', 97.1),
            ('AUG',140.2),
            ('SEP', 27.0),
            ('OCT', 89.4),
            ('NOV',128.4),
            ('DEC',142.2),
           )

# (1) Use a list comprehension to create a list of month,rainfall tuples where
# the amount of rain was greater than 100 mm.
list_great100 = [i for i in rainfall if i[1] > 100]
print(list_great100)

# (2) Use a list comprehension to create a list of just month names where the
# amount of rain was less than 50 mm. 
list_less50 = [i for i in rainfall if i[1] < 50]
print(list_less50)
# (3) Now do (1) and (2) using conventional loops (you can choose to do 
# this before 1 and 2 !). 
list_great100=[]
for i in rainfall:
    if i[1] > 100:
        list_great100.append(i[1])
print(list_great100)

list_less50=[]
for i in rainfall:
    if i[1] < 50:
        list_less50.append(i[1])
print(list_less50)


# A good example output is:
#
# Step #1:
# Months and rainfall values when the amount of rain was greater than 100mm:
# [('JAN', 111.4), ('FEB', 126.1), ('AUG', 140.2), ('NOV', 128.4), ('DEC', 142.2)]
# ... etc.

