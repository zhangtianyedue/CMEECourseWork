birds = ( ('Passerculus sandwichensis','Savannah sparrow',18.7),
          ('Delichon urbica','House martin',19),
          ('Junco phaeonotus','Yellow-eyed junco',19.5),
          ('Junco hyemalis','Dark-eyed junco',19.6),
          ('Tachycineata bicolor','Tree swallow',20.2),
         )

#(1) Write three separate list comprehensions that create three different
# lists containing the latin names, common names and mean body masses for
# each species in birds, respectively. 
Latin_names=[i[0] for i in birds]
print(Latin_names)
Common_names=[i[1] for i in birds]
print(Common_names)
Mean_body_Mass=[i[2] for i in birds  ]
print(Mean_body_Mass)
# (2) Now do the same using conventional loops (you can choose to do this 
# before 1 !). 
#For Version 
Latin_names=[]
for i in birds:
    Latin_names.append(i[0])
print(Latin_names) 
#while Version
Latin_names=[]
i=0
while i < len(birds):
    Latin_names.append(birds[i][0])
    i+=1
print(Latin_names)
#print(len(birds))

#For version 
Common_names=[]
for i in birds:
    Common_names.append(i[1])
print(Common_names)

#While Version 
Common_names=[]
i=0
while i < len(birds):
    Common_names.append(birds[i][1])
    i+=1
print(Common_names)

#For version 
Mean_body_Mass=[]
for i in birds:
    Mean_body_Mass.append(i[2])
print(Mean_body_Mass)

#While version 
Mean_body_Mass=[]
i=0
while i < len(birds):
    Mean_body_Mass.append(birds[i][2])
    i+=1
print(Mean_body_Mass)





# A nice example out out is:
# Step #1:
# Latin names:
# ['Passerculus sandwichensis', 'Delichon urbica', 'Junco phaeonotus', 'Junco hyemalis', 'Tachycineata bicolor']
# ... etc.
 