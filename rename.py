#!/usr/bin/python3

import sys
import re

file = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

# naming skemes for scientists and experiments
# Matt  groups
matt_soil = '1810 , 1840 , 1856, 1872, , NT, soil '
matt_syn = 'A0, A1, A2, B0, B1, B2, C0, C1, C2, C3 , C4' 
niel = 'a1,a2,b3,b4'





r1 = []
r2 = []
samples = []
missing =[]
r = 1
for line in file:
    line = line.strip()
    if r%2!=0: #odd
        r1.append(line)
        line = line.split('_')     # splits fastq file name at the sample name and saves this as name
        name = line[0]
        name = name.replace('N-T','NT') # replace the N-T in matt's with NT for mouse
        name=name.split('-')
        id = name[0] # identifier used to set up name skeme
        
# If statment to determine how to handle the number of split names and what group goes where


        # Matt soil

        if id  in matt_soil:
            scientist='matt'
            group = 'soil'
            cage=id
            mouse = id
            try:                               # inoclulum has inoc not number as the last entry
                day = str(26-int(name[-1]))
            except ValueError:
                day = name[-1]
       
        # Matt synthetic
        
        elif id  in matt_syn:
            scientist='matt'
            group ='syn'
            m = re.search('(\w)(\d)',id)           # The mice are in cages A,B,C followed by a 1 digit number to denote the mouse
            cage = m.group(1)
            mouse = m.group(2)
            try:                               # labeled with the date, day 1 is 10-23 so for days in nov. add 9 (31-22) to day
                if int(name[-2]) == 10:
                    day = str(int(name[-1])-22)
                if int(name[-2]) == 11:
                    day = str(int(name[-1])+9)
            
            except ValueError:
                day = name[-1]
        # Niel humanized

        elif id in niel:
            scientist = 'niel'
            if id=='a1' or 'a2':
                group = 'H3'
            if id == 'b3' or 'b4':
                group = 'C1'
            cage = id
            mouse = name[1]
            if 'd' in name[-1]:
                try: 
                    m = re.search('(?<=\w)\d+',name[-1])  # take the digits after the last letter character
                    
                    day = str(22-int(m.group(0))) 
                except ValueError:
                    day = 'WRONG'
            else:
                try:
                    day = str(22-int(name[-1]))
                except ValueError:
                    day = 'WRoNG'

        # Niels's tricky labeling a,b,c, .... dn 

        elif len(id) >4:
            if id[0]=='a':
                scientist='neil'
                group = 'H1'
                cage = '1A'
                mouse = id[1:id.find('d')]
                
                last = id[id.find('d'):len(id)]
                m = re.search('(?<=\w)\d+',last)  # take the digits after the last letter character
                day = str(22-int(m.group(0))) 

            elif id[0]=='b':
                scientist='neil'
                group = 'H2'
                cage = '2B'
                mouse = id[1:id.find('d')]
                
                last = id[id.find('d'):len(id)]
                m = re.search('(?<=\w)\d+',last)  # take the digits after the last letter character
                day = str(22-int(m.group(0))) 

            elif id[0]=='c':
                scientist='neil'
                group = 'C2'
                cage = '3C'
                mouse = id[1:id.find('d')]
                
                last = id[id.find('d'):len(id)]
                m = re.search('(?<=\w)\d+',last)  # take the digits after the last letter character
                day = str(22-int(m.group(0))) 

        
            elif id[0]=='d':
                group = 'C4'
                cage = '4D'
                mouse = id[1:id.find('d')]
                
                last = id[id.find('d'):len(id)]
                m = re.search('(?<=\w)\d+',last)  # take the digits after the last letter character
                day = str(22-int(m.group(0))) 





        # Alyx Samples


        







        newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
       
        samples.append(newname)
    else:
        missing.append(name)    
        
    if r%2==0:  #even
        r2.append(line)
    r += 1
print(len(r2))
print('You missed these you fool!')
print(missing)
for i in range(len(r1)):
    print(samples[i],r1[i],r2[i],sep='\t',end='\n',file=outfile)
    
file.close()
outfile.close()

