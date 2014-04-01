#!/usr/bin/python3

import sys
import re

file = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

# naming skemes for scientists and experiments
# Matt  groups
matt_soil = re.compile('(\d{4})-.*-(\d+)')
matt_soilNT= re.compile('(NT)-.*-(\d+)')
matt_syn = re.compile('([A-C]\d)-(\d{2})-(\d+)')
niel_2X = re.compile('([a-b]\d)-(.*)-.*')
niel_4XJT = re.compile('([1-4][A-D])-(.*)-(\d+)')
niel_no_= re.compile('([a-d])(.*)dn?(\d+)')
alyx_contam = re.compile('(\d+[AB])-(.*)-(D[-\d]\d*)')
alyx_treat = re.compile('(\w+)-(\d+)-D-(\d+)')







r1 = []
r2 = []
samples = []
missing =[]
r = 1
for line in file:
    line = line.strip()
    if r%2!=0: #odd
        fastq1=line 
        line = line.split('_')     # splits fastq file name at the sample name and saves this as name
        name1 = line[0]
        name = name1.replace('N-T','NT') # replace the N-T in matt's with NT for mouse
        
# If statement to determine how to handle the number of split names and what group goes where

       
        
               
        soil = matt_soil.match(name)
        soilNT = matt_soilNT.match(name)
        syn = matt_syn.match(name)
        a1b4=niel_2X.match(name)
        A14D=niel_4XJT.match(name)
        no_ = niel_no_.match(name) 
        contam = alyx_contam.match(name)
        treat = alyx_treat.match(name)
         # Matt soil
        if soil:
            scientist='matt'
            group = 'soil'
            cage=soil.group(1)
            mouse = soil.group(1)
            day = str(26-int(soil.group(2)))
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1
        elif soilNT: 
            scientist='matt'
            group = 'soil'
            cage=soilNT.group(1)
            mouse = soilNT.group(1)
            day = str(26-int(soilNT.group(2)))
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
       
            samples.append(newname)
            r1.append(fastq1)
            test=name1
        
        # Matt Synthetic
        elif syn:
            scientist='matt'
            group = 'syn'
            cage=syn.group(1)[0]
            mouse = syn.group(1)[1]
            if int(syn.group(2)) == 10:
                    day = str(int(syn.group(3))-22)
            if int(syn.group(2)) == 11:
                    day = str(int(syn.group(3))+9)
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
        # Niel a1a2 b3b4
        elif a1b4:
            name=name.split('-')
            if name[0]=='a1' or name[0]=='a2':
                group = 'H3'
            elif name[0]=='b3' or name[0]=='b4':
                group='C1'
            else:
                group= 'WRONG'        
            scientist='niel'
            cage = name[0]
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
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
            
            #Niel 1A-4D
        elif A14D:
            if A14D.group(1) =='1A':
                group= 'H1'
            elif A14D.group(1) =='2B':
                group = 'H2'
            elif A14D.group(1)=='3C':
                group = 'C2'
            elif A14D.group(1) == '4D':
                group = 'C3'            
            else:
                group = 'WRONG'
            scientist = 'niel'
            cage = A14D.group(1)
            mouse = A14D.group(2)
            day = str(22-int(A14D.group(3)))    
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
        #Neils tricky samples 
        elif no_:
        
            if no_.group(1) =='a':
                group= 'H1'
            elif no_.group(1) =='b':
                group = 'H2'
            elif no_.group(1)=='c':
                group = 'C2'
            elif no_.group(1) == 'd':
                group = 'C3'            
            else:
                group = 'WRONG'

            if no_.group(1) =='a':
                cage= '1A'
            elif no_.group(1) =='b':
                cage = '2B'
            elif no_.group(1)=='c':
                cage = '3C'
            elif no_.group(1) == 'd':
                cage = '4D'            
            else:
                cage = 'WRONG'

            scientist = 'niel'
            mouse = no_.group(2)
            day = str(22-int(no_.group(3)))    
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1

        #Alyx contaminated
        elif contam:
            scientist = 'alyx'
            group = contam.group(1)[0:-1]
            cage = contam.group(1)[-1]
            mouse = contam.group(2)
            day = contam.group(3)
            if '-' in day:
                day = day.split('-')
                day = day[-1]
                day = str(15-int(day))
            elif 'D0' in day:
                day = '15'
            else:
                day = day.split('D') 
                day = str(1+int(day[-1]))


            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
           
        #Alyx's treated samples (Chloroform and Carey Blair)
        elif treat:
            scientist = 'alyx'
            if 'CHL' in treat.group(1):
                group = 'Chloroform'
            elif treat.group(1)=='N' or 'Y' or 'NC':
                group = 'CairyBlair'
            cage = treat.group(1)
            mouse = treat.group(2).strip('-D')  # some samples are mouseD- others are mouse-D- this removes the D or the -D from the mouse string        
            day = str(15-int(treat.group(3)))
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
            
        else:
            missing.append(name)    
        
    if r%2==0:  #even
        fastq2=line 
        line = line.split('_')     # splits fastq file name at the sample name and saves this as name
        name2 = line[0]
        if name2==test:          # test to see if name matches that of the last successful fastq
            r2.append(fastq2)
    r += 1
print(len(r2))
print('You missed these you fool!')
print(len(missing))
print(missing)
for i in range(len(samples)):
    print(samples[i],r1[i],r2[i],sep='\t',end='\n',file=outfile)
    
file.close()
outfile.close()

