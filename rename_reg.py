#!/usr/bin/python3

import sys
import re

file = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

# naming skemes for innocula
alyx_contam_inocA = re.compile('(\d+)A-inoculum')
niel_inoc = re.compile('([CH]\D{5,6}\d+)-\d{4}')
alyx_contam_inocB = re.compile('DA(\d+)')

# naming skemes for scientists and experiments
matt_soil = re.compile('(\d{4})-.*-(\d+)')
matt_soilNT= re.compile('(NT)-.*-(\d+)')
matt_syn = re.compile('([A-C]\d)-(\d{2})-(\d+)')
niel_2X = re.compile('([a-b]\d)-(.*)-.*')
niel_4XJT = re.compile('([1-4][A-D])-(.*)-(\d+)')
niel_no_= re.compile('([a-d])([^-]*)dn?(\d+)')
alyx_contam = re.compile('(\d+[AB])-(.*)-(D[-\d]\d*)')
alyx_treat = re.compile('([NYC]\D*)-(\w+)-*D-?(\d+)')








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
        
        if name == '846dn03':  # corrects a typo which labeled one of niels samples in correctly
            name = 'd846dn03'
# If statement to determine how to handle the number and what group goes where

## INOCULA ###
        alyx_con_inocA = alyx_contam_inocA.match(name)
        n_inoc = niel_inoc.match(name)   
        alyx_con_inocB = alyx_contam_inocB.match(name)
##Samples ##       
        soil = matt_soil.match(name)
        soilNT = matt_soilNT.match(name)
        syn = matt_syn.match(name)
        a1b4=niel_2X.match(name)
        A14D=niel_4XJT.match(name)
        no_ = niel_no_.match(name) 
        contam = alyx_contam.match(name)
        treat = alyx_treat.match(name)
        
## INOCULA ###        
        if alyx_con_inocA:
            scientist='alyx'
            group = alyx_con_inocA.group(1)
            cage = 'A'
            mouse= 'INOC'
            day = mouse
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1     # save original name to be matched with other fastq later
        
        elif alyx_con_inocB:
            scientist='alyx'
            group = alyx_con_inocB.group(1)
            cage = 'B'
            mouse= 'INOC'
            day = mouse
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1  
        elif n_inoc:
            scientist = 'niel'
            if n_inoc.group(1)=='Cancer7':
                group = 'C1'
            elif n_inoc.group(1)=='Cancer16':  
                group = 'C2'
            elif n_inoc.group(1)=='Cancer21':              
                group = 'C3'
            elif n_inoc.group(1)=='Healthy2':    
                group = 'H1'
            elif n_inoc.group(1)=='Healthy11':
                group = 'H2'
            elif n_inoc.group(1)=='Healthy18':    
                group = 'H3'
            else:
            	group = 'WRONG'    
            cage = 'INOC'    
            mouse = cage
            day = cage    
                
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1  
        elif name == 'soil-Inoc':
            scientist = 'matt'
            group = 'soil'
            cage = 'INOC'
            mouse = cage
            day = cage
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1  
            
        elif name == 'Cult-inoc-10-24':
            scientist = 'matt'
            group = 'syn'
            cage = 'INOC'
            mouse = cage
            day = cage
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1  
                           
## Samples ##                
                
         # Matt soil
        elif soil:
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
            cage=syn.group(1)[0]   # All in format 'A2' were A is the cage and 2 is the mouse
            mouse = syn.group(1)[1]
            if int(syn.group(2)) == 10:         # The days were in date format with 10-23 as day 1
                    day = str(int(syn.group(3))-22)
            if int(syn.group(2)) == 11:
                    day = str(int(syn.group(3))+9)  # 11-1 was really day 10 
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
        # Niel a1a2 b3b4  Used reg exp to identify but then for some reason I split the name variable to allocate, It was just easier
        elif a1b4:
            name=name.split('-')    
            if name[0]=='a1' or name[0]=='a2':  # Had same innocula
                group = 'H3'
            elif name[0]=='b3' or name[0]=='b4':
                group='C1'
            else:
                group= 'WRONG'        
            scientist='niel'
            cage = name[0]
            mouse = name[1]
            if 'd' in name[-1]:           # Some had d in the day section
                try: 
                    m = re.search('(?<=\w)\d+',name[-1])  # take the digits after the last letter character
                    
                    day = str(22-int(m.group(0))) # d0 becomes day 22
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
            day = str(22-int(A14D.group(3)))    # again there were 21 days starting at day 0
            
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
            day = str(22-int(no_.group(3)))    # again there were 21 days starting at day 0
            
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
                day = str(15+int(day[-1]))         # 21 days starting at day 0


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
## Replace 'bad' reads with fastq from rerun samples
file.close

file = open(sys.argv[1],'r')
repeat = re.compile('^b-(.*fastq)')

r=1
replacements = []  # keep track of replacements for the sake of debugging

for line in file:
    line = line.strip()
    fastq=line 
    rep = repeat.match(fastq)
    if rep:
        searchfastq = rep.group(1)  # no b-
        searchname = searchfastq.split('_')[0]   # Everything before the _ in the fastq
        # some samples had slightly different names upon resequencing them
        if searchname == '268-B-NT-D0':
            replacements.append('b-'+searchname)
            searchname = '268B-NT-D0'
        if searchname == '518B-18-D-5':
            replacements.append('b-'+searchname)
            searchname = '581B-18-D-5'        
        if searchname == '286A-1-D-11':
            replacements.append('b-'+searchname)
            searchname = '268A-1-D-11'
        if searchname == '581B-2-d1':
            replacements.append('b-'+searchname)
            searchname = '581B-2-D1'
    
        if r%2!=0: #odd
    	    for l in range(len(r1)):         #r1 is from the odd lines and so should only include R1 reads so searching for name is fine
                if searchname +'_' in r1[l]:
                    r1[l] = fastq
                    replacename = 'b-'+searchname
                    replacements.append(replacename)

            
        if r%2==0: #even
           
            for l in range(len(r2)):         #r1 is from the odd lines and so should only include R1 reads so searching for name is fine
                if searchname +'_' in r2[l]:
                   r2[l] = fastq
                   replacename = 'b-'+searchname
                   replacements.append(replacename)
    r+=1    
print(len(r2))
print('You forgot these you fool!')
forgot = list(set(missing)-set(replacements))
print(len(forgot))
print(forgot)


for i in range(len(samples)):
    print(samples[i],r1[i],r2[i],sep='\t',end='\n',file=outfile)
    
file.close()
outfile.close()

