#!/usr/bin/python

import sys
import re

#file = open(sys.argv[1],'r')
#outfile = open(sys.argv[2],'w')
file = open('fastq.txt','r')
outfile = open('colonus.files','w')

# ------------naming skemes for innocula--------------------------------

alyx_contam_inocA = re.compile('(\d+)A-inoculum') # 268A-inoculum_S64_L001_R1_001.fastq
niel_inoc = re.compile('([CH]\D{5,6}\d+)-\d{4}')  # Cancer16-2619_S76_L001_R1_001.fastq or Healthy....
alyx_contam_inocB = re.compile('DA(\d+)')     # DA6_S380_L001_R1_001.fastq
matt_conv2 = re.compile('ConvIDneg(\d+)')    # ConvIDneg20_S137_L001_R1_001.fastq
matt_syn2_inoc = re.compile('Mix(\d)IDneg(\d+)')  # Mix1IDneg20_S138_L001_R1_001.fastq

# -------------naming skemes for scientists and experiments --------------------
matt_soil = re.compile('(\d{4})-.*-(\d+)')           #1810-D-10_S308_L001_R1_001.fastq
matt_soilNT= re.compile('(NT)-.*-(\d+)')              #N-T-dn-10_S282_L001_R1_001.fastq  The N-T is adjusted below to 'NT' 
matt_syn1 = re.compile('([A-C]\d)-(\d{2})-(\d+)')      #C4-11-5_S140_L001_R1_001.fastq
matt_syn2 = re.compile('([A-C])(\d{3})D.*?(\d{1,2})')     #C901D0_S238_L001_R1_001.fastq
niel_2X = re.compile('([a-b]\d)-(.*)-.*')               # b4-nt-d00_S155_L001_R1_001.fastq
niel_4XJT = re.compile('([1-4][A-D])-(.*)-(\d+)')         #1A-nt-5_S329_L001_R1_001.fastq
niel_no_= re.compile('([a-d])([^-]*)dn?(\d+)')            #Niel's samples that are a622d0_S219_L001_R1_001.fastq
alyx_contam = re.compile('(\d+[AB])-(.*)-(D[-\d]\d*)')      #268B-NT-D-7_S96_L001_R1_001.fastq
alyx_treat = re.compile('([NYC]\D*)-(\w+)-*D-?(\d+)')       # NCHL-106-D0_S243_L001_R1_001.fastq Also Careyblair treatments


#----------control----------------
control_matt = re.compile('C(\d)-5-(\d+)')




# -------------Set up lists---------------------------
r1 = []   # fastq1
r2 = []   #fastq2
samples = []  #Sample name
missing =[]    # those that miss the boat on the skemes above
r = 1


#------------Loop through the fastq.txt file
for line in file:
    line = line.strip()           #remove the whitespace
    if r%2!=0:                    #odd
        fastq1=line               # Save the *.fastq name for the second colum of the files file.
        line = line.split('_')     # splits fastq file name at the sample name and saves this as name
        name1 = line[0]            # Save the original name to be matched to the other fastq later
        name = name1.replace('N-T','NT') # replace the N-T in matt's with NT for mouse
        
        if name == '846dn03':  # corrects a typo which labeled one of niels samples incorrectly
            name = 'd846dn03'
        if name == 'C03Dneg18':   # corrects a typo which labeled one of matt's syn2 samples incorrectly
            name = 'C903Dneg18'    
# -------------Determining how to rename the samples-------------------------------------------------

## INOCULA ###
        alyx_con_inocA = alyx_contam_inocA.match(name)
        n_inoc = niel_inoc.match(name)   
        alyx_con_inocB = alyx_contam_inocB.match(name)
        syn2_conv = matt_conv2.match(name)
        syn2_mix = matt_syn2_inoc.match(name)
##Control##
        control = control_matt.match(name)
##Samples ##       
        soil = matt_soil.match(name)
        soilNT = matt_soilNT.match(name)
        syn1 = matt_syn1.match(name)
        syn2 = matt_syn2.match(name)
        a1b4=niel_2X.match(name)
        A14D=niel_4XJT.match(name)
        no_ = niel_no_.match(name) 
        contam = alyx_contam.match(name)
        treat = alyx_treat.match(name)
        
## INOCULA ###        
        if alyx_con_inocA:
            scientist='alyx'
            group = alyx_con_inocA.group(1)  # number
            cage = 'A'
            mouse= 'INOC'
            day = mouse
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1     # save original name to be matched with other fastq later
        
        elif alyx_con_inocB:
            scientist='alyx'
            group = alyx_con_inocB.group(1) # number 
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
            	group = 'WRONG'    # for decoding
            cage = 'INOC'    
            mouse = cage
            day = cage    
                
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1
        
        elif syn2_conv:   # conventialized mice for control in syn2 experiment 
            scientist='matt'
            group = 'syn2'
            cage = 'Conv'
            mouse = 'INOC'
            day = str(20-int(syn2_conv.group(1))) # Dneg 19 = Day 1 There were three days of inocculation Dneg 22-20
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1
        elif syn2_mix:
            scientist ='matt'
            group = 'syn2'
            if syn2_mix.group(1)=='1':
                cage = 'A'
            elif syn2_mix.group(1)=='2':
                cage = 'B'
            else:
                cage = 'WRONG'
            
            mouse = 'INOC'
            day = str(20-int(syn2_mix.group(2))) # Dneg 19 = Day 1 There were three days of inocculation Dneg 22-20
            
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
            group = 'syn1'
            cage = 'syn'
            mouse = 'INOC'
            day = mouse
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
            samples.append(newname)
            r1.append(fastq1)
            test = name1 
       
        elif name == 'Feces-inoc-10-24':
             scientist = 'matt'
             group = 'syn1'
             cage = 'Conv'
             mouse = 'INOC'
             day = mouse
            
            
             newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            
             samples.append(newname)
             r1.append(fastq1)
             test = name1      
##Control##
        elif control:
            scientist = 'matt'
            group = 'Control'
            if int(control.group(1)) <3:  # mice 1 and 2 were cohoused. 3-5 were also cohoused 
                cage = 'A'
            else:
                cage = 'B'
            mouse = control.group(1)
            day = control.group(2)  # Experiment began on May 1 so we can just take the day

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
        elif syn1:
            scientist='matt'
            group = 'syn1'
            cage=syn1.group(1)[0]   # All in format 'A2' were A is the cage and 2 is the mouse
            mouse = syn1.group(1)[1]
            if int(syn1.group(2)) == 10:         # The days were in date format with 10-23 as day -1 with 10-25 as day 1 
                    day = str(int(syn1.group(3)) - 24)
            if int(syn1.group(2)) == 11:
                    day = str(int(syn1.group(3)) + 6) # 11-1 was really day 7
            
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1
            
        elif syn2:
            #print('Matched a syn2!')
            scientist = 'matt'
            group = 'syn2'
            cage= syn2.group(1)
            mouse = syn2.group(2)
            day = str(20-int(syn2.group(3))) # there were 3 days of Inoculation they are -2,-1, and 0 day 1 corresponds to Dneg19    
            
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
            day = str(22-int(A14D.group(3)))    # again there were 22 days starting at day 0
            
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
            day = str(22-int(no_.group(3)))    # again there were 22 days starting at day 0
            
            newname = scientist + '_' + group + '_' + cage + '_' + mouse + '_' + day
            samples.append(newname)
            r1.append(fastq1)
            test=name1

        #Alyx contaminated
        elif contam:
            scientist = 'alyx'
            group = contam.group(1)[0:-1] # take cage number but save the letter for the cage cattegory in the next line
            cage = contam.group(1)[-1]
            mouse = contam.group(2)
            day = contam.group(3)
            if '-' in day:    # Negative days were denoted D_#
                day = day.split('-')
                day = day[-1]
                day = str(15-int(day))   # 14 days starting at day 0
            elif 'D0' in day:
                day = '15'
            else:
                day = day.split('D') 
                day = str(15+int(day[-1]))         # these are for the few samples that extended beyond day 0 in the original experiment


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
   #------------If the name did not match the conditions above add it to the missing list-----------------------------------                   
        else:
            missing.append(name)    
  # ------the other fastq --------even lines-------------      
    if r%2==0:  #even
        fastq2=line 
        line = line.split('_')     # splits fastq file name at the sample name and saves this as name
        name2 = line[0]
        if name2==test:          # test to see if name matches that of the last successful fastq
            r2.append(fastq2)
    r += 1
#--------------Replace 'bad' reads with fastq from rerun samples-----------------------------------------------------
file.close

file = open('fastq.txt','r')

repeat = re.compile('^b-(.*fastq)')  # all the reruns have the form 'b-name_....fastq'

r=1
r1_replacements = []  # keep track of replacements for the sake of debugging
r2_replacements=[]
for line in file:
    line = line.strip()
    fastq=line 
    rep = repeat.match(fastq)        # go through the original fastq.txt and pull out a b- 
    if rep:
        searchfastq = rep.group(1)  # no b-
        searchname = searchfastq.split('_')[0]  # Everything before the _ in the fastq
           
        # some samples had slightly different names upon resequencing them mostly due to typos
        if searchname == '268-B-NT-D0':
            searchname = '268B-NT-D0'
        elif searchname == '518B-18-D-5':
            searchname = '581B-18-D-5'        
        elif searchname == '286A-1-D-11':
            searchname = '268A-1-D-11'
        elif searchname == '581B-2-d1':
            searchname = '581B-2-D1'
        elif searchname == 'b3-1357-19':
            #print('found it')
            searchname = 'b4-1357-19'
       
        elif 'NT' in searchname:
            searchname=searchname.replace('NT-dn','N-T-dn')   # just makes sure NT from matt will match the origianl sample names
    
    
    
    
    
        if r%2!=0: #odd
    	    for l in range(len(r1)):         #r1 is from the odd lines and so should only include R1 reads so searching for name is fine
                if searchname + '_' in r1[l]:  # find the orginal fastq and replace it
                    r1[l] = fastq

            
        if r%2==0: #even
           
            for l in range(len(r2)):         #r2 is from the even lines and so should only include R2 reads so searching for name is fine
                if searchname +'_' in r2[l]:
                   r2[l] = fastq
    r+=1    
    
    
##---------------Check for missed samples--------------------------------------------#

for j in range(len(r1)):         
    replacement = repeat.match(r1[j]) # see if the fastq matches a one that was inserted above during the replacement step
    if replacement:
        replacement = replacement.group(0).split('_')[0]
        r1_replacements.append(replacement)   #Add the name to the replacement list
        
for i in range(len(r2)):        
    replacement = repeat.match(r2[i])
    if replacement:
        replacement = replacement.group(0).split('_')[0]
        r2_replacements.append(replacement)  # Does the same as the loop above but this time for the second fastq column
        
# Make sure r1 and r2 are the same
print('There were replaced differentlty between r1 and r2')
               
different1 = set(r1_replacements).difference(set(r2_replacements))
different2 = set(r2_replacements).difference(set(r1_replacements))
different = set(different1).union(set(different2))
print(different)    
    
print('You forgot these you fool!')


replacements = r1_replacements # this is good so long as there is no difference between r1 and r2

forgot = set(missing).difference(set(replacements))   #Checkes that those that were missing (which are mostly the 'b-...' reruns are now included

print(len(forgot))
print(forgot)
#-----Write to .file file------------------------------------
for i in range(len(samples)):
    outfile.write(samples[i]+'\t'+r1[i]+'\t'+r2[i]+'\n')
    
file.close()
outfile.close()