# this file reads in kepler fits headers to obtain CDPPs for each target in Dressing&Charbonneau.

from pyfits import open as pfopen
import os
import numpy as np
import pickle

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# read Kepler IDs
file = open('keplerIDs.txt','r')
keplerIDs = []
for line in file:
   keplerIDs = np.append(keplerIDs,line[:-1])
file.close()

# read files
i = 0
path = 'kepler/'
for (path, dir, file) in os.walk(path):
   i += 1
   if i == 1:
      break
files=set(file) 

# open fits files corresponding to each kepler ID and retrieve the CDPP values

len = len(keplerIDs)

CDPP = np.zeros([len,3,3])
LIVETIME = np.zeros(len)


i=0
for s in keplerIDs:
   if i % 100 == 0:
      print i,s

   matching = [ss for ss in files if s in ss]
   matching = np.array(matching)
   
   if matching.size == 0:
      print 'no files found for keplerID',s
   else:
      CDPP_temp = np.zeros([10,3])
      livet = 0e0
      ii = 0
      for j in matching:
         file = pfopen('kepler/'+j)
         
         if (file[0].header['QUARTER'] <= 6) and (file[0].header['QUARTER'] > 0):            
            
            livet = livet + file[1].header['LIVETIME']
            
            if str(file[1].header['CDPP3_0'])[0] != '<':
               CDPP_temp[ii,0] = file[1].header['CDPP3_0'] 
               CDPP_temp[ii,1] = file[1].header['CDPP6_0']
               CDPP_temp[ii,2] = file[1].header['CDPP12_0']
            else:
               CDPP_temp[ii,:] = 0e0
            
            file.close()
            
            ii = ii + 1
      
      a = CDPP_temp[:,0]
      b = CDPP_temp[:,1]
      c = CDPP_temp[:,2]
      
      if np.any(a > 0e0):
         CDPP[i,0,:] = [np.min(a[a > 0e0]),np.mean(CDPP_temp[:,0]),np.max(CDPP_temp[:,0])]
         CDPP[i,1,:] = [np.min(b[b > 0e0]),np.mean(CDPP_temp[:,1]),np.max(CDPP_temp[:,1])]
         CDPP[i,2,:] = [np.min(c[c > 0e0]),np.mean(CDPP_temp[:,2]),np.max(CDPP_temp[:,2])]
      
      LIVETIME[i] = livet
      
   i = i + 1


to_save = [keplerIDs,LIVETIME,CDPP]
save = open('CDPP.dat','w')
pickle.dump(to_save,save)
save.close()