import random

param = open('param.inp')
out = open('random_vec','w')

for line in param:
  if line.startswith('Ngrid'):
    Ngrid = int(param.next())
  else: continue

for i in range(Ngrid):
  print >> out, random.random() 
