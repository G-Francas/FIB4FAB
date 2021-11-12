# -*- coding: utf-8 -*-
"""profile.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/14rUtsFFQaJ40HthS20E-PAJuxnlPJps3
"""

import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt

import math
from IPython.display import Image as Mirage
from google.colab import drive
drive.mount('/content/drive') #This will require you enter passcode to access your drive

import sys

#Path to import extra modules from may need changing to match the path within your drive
path='/content/drive/My Drive/FIB 4 FAB/RCWA code/rcwa/examples/csvs'

def makeProf(theta, depthF, nLayers,res, material, gType):
  nm= 10**(-9)
  um=10**(-6)
  n=0
  layers=1000/(600*nm)*depthF*400*nm
  thickness=int(layers/nLayers)
  base=int(1000/6)
  profA=np.zeros((nLayers*thickness+base,res))

  place=0
  theta=90-theta
  while n < nLayers:
    if theta < 90:
        arr=np.genfromtxt(f"{path}/{gType}_{material}_slope{theta}_{depthF}_layer{n}_{res}.csv", delimiter=',')
    else:
        arr=np.genfromtxt(f"{path}/{gType}_{material}_slope{theta}_layer{n}_{res}.csv", delimiter=',')
    
    count=0
    while count < thickness:
        profA[place+count]=arr[21,:]
        count=count+1

    n=n+1
    place=n*count
  #profA=(profA/10).astype(np.uint8)
  i=0
  j=0
  x,y=np.shape(profA)
  while i < x:
      j=0
      while j <y:
          if math.isclose(profA[i][j],1,abs_tol=0.01):
              profA[i][j]=1
          else:
              profA[i][j]=0
          j=j+1
      i=i+1
  profA=profA*255
  im = Image.fromarray(profA)
  im=im.convert('L')
#   if im.mode != 'RGB':
#     im = im.convert('RGB')
  im.save(f"/content/drive/My Drive/FIB 4 FAB/RCWA code/rcwa/examples/profiles/sloped{theta}_depth{depthF}_res{res}.jpg")
  
  return im