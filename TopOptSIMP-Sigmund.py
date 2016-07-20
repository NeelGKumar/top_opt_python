# -*- coding: utf-8 -*-
"""
Topology Optimization code of Sigmund in Python
Created on Fri Apr 08 11:33:51 2011

@author: gkumar
"""
import pylab as py
import numpy as ny
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
def top(nelx,nely,volfrac,penal,rmin,problem):
    # INITIALIZE
    x = py.zeros((nely,nelx))
    x[0:nely,0:nelx] = volfrac 
    loop = 0
    change = 1.0
    dc = py.zeros((nely,nelx))
    j = 0
    #START ITERATION
    while change > 0.01  :
      j = j+1
      print j
      loop = loop + 1
      xold = x
    # FE-ANALYSIS
      U=FE(nelx,nely,x,penal,problem);         
    # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      KE = lk()
      c = 0.0
      for ely in range(nely):
          for elx in range(nelx):
              n1 = (nely+1)*elx+ely; 
              n2 = (nely+1)*(elx+1)   +ely;
              Ue = U[[2*n1,2*n1+1, 2*n2,2*n2+1, 2*n2+2,2*n2+3, 2*n1+2,2*n1+3],0]          
              c = c + (x[ely,elx]**penal) * py.inner(Ue,py.inner(KE,Ue))
              dc[ely,elx] = -penal*x[ely,elx]**(penal-1)*py.inner(Ue,py.inner(KE,Ue))
            
      # FILTERING OF SENSITIVITIES
      dc = check(nelx,nely,rmin,x,dc);    
      # DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
      x = OC(nelx,nely,x,volfrac,dc);       
      # PRINT RESULTS
      change = py.absolute(x-xold).max();
      su= py.sum(x)/(nelx*nely) 
      print 'Loop :',loop,' Vol.: ',su,' Change: ', change     
    py.imshow(x,cmap = py.cm.gray)
    py.show()      
     

    '''     title(['v=' num2str(volfrac) '; J = ' num2str(c) '; Iterations = ' num2str(loop)]);'''

''' OPTIMALITY CRITERIA UPDATE '''
def OC(nelx,nely,x,volfrac,dc)  :
    l1 = 0; l2 = 100000; move = 0.2;
    while (l2-l1 > 1e-4):
      lmid = 0.5*(l2+l1)
      xnew = py.maximum(1e-3,py.maximum(x-move,py.minimum(1.0,py.minimum(x+move,x*py.sqrt(-dc/lmid)))));
      if sum(xnew) - volfrac*nelx*nely > 0:
        l1 = lmid
      else:
        l2 = lmid
    return xnew
    
''' MESH-INDEPENDENCY FILTER '''
def check(nelx,nely,rmin,x,dc):
    dcn=py.zeros((nely,nelx));
    for i in range(1,nelx+1):
      for j in range(1,nely+1):
          sumx=0.0 
          for k in range(py.maximum(i-py.floor(rmin),1),py.minimum(i+py.floor(rmin),nelx)+1):
              
              for l in range(py.maximum(j-py.floor(rmin),1),py.minimum(j+py.floor(rmin),nely)+1):
                fac = rmin-py.sqrt((i-k)**2+(j-l)**2)
                sumx = sumx+py.maximum(0,fac)
                dcn[j-1,i-1] = dcn[j-1,i-1] + py.maximum(0,fac)*x[l-1,k-1]*dc[l-1,k-1]
                 
          dcn[j-1,i-1] = dcn[j-1,i-1]/(x[j-1,i-1]*sumx)          
    
    return dcn   
 
''' FE-ANALYSIS '''
  
def FE(nelx,nely,x,penal,problem):
    KE = lk()   
    K = py.zeros((2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)))
    #K = lil_matrix((2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)))
    #F = lil_matrix((2*(nely+1)*(nelx+1),1))
    F = py.zeros((2*(nely+1)*(nelx+1),1))
    U = py.zeros((2*(nely+1)*(nelx+1),1))
    #U = lil_matrix((2*(nely+1)*(nelx+1),1))
        
    #print 'K', K, 'F', F, 'U', U
    for elx in range(nelx):
      for ely in range(nely):
        n1 = (nely+1)*elx+ely 
        n2 = (nely+1)* (elx+1)+ely
        #print 'n', n1, n2
        edof = [2*n1,2*n1+1, 2*n2,2*n2+1, 2*n2+2,2*n2+3, 2*n1+2,2*n1+3]
        #print 'edof', edof
        E =x[ely,elx]**penal
        #print E*KE, K+.01
        for i in range(8):
            for j in range(8):
                #print edof[i],edof[j]
                   K[edof[i],edof[j]]=K[edof[i],edof[j]]+E*KE[i,j]
        #K[edof,:][:,edof] = K[edof,:][:,edof] + E*KE
        #print 'KXX',K[edof,:][:,edof]
    #py.transpose(K)    
    #print 'Final', K
    ''' DEFINE LOADS AND SUPPORTS 
    Left fixed, right tip load'''
    
    if (problem == 1):
        forcedDof = [2*(nelx+1)*(nely+1)-nely-1] # y force
        fixeddofs = py.arange(2*(nely+1))# left edge
        #print forcedDof, fixeddofs
    elif (problem == 2):
        forcedDof = 2*(nelx+1)*(nely+1) -1
        fixeddofs = py.arange(2*(nely+1))
    elif (problem == 3):
        forcedDof = py.array([2*(nelx+1)*(nely+1)-nely-1, 5*20*2+20-1]) # y force
        fixeddofs = py.arange(2*(nely+1)) # left edge
            
    F[forcedDof,0] = -1.0
    alldofs     = py.arange(2*(nely+1)*(nelx+1))
    freedofs    = list(set(alldofs) - set(fixeddofs))
    #K = K.todense
    #print K[freedofs,:][:,freedofs]
    #print F[freedofs,:]
    U[freedofs,:] = py.solve(K[freedofs,:][:,freedofs], F[freedofs,:]);
      
    U[fixeddofs,:]= 0;
    #input('abc')
    #print 'U',py.shape(U), U
    return U
''' ELEMENT STIFFNESS MATRIX '''
def lk():
    E = 1.0 
    nu = 0.3
    k=py.array([ 1./2.-nu/6.,1./8.+nu/8., -1./4.-nu/12., -1./8.+3.*nu/8., -1./4.+nu/12., -1./8.-nu/8. ,nu/6. , 1./8.-3.*nu/8.])
    KE = E/(1-nu**2)*py.array([[ k[0], k[1], k[2], k[3], k[4], k[5] ,k[6], k[7]],                [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],[k[4], k[5], k[6], k[7,], k[0], k[1], k[2], k[3]], [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],[k[6], k[3], k[4], k[1], k[2],k[7], k[0], k[5]], [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
    return KE             
''' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
def findContourValueWithArea(field,desiredVolFrac):
      #% Surround the matrix by a very low region to get closed contours
      #% But we have to go thru hoops to get the corner contours to come
      #% out correctly. Change sign of field, etc, etc.
      [nely,nelx] = py.shape(field)
      for value in py.arange(0.4,0.7,0.01):
          vf = computeAreaInContour(field,value)/(nelx*nely);
          if (vf > desiredVolFrac):
              break;
      return value
'''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
def computeAreaInContour(field,value):
    field = -field
    value = -value;
    valMin = field.min()
    [M,N] = py.shape(field);
    bufferedField = valMin*py.ones((M+2,N+2));
    #print bufferedField
    bufferedField[1:-1,1:-1] = field;
    print 'field' ,bufferedField
    [c,d] = py.shape(bufferedField)
    [X,Y] = py.mgrid[0:c,0:d]
    cs = plt.contour(X,Y,bufferedField,[value]);    
    plt.show()
    p =cs.collections[0].get_paths()
    s = py.size(p)
    a = 0.0;
    for i in range(s):
         p0 = cs.collections[0].get_paths()[i].vertices
         a = a + area(p0)
    
    #print 'p', p, 'p0', p0, 'size' , py.size(cs.collections[0].get_paths())
    #ar = area(p0);
    #print 'area', ar  
    return a;
def area(p):
    return 0.5 * abs(sum(x0*y1 - x1*y0  for ((x0, y0), (x1, y1)) in segments(p)))

def segments(p):
    return zip(p, p[1:] + [p[0]])

nelx  = 30;nely = 15; volfrac = 0.5; 
penal = 3;
rmin = 1.5;
problem = 1;
top(nelx,nely,volfrac,penal,rmin,problem)