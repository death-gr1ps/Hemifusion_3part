# packages import
import os
import numpy as np



''' Class for convenient work of bilayer circle data. '''

class point_bi:

  def __init__( self, z, n_u, n_d, H_u, H_d, M ):
    self.z = z
    self.n_u = n_u
    self.n_d = n_d
    self.H_u = H_u
    self.H_d = H_d
    self.M = M

  def __repr__(self):
    return "point_bi({0.z!r}, {0.n_u!r}, {0.n_d!r}, {0.H_u!r}, {0.H_d!r},\
      {0.M!r})".format(self)

  def __str__(self):
    return "({0.z!r}, {0.n_u!r}, {0.n_d!r}, {0.H_u!r}, {0.H_d!r}, {0.M!r})".format(self)



''' Class for convenient work of bilayer tube data. '''

class point_tub_bi:

  def __init__( self, z, n_e, n_i, R_e, R_i, R_m ):
    self.z = z
    self.n_e = n_e
    self.n_i = n_i
    self.R_e = R_e
    self.R_i = R_i
    self.R_m = R_m

  def __repr__(self):
    return "point_tub_bi({0.z!r}, {0.n_e!r}, {0.n_i!r}, {0.R_e!r}, {0.R_i!r},\
      {0.R_m!r})".format(self)

  def __str__(self):
    return "({0.z!r}, {0.n_e!r}, {0.n_i!r}, {0.R_e!r}, {0.R_i!r}, {0.R_m!r})".format(self)


''' Class for convenient work of monolayer tube data '''

class point_tub_mon:

  def __init__( self, z, n, R ):
    self.z = z
    self.n = n
    self.R = R

  def __repr__(self):
    return "point_tub_mon({0.z!r}, {0.n!r}, {0.R!r})".format(self)

  def __str__(self):
    return "{0.z!r}, {0.n!r}, {0.R!r})".format(self)



''' Function required to create a list with elements of point class.
point_type - name of class which used for making grid
result - массив результатов
z - conjugated axes vector '''

def mk_grid ( point_type, results, z ):

  # verification correctness choosing of point class
  assert_chek(point_type)

  output = list()

  if point_type == 'point_bi':
    N_points = len(results)//5

    for i in range(N_points):
      output.append( point_bi( z = z[i] , n_u = results[i*5], n_d = results[i*5 + 1],\
        H_u = results[i*5 + 2], H_d = results[i*5 + 3], M = results[i*5 + 4] ) )

  elif point_type == 'point_tub_bi':
    N_points = len(results)//5

    for i in range(N_points):
      output.append( point_tub_bi( z = z[i] , n_e = results[i*5], n_i = results[i*5 + 1],\
        R_e = results[i*5 + 2], R_i = results[i*5 + 3], R_m = results[i*5 + 4] ) )

  elif point_type == 'point_tub_mon':
    N_points = len(results)//2

    for i in range(N_points):
      output.append( point_tub_mon( z = z[i] , n = results[i*2], R = results[i*2 + 1] ) )

  return output



''' Function for converting list from elements of point_... class in vector of numpy package.
grid - source list
point_type - the class of points that form the grid '''

def grid_to_vek( point_type, grid ):
  
  row_list = [] # list formed by data that at the finish will be transformed in vector
  axis = []
  labels = define_labels(point_type) # define class attributes
  

  # traversing the original list starting from the second element
  for point in grid:

    # traversal by class attributes
    for j in labels:

      if j == 'z':
        axis.append(eval('point.' + j ))

      else:
        row_list.append(eval('point.' + j ))      

  # vector formation from list of data
  output = [ np.array(row_list), np.array(axis) ]
  
  return output



''' Function which purpose is creation of boundaries for SciPy.optimize.minimize().
sequence - np.array formed by maximum and minimum values for each point of created vector
number - elements number of ouput tuple'''

def mk_boundaries (sequence, number):

  empty_list = []

  for i in range(number):
    empty_list.extend(sequence)
    
  return tuple(empty_list)



''' Function for determination element number of list before numbers or list consisting 
from them. '''

def define_start (unknown_list):

  number = None
  n = 0

  for i in unknown_list :
  
    if type(i) == list or type(i) == tuple:
        
      for j in i:

        try:
          float(j)
          break

        except:
          number = n

    else:

      try:
        float(i)

      except:
        number = n

    n += 1

  return number



''' This function converts all numeric values in a list or tuple (may be nested once) to
 string format. '''

def convert_to_str (unknown_list):

  # required in order to avoid changes input data from list or for creation list from tuple
  data = list(unknown_list[:])

  for i in range( len(data) ):

    if type( data[i] ) == list: 

      for j in range( len(data[i]) ):

        if data[i][j] != None:        
          data[i][j] = str(data[i][j])

        else:
          data[i][j] = 'None'
    else:
      if data[i] != None:
        data[i] = str(data[i])

      else:
        data[i] = 'None'

  return data



''' This function determines file in directory with highest number and give this for saving
 another file. '''

def define_number ():

  comparator = 0
  int_set = {'0','1','2','3','4','5','6','7','8','9'}

  # traversing all files in a directory
  for file in os.listdir(os.getcwd()):
    j = len(file) - 1 
    l = None
    f = None
    counter = 0
    iterable = 0

    # traversing all symbols of file name in order to consideration of symbol _ number to not
    # consider second number in file name
    for i in file:

      if i == '_':
        counter += 1
    
    # the symbol from which the number is taken into account
    if counter == 2 and file[:5] != 'Total':
      symbol = '_' 

    # accounting for the presence of a numbered folder
    elif file[:5] == 'Total':
      symbol = True
      l = j + 1

    else:
      symbol = '.'

    while j > 1:

      # saving index value before last digit of the number
      if file[j] == symbol and symbol != True: 
        l = j 

      if l != None:

        # check that symbol is digit
        try:
          int( file[j - 1] ) 

        # saving index value of first digit number         
        except:
          f = j 
          break

      j -= 1

    # it is necessary both to avoid a null value and to avoid assigning a number in case of 
    # an empty directory
    if counter != 0:

    # slicing a vector with the same start and end produces an empty list
      if f == None or l == None:
        pass
      
      elif f == l and file[f] in int_set:
        iterable = int( file[f] )
    
      else:
        try:
          iterable = int( file[f:l] )

        except:
          print('problem with files names in data folder')

      if comparator < iterable:
        comparator = iterable

  return comparator + 1

    

''' Function required for determination class parameters '''

def define_labels(point_type):

  # checking rightness of class name
  assert_chek(point_type)

  if point_type == 'point_bi':
    point_labels = ('z', 'n_u', 'n_d', 'H_u', 'H_d', 'M')

  elif point_type == 'point_tub_bi':
    point_labels = ('z', 'n_e', 'n_i', 'R_e', 'R_i', 'R_m')

  elif point_type == 'point_tub_mon':
    point_labels = ('z', 'n', 'R')

  return point_labels



''' Function for labels determination of boundary condition for each class. '''

def define_boundary_labels(point_type):

  # checking rightness of class name
  assert_chek(point_type)

  if point_type == 'point_bi':
    left_labels = ('n_0u', 'n_0d', 'H_0u', 'H_0d', 'M_0')
    right_labels = ('n_1u', 'n_1d', 'H_1u', 'H_1d', 'M_1')

  elif point_type == 'point_tub_bi':
    left_labels = ('n_0e', 'n_0i', 'R_0e', 'R_0i', 'R_0m')
    right_labels = ('n_1e', 'n_1i', 'R_1e', 'R_1i', 'R_1m')

  elif point_type == 'point_tub_mon':
    left_labels = ('n_0', 'R_0')
    right_labels = ('n_1', 'R_1')

  return left_labels , right_labels

def assert_chek(point_type):
  assert point_type == 'point_bi' or point_type == 'point_tub_bi' or\
    point_type == 'point_tub_mon', ValueError