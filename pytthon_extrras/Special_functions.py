# packages import
import matplotlib
import numpy as np
import csv
import os
import python_extras.Generall_functions as g_f



''' The function takes the results of minimizing the input data and gives the same structure, but
    the square roots of the radii are converted to regular radii. 
x - input vector
delta_R - addition that was used for avoiding nulls
point_type - type of point_ class '''
    
def recalc_r(x, delta_R, point_type):

  if point_type == 'point_tub_bi':
    raw_n_e, raw_n_i = x[::5], x[1::5]
    raw_sqrt_R_e, raw_sqrt_R_i, raw_sqrt_R_m = x[2::5], x[3::5], x[4::5]
    raw_R_e = raw_sqrt_R_e**2 + delta_R
    raw_R_i = raw_sqrt_R_i**2 + delta_R
    raw_R_m = raw_sqrt_R_m**2 + delta_R
    x_re = (np.array([raw_n_e, raw_n_i, raw_R_e, raw_R_i, raw_R_m]).transpose()).ravel()

  elif point_type == 'point_tub_mon':
    raw_n = x[::2]
    raw_sqrt_R = x[1::2]
    raw_R = raw_sqrt_R**2 + delta_R
    x_re = (np.array([raw_n, raw_R]).transpose()).ravel()

  return x_re



''' Function reauired for saving graphs.
number - either number or False or None. In the first condition save graph, in second
doing nothing, in third only plot '''

def save_graph(number):
  
  # showing graph without saving
  if number == None:
    matplotlib.pyplot.show()

  # required for adding elements without plotting or saving.
  elif number == False: 
    pass

  # save graph
  else:
    i = 1

    if f'graph_{number}.svg' in os.listdir(os.getcwd()):

      while f'graph_{number}_{i}.svg' in os.listdir(os.getcwd()):
        i += 1
      matplotlib.pyplot.savefig(f'graph_{number}_{i}.svg')
    
    else:
      matplotlib.pyplot.savefig(f'graph_{number}.svg')



''' Function for plotting graphs. Get a list consisting of sequnce formed by elemets in the range
    from 1 to 4 elements in string format, corresponding to the parameter designations. In case of
    a TypeError, pay attention to which elements are contained in the list.
data - list formed by elemets of point_ class
printed_labels - attributes of the point_ class that need to be plotted on the chart
number - number under which file in the graph_{number}.svg format will be saved
x_label - x-axis label
reverse - responsible for swapping axles
y_label - responsible for replacing the label of the y-axis '''

def graph ( data, printed_labels, number = None, x_label = 'z, nm', y_label = None, reverse = False): 

  number_funs = len( printed_labels )
  Y = []

  # determination axis labels
  if y_label != None:
    ordin_1 = y_label
    ordin_m = y_label

  else:
    ordin_1 = str(printed_labels[0])+ ', nm'
    ordin_m = 'Shape, nm'

  # ordinate axis values generation
  for fun in printed_labels:
    point_fun_txt = f'point.{fun}'
    y = [eval(point_fun_txt) for point in data]
    Y.append(y)

  # abscissa axis
  X = [ point.z for point in data]

  # definition of axes when building the 1 dependence
  if number_funs == 1:

    if reverse == True:
      matplotlib.pyplot.plot(Y[0], X)
      matplotlib.pyplot.xlabel(ordin_1)
      matplotlib.pyplot.ylabel(x_label)

    else:
      matplotlib.pyplot.plot(X, Y[0])
      matplotlib.pyplot.xlabel(x_label)
      matplotlib.pyplot.ylabel(ordin_1)

  # definition of axes when building the several dependence
  else:
    if reverse == True:

      for i in range(number_funs):
        matplotlib.pyplot.plot(Y[i], X)     
        matplotlib.pyplot.xlabel(ordin_m)
        matplotlib.pyplot.ylabel(x_label)

    else:
      for i in range(number_funs):
        matplotlib.pyplot.plot(X, Y[i])     
        matplotlib.pyplot.xlabel(x_label)
        matplotlib.pyplot.ylabel(ordin_m)
        
  # grid plotting
  matplotlib.pyplot.grid()

  # graph saving
  save_graph(number)



''' Energy density graph plotting.
ener_den - list or tuple with energy density values for each point
number - number under which the graph will be saved, if None, then it will only be plotted
axis_labes - list with both axes labels
reverse -responsible for the arrangement of the axles between themselves '''

def en_den_graph( ener_den, axis, number = None, axis_lables = ('r, nm','Energy density'),\
  reverse = None ):

  # checking for Multiple Curves in a Graph
  if type( ener_den[0] ) == list or type( ener_den[0] ) == tuple:

    for i in range( len(ener_den) ):

      y = []

      # traversal of energy densities and their recalculation
      for j in range(len (ener_den[i])):
        y.append(ener_den[i][j]/2/axis[i][j])
      
      x = axis[i][:-1]

      if reverse == None:
        matplotlib.pyplot.plot(x, y)
        matplotlib.pyplot.xlabel(axis_lables[0])
        matplotlib.pyplot.ylabel(axis_lables[1])

      elif reverse == True:
        matplotlib.pyplot.plot(y, x)
        matplotlib.pyplot.xlabel(axis_lables[1])
        matplotlib.pyplot.ylabel(axis_lables[0])

  else:
    
    y = []

    for j in range(len (ener_den)):
      y.append(ener_den[j]/2/axis[j])
      
    x = axis[:-1]
    
    if reverse == None:
      matplotlib.pyplot.plot(x, y)
      matplotlib.pyplot.xlabel(axis_lables[0])
      matplotlib.pyplot.ylabel(axis_lables[1])

    elif reverse == True:
      matplotlib.pyplot.plot(y, x)
      matplotlib.pyplot.xlabel(axis_lables[1])
      matplotlib.pyplot.ylabel(axis_lables[0])

  matplotlib.pyplot.grid()

  # graph saving
  save_graph(number)



''' Saving results in csv file.
point_type - type of point_ class 
grid - source list with data
number - csv file number
energy_res - calculated energy
left_values и right_values - initial values on the left and right, respectively
point_labels - point_ class parameters index
left_labels и right_labels - designation of the left and right boundary values
energy_density - energy in point '''

def save_grid ( point_type, grid, number, energy_res = None, left_values = None,\
  right_values = None, energy_density = None):

  # class parameters determination
  point_labels = g_f.define_labels(point_type)
  left_labels, right_labels = g_f.define_boundary_labels(point_type)

  if left_values != None:
    str_left_values = g_f.convert_to_str(left_values)

  if right_values != None:  
    str_right_values = g_f.convert_to_str(right_values)

  # create file header
  if left_values != None or right_values != None:
    header = [ ['', '', '', 'Results', '', '', '', '', '', '', '', '', 'Init values'],\
      ['numb'] ]
    header[1].extend(point_labels)

    if energy_density != None:
      header[1].append('Energy')

    if energy_res != None:
      header[1].append('Full_energy')

    if left_values != None:
      header[1].extend(left_labels)

    if right_values != None:
      header[1].extend(right_labels)
  else:
    header = ['numb']
    header.extend(point_labels)

    if energy_density != None:
      header.append('Energy')

    if energy_res != None:
      header.append('Full_energy')
      
  data = [] # data table row

  # сreating a nested list from string data that is subsequently stored
  for i in range(len(grid)):
    empty_list = [str(i)] # empty data list 

    # adding in empty data list point_ class value in string format
    for var_lab in point_labels:
      dot_value = eval( 'grid[i].' + var_lab )  # calling the values of the point_ class
      empty_list.append( str(dot_value) ) 

    # Добавление в строку со значениями значений энергии в точке
    if energy_density != None and i < len(grid) - 1:
      empty_list.append( str(energy_density[i]))

    # adding in list which will be used for results saving
    data.append(empty_list) 

    # сhecking for additional data to store it in the first row of values (excluding the header)
    if i == 0:

      if energy_res != None:
        data[0].append( str(energy_res) )

      if left_values != None:
        data[0] += str_left_values

      if right_values != None:
        data[0] += str_right_values

    i += 1

  # select the name under which the results will be saved
  if number == None:
    file_name = 'grid_.csv'

  else: 
    file_name = f'grid_{number}.csv'

  # сreating a file and saving the results in it
  with open(file_name, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    if left_values != None or right_values != None:
      writer.writerows(header)

    else:
      writer.writerow(header)

    writer.writerows(data)



''' Reads the results from a csv file, resulting in either a list of elements of the point
    class, or of the energy of the left and right edge values, as well as a list of elements
    of the point class.
point_type - type of point_ class
file_name - the name of the file from which the data will be loaded '''

def load_grid ( point_type, file_name):

  # class parameters determination
  point_labels = g_f.define_labels(point_type)
  left_labels, right_labels = g_f.define_boundary_labels(point_type)

  left_values = []
  right_values = []
  energy_res = None
  energy_density = []

  # opening a file for reading and transferring the data in the content list
  with open(file_name, newline='') as csvfile:
    array = csv.reader(csvfile, delimiter=';' )
    content=[]
    for row in array:
      if len( row ) > 0:
        content.append( row[0].split(',') )

  # determination last row with string elements.
  j = g_f.define_start( content )

  # assigning the first line of the file to the header variable, which will be
  # responsible for assigning values to the point class
  header = content[j]
  grid = list()
  major_part =  content[j+1:]

  # Traversing rows of a file
  for n in range( len( major_part ) ):
    dot = str()

    # Traversing elements of rows (columns)
    for i in range( len( major_part[n] ) ):

      # Checking for the presence of boundary conditions, or energy values
      if n == 0:
        if header[i] in left_labels:

          if major_part[n][i] == 'None':
            left_values.append( None )

          else:
            left_values.append( float( major_part[n][i] ) )
        
        elif header[i] in right_labels:

          if major_part[n][i] == 'None':
            right_values.append( None )

          else:
            right_values.append( float( major_part[n][i] ) )

        elif header[i] == 'Full_energy':
          energy_res = float(major_part[n][i])

      # Creating a string with results for a single element of point_ class
      if header[i] in point_labels:
        dot += f'{header[i]} = {major_part[n][i]},' # for instance, n_u = 64

      # adding elemets in energy density list
      elif header[i] == 'Energy':
        energy_density.append(float(major_part[n][i]))

    # creation complete point of point class with with followed adding to output list
    dot = dot[:-1]
    exec_code = f'g_f.{point_type}({dot})' # point( r = 2, ...
    grid.append( eval( exec_code )  )

  # checking for presence boundary conditions or allocated energy
  if len(left_values) == 0 and len(right_values) == 0 and energy_res == None and\
    len(energy_density) == 0:
    return grid

  # formation of a list with parameters to be output
  else:
    output = [grid]

    if energy_density != None:
      output.append(energy_density)

    if energy_res != None:
      output.append(energy_res)
    
    if len(left_values) != 0:
      output.append(left_values)

    if len(right_values) != 0:
      output.append(right_values)

    return output