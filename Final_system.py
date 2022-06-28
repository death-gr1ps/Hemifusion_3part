# Module for finding the shape of the lipid bilayer with the initial direction of the director

# packages import
import matplotlib
import numpy as np
import scipy
import os
import python_extras.Generall_functions as g_f
import python_extras.Special_functions as s_f
import python_extras.Energy_functions as e_f



''' Definition of a slice along the x axis. Returns a slice of type [A:B].
x - slice end
z - two bonded axes
numb_list - a list of point numbers separating different task areas '''

def x_ax (x, z, numb_list):

  for i in numb_list:

    if i == x:
      ind = numb_list.index(i)
        
  if ind == 0:
    output = z[ : x]
    
  else:
    output = z[numb_list[ind-1] : x]

  return output



''' Definition of a slice along the y-axis. Returns a slice of type [A:B].
y - slice end
z - two bonded axes
zero - number from which start new axis in bounded axes
numb_list - a list of point numbers separating different task areas '''

def y_ax (y, z, zero, numb_list):

  for i in numb_list:

    if i == y:
      ind = numb_list.index(i)

  if ind == 0:
    output = z[zero : y]
    
  else:
    output = z[numb_list[ind-1] : y]

  return output



''' Definition of the parameter vector for the x-axis. Returns a slice of type [A*N:B*N].
x - slice end
x0 - parameters vector
numb_list - a list consisting of point numbers separating different task areas '''

def x_vec (x, x0, numb_list):

  for i in numb_list:

    if i == x:
      ind = numb_list.index(i)
        
  if ind == 0:
    output = x0[ : x*5]
    
  else:
    output = x0[numb_list[ind-1]*5: x*5]

  return output



''' Determination of the parameter vector for the y-axis. Important: due to the difference
    in the length of the vector for a single iteration for a bilayer and monolayer tube,
    we have a more complex function. Returns a slice of type [A*N:B*N].
y - slice end
x0 - parameters vector
zero - number from which start new axis in bounded axes
numb_list - a list consisting of point numbers separating different task areas '''

def y_vec (y, x0, zero, vec_start, numb_list):

  # sum of vector elements
  sum = 0

  # cycle counter
  j = 0

  for i in numb_list:

    if i == y:
      ind = numb_list.index(i) 
      break

    else:
      
      if j == 0:
        sum += (i - zero)*2 
      else:
        if j % 2 == 0:
          sum += (i - numb_list[j-1]) *2

        elif j % 2 == 1:
          sum += (i - numb_list[j-1]) *5

    j += 1

  if j == 0:
    # y-c since indexes are defined before expanding from C to C1
    output = x0[vec_start*5 : vec_start*5 + (y - zero)*2] 

  else:

    if j % 2 == 0:
      output = x0[vec_start*5 + sum: vec_start*5 + sum + (y - numb_list[ind-1])*2]

    elif j % 2 == 1:
      output = x0[vec_start*5 + sum: vec_start*5 + sum + (y - numb_list[ind-1])*5]

  return output



''' Energy function of a composite system consisting of various tasks: three bilayer circles
    (0A, BC, CC1), two bilayer tubes (DE, FG) and three monolayer tubes (GH, 0D, EF). For a
    better understanding, you need to view the attached figure with segments of hemifusion
    structure. Important note # 1: there is no point C1 in the figure. The fact is that the
    areas of the vector BC and CC1 are summed over the same slice of the x-axis BC. Important
    note #2: The received energy is reduced by np.pi. The correct answer requires multiplication.
x0 - vector of parameters to be optimized with a predetermined length
z - values of the R and Z axes connected together
gen_param - all the constants needed to determine the energy. Important, the position of the upper
and lower layers is determined by the sections with a bilayer tube
gen_left_values - tuple of boundary conditions for the left and right sides of the bilayer circle
from 0 to A
gen_right_values - tuple of two tuples of boundary conditions for the upper and lower bilayer circle
from B to C
sep_numbers - a tuple of tuples with position numbers in the vector responsible for dividing task
areas in accordance with the above scheme '''

def HF_tub_en( x0, z, gen_param, gen_left_values, gen_right_values, sep_numbers, C1 = None,
  tot_en_den = None ):

  if type(tot_en_den) == list:

    if len(tot_en_den) != 8:
      tot_en_den = [ None, None, None, None, None, None, None, None]
  
  else:
    tot_en_den = [ None, None, None, None, None, None, None, None]

  # determination of parameters for calculations
  if type(gen_param[0]) != tuple:
    (sigma, K_t, K_a, K_g, B, h_0e, h_0i, J_Se, J_Si) = gen_param
    h2_0i = h1_0i = h_0i
    J2_Si = J1_Si = J_Si

  elif len(gen_param) == 4:
    (sigma, K_t, K_a, K_g, B) = gen_param[0]
    h_0e, J_Se = gen_param[1] # parameters for outer layer
    h1_0i, J1_Si = gen_param[2] # parameters for inner layer below E
    h2_0i, J2_Si = gen_param[3] # parameters for inner layer above E

  else:
    raise TypeError 

  # definition of constants for each layer
  hor_param_0 = (sigma, K_t, K_a, B, h2_0i, h1_0i, J2_Si, J1_Si) 
  hor_param_1 = (sigma, K_t, K_a, B, h_0e, h1_0i, J_Se, J1_Si)
  hor_param_2 = (sigma, K_t, K_a, B, h2_0i, h_0e, J2_Si, J_Se)
  tub_mon_param_0 = (sigma, K_t, K_g, B, J1_Si)
  tub_mon_param_1 = (sigma, K_t, K_g, B, J_Se)
  tub_mon_param_2 = (sigma, K_t, K_g, B, J2_Si)
  tub_bi_param_0 = (sigma, K_t, K_a, K_g, B, h_0e, h1_0i, J_Se, J1_Si)
  tub_bi_param_1 = (sigma, K_t, K_a, K_g, B, h_0e, h2_0i, J_Se, J2_Si)

  # determination of point numbers, counting from 1, dividing the axis into areas of different tasks
  A, B, C = sep_numbers[0]
  D, E, F, G, H = sep_numbers[1]

  if C1 == None:
    C1 = C*2 - B

  vec_delim = sep_numbers[0] + tuple([C1]) 
 
  # determining the type (left or right) of boundary conditions returned by functions
  x_values = [ ['left'] ]*2
  y_values = [ ['right'], ['left'], ['right'], ['left'] ]

  # determining expressions for director matching in horizontal and vertically oriented layers
  
  # matching conditions n_z for monolayer tube and bilayer circle from above
  dir_low_mon = lambda x: 1 - x
  # matching conditions n_z for bilayer tube and bilayer circle from above
  dir_low_bi =  lambda x: -1 -x
  # matching conditions n_z for monolayer tube and bilayer circle from below
  dir_up_mon = lambda x: x - 1
  # matching conditions n_z for bilayer tube and bilayer circle from below
  dir_up_bi = lambda x: 1 + x
  # matching conditions n_r for the upper layer of the bilayer circle on the left and the bilayer tube
  dir_up_left = lambda x: 1 + x
  # matching conditions n_r for the lower layer of the bilayer circle on the left and the bilayer tube
  dir_low_left = lambda x: 1 - x 

  # lowermost horizontal bilayer circle
  out_1 = e_f.horiz_bi( x_vec(C, x0, vec_delim), x_ax(C, z, sep_numbers[0]), hor_param_1,
    left_values = gen_right_values[0][0], right_values = gen_right_values[0][1],
      output_values = x_values[0], energy_density = tot_en_den[0] )

  # uppermost horizontal bilayer circle
  out_2 = e_f.horiz_bi( x_vec(C1, x0, vec_delim), x_ax(C, z, sep_numbers[0]), hor_param_2,
    left_values = gen_right_values[1][0], right_values = gen_right_values[1][1],
      output_values = x_values[1], energy_density = tot_en_den[1] )
  
  # monolayer tube from 0 to D
  input_values_3 = [ dir_low_mon(x_values[0][1]), z[B] ] # n, R
  
  out_3 = e_f.tub_in_mon( y_vec(D, x0, C, C1, sep_numbers[1]), y_ax(D, z, C, sep_numbers[1]),
    tub_mon_param_0, left_values = input_values_3, right_values = None,
      output_values = y_values[0],energy_density = tot_en_den[2] )
  
  # monolayer tube from G to H
  input_values_4 = [ dir_up_mon(x_values[1][0]), z[B] ] # n, R
  
  out_4 = e_f.tub_in_mon( y_vec(H, x0, C, C1, sep_numbers[1]), y_ax(H, z, C, sep_numbers[1]),
    tub_mon_param_2, left_values = None , right_values = input_values_4,
      output_values = y_values[1], energy_density = tot_en_den[3] )
  
  # bilayer tube from D to E
  input_values_5l = [ dir_low_bi(x_values[0][0]), y_values[0][0], z[B], y_values[0][1], None ]
  input_values_5r = [None, None, None, z[A], None] # n_e, n_i, R_e, R_i, R_m
  
  out_5 = e_f.tub_bi( y_vec(E, x0, C, C1, sep_numbers[1]), y_ax(E, z, C, sep_numbers[1]),
    tub_bi_param_0, left_values = input_values_5l, right_values = input_values_5r,
      output_values = y_values[2], energy_density = tot_en_den[4] )

  # bilayer tube from F to G
  input_values_6r = [ dir_up_bi(x_values[1][1]), y_values[1][0], z[B], y_values[1][1], None ] 
  input_values_6l = [None, None, None, z[A], None] # n_e, n_i, R_e, R_i, R_m

  out_6 = e_f.tub_bi( y_vec(G, x0, C, C1, sep_numbers[1]), y_ax(G, z, C, sep_numbers[1]),
    tub_bi_param_1, left_values = input_values_6l , right_values = input_values_6r,
      output_values = y_values[3], energy_density = tot_en_den[5] )
  
  # central monolayer tube from E to F
  input_values_7l = [ y_values[2][0], y_values[2][2] ] # n, R -- left
  input_values_7r = [ y_values[3][0], y_values[3][2] ] # n, R -- right

  out_7 = e_f.tub_ex_mon( y_vec(F, x0, C, C1, sep_numbers[1]), y_ax(F, z, C, sep_numbers[1]),
    tub_mon_param_1, left_values = input_values_7l, right_values = input_values_7r,
      output_values = None, energy_density = tot_en_den[6] )
  
  # leftmost horizontal bilayer circle
  input_values_8 = [dir_up_left(y_values[3][1]), dir_low_left(y_values[2][1]), z[F], z[E], None]

  out_8 = e_f.horiz_bi( x_vec(A, x0, vec_delim), x_ax(A, z, sep_numbers[0]), hor_param_0,
    left_values = gen_left_values, right_values = input_values_8, output_values = None,
      energy_density = tot_en_den[7] )
  
  return out_1 + out_2 + out_3 + out_4 + out_5 + out_6 + out_7 + out_8



''' A function that generates the parameter vector and the ordinate axis necessary to optimize the
    overall function.
sep_numbers - a tuple of tuples with position numbers in the vector responsible for dividing task
lengths - a tuple of tuples with region lengths for different task regions
width - tuple of layer thicknesses in a bilayer tube
    Importantly, the following parameters are entered for the problem when the H_m of the bilayer
    circle is at zero on the z-axis.
gen_left_values - tuple of boundary conditions for the left and right sides of the bilayer circle
from 0 to A
gen_right_values - tuple of two tuples of boundary conditions for the upper and lower bilayer circle
from B to C.'''

def input_parameters (sep_numbers, lengths, width, gen_right_values, C1 = None):

  # determination of layer thicknesses in a bilayer tube
  h_0e, h_0i = width

  # determination of the points numbers (starting from 1) dividing the axis into areas of
  # different tasks
  A, B, C = sep_numbers[0]
  D, E, F, G, H = sep_numbers[1]

  if C1 == None:
    C1 = C*2 - B

  # determination of the region lengths of various tasks
  A0, AB, BC = lengths[0]
  D0, DE, EF, FG, GH = lengths[1]

  # generation of boundary conditions for bilayer tubes B to C
  out_right_val = [ [[], []], [[], []] ]
  out_right_val[0][0] = [None, None, D0, 0, None]
  out_right_val[0][1] = gen_right_values[0]
  length = D0 + DE + EF + FG
  out_right_val[1][0] = [None, None, length + GH, length, None]
  out_right_val[1][1] = gen_right_values[1]

  # determination of the ordinate vector by which the optimization will be carried out
  x = [' ']*9
  z = np.linspace(0.1, A0//1.333, num = A//2)
  x[0] = np.linspace(A0//1.333 + z[-1] - z[-2], A0, num = A - A//2)
  x[1] = np.linspace(A0, A0+AB, num = B - A)
  x[2] = np.linspace(A0+AB, A0+AB+BC//4, num = (C - B)//2)
  x[3] = np.linspace(A0+AB+BC//4 + x[2][-1] - x[2][-2], A0+AB+BC, num = C - B - (C - B)//2)
  x[4] = np.linspace(0, D0, num = D - C)
  x[5] = np.linspace(D0, D0+DE, num = E - D)
  x[6] = np.linspace(D0+DE, D0+DE+EF, num = F - E)
  x[7] = np.linspace(D0+DE+EF, D0+DE+EF+FG, num = G - F)
  x[8] = np.linspace(D0+DE+EF+FG, D0+DE+EF+FG+GH, num = H - G)

  # conjunction of resulting sequences
  for i in x:
    z = np.append(z, i)

  # Determination of a vector of initial values. For the best work it is necessary to do close
  # to the possible solution. Note: under area AB, for this no element in the vector area of
  # values is allocated.

  x = [' ']*8

  # bilipid circle from 0 to A
  length = D0+DE+EF/2
  N = A
  x0 = ((np.array([np.linspace(0.0, 0.5, num = N), np.linspace(0.0, 0.5, num = N),
    np.linspace(h_0i + length, length + EF/2, num = N),
      np.linspace(length - h_0i, length - EF/2, num = N),
        np.linspace(length, length, num = N)])).transpose()).ravel() # n_u, n_d, R_u, R_d, R_m

  # a vector piece for an element from A to B where there is no task
  length = D0/2
  N = B - A
  x[0] = ((np.array([np.linspace(0.0, 0.0, num = N), np.linspace(0.0, 0.0, num = N),
    np.linspace(0.0, 0.0, num = N), np.linspace(0.0, 0.0, num = N),
      np.linspace(0.0, 0.0, num = N)])).transpose()).ravel() # n_u, n_d, R_u, R_d, R_m

  # lowermost horizontal bilayer
  length = D0/2
  N = C - B
  x[1] = ((np.array([np.linspace(-0.5, 0.0, num = N), np.linspace(0.2, 0.0, num = N),
    np.linspace(h_0e + length, h_0e + length, num = N),
      np.linspace(length - h_0i, length - h_0i, num = N),
        np.linspace(length, length, num = N)])).transpose()).ravel() # n_u, n_d, R_u, R_d, R_m

  # uppermost horizontal bilayer
  length = D0+DE+EF+FG+GH/2
  N = C1 - C
  x[2] = ((np.array([np.linspace(0.2, 0.0, num = N), np.linspace(-0.5, 0.0, num = N),
    np.linspace(h_0i + length, h_0i + length, num = N),
      np.linspace(length - h_0e, length - h_0e, num = N),
        np.linspace(length, length, num = N)])).transpose()).ravel() # n_u, n_d, R_u, R_d, R_m

  # monolayer tube from 0 to D
  length = A0+AB
  N = D - C # not from point C1 so we need difference between points
  x[3] = ((np.array([ np.linspace(0.8, 0.2, num = N),
    np.linspace(length, length - h_0i - h_0e, num = N) ])).transpose()).ravel() # n, R

  # bilayer tube from D to E
  length = A0+AB/2
  N = E - D
  x[4] = ((np.array([np.linspace(-0.5, 0.0, num = N), np.linspace(0.2, 0.5, num = N),
    np.linspace(length + h_0e, length + h_0e, num = N),
      np.linspace(length - h_0i, length - h_0i, num = N),
        np.linspace(length, length, num = N)])).transpose()).ravel() # n_e, n_i, R_e, R_i, R_m 

  # monolayer tube from E to F
  length = A0+AB/2 + h_0e
  N = F - E
  x[5] = ((np.array([ np.linspace(0.0, 0.0, num = N),
    np.linspace(length, length , num = N) ])).transpose()).ravel() # n, R

  # bilayer tube from F to G
  length = A0+AB/2
  N = G - F
  x[6] = ((np.array([np.linspace(0.0, 0.5, num = N), np.linspace(-0.5, -0.2, num = N),
    np.linspace(length + h_0e, length + h_0e, num = N),
      np.linspace(length - h_0i, length - h_0i, num = N),
        np.linspace(length, length, num = N)])).transpose()).ravel() # n_e, n_i, R_e, R_i, R_m 

  # monolayer tube from G to H
  length = A0+AB
  N = H - G
  x[7] = ((np.array([ np.linspace(-0.2, -0.8, num = N),
    np.linspace(length - h_0i - h_0e, length , num = N) ])).transpose()).ravel() # n, R

  # conjunction of results
  for i in x:
    x0 = np.append(x0, i)

  return z, x0, out_right_val



''' The function needed to translate the resulting vector into a list of lists with elements of
    the grid class for different areas of the task.
z - values of the R and Z axes connected together
x0 - parameters vector
sep_numbers - a tuple of tuples with position numbers in the vector responsible for dividing task '''

def total_grid (z, x0, sep_numbers):

  # results list
  out = [[]] *8

  # determination of the points numbers (starting from 1) dividing the axis into areas of
  # different tasks
  A, B, C = sep_numbers[0]
  D, E, F, G, H = sep_numbers[1]
  C1 = C*2 - B
  vec_delim = sep_numbers[0] + tuple([C1])

  # results for bilayer circle from 0 to A
  out[0] = g_f.mk_grid( 'point_bi', x_vec(A, x0, vec_delim), x_ax(A, z, sep_numbers[0]) )

  # results for downmost bilayer circle from B to C
  out[1] = g_f.mk_grid( 'point_bi', x_vec(C, x0, vec_delim), x_ax(C, z, sep_numbers[0]) )

  # results for uppermost bilayer circle from B to C
  out[2] = g_f.mk_grid( 'point_bi', x_vec(C1, x0, vec_delim), x_ax(C, z, sep_numbers[0]) )

  # results for monolayer tube from 0 to D
  out[3] = g_f.mk_grid( 'point_tub_mon', y_vec(D, x0, C, C1, sep_numbers[1]),
    y_ax(D, z, C, sep_numbers[1]) )

  # results for bilayer tube from D to E
  out[4] = g_f.mk_grid( 'point_tub_bi', y_vec(E, x0, C, C1, sep_numbers[1]),
    y_ax(E, z, C, sep_numbers[1]) )
  
  # results for monolayer tube from E to F
  out[5] = g_f.mk_grid( 'point_tub_mon', y_vec(F, x0, C, C1, sep_numbers[1]),
    y_ax(F, z, C, sep_numbers[1]) )

  # reults for monolayer tube from F to G
  out[6] = g_f.mk_grid( 'point_tub_bi', y_vec(G, x0, C, C1, sep_numbers[1]),
    y_ax(G, z, C, sep_numbers[1]) )

  # reults for monolayer tube from G to H
  out[7] = g_f.mk_grid( 'point_tub_mon', y_vec(H, x0, C, C1, sep_numbers[1]),
    y_ax(H, z, C, sep_numbers[1]) )

  return out



''' A function that builds a graph of the location of lipid layers for the task represented
    by the HF_tub_en function (separately for each task piece)
data - data list which was generated by total_grid
layer - string variable (all, upper, lower) or list of string responsible for selecting layer
or layers which will be plotted
number - the number under which the graph will be saved. If None, it will only be built '''

def total_graph (data, layer = 'all', number = None, x_label = 'r, nm', y_label = 'z, nm'):
  
  pr_r = [None, [None, None], [None, None, None, None, None]]
  pr_z = [None, [None, None], [None, None, None, None, None]]

  if len(data) == 10:
    pr_r[2].extend([None, None])
    pr_z[2].extend([None, None])

  elif len(data) == 9:
    pr_r[2].append(None)
    pr_z[2].append(None)

  # adding upper layer in the graph
  if 'all' in layer or 'upper' in layer:

    # replacing None with a list that will expand after
    pr_r[0] = [data[1][-1].z]
    pr_z[0] = [data[1][-1].H_u]
    j = len(data[1]) - 2

    while j >= 0:
      pr_r[0].append(data[1][j].z)
      pr_z[0].append(data[1][j].H_u)
      j -= 1

    for i in data[4]:
      pr_r[0].append(i.R_e)
      pr_z[0].append(i.z)

    for i in data[5]:
      pr_r[0].append(i.R)
      pr_z[0].append(i.z)

    for i in data[6]:
      pr_r[0].append(i.R_e)
      pr_z[0].append(i.z)

    for i in data[2]:
      pr_r[0].append(i.z)
      pr_z[0].append(i.H_d)

  # adding lower layer in the graph
  if 'all' in layer or 'lower' in layer:

    # replacing None with a list that will expand after
    pr_r[1][0] = [data[1][-1].z]
    pr_z[1][0] = [data[1][-1].H_d]
    j = len(data[1]) - 2

    while j >= 0:
      pr_r[1][0].append(data[1][j].z)
      pr_z[1][0].append(data[1][j].H_d)
      j -= 1

    for i in data[3]:
      pr_r[1][0].append(i.R)
      pr_z[1][0].append(i.z)

    for i in data[4]:
      pr_r[1][0].append(i.R_i)
      pr_z[1][0].append(i.z)

    # reverse list elements traversing 
    j = len(data[0]) - 1
    
    while j >= 0:
      pr_r[1][0].append(data[0][j].z)
      pr_z[1][0].append(data[0][j].H_d)
      j -= 1

    # replacing None with a list that will expand after
    pr_r[1][1] = [data[0][0].z]
    pr_z[1][1] = [data[0][0].H_u]
    j = len(data[0])
    i = 1

    while i < j:
      pr_r[1][1].append(data[0][i].z)
      pr_z[1][1].append(data[0][i].H_u)
      i += 1

    for i in data[6]:
      pr_r[1][1].append(i.R_i)
      pr_z[1][1].append(i.z)

    for i in data[7]:
      pr_r[1][1].append(i.R)
      pr_z[1][1].append(i.z)

    for i in data[2]:
      pr_r[1][1].append(i.z)
      pr_z[1][1].append(i.H_u)

  # adding lower layer in the graph (boundary between two lipid layers)
  if 'all' in layer or 'midle' in layer:

    # replacing None with a list that will expand after
    pr_r[2][0] = [data[0][0].z]
    pr_z[2][0] = [data[0][0].M]
    j = len(data[0])
    i = 1

    while i < j:
      pr_r[2][0].append(data[0][i].z)
      pr_z[2][0].append(data[0][i].M)
      i += 1

    # replacing None with a list that will expand after
    pr_r[2][1] = [data[1][0].z]
    pr_z[2][1] = [data[1][0].M]
    j = len(data[1])
    i = 1

    while i < j:
      pr_r[2][1].append(data[1][i].z)
      pr_z[2][1].append(data[1][i].M)
      i += 1

    # replacing None with a list that will expand after
    pr_r[2][2] = [data[2][0].z]
    pr_z[2][2] = [data[2][0].M]
    j = len(data[2])
    i = 1

    while i < j:
      pr_r[2][2].append(data[2][i].z)
      pr_z[2][2].append(data[2][i].M)
      i += 1

    # replacing None with a list that will expand after
    pr_r[2][3] = [data[4][0].R_m]
    pr_z[2][3] = [data[4][0].z]
    j = len(data[4])
    i = 1

    while i < j:
      pr_r[2][3].append(data[4][i].R_m)
      pr_z[2][3].append(data[4][i].z)
      i += 1

    # replacing None with a list that will expand after
    pr_r[2][4] = [data[6][0].R_m]
    pr_z[2][4] = [data[6][0].z]
    j = len(data[6])
    i = 1

    while i < j:
      pr_r[2][4].append(data[6][i].R_m)
      pr_z[2][4].append(data[6][i].z)
      i += 1

  # axes assembling for graph plotting
  i = 0

  while i < 3:

    if type(pr_r[i]) == list:
      j = 0

      while j < len(pr_r[i]):

        if type(pr_r[i][j]) == list:
          matplotlib.pyplot.plot(pr_r[i][j], pr_z[i][j])

        elif pr_r[i][j] != None:
          matplotlib.pyplot.plot(pr_r[i], pr_z[i])
          break
        j +=1

    i += 1
  matplotlib.pyplot.xlabel(x_label)
  matplotlib.pyplot.ylabel(y_label)
  matplotlib.pyplot.grid()

  # graph saving
  s_f.save_graph(number)



''' Function that builds a graph of the location of lipid layers for the task represented
    by the HF_tub_en function.
data - data list which was generated by total_grid
layer - string variable (all, upper, lower) or list of string responsible for selecting layer
or layers which will be plotted
number - the number under which the graph will be saved. If None, it will only be built '''

def sep_graph (data,  layer = 'all', number = None):

  if 'upper' in layer or 'all' in layer:
    s_f.graph( data[1], ['H_u'], False, 'r, nm', 'z, nm')
    s_f.graph( data[2], ['H_d'], False, 'r, nm', 'z, nm')
    s_f.graph( data[4], ['R_e'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[5], ['R'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[6], ['R_e'], False, 'z, nm', 'r, nm', True)

  if 'lower' in layer or 'all' in layer:
    s_f.graph( data[0], ['H_u', 'H_d'], False, 'r, nm', 'z, nm')
    s_f.graph( data[1], ['H_d'], False, 'r, nm', 'z, nm')
    s_f.graph( data[2], ['H_u'], False, 'r, nm', 'z, nm')
    s_f.graph( data[3], ['R'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[4], ['R_i'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[6], ['R_i'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[7], ['R'], False, 'z, nm', 'r, nm', True)

    if len(data) >= 9:
      s_f.graph( data[8], ['R'], False, 'z, nm', 'r, nm', True)

    if len(data) == 10:
      s_f.graph( data[9], ['R'], False, 'z, nm', 'r, nm', True)

  if 'midle' in layer or 'all' in layer:
    s_f.graph( data[0], ['M'], False, 'r, nm', 'z, nm')
    s_f.graph( data[1], ['M'], False, 'r, nm', 'z, nm')
    s_f.graph( data[2], ['M'], False, 'r, nm', 'z, nm')
    s_f.graph( data[4], ['R_m'], False, 'z, nm', 'r, nm', True)
    s_f.graph( data[6], ['R_m'], False, 'z, nm', 'r, nm', True)

  # grpah saving
  matplotlib.pyplot.grid()
  s_f.save_graph(number)



''' The function needed to save the results as a file folder for each task area
data - data list which was generated by total_grid
number - the number under which the folder will be saved. If None, it will be saved
without number.
energy_res - calculated energy of full system
gen_left_values - left boundary conditions for far left lipid circle 
gen_right_values - list of 2 lists of boundary conditions for far right bilipid circles
path - path where the file will be saved 
tot_en_den - empty list where energy density will be sorted'''

def save_total_grid (data, number = None, energy_res = None, gen_left_values = None,
  gen_right_values = [ [None, None], [None, None] ], path = None,
    tot_en_den = [None, None, None, None, None, None, None, None] ):

  # saving path to old directory
  old_dir = os.getcwd()

  # selecting the directory in which the folder will be created
  if path == None:
    path = old_dir

  # file name selecting
  if number == None:
    folder  = 'Total_grid'

  else:
    folder  = f'Total_grid_{number}'

  direct = os.path.join(path, folder)

  if not folder in os.listdir(os.getcwd()): 
    os.mkdir(direct)

  os.chdir(direct)

  # results saving in file of csv format
  s_f.save_grid('point_bi', data[0], 0, energy_res = energy_res,
    left_values = gen_left_values, energy_density = tot_en_den[7])
  s_f.save_grid('point_bi', data[1], 1, energy_res =energy_res, left_values = gen_right_values[0][0],
    right_values = gen_right_values[0][1], energy_density = tot_en_den[0])
  s_f.save_grid('point_bi', data[2], 2, energy_res = energy_res, left_values = gen_right_values[1][0],
    right_values = gen_right_values[1][1], energy_density = tot_en_den[1])
  s_f.save_grid('point_tub_mon', data[3], 3, energy_res =energy_res, energy_density = tot_en_den[2])
  s_f.save_grid('point_tub_bi', data[4], 4, energy_res = energy_res, energy_density = tot_en_den[4])
  s_f.save_grid('point_tub_mon', data[5], 5, energy_res = energy_res, energy_density = tot_en_den[6])
  s_f.save_grid('point_tub_bi', data[6], 6, energy_res = energy_res, energy_density = tot_en_den[5])
  s_f.save_grid('point_tub_mon', data[7], 7, energy_res = energy_res, energy_density = tot_en_den[3])

  # return to the old directory
  os.chdir(old_dir)



''' Function that reads the results obtained from the minimization of HF_tub_en function.
folder_name - name of folder that contains results'''

def load_total_grid (folder_name):

  # saving path to old directory
  old_dir = os.getcwd()

  # directory changing
  direct = os.path.join(old_dir, folder_name)
  os.chdir(direct)

  # reading results from file of csv format
  grid = [ s_f.load_grid ('point_bi', 'grid_0.csv') ]
  grid.append( s_f.load_grid ('point_bi', 'grid_1.csv') )
  grid.append( s_f.load_grid ('point_bi', 'grid_2.csv') )
  grid.append( s_f.load_grid ('point_tub_mon', 'grid_3.csv') )
  grid.append( s_f.load_grid ('point_tub_bi', 'grid_4.csv') )
  grid.append( s_f.load_grid ('point_tub_mon', 'grid_5.csv') )
  grid.append( s_f.load_grid ('point_tub_bi', 'grid_6.csv') )
  grid.append( s_f.load_grid ('point_tub_mon', 'grid_7.csv') )

  output = [ [], [], [], []]
  output[0] = [ grid[0][0] ]
  output[0].append(grid[1][0])
  output[0].append(grid[2][0])
  output[0].append(grid[3][0])
  output[0].append(grid[4][0])
  output[0].append(grid[5][0])
  output[0].append(grid[6][0])
  output[0].append(grid[7][0])

  # checking presence of energy_density
  if len(grid[7]) == 3:
    output[1] = [ grid[0][1] ]
    output[1].append(grid[1][1])
    output[1].append(grid[2][1])
    output[1].append(grid[3][1])
    output[1].append(grid[4][1])
    output[1].append(grid[5][1])
    output[1].append(grid[6][1])
    output[1].append(grid[7][1])
    output[2] = grid[0][2]
    output[3] = grid[0][3]
    output.append( [ [ grid[1][3], grid[1][4] ], [ grid[2][3], grid[2][4] ] ] )
  
  elif len(grid[7]) == 2:
    output[1] = grid[0][1]
    output[2] = grid[0][2]
    output[3] = [ [ grid[1][2], grid[1][3] ], [ grid[2][2], grid[2][3] ] ]

  # сhange directory back to original
  os.chdir(old_dir)

  return output



''' Function for determining the optimal lengths of task regions.
x - task regions lengths vector
sep_numbers - a tuple of tuples with position numbers in the vector responsible for dividing task
gen_left_values - left boundary conditions for the left bilayer circle
right_values - right boundary conditions for bilayer circles on the right
gen_param - all the constants needed to determine the energy. Important, the position of the upper
width - tuple of layer thicknesses in a bilayer tube '''

def define_lengths_as (x, sep_numbers, gen_left_values, right_values, gen_param, total_lengths,
  save_data = None):

  # determining the thickness of the layers required at the stage of generating the vector of
  # parameters and axes
  if len(gen_param) == 9:
    width = gen_param[5], gen_param[6]

  elif len(gen_param) == 4:
     width = gen_param[1][0], (gen_param[2][0] + gen_param[3][0])/2

  else:
    raise TypeError

  # determination of full lengths for two dimensions
  r_len, z_len = total_lengths

  # determination of thicknesses for different regions of the task
  A0, AB, D0, DE, EF, FG = x[0], x[1], x[2], x[3], x[4], x[5]
  BC = r_len - A0 - AB
  GH = z_len - D0 - DE - EF - FG
  lengths = ((A0, AB, BC), (D0, DE, EF, FG, GH))

  # generation of axes and parameter vectors depending on initial conditions
  z, x0, gen_right_values = input_parameters(sep_numbers, lengths, width, right_values)

  # calculation by the minimizer of system parameters
  res = scipy.optimize.minimize( HF_tub_en, x0, args = ( z, gen_param, gen_left_values, gen_right_values, sep_numbers ) )
  
  # calculation of energy for given values
  en = HF_tub_en( res.x, z, gen_param, gen_left_values, gen_right_values, sep_numbers )

  # saving results at each iteration
  if save_data == True:
    x_min = res.x
    data_number = g_f.define_number()
    results = total_grid (z, x_min, sep_numbers)
    save_total_grid(results, number = data_number, energy_res = np.pi*en, gen_left_values = gen_left_values, gen_right_values = gen_right_values)
  
  return en



''' Function for constructing the energy densities of the central circle of the hemifusion structure.
numbers - list consisting of numbers under which grids are saved
number - the number under which the graph will be saved '''

def en_den_com_graph (numbers, number = None):

  data = []
  axis_y = []
  axis_x = []

  for one in numbers:

    data = ( load_total_grid (f'Total_grid_{one}') )

    axis_y.append( data[1][0] )
    axis_x.append( g_f.grid_to_vek( 'point_bi',data[0][0] )[1] )
    
  s_f.en_den_graph(axis_y, axis_x, number)