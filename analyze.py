#!/usr/bin/python
import random
import math

hbarc_GeV_fm = 0.197326972 #according to Google

#returns one column from a list of rows of data
def select_column(data, column):
  return [line[column] for line in data]

#returns a list of columns consisting of the rows of the set of input columns
def transpose(columns):
  return [[columns[i][j] for i in range(len(columns))] for j in range(len(columns[0]))]

def flatten(list_of_lists):
  return [x for li in list_of_lists for x in li]

#reads a simple data file consisting of whitespace delimited columns and 
#returns a list of columns of data. will throw away the first exclude_lines rows
def read_columns(filename, exclude_lines=0):
  with open(filename) as f:
    lines = f.readlines()
    lines = lines[exclude_lines:]
    data = [[float(s) for s in line.split()] for line in lines if not line.startswith("#")]
  columns = []
  for c in range(len(data[0])):
    columns.append(select_column(data, c))
  return columns


#reads in a file written by AlgPlaqT, which records the average
#plaquette on each time slice for each configuration. 
#returns (spacelike, timelike), which are each a list of columns of
#average plaquettes, one column per time slice 
def read_plaq_t(filename, exclude_lines=0):
  spacelike_rows = []
  timelike_rows = []
  with open(filename) as f:
    for line in f:
      if "time_slice" in line:
        #start of a new configuration
        spacelike_rows.append([])
        timelike_rows.append([])
      else:
        #record data for one time slice
        data = [float(s) for s in line.split()]
        spacelike_rows[-1].append(data[1])
        timelike_rows[-1].append(data[2])
  
  spacelike_rows = spacelike_rows[exclude_lines:]
  timelike_rows = timelike_rows[exclude_lines:]

  spacelike_columns = []
  timelike_columns = []
  for c in range(len(spacelike_rows[0])):
    spacelike_columns.append(select_column(spacelike_rows, c))
    timelike_columns.append(select_column(timelike_rows, c))
  
  return (spacelike_columns, timelike_columns)


#from de Forcrand et al.
def improved_Q(Q_11, Q_12, Q_22, Q_33, Q_13, c5):
  c_11 = (19. - 55. * c5) / 9.
  c_22 = (1. - 64. * c5) / 9.
  c_12 = (-64. + 640. * c5) / 45.
  c_13 = 1./5. - 2. * c5
  c_33 = c5
  return c_11 * Q_11 + c_12 * Q_12 + c_13 * Q_13 + c_22 * Q_22 + c_33 * Q_33 


#reads in a file written by AlgTcharge, and returns a column of global
#topological charge measurements. 
#type is '5Li', '4Li', '3Li', or '1x1'
def read_tcharge_global(filename, type):
  if type == '5Li': return read_tcharge_global_NLi(filename, 1./20.)
  elif type == '4Li': return read_tcharge_global_NLi(filename, 1./10.)
  elif type == '3Li': return read_tcharge_global_NLi(filename, 0.)
  elif type == '1x1': return read_tcharge_global_1x1(filename)
  else: raise Exception("READ_TCHARGE_GLOBAL GOT UNRECOGNIZED TYPE: " + str(type))

def read_tcharge_global_NLi(filename, c5):
  with open(filename) as f:
    lines = f.readlines()

  impr = []

  line = 0
  waiting_for_new = True

  while line < len(lines):
    if "AlgTcharge:" in lines[line]:
      Q_11 = float(lines[line + 7].split()[3])
      Q_12 = float(lines[line + 12].split()[3])
      Q_22 = float(lines[line + 16].split()[3])
      Q_33 = float(lines[line + 19].split()[3])
      Q_13 = float(lines[line + 21].split()[3])
      line += 22

      impr.append(improved_Q(Q_11, Q_12, Q_22, Q_33, Q_13, c5)) 
    else:
      line += 1
  return impr

def read_tcharge_global_1x1(filename):
  with open(filename) as f:
    lines = f.readlines()

  impr = []

  line = 0

  while line < len(lines):
    if "AlgApeSmear" in lines[line] or "AlgWilsonFlow" in lines[line]:
      line += 1
      continue
    else:
      Q_11 = float(lines[line + 7].split()[3])
      line += 27

      impr.append(Q_11)

  return impr



#Returns a list of Nt lists, each representing the MD-time
#history of topological charge measurements on a single time slice.
#type is '5Li', '4Li', '3Li', or '1x1'
def read_tcharge_t(filename, type):
  if type == '5Li': return read_tcharge_t_NLi(filename, 1./20.)
  elif type == '4Li': return read_tcharge_t_NLi(filename, 1./10.)
  elif type == '3Li': return read_tcharge_t_NLi(filename, 0.)
  elif type == '1x1': return read_tcharge_t_1x1(filename)
  else: raise Exception("READ_TCHARGE_T GOT UNRECOGNIZED TYPE: " + str(type))

def read_tcharge_t_NLi(filename, c5): 
  with open(filename) as f:
    lines = f.readlines()

  impr_t = []

  line = 0

  while line < len(lines):
    if "AlgApeSmear" in lines[line] or "AlgWilsonFlow" in lines[line]:
      line += 1
      continue
    else:
      Q_11_t = [float(s) for s in lines[line+22].split()[2:]]
      Q_12_t = [float(s) for s in lines[line+23].split()[2:]]
      Q_22_t = [float(s) for s in lines[line+24].split()[2:]]
      Q_33_t = [float(s) for s in lines[line+25].split()[2:]]
      Q_13_t = [float(s) for s in lines[line+26].split()[2:]]
      line += 27

      impr_t.append([improved_Q(Q_11_t[t], Q_12_t[t], Q_22_t[t], Q_33_t[t], Q_13_t[t], c5) for t in range(len(Q_11_t))])

  return [select_column(impr_t, t) for t in range(len(impr_t[0]))]

#Reads the topological charge computed from the 1x1 clover on each time slice
def read_tcharge_t_1x1(filename):
  with open(filename) as f:
    lines = f.readlines()

  impr_t = []

  line = 0

  while line < len(lines):
    if "AlgApeSmear" in lines[line] or "AlgWilsonFlow" in lines[line]:
      line += 1
      continue
    else:
      Q_11_t = [float(s) for s in lines[line+22].split()[2:]]
      line += 27

      impr_t.append(Q_11_t)

  return [select_column(impr_t, t) for t in range(len(impr_t[0]))]
  


def read_moments(filename):
  columns = []
  Q = []
  with open(filename) as f:
    for line in f:
      if(len(line.split()) != 2):
	continue
      if line.split()[0] == "Q":
        Q.append(float(line.split()[1]))
      else:
        power = int(line.split()[0])
	moment = float(line.split()[1])
        while len(columns) < power:
 	  columns.append([])
        columns[power-1].append(moment)
  return (columns, Q)


#write out a data file consisting of the given columns in a form readable by gnuplot
def write_gnuplot_file(filename, columns):
  with open(filename, 'w') as f:
    for i in range(len(columns[0])):
      for column in columns:
	f.write(str(column[i]) + '\t')
      f.write('\n')

def write_file(filename, text):
  with open(filename, 'w') as f:
    f.write(text)

#only keep each sample_time'th data point
def thin_data(data, sample_time, start_offset = 0):
  return [data[i] for i in range(start_offset, len(data), sample_time)]

#separates an interleaved column into num_columns separate columns
def split_column(data, num_columns):
   truncated_data = data[:((len(data)/num_columns) * num_columns)]
   return [thin_data(truncated_data, num_columns, c) for c in range(num_columns)]

#compresses a data column into bins by replacing each bin of data with its mean
def bin_column(data, bin_size):
  return [mean(data[i:(i+bin_size)]) for i in range(0, len(data), bin_size)]

#selects a random sample from rows with replacement, with length equal to len(rows)
def bootstrap_sample(rows):
  sample = []
  for i in range(len(rows)):
    sample.append(rows[random.randint(0, len(rows)-1)])
  return sample

#mean of a list of numbers
def mean(data):
  return sum(data)/float(len(data))

#given a list of lists, all of the same length,
#averages all the lists together
def list_mean(data):
  N = len(data)
  L = len(data[0])
  return [mean([data[n][i] for n in range(N)]) for i in range(L)]

#given a list of lists, all of the same length,
#sums all the lists together
def list_sum(data):
  N = len(data)
  L = len(data[0])
  return [sum([data[n][i] for n in range(N)]) for i in range(L)]

#variance of a list of numbers
def var(data):
  N = len(data)
  # abs protects from negative numbers due to floating point precision errors:
  return N/(N-1.0) * (mean([x*x for x in data]) - mean(data)**2)

#standard deviation of a list of numbers
def std(data):
  # abs protects from negative numbers due to floating point precision errors:
  return math.sqrt(abs(var(data)))

#root-mean-square of a set of numbers
def rms(data):
  return math.sqrt(mean([x*x for x in data]))

#statistical error of the measurement of an expectation value. Generates a number
#of bootstrap samples and finds the standard deviation of the means of the samples
def statistical_error_bootstrap(data, bin_size, num_bootstraps):
  binned = bin_column(data, bin_size)
  means = []
  for i in range(num_bootstraps):
    means.append(mean(bootstrap_sample(binned)))
  return std(means)

#compute the statistical error of an rms measurement
def statistical_error_rms_bootstrap(data, num_bootstraps):
  rmses = []
  for i in range(num_bootstraps):
    rmses.append(rms(bootstrap_sample(data)))
  return std(rmses)

#autocorrelation, normalized so that autocorrelation(delay = 0) = 1
def autocorrelation(data, delay):
  avg = mean(data)
  return autocorrelation_known_mean(data, delay, avg)

def autocorrelation_known_mean(data, delay, avg):
  covar = mean([(data[i] - avg) * (data[i+delay] - avg) for i in range(len(data)-delay)])
  covar0 = mean([(x - avg) * (x - avg) for x in data])
  if covar==0 or covar0==0:
    return 1
  else:
    return covar / covar0

def autocorrelation_known_mean_replicas(datas, delay, avg):
  return mean([autocorrelation_known_mean(data, delay, avg) for data in datas])

#an estimate of the statistical error on autocorrelation estimator above
# delay is the delay of the autocorrelation whose error we want to estimate
# N is the total number of data points
# autocorrs is a previously computed list of autocorrelation estimates, indexed by delay
# cut is a cutoff for the sum: we effectively set rho(k) to zero for |k| > cut
def autocorrelation_error(delay, N, autocors, cut):
  assert delay >= 0 # code doesn't handle negative delays properly
  assert cut < len(autocors) # make sure we don't need rho beyond where it has been computed
  t = delay
  var = 0
  for k in range(1, cut+t+1):
    summand = 0
    if k+t <= cut: summand += autocors[k+t]
    if abs(k-t) <= cut: summand += autocors[abs(k-t)]
    if k <= cut and t <= cut: summand += -2*autocors[k]*autocors[t]
    summand **= 2
    var += summand
  var /= N
  return math.sqrt(var)

  #return math.sqrt(1./N * sum([(autocors[abs(k+t)] + autocors[abs(k-t)] - 2*autocors[k]*autocors[t])**2 for k in range(1, min(t+cut, N-t))]))

def autocorrelation_error_replicas(delay, N, autocors, cut, Nreplicas):
  return autocorrelation_error(delay, N, autocors, cut) / math.sqrt(Nreplicas)

#sum of a autocorrelation over all delays
def integrated_autocorrelation_time(data, max_delay):
  return 0.5 + sum([autocorrelation(data, delay) for delay in range(1, max_delay+1)])

def integrated_autocorrelation_time_known_mean(data, max_delay, avg):
  return 0.5 + sum([autocorrelation_known_mean(data, delay, avg) for delay in range(1, max_delay+1)])
 
#estimate of the statistical error of the above
def iat_error(iat, max_delay, num_data_points):
  return math.sqrt(2.0 * (2.0 * max_delay + 1.0) / num_data_points) * iat #the "Madras-Sokal formula"

def iat_error_replicas(iat, max_delay, num_data_points, Nreplicas):
  return iat_error(iat, max_delay, num_data_points) / math.sqrt(Nreplicas)

#returns an analysis of the IATs and errors on the IATs for a given range of delays
def do_iats_and_errors(data, max_delay, time_step):
  max_delays = [i for i in range(max_delay+1)]
  iats = [integrated_autocorrelation_time(data, m) for m in max_delays]
  iat_errors = [iat_error(iats[i], max_delays[i], len(data)) for i in range(len(max_delays))]
  max_delays = [time_step * x for x in max_delays]
  iats = [time_step * x for x in iats]
  iat_errors = [time_step * x for x in iat_errors]
  return [max_delays, iats, iat_errors]

def do_iats_and_errors_known_mean(data, max_delay, time_step, avg):
  max_delays = [i for i in range(max_delay+1)]
  iats = [integrated_autocorrelation_time_known_mean(data, m, avg) for m in max_delays]
  iat_errors = [iat_error(iats[i], max_delays[i], len(data)) for i in range(len(max_delays))]
  max_delays = [time_step * x for x in max_delays]
  iats = [time_step * x for x in iats]
  iat_errors = [time_step * x for x in iat_errors]
  return [max_delays, iats, iat_errors]


#Suppose we have a quantity measured on each time slice on each configuration. 
#Denote it X(t, i) where t is the time slice and i is the configuration.
#This function computes < X(t,i) X(t+delta_slice,i+delay) >, where 
#the lattice is assumed to be periodic in the time direction.
#
#data should be a list of lists; the first index should be the time slice
#and the second index should be the configuration
def autocorrelation_time_slice_avg(data, delta_slice, delay, known_mean):
  T = len(data)
  N = len(data[0])
  sample_sum = 0
  Nsamples = 0
  for t in range(T):
    for i in range(N-delay):
      sample_sum += (data[t][i] - known_mean) * (data[(t+delta_slice)%T][i+delay] - known_mean)
      Nsamples += 1
  return sample_sum / Nsamples
  
def normalized_autocorrelation_time_slice_avg(data, delta_slice, delay, known_mean, variance):
  return autocorrelation_time_slice_avg(data, delta_slice, delay, known_mean) / variance

#Computes the generalized autocorrelation functions for a set of Euclidean time translations
#of MD time histories
def build_rho(data, max_delay, known_mean):
  T = len(data)
  N = len(data[0])
  variance = autocorrelation_time_slice_avg(data, 0, 0, known_mean)

  rho = []
  for t in range(T):
    rho.append([0]*(N+max_delay))

  for t in range(0, T):
    for i in range(0, max_delay+1):
      rho[t][i] = normalized_autocorrelation_time_slice_avg(data, t, i, known_mean, variance)

  return rho

def build_rho_replicas(datas, max_delay, known_mean):
  rhos = [build_rho(data, max_delay, known_mean) for data in datas]
  ret = []
  T = len(datas[0])
  N = min([len(data[0]) for rep in data])
  for t in range(T):
    ret.append([0]*(N+max_delay))
  for t in range(0, T):
    for i in range(0, max_delay+1):
      ret[t][i] = mean([rho[t][i] for rho in rhos])
  return ret 
  

def autocorrelation_time_slice_avg_error(data, rho, delta_slice, delay, cut):
  s = delta_slice
  m = delay
  T = len(data)
  N = len(data[0])

  assert m >= 0
  assert cut < len(rho[0]) # make sure we don't need rho beyond where it has been computed

  var = 0
  for v in range(T):
    for k in range(-(cut+m), cut+m+1):
      summand = 0
      if abs(k+m) <= cut: summand += rho[(v+s)%T][abs(k+m)]
      if abs(k-m) <= cut: summand += rho[(v-s+T)%T][abs(k-m)]
      if m <= cut and abs(k) <= cut: summand += -2 * rho[s][m] * rho[v][abs(k)]
      summand **= 2
      var += summand
  var /= (2*N*T)
  return math.sqrt(var)

  #summands = [[(rho[(v+s)%T][abs(k+m)] + rho[(v-s+T)%T][abs(k-m)] - 2*rho[s][m]*rho[v][abs(k)])**2 \
  #                for k in range(-max_delay-m, max_delay+m+1)] for v in range(T)]
  #return math.sqrt(1./(2*N*T) * sum([sum(summands[t]) for t in range(T)]))

#FIXME!!!!!
def autocorrelation_time_slice_avg_error_replicas(datas, rho, delta_slice, delay, max_delay, Nreplicas):
  return autocorrelation_time_slice_avg_error(datas[0], rho, delta_slice, delay, max_delay) / math.sqrt(Nreplicas)

def time_slice_avg_iat(rho, delta_slice, max_delay):
  v = delta_slice
  return 0.5 * sum([rho[v][abs(k)] for k in range(-max_delay, max_delay+1)])

def time_slice_avg_iat_error(generalized_iats, max_delay, num_confs):
  M = max_delay
  N = num_confs
  return math.sqrt(2.0*(2.0*M + 1)/N * mean([x*x for x in generalized_iats]))

#FIXME!!!!!
def time_slice_avg_iat_error_replicas(generalized_iats, max_delay, num_confs, Nreplicas):
  return time_slice_avg_iat_error(generalized_iats, max_delay, num_confs) / math.sqrt(Nreplicas)



#lattice spacing parametrizations:

#computes the iwasaki lattice spacing in units of r0
#error should be < 0.8% in the range 2.16 < beta < 2.71 according to hep-lat/0306005
def iwasaki_spacing_r0(beta):
  if beta < 2.16 or beta > 2.71:
    print "Warning: iwasaki_spacing_r0 is being used outside the domain of its parametrization."
  c1 = -2.1281
  c2 = -1.0056
  c3 = 0.6041
  return math.exp(c1  +  c2 * (beta - 3) +   c3 * (beta - 3)**2)

#computes the wilson lattice spacing in units of r0
#parametrization is from hep-lat/0306005 and is claimed valid for 5.7 < beta < 6.92
#with errors < 1.0%
def wilson_spacing_r0(beta):
  if beta < 5.7 or beta > 6.92:
    print "Warning: wilson_spacing_r0 is being used outside the domain of its parametrization."
  a0 = -1.6804
  a1 = -1.7331
  a2 = 0.7849
  a3 = -0.4428
  return math.exp( a0  +  a1 * (beta - 6)  +  a2 * (beta - 6)**2  +  a3 * (beta - 6)**3 )

#computes the DBW2 lattice spacing in units of r0
#parametrization is from hep-lat/0306005 and is claimed valid for 0.75696 < beta < 1.04
#with errors < 0.8%
#this paper defines the DBW2 action to use c1 = -1.4088
def dbw2_spacing_r0(beta):
  if beta < 0.75696 or beta > 1.04:
    print "Warning: dbw2_spacing_r0 is being used outside the domain of its parametrization"
  d1 = -1.6007
  d2 = -2.3179
  d3 = -0.8020
  d4 = -19.8509
  return math.exp( d1 + d2 * (beta - 1) + d3 * (beta - 1)**2 + d4 * (beta - 1)**3 )
  

def special_f(beta):
  pi = 3.14159265358979
  b0 = 11.0 / (4 * pi)**2
  b1 = 102.0 / (4 * pi)**4
  g = math.sqrt( 6.0 / beta )
  return (b0 * g**2) ** (-b1 / (2 * b0**2)) * math.exp(- 1.0 / (2 * b0 * g**2))

#from hep-lat/9905005:
def iwasaki_spacing_times_root_sigma(beta):
  if beta < 2.15 or beta > 3.2:
    print "Warning: iwasaki_spacing_times_root_sigma is being used outside the domain of its parametrization"
  ahat = special_f(beta) / special_f(2.4)
  c0 = 0.524
  c2 = 0.274
  c4 = 0.105
  return special_f(beta) * (1 + c2 * ahat**2 + c4 * ahat**4) / c0

#iwasaki spacing in physical units using the convention r0 = 0.5 fm
def iwasaki_inverse_spacing_GeV(beta):
  hbarc = 0.197326938 #in GeV*fm
  r0 = 0.5 #by convention we take r0 = 0.5 fm
  r0rootsigma = 1.1525 #dimensionless, implies sqrt(sigma) = 455 MeV, chosen to match spacing from r0 and spacing from sigma  

  print "Using convention r0 = 0.5 fm"

  from_r0 = hbarc / (iwasaki_spacing_r0(beta) * r0)
  print "1/a from r0 =", from_r0, "GeV"
  from_sigma = hbarc / (iwasaki_spacing_times_root_sigma(beta) * r0 / r0rootsigma)
  print "1/a from sigma =", from_sigma, "GeV"

  print "difference =", 100*(from_sigma - from_r0)/from_sigma, "%"

def wilson_inverse_spacing_GeV(beta):
  hbarc = 0.197326938 #in GeV*fm
  r0 = 0.5 #by convention we take r0 = 0.5 fm

  print "Using convention r0 = 0.5 fm"

  from_r0 = hbarc / (wilson_spacing_r0(beta) * r0)
  print "1/a from r0 =", from_r0, "GeV"

def dbw2_inverse_spacing_GeV(beta):
  hbarc = 0.197326938 #in GeV*fm
  r0 = 0.5 #by convention we take r0 = 0.5 fm
  
  print "Using convention r0 = 0.5 fm"

  from_r0 = hbarc / (dbw2_spacing_r0(beta) * r0)
  print "1/a from r0 =", from_r0, "GeV"
  


#put a value and error into a string of the form 1.23456(78)
def format_error(value, error):
  value_place_value = int(math.floor(math.log10(abs(value)))) #the place in which the value has its first nonzero digit
  error_place_value = int(math.floor(math.log10(error))) #the place in which the error has its first nonzero digit
  if error_place_value >= 0:
    print "format_error is not designed for errors >= 1!!!!!!!!"
  value_digits_past_decimal = 1 - error_place_value
  value_string = ("%." + str(value_digits_past_decimal) + "f") % value
  truncated_error = int(error * 10**(1 - error_place_value)) #should be a 2 digit number
  error_string = str(truncated_error)
  return value_string + "(" + error_string + ")"

  
#t are value of t/a^2
#E is a^4 E = a^4 / 4 * tr(G_mu_nu G_mu_nu)
#Finds the value of sqrt(t) such that t^2 E = 0.3
#Result is a physical length scale in expressed in units of a
def flow_scale_luscher(t, E, threshold = 0.3):
  t2E = [t[i]**2 * E[i] for i in range(len(t))]
  if max(t2E) < threshold:
    raise Exception("Didn't flow long enough!")
    #(lower_i, upper_i) = (len(t)-2, len(t)-1)
  else:
    (lower_i, upper_i) = next((i, i+1) for i in range(len(t)-1) if t2E[i+1] > threshold)
  lower_t2E = t2E[lower_i]
  upper_t2E = t2E[upper_i]
  
  scale2 = (t[lower_i] * (upper_t2E - threshold) + t[upper_i] * (threshold - lower_t2E)) / (upper_t2E - lower_t2E)
  return math.sqrt(max(0, scale2))

#t are values of t/a^2
#E is a^4 E = a^4 / 4 * tr(G_mu_nu G_mu_nu)
#Finds the value of sqrt(t) such that t (d/dt) (t^2 E)
#Result is a physical length scale expressed in units of a
def flow_scale_BMW(t, E, threshold = 0.3):
  if any([math.isnan(x) for x in E]): return 0.0

  t2E = [t[i]**2 * E[i] for i in range(len(t))]
  
  #ddt[i] is (d/dt)(t^2 E) at t = (t[i] + t[i+1])/2
  ddt = [(t2E[i+1] - t2E[i]) / (t[i+1] - t[i]) for i in range(len(t)-1)]

  #tddt[i] is t (d/dt)(t^2 E) at t = (t[i] = t[i+1])/2
  tddt = [(t[i] + t[i+1])/2 * ddt[i] for i in range(len(ddt))]

  if max(tddt) < threshold:
    #(lower_i, upper_i) = (len(tddt)-2, len(tddt)-1)
    raise Exception("Didn't flow long enough!")
  else:
    (lower_i, upper_i) = next((i, i+1) for i in range(len(tddt)-1) if tddt[i+1] > threshold)

  lower_t = (t[lower_i] + t[lower_i+1])/2
  upper_t = (t[upper_i] + t[upper_i+1])/2

  lower_tddt = tddt[lower_i]
  upper_tddt = tddt[upper_i]

  scale2 = (lower_t * (upper_tddt - threshold) + upper_t * (threshold - lower_tddt)) / (upper_tddt - lower_tddt)
  return math.sqrt(max(0, scale2))


def read_scale_setting(filename):
  ts = []
  Es = []
  t = 0
  with open(filename, 'r') as f:
    for line in f:
      if 'AlgWilsonFlow' in line:
        #lines have the form: "AlgWilsonFlow: dt = 1.000000e-02"
	dt = float(line.split()[3])
	t += dt
      else:
	E = float(line.strip())
	ts.append(t)
	Es.append(E)
  return (ts, Es)

def read_scale_setting_ape(filename):
  ts = []
  Es = []
  t = 0
  with open(filename, 'r') as f:
    for line in f:
      if 'AlgApeSmear' in line: 
        #lines have the form: "AlgApeSmear: coef = 4.500000e-01"
	coeff = float(line.split()[3])
        dt = coeff / 6.0 #comes from matching leading terms of APE and WFlow kernels
	t += dt
      else:
	E = float(line.strip())
	ts.append(t)
	Es.append(E)
  return (ts, Es)

def mean_columns(columns):
  N = len(columns[0])
  return [mean([columns[c][i] for c in range(len(columns))]) for i in range(N)]

def differentiate(x, y):
  diffx = [(x[i] + x[i+1])/2 for i in range(len(x)-1)]
  derivative = [(y[i+1]-y[i])/(x[i+1]-x[i]) for i in range(len(x)-1)]
  return (diffx, derivative)

# Given a list of values, produces a set of jackknife samples, where
# each sample consists of the list of values with a set of block_size
# adjacent entries deleted.
def make_jackknife_samples(values, block_size):
  Njack = len(values) / block_size
  ret = []
  for i in range(Njack):
    ret.append(values[:(block_size*i)] + values[(block_size*(i+1)):])
  return ret

# Given a list of lists of values, produces a set of superjackknife samples
# E.g., values = [[a, b], [c, d]], block_sizes = [1, 1]
# returns [ [[b], [c, d]], 
#           [[a], [c, d]],
#           [[a, b], [d]],
#           [[a, b], [c]] ]
def make_superjackknife_samples(values, block_sizes):
  Nens = len(values) 
  ensemble_jackknifes = [make_jackknife_samples(values[n], block_sizes[n]) for n in range(Nens)]

  superjack_samples = []
  for n in range(Nens):
    ens_samples = ensemble_jackknifes[n]
    for j in range(len(ens_samples)):
      superjack_sample = values[:n] + [ens_samples[j]] + values[n+1:]
      superjack_samples.append(superjack_sample)

  return superjack_samples

def histogram(data, box_min, box_width, n_boxes):
  boxes = [0]*n_boxes
  for x in data:
    b = int(math.floor((x - box_min)/box_width))
    if b < 0 or b >= n_boxes: 
      raise Exception('Data value %f is outside the range of the histogram' % x)
    boxes[b] += 1

  return boxes
    

# Given a list of numbers computed from jackknife samples, and a central numbers,
# applies the jackknife error formula.
def jackknife_error(jackknife_values, central_value):
  N = len(jackknife_values)
  return math.sqrt((N-1.0)/N * sum([(x - central_value)**2 for x in jackknife_values]))

# Computes the jackknife error on a number, or a list of numbers, 
# or a list of lists of numbers, or a tuple of lists of tuples of numbers...
def jackknife_error_recursive(jackknife_values, central_value):
  Njack = len(jackknife_values)
  if isinstance(central_value, (list, tuple)):
    # Get the jackknife error on each element of the list or tuple:
    Nelem = len(central_value)
    element_jackknife_values = [[jackknife_values[j][i] for j in range(Njack)] for i in range(Nelem)]
    element_jackknife_errors = [jackknife_error_recursive(element_jackknife_values[i], central_value[i]) for i in range(Nelem)]
    if isinstance(central_value, list): return element_jackknife_errors
    else: return tuple(element_jackknife_errors)
  else:
    return jackknife_error(jackknife_values, central_value)

# Accepts a function, "computation", that accepts a list of values and 
# returns a number, or list of numbers, etc. Returns both the central
# value and jackknife error on the output of "computation"
def jackknife_central_and_error(computation, values, block_size):
  central_value = computation(values)
  
  jackknife_samples = make_jackknife_samples(values, block_size)
  jackknife_values = [computation(sample) for sample in jackknife_samples]

  error = jackknife_error_recursive(jackknife_values, central_value)
  return (central_value, error)

# Accepts a function, "computation", that accepts a list of lists of values and 
# returns a number, or list of numbers, etc. Returns both the central
# value and jackknife error on the output of "computation"
def superjackknife_central_and_error(computation, values, block_size):
  central_value = computation(values)
  
  superjackknife_samples = make_superjackknife_samples(values, block_size)
  jackknife_values = [computation(sample) for sample in superjackknife_samples]

  error = jackknife_error_recursive(jackknife_values, central_value)
  return (central_value, error)
  
def jackknife_mean_and_error(values, block_size):
  return jackknife_central_and_error(mean, values, block_size)



def format_error(value, error):
  assert value >= 0

  if error == 0: return "%0.16f" % value

  if error >= 10:
    return "%d(%d)" % (int(round(value)), int(round(error)))

  if error >= 1:
    return "%0.1f(%0.1f)" % (value, error)

  d = int(1.999999999-math.log10(error))
  format = "%%0.%df(%%d)" % d
  return format % (value, int(round(10**d * error)))








