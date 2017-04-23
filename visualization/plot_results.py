import pandas as pd
import matplotlib.pyplot as plt

datafile = pd.read_csv('../build/dataset-2-output.txt', sep='\t')

x_true = datafile.px_true
y_true = datafile.py_true
x_estimated = datafile.px
y_estimated = datafile.py
x_measured = datafile.px_measured
y_measured = datafile.py_measured

true_position, = plt.plot(x_true, y_true, 'b', label="True position")
estimated_position, = plt.plot(x_estimated, y_estimated, 'g', label="Estimated position")
measured_position, = plt.plot(x_measured, y_measured, 'r', label="Measured position")
plt.legend(handles=[true_position, estimated_position, measured_position])
plt.show()