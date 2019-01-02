from statistics import stdev
filename = "config/log/Graph1.txt"

gps_x_data = []
acc_x_data = []

with open(filename , 'r') as f:
    first_line = f.readline()
    for line in f:
        gps_x_data.append(float(line.split(',')[1]))

filename = "config/log/Graph2.txt"
with open(filename , 'r') as f:
    first_line = f.readline()
    for line in f:
        acc_x_data.append(float(line.split(',')[1]))

print(len(gps_x_data))
print(gps_x_data)

print(len(acc_x_data))
print(acc_x_data)

print(stdev(gps_x_data))
print(stdev(acc_x_data))
