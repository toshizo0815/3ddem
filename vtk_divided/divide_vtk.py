import numpy as np
import os
import matplotlib

matplotlib.use("Agg")

output_folder_master = "data_divided"
number = 2000

output_folder1 = output_folder_master + "_small"
output_folder2 = output_folder_master + "_large"

os.makedirs(output_folder1, exist_ok=True)
os.makedirs(output_folder2, exist_ok=True)

count_out = 0
for i in range(0, 500000, 500):
    count_out += 1
    print(i)
    
    datafile = "./../vtk_data_mb/data%07d.vtk" % (i)
    f1 = open(output_folder1 + "/data%07d.vtk" % (count_out-1), "w")
    f2 = open(output_folder2 + "/data%07d.vtk" % (count_out-1), "w")
    x_data = np.genfromtxt(datafile, skip_header=5, skip_footer=2*number+6)
    d_data = np.genfromtxt(datafile, skip_header=5+number+3, skip_footer=number+2)
    
    count1 = 0
    count2 = 0
    x1_data = []
    x2_data = []
    
    for j in range(number):
        if (d_data[j] == 0.5):
            count1 += 1
            x1_data.append(x_data[j])
        elif (d_data[j] == 0.4):
            count2 += 1
            print(j)
            x2_data.append(x_data[j])
            
    f1.write("# vtk DataFile Version 3.0\nx-t output\nASCII \nDATASET UNSTRUCTURED_GRID\nPOINTS %9d float" % count1)
    f2.write("# vtk DataFile Version 3.0\nx-t output\nASCII \nDATASET UNSTRUCTURED_GRID\nPOINTS %9d float" % count2)
    
    for j in range(count1):
        f1.write("\n%8.2f%8.2f%8.2f" % (x1_data[j][0], x1_data[j][1], x1_data[j][2]))
    for j in range(count2):
        f2.write("\n%8.2f%8.2f%8.2f" % (x2_data[j][0], x2_data[j][1], x2_data[j][2]))
    
    f1.write("\nPOINT_DATA %9d\nSCALARS radius float\nLOOKUP_TABLE default" % count1)
    f2.write("\nPOINT_DATA %9d\nSCALARS radius float\nLOOKUP_TABLE default" % count2)
    
    for j in range(count1):
        f1.write("\n1.00")
    for j in range(count2):
        f2.write("\n1.80")
    
    f1.close()
    f2.close()
    
