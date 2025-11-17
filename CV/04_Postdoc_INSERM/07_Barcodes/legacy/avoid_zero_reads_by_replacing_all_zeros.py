

import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_all_replaced.csv", 'w') as file_write:
    with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered.csv") as file_read:
        data = file_read.readlines()
        file_write.write(data[0])
        
        ## write a new matrix converting the zero reads into 0.001
        for line in data[1:] :
            line = line.replace("\n","").split(";")
            line_fill = [i for i in line]
            for i in range(1, len(line)) :
                if float(line[i]) == float(0) :
                    line_fill[i] = str(0.01)

            # fill the new line in new file
            file_write.write(';'.join(line_fill))
            file_write.write('\n')




stop = timeit.default_timer()
print(stop - start)  