

import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads.csv", 'a+') as file_write:
    with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered.csv") as file_read:
        data = file_read.readlines()
        file_write.write(data[0])
        
        for line in data[:1] : # we parse only the column names, 1st line of the matrix
            line = line.replace("\n","").split(";")

            # create a list of all experiment names
            exp_list = []
            for full_id in line[1:] : # we start from 1 to avoid the empty first cell
                #(print(full_id))

                experiment_name = (full_id.split("_"))[2] # we get the experiment name
                exp_list.append(experiment_name) # and append it to exp_list
            exp_list = list(set(exp_list))
            # print(len(exp_list))
            # for exp in exp_list :
            # print(exp)

            # get indexes of the treatment replicates for each drug/dose/exp combination and store in a dict
            drug_dose_exp_dict = {}
            for i in range(1,len(line)): # we start from 1 to avoid the empty first cell
                #print(line[i])
                drug = line[i].split("_")[0][:-1]
                dose = line[i].split("_")[1]
                exp = line[i].split("_")[2]
                drug_dose_exp = "_".join([drug, dose, exp])
                if drug_dose_exp not in drug_dose_exp_dict :
                    drug_dose_exp_dict[drug_dose_exp] = [i]
                else :
                    drug_dose_exp_dict[drug_dose_exp].append(i)
            # print(drug_dose_exp_dict)
            for drug_dose_exp in drug_dose_exp_dict :
                print("%s: %s" % (drug_dose_exp, drug_dose_exp_dict[drug_dose_exp]))

            # get indexes of the controls for each exp and store in a dict
            control_dict = {}
            for k in range(1,len(line)): # we start from 1 to avoid the empty first cell
                #print(line[k])
                exp = line[k].split("_")[2]
                if "Contro" in line[k] :
                    if exp not in control_dict :
                        control_dict[exp] = [k]
                    else :
                        control_dict[exp].append(k)
            # print(control_dict)
            # for exp in control_dict :
            #     print("%s: %s" % (exp, control_dict[exp]))

        # in each experiment, if sum reads for a treatment is zero and sum controls is non-zero, store the barcode and the treatment id in a dict
        # and write a new matrix converting the zero reads into 0.001
        zero_reads_dict = {}
        for line in data[1:] :
            line = line.replace("\n","").split(";")
            line_fill = [i for i in line]
            for drug_dose_exp in drug_dose_exp_dict :
                #print(drug_dose_exp)
                
                # get experiment id
                experiment = drug_dose_exp.split("_")[2]
                
                # get control reads indexes from that experiment
                control_reads_indexes = control_dict[experiment]
                
                # make sum of controls
                sum_controls = sum([float(line[j]) for j in control_reads_indexes])
                
                # if sum controls is non-zero, then check sum of treatments
                if (sum_controls != 0) :

                    treatment_indexes = drug_dose_exp_dict[drug_dose_exp]
                    reads = [float(line[j]) for j in treatment_indexes]
                    sum_reads = sum(reads)
                    #print(sum_reads)
                    if (float(sum_reads) == float(0)) :
                        barcode = line[0]
                        #print(barcode, drug_dose_exp)
                        if (barcode not in zero_reads_dict) :
                            zero_reads_dict[barcode] = [drug_dose_exp]
                        else :
                            zero_reads_dict[barcode].append(drug_dose_exp)
                        # replace in line_fill the zero values by 0.001
                        for h in treatment_indexes :
                            line_fill[h] = '0.01'
            
            # fill the new line in new file
            file_write.write(';'.join(line_fill))
            file_write.write('\n')


        # for barcode in zero_reads_dict :
        #     print(barcode, zero_reads_dict[barcode])

        #print(len(zero_reads_dict)) # 11777 (out of 12305... seems quite high)

stop = timeit.default_timer()
print(stop - start)  