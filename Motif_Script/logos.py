# useful imports
import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import logomaker as lm
import csv
import sys




def main():
    
    #tables to hold motif bases
    data2 = []
    gibbs_data = []

    #count files procuced from PSSM and GIBBS scripts
    file_path = "PSSMCounts.txt"
    gibbs_file = "GIBBSCounts.txt"

    #open the gibbs file and add it to gibbs_data
    with open(gibbs_file, 'r') as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            if row[-1] == '':
                row.pop()
            row = [float(value) for value in row]
            gibbs_data.append(row)

    #open the pssm file and add it to data2
    with open(file_path, 'r') as file2:
        reader2 = csv.reader(file2, delimiter="\t")
        for row2 in reader2:
            if row2[-1] == '':
                row2.pop()
            row2 = [float(value) for value in row2]
            data2.append(row2)

    #convert each of the tables to pandas data frames (it is what the logo function requires)
    df2 = pd.DataFrame(data2)
    gibbs_df = pd.DataFrame(gibbs_data)

    #add headers to each pandas data frame
    new_column_names = ['A','C','G','T']
    df2.columns = new_column_names
    gibbs_df.columns = new_column_names

    #print(df2)

    #create logo for PSSM
    lm.Logo(df2)
    logo = lm.Logo(df2)
    plt.savefig('PSSM_logo_image1.png')
    plt.show()

    #create logo for Gibbs
    lm.Logo(gibbs_df)
    logoGibbs = lm.Logo(gibbs_df)
    plt.savefig('Gibbs_logo_image.png')
    plt.show()

    #print(data2)

    print("\nLogos Creation Program Ran")
if __name__ == '__main__':
    main()