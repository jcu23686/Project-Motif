import sys
import math

def main():
    #open motif file with motif information, background fasta file, and list size
    motif_file = open(sys.argv[1], 'r')
    background_file = open(sys.argv[2], 'r')
    min_score = int(sys.argv[3])

    motif_file_path = motif_file.name
    background_file_path = background_file.name

    #get base statistics regarding the background file
    char_count = 0
    count = {'A':0, 'C':0, 'G':0, 'T':0}
    motif_length = 0
    number_motifs = 0
    background_array = []

    other_count = 0
    header = False

    #find statistics of the background file
    for line in background_file:
        for char in line:
            if char == '>':
                header = True
            elif header and char == '\n':
                header = False
            elif char == '\n':
                continue
            elif not header:
                char_count += 1
                if char in count:
                    count[char] += 1
                    background_array.append(char)
                    if background_file == 'A':
                        count['A']+=1
                    elif background_file == 'C':
                        count['C']+=1
                    elif background_file == 'G':
                        count['G']+=1
                    elif background_file == 'T':
                        count['T']+=1    
                else:
                    other_count += 1   

    percent_counts = dict(count)

    for key in percent_counts:
        percent_counts[key] = count[key]/char_count


    for line in motif_file:
        for char in line:
            if char !='\n':
                motif_length+=1
        number_motifs +=1


    motif_length = int(motif_length/number_motifs)
    ##display statistics
    #print(f"The length of the background file is {char_count}\nIn the background file there are {other_count} other characters")
    #print(f"In the in the background file there are:\n{count['A']} Adenines ({percent_counts['A']:.6f}%)\n{count['C']} Cytosines ({percent_counts['C']:.6f}%)\n{count['G']} Guanines ({percent_counts['G']:.6f}%)\n{count['T']} Thymines ({percent_counts['T']:.6f}%)")
    #print(f"There are {number_motifs} motifs in the motif file provided\nEach motif is ({motif_length}) base pairs long")

    #write to file
    file_path = "PSSM_OUTPUT.txt"
    with open(file_path, "w") as file:
        file.write(f"Position Specific Scoring Matrix (PSSM) for motif finding\n_________________________________________________\nBackgournd statistics:\nThe Background file is {background_file_path}\nThe Motifs file is {motif_file_path}\n")
        file.write(f"The length of the background file is {char_count}\nIn the Background file there are {other_count} other characters\n")
        file.write(f"In the Background file there are:\n{count['A']}\tAdenine\t({percent_counts['A']:.6f}%)\n{count['C']}\tCytosine({percent_counts['C']:.6f}%)\n{count['G']}\tGuanine\t({percent_counts['G']:.6f}%)\n{count['T']}\tThymine\t({percent_counts['T']:.6f}%)")
        file.write(f"\n\n_________________________________________________\n\nMotif Statistics:\nThere are {number_motifs} motifs\nEach motif is {motif_length} base pairs long")

    background_file.seek(0)
    
    motifs_table = [['0' for _ in range(motif_length)] for _ in range(number_motifs)]
    
    motif_file.seek(0)
    
    i =0
    for line in motif_file:
        j=0

        with open(file_path, "a") as file:
            file.write(f"\nMotif #{i} is: ")

        for char in line:
            if char!='\n':
                motifs_table[i][j] = char.upper()
                #print(f"{motifs_table[i][j]} ", end="")
                with open(file_path, "a") as file:
                    file.write(f"{motifs_table[i][j]}")
            j+=1
        #print()
        i+=1

    flipped_scores = frequency_table = scores_matrix = [[0 for _ in range(4)] for _ in range(motif_length)]

    #print("\nNucleotides at each Position\n")
    
    #find nucletides at each position
    for i in range(number_motifs):
        for j in range(motif_length):
            if motifs_table[i][j] == 'A':
                frequency_table[j][0] += 1
            elif motifs_table[i][j] == 'C':
                frequency_table[j][1] += 1
            elif motifs_table[i][j] == 'G':
                frequency_table[j][2] += 1
            elif motifs_table[i][j] == 'T':
                frequency_table[j][3] += 1

    #add .25 to all data points for smoothing (we don't want issues with probabilities of 0)
    for i in range(motif_length):
        for j in range(4):
            frequency_table[i][j] += 0.25

    #find counts at each position and write to file
    count_file = "PSSMCounts.txt"
    with open(count_file, "w") as write_count:
        for i in range(motif_length):
            for j in range(4):
                    write_count.write(f"{frequency_table[i][j]}\t")
            write_count.write("\n")
    
    #find frequencies
    for i in range(motif_length):
        for j in range(4):
            frequency_table[i][j] = ((frequency_table[i][j])/(number_motifs+1))

    #display nucletides at each position 
    with open(file_path, "a") as file:
        file.write(f"\n\n_________________________________________________\nNucleotides at Each Position\nBase#\n  ")
        for i in range(motif_length):
            file.write(f"{i+1}\t")

    for i in range(4):
        #print(f"\nBase #{i}\t", end="")
        with open(file_path, "a") as file:
            if i == 0:
                file.write(f"\nA ")
            if i == 1:
                file.write(f"\nC ")
            if i == 2:
                file.write(f"\nG ")
            if i == 3:
                file.write(f"\nT ")


        for j in range(motif_length):
            #print(f"{frequency_table[i][j]}\t", end="")
            with open(file_path, "a") as file:
                file.write(f"{frequency_table[j][i]:.3f}\t")
        #print()
    
    #display nucletides at each position    
    #for i in range(motif_length):
        #print(f"\nBase #{i}\t", end="")
        #for j in range(4):
            #print(f"{frequency_table[i][j]:.6f}\t", end="")
        #print()

    #scores assigned including the background model
    for i in range(motif_length):
        for j in range(4):
            if j == 0:
                scores_matrix[i][j] =math.log2((frequency_table[i][j])/(percent_counts['A']))
            if j == 1:
                scores_matrix[i][j]=math.log2((frequency_table[i][j])/(percent_counts['C']))
            if j == 2:
                scores_matrix[i][j]=math.log2((frequency_table[i][j])/(percent_counts['G']))
            if j == 3:
                scores_matrix[i][j]=math.log2((frequency_table[i][j])/(percent_counts['T']))   

    #Display the scores matrix
    #print("Scores including the background model")
    #for i in range(motif_length):
     #   print(f"\nBase #{i}\t", end="")
     #   for j in range(4):
      #      print(f"{scores_matrix[i][j]:.6f}\t", end="")
      #  print() 

    #reverse strand
    flipped_scores = [row[::-1] for row in scores_matrix[::-1]]

    #Display the flipped scores matrix
    #print("Scores including the background model")
    #for i in range(motif_length):
        #print(f"\nBase #{i}\t", end="")
        #for j in range(4):
            #print(f"{flipped_scores[i][j]:.6f}\t", end="")
        #print()

    found_motif = ['']*motif_length

    #print("Forward Strand Motifs")

    forward_strand_count =1
    reverse_strand_count =1
    
    with open(file_path, "a") as file:
        file.write(f"\n\nThe minimum score provided for a motif is: {min_score} \n_________________________________________________\nFordward Strand Motifs\n")
    
    #function is for searching through the background file and trying to find motifs
    for a in range(char_count):
        motif_score = 0.00

        if a + motif_length > len(background_array):  # Check if accessing beyond the list
            break

        i = 0
        for b in range(a, a + motif_length):
            if background_array[b] == 'A':
                found_motif[i] = 'A'
                motif_score += scores_matrix[i][0]
                i += 1
            elif background_array[b] == 'C':
                found_motif[i] = 'C'
                motif_score += scores_matrix[i][1]
                i += 1
            elif background_array[b] == 'G':
                found_motif[i] = 'G'
                motif_score += scores_matrix[i][2]
                i += 1
            elif background_array[b] == 'T':
                found_motif[i] = 'T'
                motif_score += scores_matrix[i][3]
                i += 1

        if motif_score > min_score:
            
            with open(file_path, "a") as file:
                file.write(f"\nForward Strand Score #{forward_strand_count}: {motif_score:.6f}\t",)
                file.write(f"\tMotif Start is {a + 1}\tMotif End is {a + motif_length}\t")


            #print(f"\nForward Strand Score {forward_strand_count}: {motif_score:.6f}", end='\t')
            #print(f"Motif Start is {a + 1}\tMotif End is {a + motif_length}", end='\t')
            forward_strand_count +=1
            for p in range(motif_length):
                
                with open(file_path, "a") as file:
                    file.write(f"{found_motif[p]}")


                #print(found_motif[p], end='')

    reverse_found_motif = ['']*motif_length

    #print("\nReverse Strand Motifs")
    with open(file_path, "a") as file:
        file.write("\n_________________________________________________\nReverse Strand motifs\n")

    
    for a in range(char_count):
        motif_score = 0.00

        if a + motif_length > len(background_array):  # Check if accessing beyond the list
            break

        i = 0
        for b in range(a, a + motif_length):
            if background_array[b] == 'A':
                reverse_found_motif[i] = 'A'
                motif_score += flipped_scores[i][0]
                i += 1
            elif background_array[b] == 'C':
                reverse_found_motif[i] = 'C'
                motif_score += flipped_scores[i][1]
                i += 1
            elif background_array[b] == 'G':
                reverse_found_motif[i] = 'G'
                motif_score += flipped_scores[i][2]
                i += 1
            elif background_array[b] == 'T':
                reverse_found_motif[i] = 'T'
                motif_score += flipped_scores[i][3]
                i += 1

        if motif_score > min_score:
            with open(file_path, "a") as file:
                file.write(f"\nReverse Strand Score #{reverse_strand_count}: {motif_score:.6f}\t",)
                file.write(f"\tMotif Start is {a + 1}\tMotif End is {a + motif_length}\t")
            reverse_strand_count +=1

            #print(f"\nReverse Strand Score: {motif_score:.6f}", end='\t')
            #print(f"Motif Start is {a + 1}\tMotif End is {a + motif_length}", end ='\t')
            for p in range(motif_length):
                with open(file_path, "a") as file:
                    file.write(f"{found_motif[p]}")

                #print(reverse_found_motif[p], end='')

    print("\nPSSM Program Ran")

if __name__ == '__main__':
    main()