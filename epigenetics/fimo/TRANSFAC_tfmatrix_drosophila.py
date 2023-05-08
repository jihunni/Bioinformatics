# to connect the TRANSFAC DB on PostgreSQL
import psycopg2
try:
    conn = psycopg2.connect("dbname=transfac_db user=postgres password=sandiaadmin")
except psycopg2.DatabaseError as db_err:
    print("Not connected!")
    print(db_err)

try:
    cur = conn.cursor() #create a new cursor
except psycopg2.DatabaseError as db_err:
    print("Not connected!")
    print(db_err)
    
# query
query = "SELECT accession, identifier, name, tfmatrix, organism FROM transfac.matrix WHERE organism = 'fruit fly, Drosophila melanogaster' or organism = 'Vertebrata';"

# query
query = "SELECT accession, identifier, name, tfmatrix, organism FROM transfac.matrix WHERE organism = 'fruit fly, Drosophila melanogaster' or organism = 'Insecta' or organism = 'Vertebrata';"

# to execute query
cur.execute(query)
rows=cur.fetchall()
print(rows)




##########################
# query all TF
#####################

# load query data
query_result = rows
print("The number of TF motif : ", len(query_result))

# environment parameter to save file
working_directory = "/home/jihun/python/tfmotif_drosophila_new/"
file_name = "TRANSFAC_tfmatrix_drosophila_MEME.txt"
file = open(working_directory + file_name, 'w')

# write header for MEME Motif format
write_contents_0 = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + - \n\nBackground letter frequencies (from uniform)\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
file.write(write_contents_0)

for i in range(len(query_result)): 
    motif_length = len(rows[i][3]['matrixPositions']) # length of tfmatrix (motif length)
   
    write_contents_1 = "MOTIF " + query_result[i][1] + " " + query_result[i][2] + "\n" + \
                        "letter-probability matrix: alength= 4 " + "w= " + str(motif_length) + "\n"
    write_contents_2 = ""
    for k in range(motif_length):
        tf_matrix_oneLine = query_result[i][3]['matrixPositions'][k]['bases']
        sum = tf_matrix_oneLine[0] + tf_matrix_oneLine[1] + tf_matrix_oneLine[2] + tf_matrix_oneLine[3]
        write_contents_2 = write_contents_2 + str(tf_matrix_oneLine[0]/sum) + " " + str(tf_matrix_oneLine[1]/sum) + \
        " " + str(tf_matrix_oneLine[2]/sum)+ " " + str(tf_matrix_oneLine[3]/sum) + "\n"
    write_contents = write_contents_1 + write_contents_2 + "\n"
    file.write(write_contents)
    if (i%10==0):
        print("current location: ", i)
file.close()

print("finsih")


#############################################
# query all TF and split into n files
##############################################
# load query data
query_result = rows
print("The number of TF motif : ", len(query_result))

# environment parameter to save file
num_motif_per_file = 60 #for 1 day job schedule ; 10 motif takes 4 hours.
if len(query_result) % num_motif_per_file == 0 :
    num_files = len(query_result) // num_motif_per_file
else : 
    num_files = len(query_result) // num_motif_per_file + 1
working_directory = "/home/jihun/python/tfmotif_drosophila_new_split/"

for file_iteration in range(num_files): 
    if (file_iteration < 10) : 
        file_name = "TRANSFAC_tfmatrix_drosophila_MEME_" + "0" + str(file_iteration) + ".txt"
    else : 
        file_name = "TRANSFAC_tfmatrix_drosophila_MEME_" + str(file_iteration) + ".txt"
    file = open(working_directory + file_name, 'w')

    # write header for MEME Motif format
    write_contents_0 = "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + - \n\nBackground letter frequencies (from uniform)\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
    file.write(write_contents_0)

    for i in range(60*file_iteration, 60*(file_iteration+1)): #fix the last col
        motif_length = len(rows[i][3]['matrixPositions']) # length of tfmatrix (motif length)

        write_contents_1 = "MOTIF " + query_result[i][1] + " " + query_result[i][2] + "\n" + \
                            "letter-probability matrix: alength= 4 " + "w= " + str(motif_length) + "\n"
        write_contents_2 = ""
        for k in range(motif_length):
            tf_matrix_oneLine = query_result[i][3]['matrixPositions'][k]['bases']
            sum = tf_matrix_oneLine[0] + tf_matrix_oneLine[1] + tf_matrix_oneLine[2] + tf_matrix_oneLine[3]
            write_contents_2 = write_contents_2 + str(tf_matrix_oneLine[0]/sum) + " " + str(tf_matrix_oneLine[1]/sum) + \
            " " + str(tf_matrix_oneLine[2]/sum)+ " " + str(tf_matrix_oneLine[3]/sum) + "\n"
        write_contents = write_contents_1 + write_contents_2 + "\n"
        file.write(write_contents)
        if (i%10==0):
            print("file_number : ", file_iteration, "\tcurrent location: ", i)
    file.close()

print("finsih")
