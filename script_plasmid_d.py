# -*- coding: utf-8 -*-

"""
Created on Thu Apr 23 17:23:43 2020

@author: dalonsom and davidcalfran
"""
# Load the Pandas libraries with alias 'pd' 
import pandas as pd 
import glob

path = r"/Volumes/TOSHIIBA/TFM BIOINFORMATICA/bin3c_clust/fasta" # use your path
all_files = glob.glob(path + "/*.fna")
#

pd_list = [] #una lista de ficheros donde voy a guardar los ficheros filtrados. Se genera vacia al principio.
for f in all_files:  #all_files es la lista y la f es cada uno de los elementos
    print(f)
    doc2_temp = pd.read_csv(f, sep='\s+', header=None) #leo el fichero 1
    #Doc2 - Filter Rows with its first column not ">.*"
    doc2_temp = doc2_temp[doc2_temp[0].str.contains(">")].ix[:,[0,1]] #tratamiento que haciamos antes: filtrar del mayor que y quedarse solo con la columna 0 y 1.
    #Doc2 - Filter First and Second Column
    doc2_temp.ix[:,[0,1]]
    pd_list.append(doc2_temp)

doc2 = pd.concat(pd_list, axis=0, ignore_index=True)#coge la lista y lo concatena, para tener 1 solo dataframe en lugar de una lista de data.frames. El axis 0 es que lo juntas por el eje 0, que me ponga uno debajo del otro y au. Esto me genera un fichero/documento.
    
# Read data from file 'filename.csv' 
# (in the same directory that your python process is based)
# Control delimiters, rows, column names with read_csv (see later) 
doc1 = pd.read_csv("plasmid_top_hits.txt", sep='\s+', header=None)
#doc2 = pd.read_csv("CL0003.txt", sep='\s+', header=None) 
doc3 = pd.read_csv("contig_links.txt", sep='\s+', header=None) 

# Data cleaning 

#Doc1 - Keep column 1 delete others
doc1 = doc1.ix[:,[0]]

#Doc 3 -Delete rows with equal sign in any column
doc3 = doc3[~doc3[0].str.contains("=")]
doc3 = doc3[~doc3[1].str.contains("=")]
#if k99 in doc3 is not in doc2 - delete interaction


# Data preparation
doc2[0] = doc2[0].str.replace(">","")
doc2[1] = doc2[1].str.replace("contig:","")
k99_list=doc2[1].drop_duplicates().values.tolist()
doc3 = doc3[doc3[0].isin(k99_list)]
doc3 = doc3[doc3[1].isin(k99_list)]


#Join doc1 % doc2
doc12=doc1.merge(doc2, on=[0], how='left')

#Join doc12 % doc3
doc12_doc3_c1=doc12.merge(doc3, left_on=[1], right_on=[0], how='left')
doc12_doc3_c1=doc12_doc3_c1.ix[:,['0_x','1_x','1_y',2]]
doc12_doc3_c2=doc12.merge(doc3, left_on=[1], right_on=[1], how='left')

#Rename columns
header_list = ["CL-doc1", "k99-doc2", "k99-interac","value"]
doc12_doc3_c1.columns = header_list
doc12_doc3_c2.columns = header_list

#Append
resultado=doc12_doc3_c1.append(doc12_doc3_c2).sort_values(by=['CL-doc1'])

#Write Result
resultado.to_csv("result_plasmid.txt", index=False)

#GroupBy["CL-doc1", "k99-doc2"] and sum
resultado_suma=resultado.groupby(["CL-doc1", "k99-doc2"]).sum()
resultado_suma.to_csv("result-sum_plasmid.txt")

#Calculate normalization file by number of contigs per cluster
doc1_norm = pd.read_csv("plasmid_top_hits.txt", sep='\s+', header=None)
doc1_norm = doc1_norm.ix[:,[0]]
doc1_norm = doc1_norm[0].str.slice(0, 6).to_frame().drop_duplicates(keep='first').dropna().reset_index(drop=True)
doc2_filter= doc2[0].str.slice(0, 6).to_frame()
doc2_filter[1]= doc2_filter[0]
doc2_filter=doc2_filter.groupby([1], as_index=False).count()
join_norm=doc1_norm.merge(doc2_filter, left_on=[0], right_on=[1], how='left').ix[:,['key_0','0_y']].fillna(0)
header_list = ["CLXXXX", "count-doc2"]
join_norm.columns = header_list
join_norm["count-doc2"]=join_norm["count-doc2"].astype(int)
join_norm.to_csv("result_CLXXXXOcurrences_plasmid.txt", index=False)
