import math
import sys
import os
import argparse



def create_dico() :                 #create the nested dictionnary and fill it with 0
    pair_r = ["AA" ,"AU","AG","AC","GC","GU","GG","UC","UU","CC","Nxx"]
    dico_pair_r = dict()
    
    for i in pair_r:
        dico_pair_r[i] = {}
        for j in range(21):
            dico_pair_r[i][j] = int()
        dico_pair_r[i]['Nij'] = int()

    return dico_pair_r



def sep_line (l_f):                 #function that transform each line into a list
    ligne_sep = [ l_f[:6],  l_f[6:11],      l_f[12:16],     l_f[19:20],    l_f[21:22],     l_f[22:26],     l_f[30:38],     l_f[38:46],     l_f[46:54],     l_f[54:60]  ]
                  #ATOM    atom_serial_nb    atom_name      residue_name    chain_id        residue_seq_nb   coord_x        coord_y         coord_z         occupancy
                  
                  
    return ligne_sep


def calcul_distance(l1, l2) :       #function that compute the distance
    
    #attribuate each coordinate to variables
    x_a = float(l1[6])
    y_a = float(l1[7])
    z_a = float(l1[8])
    x_b = float(l2[6])
    y_b = float(l2[7])
    z_b = float(l2[8])

    d = ( ((x_b - x_a)**2) + ((y_b - y_a)**2) + ((z_b - z_a)**2) )**(1/2)   #compute the distance
    return d


def fill_dico(dico, d, l1,l2):      #fill the nested dictionnary with the count for each distances

    couple_r = l1[3] + l2[3]
    if couple_r not in dico :       #verify if the pair (or its opposite) is one of the dictionnary's key 
        couple_r = l2[3]+l1[3]      # ex: if "UA" not in the keys we trandform it as "AU
    dico[couple_r][d] = dico[couple_r][d]+1
    #print("pair :",couple_r, " dico_coupler :",dico[couple_r], " count  :",dico[couple_r][d])

    return dico

def fill_sum(dico):                  #make the total of each line and each column, fill Nxx and Nij
    for i in dico:
        for j in dico[i]:
            if j == "Nij" : 
                dico[i][j] = sum(dico[i].values() )    
    for i in dico :
        for j in dico[i] :
            sumValue_col = sum(d[j] for d in list(dico.values() )[:-1] if d)
            dico["Nxx"][j] = sumValue_col
    return dico


def freq_obs(dico) :                # compute the observed frequence 

    
    for i in list(dico)[:-1]:
        for j in list(dico[i])[:-1]:
            if dico[i]["Nij"] == 0 :
                    dico[i]["Nij"] = 1
            dico[i][j] =  (dico[i][j] / dico[i]["Nij"])
    return dico

def freq_ref(dico) :                # compute the reference frequence
    for i in list(dico["Nxx"])[:-1]:
        if dico["Nxx"][i] == 0 :
                    dico["Nxx"][i] = 1
        dico["Nxx"][i] =   ( dico["Nxx"][i] / dico["Nxx"]["Nij"])
    return dico

def comput_score(dico) :            # compute the score
    for i in list(dico)[:-1]:
        for j in list(dico[i])[:-1]:
            if dico["Nxx"][j] == 0 or dico[i][j] == 0:  #we replace 0 values by 10 to avoid errors
                dico[i][j] = 10
            else :
                dico[i][j] =   -(math.log( dico[i][j] / dico["Nxx"][j]))
    return dico




def dir_path(string):               #check if the given argument is a path
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
    

def main():    
    dico_pair_r = create_dico()    
    parser = argparse.ArgumentParser()
    parser.add_argument("-path",type=dir_path,  help="file containning the trainning set of pdb")
    #parser.add_argument('-filename',type =str, help="pdb file")
    args = parser.parse_args()

    dossier = os.listdir(args.path)           
               #create the nested dictionnary

    
    for pdb in dossier:                                  #search for pdb file
        print(args.path+pdb)
        with open(args.path+pdb,"r") as f1_in:       #open the pdb file
            for ligne_f1 in f1_in:
                ligne_sep1 = sep_line(ligne_f1)                                                                                             #separate the line
                if ligne_f1.startswith("ATOM") and  ligne_sep1[2]== " C3'" :                                                                # check if the line start with ATOM and if the atom is the C3' carbon
                    if ligne_f1[17:20] == "  A" or ligne_f1[17:20] ==  "  G" or ligne_f1[17:20] == "  U" or ligne_f1[17:20] == "  C" :      #verify that we have acid nucleic bases and not amino acids
                        
                        with open(args.path+pdb,"r") as f2_in :                                                                          #open the file a second time and apply the same conditions
                            for ligne_f2 in f2_in: 
                                ligne_sep2 = sep_line(ligne_f2)
                                if ligne_f2.startswith("ATOM") and ligne_sep2[2]== " C3'" :
                                    if ligne_f2[17:20] == "  A" or ligne_f2[17:20] ==  "  G" or ligne_f2[17:20] == "  U" or ligne_f2[17:20] == "  C" : 
                                        
                                        if int(ligne_sep2[5]) >= int(ligne_sep1[5])+4 and ligne_sep1[4] == ligne_sep2[4]:                   #check if the pair of bases are separated by at least 3 bases and if they are on the same chain
                                            d = calcul_distance(ligne_sep1, ligne_sep2)
                                            if d <=20 :                                                                                     #if the computed distance is inferior or equal to 0 then we can fill the dictionnary
                                                dico_pair_r = fill_dico(dico_pair_r, math.ceil(d), ligne_sep1, ligne_sep2)

    
    dico_pair_r = fill_sum(dico_pair_r)                 #fill the sum
    dico_pair_r = freq_obs(dico_pair_r)                 #fill the observed frequencies
    dico_pair_r = freq_ref(dico_pair_r)                 #fill the reference frequencies
    dico_pair_r = comput_score(dico_pair_r)             #fill with the computed score




    for i in list(dico_pair_r)[:-1]:                    #write the file with the scores for each pair of bases
            with open("liste_score/score_"+i+".txt", "w") as f_out : 

                for j in list(dico_pair_r[i])[:-1]:
                    f_out.write(str(j)+ " \t"+str(dico_pair_r[i][j] )+"\n" ) 
                    
                    
if __name__ == '__main__':
    main()