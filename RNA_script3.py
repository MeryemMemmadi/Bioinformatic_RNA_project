from  RNA_script1 import calcul_distance, sep_line
import argparse
import math
import sys
import os

def interpolation(d, couple_r, eg):
    e=0
    while not os.path.exists("score_"+couple_r+".txt") :
        couple_r = "".join(reversed(couple_r))

    #get the file of the corresponding couple_r
    
    with open("score_"+couple_r+".txt", "r") as f1 :  
        lst = []
        for line in f1:
            lst.append([ x for x in line.split()])

        dist = [ x[0] for x in lst]
        scores = [ x[1] for x in lst]

        for i in range(len(dist)) : 
            
            
            if int(dist[i]) == int(math.floor(d)) :
                
                x1 = float(math.floor(d))
                x2 = float(math.ceil(d))
                y1 = float(scores[i])
                y2 = float(scores[i+1])
                
                e = y1 + (d - x1) * (y2 - y1)/(x2 - x1) 
        
                eg = eg+e
    return eg


   
parser = argparse.ArgumentParser()
parser.add_argument('-f', type=str, help="pdb file")
args = parser.parse_args()

pdb = args.f
energie_G=0

#do the same as in th script 1, open the pdb 2 times and go through it to compute the score for each C3' base pair
with open(pdb,"r") as f1_in:                        
    for ligne_f1 in f1_in:
        ligne_sep1 = sep_line(ligne_f1)            
        
        if ligne_f1.startswith("ATOM") and  ligne_sep1[2]== " C3'" :
            if ligne_f1[17:20] == "  A" or ligne_f1[17:20] ==  "  G" or ligne_f1[17:20] == "  U" or ligne_f1[17:20] == "  C" :
                
                with open(pdb,"r") as f2_in : 
                    for ligne_f2 in f2_in: 
                        ligne_sep2 = sep_line(ligne_f2)
                        
                        if ligne_f2.startswith("ATOM") and ligne_sep2[2]== " C3'" :
                            if ligne_f2[17:20] == "  A" or ligne_f2[17:20] ==  "  G" or ligne_f2[17:20] == "  U" or ligne_f2[17:20] == "  C" : 
                                if int(ligne_sep2[5]) >= int(ligne_sep1[5])+4 and ligne_sep1[4] == ligne_sep2[4]: 
                        
                                    d = calcul_distance(ligne_sep1, ligne_sep2)
                                    if d <=20 : 
                                        couple_r = ligne_sep1[3] + ligne_sep2[3]
                                        energie_G = interpolation(d, couple_r, energie_G)
                                        
print("Gibbs energy : ",energie_G)                                                              
