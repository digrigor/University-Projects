##experiment with two simple data classification approaches on the 
#Wisconsin Diagnostic Breast Cancer Dataset. The dataset contains morphological 
#data on the cell nuclei from biopsies of 569 patients with benign or malign breast cancer.

###File I/O - reading the wdbc data###

#Function that take as its argument the name of the data file, 
#and return a tuple of two dictionaries. For both dictionaries, 
#the patient IDs serve as key. The first dictionary contains the 
#diagnosis for each patient, and the second the feature vector.


def loadData(filename):
    '''
    ***Remember to insert the name of the data file followed by the file extension (e.g. wdbc.data) as the argument of 
    the function***\n
    This function takes the name of the  Wisconsin Diagnostic Breast Cancer Dataset data file (.data) as its argument, 
    reads the corresponding file and returns a tuple of two dictionaries. The first dictionary contains the 
    diagnosis for each patient (keys->Patient ID, values->diagnosis) and the second contains the feature vector
    of each patient (keys->Patient ID, values->List of 30 float numerical attributesdescribing the morphology of the nuclei).''' 
    
    wdbc_file=open(str(filename),"r")
    id_no=""
    diagnosis_dict={} #First dictionary which will contain the diagnosis for each patient's ID no.
    fvector_dict={} #Second dictionary which will contain the feature vector of each patient
    diag="" 
    flfv=[]
    float_feats=[] #List of float representing feature vectors.This is the value for each patient's ID in the features dictionary
    
    for ll in wdbc_file: 
        for i in range(0,len(ll)): #Reads the text file as a list. 
            if ll[i]!=",":
                id_no+=str(ll[i]) #First numbers of each list element of the file
                #until the first comma represent the patient id.
            else:
                diag=ll[i+1] #The letter after the first comma of each list element (M or B) represent the diagnosis 
                #of a listed patient (M or B).
                fvector=ll[i+3:len(ll)] #The remaining string contents of each list element represent the feature vector.
                break
        
        diagnosis_dict[id_no]=diag 
        lfv=fvector.split(",") #Converts fvector str into a list 
        
        for l in lfv:
            flfv.append(float(l)) #Converts every attribute of the list to a float number
        fvector_dict[id_no]=flfv
        id_no=""
        diag=""
        fvector=""
        lfv=[]
        flfv=[] #Reset of all variables
    wdbc_file.close()
    return (diagnosis_dict, fvector_dict)
	
# this cell should run without any errors (assuming wdbc.data is the name of the dataset).
# it performs a couple of checks to rule out common errors in reading the dataset.
# you should use the global variables initialized here in the rest of your code

# Use truelabels and features as the label and feature vector dictionaries respectively
truelabels, features=loadData('wdbc.data')
assert len(truelabels)==569, "Wrong data count - have all records been read?"
dimensions=30 # You are allowed to use this variable for the length of each feature vector
assert len(features['852552'])==dimensions, "Feature vectors have the wrong number of components"

###Computing averages###

def averageVec(vec_list):
    '''This function takes as its argument a list of vectors of equal lengths and returns the average vector 
    (ie, the vector formed by the average of the input vectors, computed component by component)'''
    
    average_vec_list=[]
    sum_vec=0
    z=0
    
    try:
        while z<len(max(vec_list, key=len)):    #z represents one number of a vector at a time 
            #until the last number of the largest vector (all vectors should be equal though). see below (except)
            for i in range(0,len(vec_list)):    #i represents one vector at a time 
                sum_vec+=float(vec_list[i][z])  #sums the nth attribute(z) of all vectors(i)
            av_vec=float(sum_vec)/len(vec_list) 
            average_vec_list.append(av_vec) #each time that a sum is calculated it is added to the final list
            sum_vec=0
            av_vec=0 #reset of the variables for the next loop
            z+=1 #loop will start again for the next attribute of the vectors
        return average_vec_list
    
    except IndexError: #if at least one vector is not equal with the other vectors, function returns an error"
        print "Your vectors are not equal. Please try again"   
	
#Use the averageVec function to compute the average feature vector 
#for patients in the Benign and in the Malign class
	
M_vec_av=[]
B_vec_av=[]

for j in truelabels: #Searches in truelabels for patient ids of M or B patients. Then uses this
    #patient ids to search in features for the patient's vector. Then adds these vectors into a list. 
    if truelabels[j]=="M":
        M_vec_av.append(features[j])
    elif truelabels[j]=="B":
        B_vec_av.append(features[j])

mavg=averageVec(M_vec_av) 
bavg=averageVec(B_vec_av)

print "The average feature vector for patients in the Benign class is:"
print bavg
print ""
print "The average feature vector for patients in the Malign class is:"
print mavg

###Computing Distances###
#Two functions named euclid and manhattan, each taking two vectors as their argument, 
#that return the Euclidean or, respectively, the Manhattan distance between the argument vectors.
import math
def euclid(veclist_1,veclist_2):
    '''Takes two vectors as an argument and returns the Euclidean distance
    between these argument vectors'''
    
    eu_sums=0
    
    for nu in range(0,len(veclist_1)): 
        eu_sums+=(float(veclist_1[nu])-float(veclist_2[nu]))**2 #Calculates the squares 
        #of the differences of the nth number of both given lists
    eu_dist=math.sqrt(eu_sums) 
    
    return eu_dist

def manhattan(veclist_1,veclist_2):
    '''Takes two vectors as an argument and returns the the Manhattan distance
    between these argument vectors'''
    
    manh_dist=0
    
    for nu2 in range(0,len(veclist_1)):
        manh_dist+=math.fabs(veclist_1[nu2]-veclist_2[nu2])
    
    return manh_dist

#Testing the functions
print euclid(features["8611555"],features["8611792"])
print manhattan(features["8611555"],features["8611792"])

###Partitioning the data###

def partition(centres,feat_vecs,e_or_m):
    '''This function that takes as its argument a list of two vectors (the "centres") and a dictionary 
       with patient IDs as keys and the feature vectors as values and an e_or_m variable (see below). The function 
       returns a newdictionary with the patient IDs as its keys. For each ID, the value is `0` if the feature vector 
       for the patient is closest to the first centre, and `1` otherwise. Distance between a patient feature
       vector and each center was computed using euclid distance if e_or_m is "Euclidean" or using manhattan
       distance if e_or_m is "Manhattan".'''
    
    part_dict1={}
    
    for i in feat_vecs:
        
        if e_or_m=="Euclidean":
            dist1=euclid(feat_vecs[i],centres[0]) #euclidean distances are computed
            dist2=euclid(feat_vecs[i],centres[1])
        elif e_or_m=="Manhattan":
            dist1=manhattan(feat_vecs[i],centres[0]) #Manhattan distances are computed
            dist2=manhattan(feat_vecs[i],centres[1])
        
        if dist1<dist2: #if the vector is closer to the first centre
            part_dict1[i]="0"
        elif dist1>dist2: #if the vector is closer to the second centre
            part_dict1[i]="1"
    
    return part_dict1
	
	
#Data partitioning using a list with the two centres bavg, mavg as its first argument 
#and the features dictionary as its second argument. Output is stored in variable 
#computed_labels. This dictionary will contain a 0 for each patient ID that is closer to 
#the average of the Benign class, or respectively a 1 for each patient closer to the average 
#of the Malign class. This is a simple classification of the data based on the distance 
#from a template.
	
m_bavg=[bavg,mavg]
computed_labels_e=partition(m_bavg,features,"Euclidean") #computed labels using euclidean distance for partitioning
computed_labels_m=partition(m_bavg,features,"Manhattan") #computed labels using manhattan distance for partitioning

###Evaluating the clustering algorithm###

from decimal import *
def normalization(a,b):
    '''This function takes 2 numbers as its arguments (a and b) and
       normalizes the "a" value in order a+b=1. The return value precision
       is restricted to 4 places maximum.'''
    
    getcontext().prec = 4
    
    if a!=0 and b!=0: 
        nor_a=Decimal(a)/Decimal(a+b)
        return Decimal(nor_a)
    elif a==0 and b!=0: 
        nor_a=0
        return Decimal(nor_a)
    elif a!=0 and b==0:
        nor_a=1
        return Decimal(nor_a)
    elif a==0 and b==0:
        nor_a=0
        return nor_a
    

def confusionMatrix(dict_BorM,dict_0or1):
    '''This function takes  dictionary with values B, M as its first argument (dict_BorM) 
       and a dictionary with values 0, 1 as its second argument (dict_0or1). Both dictionaries 
       should have the same keys. The function returns a list of values:a,b,c,d. a is the number
       of keys that have a value of B in the first dictionary and a value of 0 in the second 
       dictionary, divided by the total number of keys with a B label (b is the opposite: B and 1). Likewise, d 
       is the number of keys with values M in the first dictionary and 0 in the second divided by the total 
       number of keys with an M label (c is the opposite: M and 1). The values are normalized, meaning that a+b=1 
       and c+d=1.
       Use confusionMatrix_table_print function (See documentation) which takes as input the output of this function,  
       to print a table for these values.'''
    
    Bn0=0 #no. of patients with "B" value in the 1st dict and "0" in the 2nd dict 
    Bn1=0 #no. of patients with "B" value in the 1st dict and "1" in the 2nd dict  
    Mn0=0 #no. of patients with "M" value in the 1st dict and "0" in the 2nd dict
    Mn1=0 #no. of patients with "M" value in the 1st dict and "1" in the 2nd dict
    
    for ids in dict_BorM: #for every patient ID in the 1st dict
        
        if ids in dict_0or1: 
            if dict_0or1[ids]=="0":
                if dict_BorM[ids]=="B":
                    Bn0+=1
                elif dict_BorM[ids]=="M":
                    Mn0+=1                  #Numbers of patients which have the same values in both directories or not 
            elif dict_0or1[ids]=="1":       #are stored in the above defined values (Bn0,etc.)
                if dict_BorM[ids]=="M":
                    Mn1+=1
                elif dict_BorM[ids]=="B": 
                    Bn1+=1
        else:
            continue
    return [normalization(Bn0,Bn1),normalization(Bn1,Bn0),normalization(Mn0,Mn1),normalization(Mn1,Mn0)]
    #only returns the 4 requested values in the following order: Bn0, Bn1, Mn0, Mn1 to use them in following questions. 
    #confusionMatrix_table_print function is used for the table printing. 

def confusionMatrix_table_print(out_list):
    '''This function takes the output values of the confusionMatrix function as its arguments 
       and returns a table with these values. The function prints a table as the following:
       ------------------
       Counts    0  1
            B    a  b
            M    c  d
       ------------------ 
       see confusionMatrix function documentation for values definition'''
    
    #Just pretty printing below:
    return "Count"+"    0 "+"     1  "+'\n'+"    B  "+\
    str(out_list[0])+"   "+str(out_list[1])+\
    '\n'+"    M  "+str(out_list[2])+"   "+str(out_list[3])    
	
#Application of the confusionMatrix function to truelabels and computed_labels.
euclidean_output=confusionMatrix(truelabels,computed_labels_e)
manhattan_output=confusionMatrix(truelabels,computed_labels_m)
print "Confusion matrix using Euclidean distance for partitioning:"
print confusionMatrix_table_print(euclidean_output)
print "----------------------"
print "Confusion matrix using Manhattan distance for partitioning:"
print confusionMatrix_table_print(manhattan_output)

###Computing the average confusion matrix###

#A function avgConfusionMatrix that works like confusionMatrix,
#except that it takes a list of dictionaries with values 0 or 1 as its second argument.

def avgConfusionMatrix(dict_BorM, dict_list):
    
    '''This function works like confusionMatrix, except that it takes a list of dictionaries with values 0 or 1 
    as its second argument (see confusionMatrix documentation). This Function avgConfusionMatrix computes a confusion 
    matrix between its first argument and each of the dictionaries in its second argument calculating a single average 
    confusion Matrix. This function returns a list of values like confusionMatrix function (See confusionMatrix documentation 
    more info. The confusionMatrix function which is called inside this function has the ability to check which keys are 
    in common between the two under-comparison dictionaries. Therefore, only the common keys between the two dictionaries 
    are considered for the computation of the corresponding confusion matrix.'''
    
    Bn0_sum=0
    Mn0_sum=0
    
    for dn in range(0,len(dict_list)): 
        output_list=confusionMatrix(dict_BorM, dict_list[dn]) #Computes the confusion matrix for each two dictionaries
        Bn0_sum+=(output_list[0])                             #each time and...
        Mn0_sum+=(output_list[2])
    Bn0_av=Bn0_sum/len(dict_list)
    Mn0_av=Mn0_sum/len(dict_list)                             #...calculates the average for all of them
    Mn1_av=1-Mn0_av
    Bn1_av=1-Bn0_av
    
    return [Bn0_av,Bn1_av,Mn0_av,Mn1_av]                      #same return as confusionMatrix function

	
	
	
	
	import random
import itertools

#Code that implements 10-fold cross-validation by randomly splitting the set 
#of patients into ten subsets of (as much as possible) equal size
def ten_fold_validation(e_or_m):
    '''This function takes "Euclidean" or "Manhattan" as its argument which define which of these 
    distances will be used in the partition function which is called inside the function. The function 
    randomly splits the set of patients included in the truelabels dictionary into ten subsets of 56 patients
    Each one of the subsets is called test set and the others constitute the training set. Average feature vectors
    for the Benign and Malign classes are computed over the training set and used to call partition on 
    the test set, thus yielding a dictionary of computed_labels for it.  The function rotates the test set over the
    ten possible choices, each time adding the corresponding dictionary of computed_labels to a list (tenfold_val_e)
    if Euclidean distance was used or (tenfold_val_m) if Manhattan distance was used for partition'''
    
    #definition of used variables:
    random_ID_list=[] #a list of random patients IDs
    test_set_features={} #a dictionary test_set IDs and their feature vectors
    training_vectors=[] #a list of features for every training_set ID
    training_M=[]
    training_B=[]
    tenfold_val_e=[]
    tenfold_val_m=[]
    true_keys_list=list(truelabels.keys())
    
    for countl in range(0,10):
        random_sample=random.sample(true_keys_list,56)
        random_ID_list.append(random_sample) #This list contains 56 random IDs derived from the truelabels dictionary.
        for ele in random_sample: #Each time that an ID is added to the random_ID_list is removed from the list of the 
            #true labels IDs
            true_keys_list.remove(ele)
    
    for rans in random_ID_list: 
        training_IDs_splitted=list(random_ID_list)
        training_IDs_splitted.remove(rans)
        training_IDs=list(itertools.chain.from_iterable(training_IDs_splitted)) #this list contains all the IDs of the 
        #random_ID list appart from the IDs of the test_set (the current 'rans' of the random_ID_list). So this list  
        #contains all the IDs of the training set
        for tID in training_IDs: 
            if truelabels[tID]=="M":
                training_M.append(features[tID]) #Creates a dictionary with the Malign-patients' IDs of the 
                #training set as keys and their corresponding feature vectors as values.
            if truelabels[tID]=="B": 
                training_B.append(features[tID]) #Creates a dictionary with the Benign-patients' IDs of the 
                #training set as keys and their corresponding feature vectors as values.
        averM=averageVec(training_M) 
        averB=averageVec(training_B) #average feature vectors for the Benign and Malign classes are computed
        M_Bavg=[averB,averM] #and they are concatenated into one list
    
        for ent in rans: #For every ID of the test_set 
            test_set_features[ent]=features[ent] #Creates a dictionary with the test_set patient's IDs as keys
            #and their corresponding feature vectors as values
            t_computed_labels_e=partition(M_Bavg,test_set_features,"Euclidean")
            t_computed_labels_m=partition(M_Bavg,test_set_features,"Manhattan") #Partitioning: test_set vs. training_set
            #based averages. t_computed_labels --> dictionaries {patient's ID: 0 or 1}
        
        tenfold_val_e.append(t_computed_labels_e)
        tenfold_val_m.append(t_computed_labels_m) #Each corresponding dictionary of computed_labels produced is added 
        #to the final list
        training_M=[]
        training_B=[]
        test_set_features={}
    
    #Euclidean or Manhattan distance for partitioning?
    if e_or_m=="Euclidean":
        return tenfold_val_e
    elif e_or_m=="Manhattan":
        return tenfold_val_e
		
#The ground truth labels truelabels are being passed together with the list of computed label 
#dictionaries from the cross-validation to avgConfusionMatrix.
#The calculation is repeated a few times, each time with a different random partitioning 
#of the data in ten sets.
		
b0_val_e=m0_val_e=b0_val_m=m0_val_m=0
times_count=0

for ki in range(1,11): #calculation is repeated 10 times.
    #Calculation and printing of each average computation matrix for each random partitioning of the data 
    tenfold_random_e=ten_fold_validation("Euclidean") #random partitioning using the Euclidean distance
    tenfold_random_m=ten_fold_validation("Manhattan") #random partitioning using the Manhattan distance
    b0_e,b1_e,m0_e,m1_e=avgConfusionMatrix(truelabels, tenfold_random_e) 
    
    print "---Random Partioning:",ki,"---"
    print "Average Confusion matrix using Euclidean distance for partitioning:"
    print confusionMatrix_table_print([b0_e,b1_e,m0_e,m1_e]) #printing function defined earlier is used
    b0_m,b1_m,m0_m,m1_m=avgConfusionMatrix(truelabels, tenfold_random_m)
    print "Average Confusion matrix using Manhattan distance for partitioning:"
    print confusionMatrix_table_print([b0_m,b1_m,m0_m,m1_m])
    print "---------------------------"
    
    #Calculation and printing of the average confusion Matrix of all the confusion matrixes produced by multiple
    #random partitioning
    b0_val_e+=b0_e
    m0_val_e+=m0_e
    b0_val_m+=b0_m    
    m0_val_m+=m0_m    
    times_count+=1
b0_avg_e=b0_val_e/times_count
m0_avg_e=m0_val_e/times_count
b0_avg_m=b0_val_m/times_count
m0_avg_m=m0_val_m/times_count


print ""
print ""
print "---Average confusion matrices for all the random partitioning data---"
print "Using Euclidean distance:"
print confusionMatrix_table_print([b0_avg_e,1-b0_avg_e,m0_avg_e,1-m0_avg_e]) #printing function defined earlier is used
print "Using Manhattan distance"
print confusionMatrix_table_print([b0_avg_m,1-b0_avg_m,m0_avg_m,1-m0_avg_m])

###Unsupervised classification###
###Sampling the Dataset###

def randomCentres(inp_dict):
    '''This function takes as its argument a dictionary. The function returns 
    a list with two values chosen at random from those in the dictionary. '''
    out_val_list=[]
    out_val_list.append(random.sample(inp_dict.values(),2))
    return out_val_list[0]
	
	
#randomCentres was used on the features dictionary to obtain 
#a list of two random feature vectors; this is stored in the centres variable.
#dio_toy dictionary was used for testing the code.
#dio_toy={'1267':[4.0,3.1,3.9,4.2],'1268':[3.1,3.7,4.12,4],'1269':[3.4,4.5,2.9,3],'1270':[8.3,9.1,8.8,8.2],'1271':[9,9.2,8.3,7.9],'1272':[10,8.6,8.8,8],'1273':[4.0,3.1,3.9,4.2],'1274':[3.1,3.7,4.12,4],'1275':[3.4,4.5,2.9,3],'1276':[8.3,9.1,8.8,8.2],'1277':[9,9.2,8.3,7.9],'1278':[10,8.6,8.8,8]}
centres=randomCentres(features)

#Code that implements the following loop:
#call partitions with the (initially) random centres and features to obtain computed_labels
#use averageVec on the list of feature vectors in the 0 class, respectively in the 1 class, to compute centres for the two classes
#compare the new centres to the old ones using the distance function: if the distance is zero (they have not changed), exit
#replace the previous centres with the newly computed centres of the 0 and 1 classes, iterate

iterats=0 #number of iterations
computed_labels=[]
given_centres=centres

#If the code runs without randomizing the data will throw an IndexError. Let's prevent this from happening:
try:
    while True:
        class1_list=[]
        class0_list=[]
        computed_labels=partition(given_centres,features,"Euclidean")

        for ID_c in features:
            if computed_labels[ID_c]=='1':
                class1_list.append(features[ID_c])
            elif computed_labels[ID_c]=='0':
                class0_list.append(features[ID_c])
        avg_class1=averageVec(class1_list)
        avg_class0=averageVec(class0_list)
        
        if euclid(given_centres[1],avg_class1)==0 and euclid(given_centres[0],avg_class0)==0:
            print "Done!"
            centres=[]
            break
        else:
            iterats+=1
            given_centres=[avg_class0,avg_class1]

    print "Iterations required for convergence:",iterats

except IndexError:
    print"Go on the previous cell and randomize your Data!"
	
		
#dio_labels dictionary was used for testing the code.
#dio_labels={'1267':'B','1268':'B','1269':'B','1270':'M','1271':'M','1272':'M','1273':'B','1274':'B','1275':'B','1276':'M','1277':'M','1277':'M'}
out_val_list=confusionMatrix(truelabels,computed_labels)
print confusionMatrix_table_print(out_val_list)

#the distance between the centre of the Benign (resp. Malign) class and the centre 
#of the 0 or 1 class  (whichever is closer) is computed:
if out_val_list[0]>out_val_list[1]:
    B_class_dist=euclid(bavg,avg_class0)
    M_class_dist=euclid(mavg,avg_class1)
elif out_val_list[0]<out_val_list[1]:
    B_class_dist=euclid(bavg,avg_class1)
    M_class_dist=euclid(mavg,avg_class0)

#Print the distances and divide each one of them by the distance between the centres of the Benign and Malign classes
print "The distance between the Bening centre and its closest centre between 0 and 1 class centres (B_class_distance) is:"
print B_class_dist
print ""

print "The distance between the Malign centre and its closest centre between 0 and 1 class centres (M_class_distance) is:"
print M_class_dist
print ""

print "B_class_distance divided by the distance between the centres of the Benign and Malign classes is:"
print B_class_dist/euclid(mavg,bavg)
print ""

print "M_class_distance divided by the distance between the centres of the Benign and Malign classes is:"
print M_class_dist/euclid(mavg,bavg)

