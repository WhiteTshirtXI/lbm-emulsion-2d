def find(A,val):
    '''finds the position of a number in an array'''
    import numpy as np
    counter=0
    counter2=0
    result=[]
    b=A.shape[0]
    for i in range(A.shape[0]):
        if len(A.shape) is not 1:
            for j in range(A.shape[1]):
               if A[i,j] == val:
                    result.append(counter)
                    counter = counter+1
               elif A[i,j] is not val:
                    counter = counter+1
        elif (len(A.shape) is 1) and (A[i] == val):
            result.append(counter)
            counter = counter+1
            return result
        elif (len(A.shape) is 1) and (A[i] is not val):
            counter=counter+1

    resultf=np.zeros([len(result),2])
    for i in range(len(result)):
        if np.floor(result[i]/A.shape[1]) is 0:
            resultf[i,0]= np.floor(result[i]/A.shape[1])
            resultf[i,1]= result[i]
        elif np.floor(result[i]/A.shape[1])is not 0:
            resultf[i,0]= np.floor(result[i]/A.shape[1])
            resultf[i,1]= result[i]-A.shape[1]*resultf[i,0]
        
    return resultf
