def circshift(A,fac1,fac2):
        '''shifts matrix circularly'''
        import numpy as np
        if fac1>0:
                temp1=np.array([A[0]])
                temp2=np.array([A[1]])
                for i in range(A.shape[0]):
                        A[i+1,:]=temp1
                        temp1[:]=temp2[:]
                        if (i+2)>(A.shape[0]-1):
                                A[0,:]=temp1
                                break
                        elif ((i+2)<(A.shape[0]-1)) or ((i+2) == (A.shape[0]-1)):
                                temp2[:]=A[i+2,:]

        elif fac1<0:
                temp1=np.array([A[0]])
                temp2=np.array([A[-1]])
                for i in range(A.shape[0]):
                        A[-i-1,:]=temp1
                        temp1[:]=temp2[:]
                        if A.shape[0]-2-i == 0:
                                A[0,:]=temp1
                                break
                        elif A.shape[0]-2-i > 0:
                                temp2[:]=A[-i-2,:]
        elif fac1==0:
                A=A

        if fac2>0:
                temp3=np.array([A[:,0]])
                temp4=np.array([A[:,1]])
                for i in range(A.shape[1]):
                        A[:,i+1]=temp3
                        temp3[:]=temp4[:]
                        if (i+2)>(A.shape[1]-1):
                                A[:,0]=temp3
                                break
                        elif ((i+2)<(A.shape[1]-1)) or ((i+2) == (A.shape[1]-1)):
                                temp4[:]=A[:,i+2]
        elif fac2<0:
                temp3=np.array([A[:,0]])
                temp4=np.array([A[:,-1]])
                for i in range(A.shape[1]):
                        A[:,-i-1]=temp3
                        temp3[:]=temp4[:]
                        if A.shape[1]-2-i == 0:
                                A[:,0]=temp3
                                break
                        elif A.shape[1]-2-i > 0:
                                temp4[:]=A[:,-i-2]
        elif fac2==0:
                A=A

        return A
