
import numpy as np

def traj_evaluate(t0,tf,currT,P0,PF):
        Amat = np.array([[1, t0, t0**2, t0**3, t0**4, t0**5],
        [0, 1, 2*t0, 3*t0**2, 4*t0**3, 5*t0**4],
        [0, 0, 2, 6*t0, 12*t0**2, 20*t0**3],
        [1, tf, tf**2, tf**3, tf**4, tf**5],
        [0, 1, 2*tf, 3*tf**2, 4*tf**3, 5*tf**4],
        [0, 0, 2, 6*tf, 12*tf**2, 20*tf**3]])

        xCoeffs = np.matmul(np.linalg.inv(Amat), np.array([[P0[0]],[0],[0],[PF[0]],[0],[0]]))
        yCoeffs = np.matmul(np.linalg.inv(Amat), np.array([[P0[1]],[0],[0],[PF[1]],[0],[0]]))
        zCoeffs = np.matmul(np.linalg.inv(Amat), np.array([[P0[2]],[0],[0],[PF[2]],[0],[0]]))
        coEffs = np.concatenate((xCoeffs,yCoeffs,zCoeffs),axis=1)

        A = np.array([[1, currT, currT**2, currT**3, currT**4, currT**5],
            [0, 1, 2*currT, 3*currT**2, 4*currT**3, 5*currT**4],
            [0, 0, 2, 6*currT, 12*currT**2, 20*currT**3]])

        traj = np.zeros([3,3])
        traj[:,0] = np.matmul(A,coEffs[:,0])
        traj[:,1] = np.matmul(A,coEffs[:,1])
        traj[:,2] = np.matmul(A,coEffs[:,2])
        return traj


P0 = np.array([0,0,0])
P1 = np.array([0,0,1])

traj = traj_evaluate(0,5,5,P0,P1)

print(traj)