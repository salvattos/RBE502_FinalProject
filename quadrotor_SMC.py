#!/usr/bin/env python3
from math import pi, sqrt, atan2, cos, sin
from turtle import position
import numpy as np
from numpy import NaN
import rospy
import tf
from std_msgs.msg import Empty, Float32
from nav_msgs.msg import Odometry
from mav_msgs.msg import Actuators
from geometry_msgs.msg import Twist, Pose2D
import pickle
import os

class Quadrotor():
    def __init__(self):
        # publisher for rotor speeds
        self.motor_speed_pub = rospy.Publisher("/crazyflie2/command/motor_speed", Actuators, queue_size=10)
        # subscribe to Odometry topic
        self.odom_sub = rospy.Subscriber("/crazyflie2/ground_truth/odometry",Odometry, self.odom_callback)
        self.t0 = None
        self.t = None
        self.t_series = []
        self.x_series = []
        self.y_series = []
        self.z_series = []
        self.mutex_lock_on = False
        self.inTransit = True
        self.U = np.zeros((4,1))
        rospy.on_shutdown(self.save_data)
        # TODO: include initialization codes if needed

    def traj_evaluate(self,t0,tf,currT,P0,PF):
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

    def sat(self,s,boundary):
        return min(max(s/boundary,-1),1)


    def smc_control(self, xyz, xyz_dot, rpy, rpy_dot):
        #Known Variables
        m,L,Ix,Iy,Iz = .027, .046, 16.571710E-6, 16.571710E-6, 29.261652E-6
        Ip,Kf,Km,Wmax,Wmin,g = 12.65625E-8, 1.28192E-8, 5.964552E-3, 2618, 0, 9.81
        x,y,z = 0,1,2
        phi,theta,psi = 0,1,2

        #Set Gains
        kp = 100
        kd = 10
        K = np.array([5,140,140,25])
        Lam = np.array([2,13,13,5])
        boundary = np.array([.1,.1,.1,.1])

        # obtain the desired values by evaluating the corresponding trajectories
        P0 = np.array([0,0,0])
        P1 = np.array([0,0,1])
        if(self.inTransit):
            self.inTransit = False
            self.endTime = self.t0 + 10
            self.startTime = self.t0
        currT = self.t + self.t0

        desiredPts = np.zeros((3,6))
        desiredPts[:,0:3] = self.traj_evaluate(self.startTime,self.endTime,currT,P0,P1)

        #Calculate Omega based off of previous U input
        #Allocation Matrix
        allocMat = np.array([[1/(4*Kf*L), -sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), -1/(4*Km*Kf)],
                    [1/(4*Kf*L), -sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), 1/(4*Km*Kf)],
                    [1/(4*Kf*L), sqrt(2)/(4*Kf*L), sqrt(2)/(4*Kf*L), -1/(4*Km*Kf)],
                    [1/(4*Kf*L), sqrt(2)/(4*Kf*L), -sqrt(2)/(4*Kf*L), 1/(4*Km*Kf)]])

        # TODO: maintain the rotor velocities within the valid range of [0 to 2618]
        Wdesired = np.matmul(allocMat,self.U)
        Wdesired = np.clip(Wdesired,Wmin,Wmax)
        Omega = Wdesired[0] - Wdesired[1] + Wdesired[2] - Wdesired[3]

        #Z Control Law
        eZ = np.array([[desiredPts[0,z]-xyz[z]],
                        [desiredPts[1,z]-xyz_dot[z]]])
        #print(desiredPts)
        #print(eZ[1][0])
        sZ = eZ[1][0] + Lam[0]*eZ[0][0] 
        satZ = self.sat(sZ,boundary[0])
        UrZ = K[0] * satZ
        self.U[0] = (m/(np.cos(phi)*np.cos(theta)))*(desiredPts[2,z]+g+Lam[0]*eZ[1][0] + UrZ)


        #Calculate X and Y forces
        Fx = m*(-kp*(xyz[x]-desiredPts[0,0]) - kd*(xyz_dot[x]-desiredPts[1,1]) + desiredPts[2,1])
        Fy = m*(-kp*(xyz[y]-desiredPts[0,1]) - kd*(xyz_dot[y]-desiredPts[1,2]) + desiredPts[2,2])
        print(self.U[0])

        #Update desired points
        thetaDesired = np.arcsin(Fx/self.U[0])
        phiDesired = np.arcsin(-Fy/self.U[0])
        psiDesired = 0
        print(desiredPts[0,3:6])
        desiredPts[0,3:6] = [phiDesired,thetaDesired,psiDesired]
        print(desiredPts[0,3:6])
        # TODO: implement the Sliding Mode Control laws designed in Part 2 to calculate the control inputs "u"
        # REMARK: wrap the roll-pitch-yaw angle errors to [-pi to pi]
        # TODO: convert the desired control inputs "u" to desired rotor velocities "motor_vel" by using the "allocation matrix"

        # publish the motor velocities to the associated ROS topic
        motor_speed = Actuators()
        motor_speed.angular_velocities = [Wdesired[0], Wdesired[1], Wdesired[2], Wdesired[3]]

        rospy.loginfo(self.t)
        #rospy.loginfo(self.startTime)
        #rospy.loginfo(self.endTime)
        #print(xyz)
        self.motor_speed_pub.publish(motor_speed)

    # odometry callback function (DO NOT MODIFY)

    def odom_callback(self, msg):
        if self.t0 == None:
            self.t0 = msg.header.stamp.to_sec()
        self.t = msg.header.stamp.to_sec() - self.t0
        # convert odometry data to xyz, xyz_dot, rpy, and rpy_dot
        w_b = np.asarray([[msg.twist.twist.angular.x], [msg.twist.twist.angular.y], [msg.twist.twist.angular.z]])
        v_b = np.asarray([[msg.twist.twist.linear.x], [msg.twist.twist.linear.y], [msg.twist.twist.linear.z]])
        xyz = np.asarray([[msg.pose.pose.position.x], [msg.pose.pose.position.y], [msg.pose.pose.position.z]])
        q = msg.pose.pose.orientation
        T = tf.transformations.quaternion_matrix([q.x, q.y, q.z, q.w])
        T[0:3, 3] = xyz[0:3, 0]
        R = T[0:3, 0:3]
        xyz_dot = np.dot(R, v_b)
        rpy = tf.transformations.euler_from_matrix(R, 'sxyz')
        rpy_dot = np.dot(np.asarray([
                [1, np.sin(rpy[0])*np.tan(rpy[1]), np.cos(rpy[0])*np.tan(rpy[1])],
                [0, np.cos(rpy[0]), -np.sin(rpy[0])],
                [0, np.sin(rpy[0])/np.cos(rpy[1]), np.cos(rpy[0])/np.cos(rpy[1])]]), w_b)
        rpy = np.expand_dims(rpy, axis=1)
        # store the actual trajectory to be visualized later
        if (self.mutex_lock_on is not True):
            self.t_series.append(self.t)
            self.x_series.append(xyz[0, 0])
            self.y_series.append(xyz[1, 0])
            self.z_series.append(xyz[2, 0])
        # call the controller with the current states
        self.smc_control(xyz, xyz_dot, rpy, rpy_dot)

    # save the actual trajectory data
    def save_data(self):
        # TODO: update the path below with the correct path
        with open("/home/salvattos/rbe502_project/src/project/scripts/log.pkl","wb") as fp:
            self.mutex_lock_on = True
            pickle.dump([self.t_series,self.x_series,self.y_series,self.z_series], fp)

if __name__ == '__main__':
    rospy.init_node("quadrotor_control")
    rospy.loginfo("Press Ctrl + C to terminate")
    rospy.Rate(100)
    whatever = Quadrotor()
    whatever.inTransit = True
    try:
        rospy.spin()
    except KeyboardInterrupt:
        rospy.loginfo("Shutting down")