Main240305 file is the main program, please open this file for simulation
Author: Ruixin Miao
Data : 24.3.5
The following is a brief explanation of the functions of each part of the code. The specific IO has complete comments in the program. If there is anything that is not detailed, you can consult me.
email: mrxaiwu@gmail.com
The following is a description of the inspection of each part:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Changes in magnetization vector under on-resonance
%%%%%%%%%%%%  Initialization parameters
ni： the ordinal number of the magnetization vector of each component
gamma ： the magnetic gyroscopic ratio of hydrogen protons
forwarding.B0.Hsum ： the magnetic field strength of the static magnetic field |B0|
B1：RF field strength
T ： represents the time interval of the calculation
%%%%%%%%%%% IntegralM0_GM : the magnetization vector calculated by the GM method after the pulse excitation time t
[M0] = IntegralM0_GM(forwarding,B1,ni,t)
Input： 
forwarding：Forward structure variables
forwarding.B0：static magnetic field
forwarding.TLoop ：RF Coil
B1 ：RF field strength
ni : the ordinal number of the magnetization vector of each component
t：Pulse excitation time
M0:Magnetization vector structure
M0.rx：x-direction magnetization vector
M0.ry：y-direction magnetization vector
M0.rz：z-direction magnetization vector
M0.rxy：Transverse magnetization vector
M0.rsum: The magnitude of the magnetization vector
%%%%%%%%%%% IntegralM0_RM: the magnetization vector calculated by the RM method after the pulse excitation time t
[M0] = IntegralM0_GM(forwarding,B1,ni,t)
Input： 
forwarding：Forward structure variables
forwarding.B0：static magnetic field
forwarding.TLoop ：RF Coil
B1 ：RF magnetic field
ni : the ordinal number of the magnetization vector of each component
t：Pulse excitation time
M0:Magnetization vector structure
M0.rx：x-direction magnetization vector
M0.ry：y-direction magnetization vector
M0.rz：z-direction magnetization vector
M0.rxy：Transverse magnetization vector
M0.rsum: The magnitude of the magnetization vector
%%%%%%%%%%% Calculation results
MAPE_X：MAPE of magnetization vector X
MAPE_Y：MAPE of magnetization vector Y
MAPE_Z：MAPE of magnetization vector Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.Changes of each magnetization vector under the condition of multi-component Larmor frequency
%%%%%%%%%%%%  Initialization parameters
tdd : Pulse excitation time
gamma : The magnetic gyroscopic ratio of hydrogen protons
B0Hsum : Magnetic field strength at the center of static magnetic field
B1 : RF field strength
NM : Number of Larmor frequencies
N : Avogadro constant 
hq：Planck's constant/2*pi [J.s]
K：Boltzmann's constant  [J/K]
T0：absolute temperature  [K]
%%%%%%%%%%% for loop
1. Calculate the magnetization vector at different Larmor frequencies and compare the RM and GM algorithms;
2. In the loop process, first set the range of alpha between [-pi/3, pi/3], because excessive off-resonance contributes little to the calculation of the signal. If the reader is interested, this value can be changed;
3. Calculate each component B0 through different alpha pairs to obtain NM(i) Larmor frequencies;
4. In the calculation process, save the calculation results and time consumption of the two algorithms for different Larmor frequencies;
5. mape_x,mape_y,mape_z：MAPE of the three components of the magnetization vector
6. t1，t2 : the Time Consuming of GM and RM
%%%%%%%%%%% IntegralM0_GM_V1 : the magnetization vector calculated by the GM method after the pulse excitation time t
The programming method of matrix calculation is used to reduce the calculation time,
 but it increases the memory usage. The principle is the same as IntegralM0_GM, so I will not go into details.
%%%%%%%%%%% IntegralM0_RM_V1 : the magnetization vector calculated by the RM method after the pulse excitation time t
The programming method of matrix calculation is used to reduce the calculation time, 
but it increases the memory usage. The principle is the same as IntegralM0_GM, so I will not go into details.

