      module  constants
        implicit none

      real, parameter :: N=1E7 ! total number of bacteria in largeintestine

        character(22) :: filename,abc,ab,filename1

        !real,parameter  ::l=250.0  !Total length of large intestine (in
        !cm) 
        real,parameter  ::l=1.0
        real,parameter::c=0.008 !drift of protein in cm^2 per hour  
        integer, parameter :: t=100*24*3600 !Total time in seconds
        real, parameter ::dt=1.0
        !real, parameter :: dx=0.04
        real, parameter :: dx=1.0
        integer,parameter :: sp=int(l/dx)
        integer, parameter ::tp=int(t/dt)
        real, parameter :: bt=(N)/float(sp)
                !!!!!!!Giving controlled diet 86g of protein, 81g of fat
                !and
                !373g of carbohydrates !!!!!!!!!!!!
        real,parameter :: p0=1.0! in hecto gram
        real,parameter :: f0=0.45 ! in hecto gram 
        real,parameter :: c0=1.0 ! in hectogram
        real,parameter :: fi0=2.0

        real:: D_b=0.0 !Difusion coefficient of bacteri

        integer,parameter :: ts1=int(3600.0/1.0)

                real,parameter :: delta2= 0.2E-12*0.25 !Nutrientconsumption for host
                real,parameter :: delta= 0.2E-15*0.25  !Nutrientconsumption per bacteria per second 
                real,parameter :: delta1=1.0E-12 ! mass of a bacteria
        real,parameter :: death=0.0006*0.005*0.05
        real,parameter :: K_s = 0.1 ! concentration of nutrient when the growth rate is half of its maximum value in mili molar.
        real ::cfl
        integer,parameter :: tnob=12,tnon=4
        real,parameter :: k = 1.0E08, k1=1.0E07, k3=1.0E06,k9=5.0E06
        real,parameter :: ae=1000
        integer, parameter :: ub=12, lb=1,ub1=3,lb1=1
        real, parameter :: sc=6.0
        !real, parameter :: sc=1.0
 
        real, parameter  :: r1=0.05
        real,parameter :: st = 0.1 
        real,parameter :: ctm = 0.06
        real,parameter :: wm=0.0
        integer,parameter :: tnos = 120000
      end module constants

