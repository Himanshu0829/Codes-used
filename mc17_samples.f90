        implicit none
        real, parameter :: N=1E6 ! total number of bacteria in large intestine
        
        character(22) :: filename,abc,ab,filename1
        
        !real,parameter  ::l=250.0  !Total length of large intestine (in cm) 
        real,parameter  ::l=1.0 
        real,parameter::c=0.008 !drift of protein in cm^2 per hour  
        integer, parameter :: t=50*24*3600 !Total time in seconds
        real, parameter ::dt=1.0
        !real, parameter :: dx=0.04
        real, parameter :: dx=1.0
        integer,parameter :: sp=int(l/dx)
        integer, parameter ::tp=int(t/dt)
        real, parameter :: bt=(N)/float(sp)
                !!!!!!!Giving controlled diet 86g of protein, 81g of fat and
                !373g of carbohydrates !!!!!!!!!!!!
        real,parameter :: p0=1.0! in hecto gram
        real,parameter :: f0=0.45 ! in hecto gram 
        real,parameter :: c0=1.0 ! in hectogram
        real,parameter :: fi0=2.0
 
        real:: D_b=0.0 !Difusion coefficient of bacteri
      
        integer,parameter :: ts1=int(2160000.0/1.0)
                
                real,parameter :: delta2= 0.2E-12*0.25 !Nutrient consumption for host
                real,parameter :: delta= 0.2E-15*0.25  !Nutrient consumption per bacteria per second 
                real,parameter :: delta1=1.0E-12 ! mass of a bacteria
        real,parameter :: death=0.0006*0.005*0.05
        real,parameter :: K_s = 0.1 ! concentration of nutrient when the growth rate is half of its maximum value in mili molar.
        real ::cfl
        integer,parameter :: tnob=12,tnon=4
        real,parameter :: k = 1.0E08, k1=1.0E07, k3=1.0E06,k9=5.0E06
        real,parameter :: ae=1000
        real, dimension(tnob,tnon) ::ds,pb,cb,fb,fib,ph,ch,fh,fih,pt,ct,ft,fit
        real,dimension(tnob,tnon) :: growth_p,growth_c,growth_f,growth_fi
        real,dimension(tnob,tnon) :: growth_pn,growth_cn,growth_fn,growth_fin
        real, dimension(tnob,tnob) :: bs,ec,bs_new,bs1,bs_ex
        real,dimension(tnob,sp):: b,b_new,br, br_new,br_newn,b_newn,b_com1,b_com2
        integer :: i, j,ss,ts,xp1,xm1,un,m,um,q,mm
        real, dimension(sp) :: x, c_p, c_f, c_c,c_fi,rn,ln, c_p_new, c_f_new,c_c_new,c_fi_new,&
         & c_p_newn, c_f_newn,c_c_newn,c_fi_newn
        real ::asso,p_con,c_con,f_con,fi_con,a,asso1,asso2,asso3,asso4
        real,dimension(tnob) ::alpha_p,alpha_c,alpha_f,alpha_fi,deathb,growth,growthf,k0,b_av,b_fi
        real,dimension(tnob) ::alpha_pn,alpha_cn,alpha_fn,alpha_fin,deathbn,growthn,growthfn,re,re1
        real :: bsum,t1,rmsd1,rmsd2
        integer, parameter :: ub=12, lb=1,ub1=3,lb1=1
        real :: ri3,rmsd   
        integer :: ri1,ri2,ri4,ri5,f
        real, parameter :: sc=6.0
        real, parameter  :: r1=0.05
     ! rr4=sqrt(rr1)/4.0
     ! rr5=sqrt(rr2)/4.0
        real ::rr0, rrr0, rr1 , rr2, rr3,rr4,rr5,rr6,re2,re3,btav,ree
        integer :: cim, acim,d,tss
        real:: e


        cim=0
        btav=0
        acim=0
        un=7   ! Initial value is 7 is because 6 is used to write on screen 
        um=2
         cfl=c*(dt/dx)
        open(101, file='bd.dat',action='read',status='old')
        !open(102, file='test_inter.dat',action='read',status='old')

        open(102, file='in_ma.dat',action='read',status='old')
        open(14,file='co_ma_ex.dat',action='read',status ='old')
        open(20,file='av_bi.dat',action='read',status ='old')
      !  open(12,file='pop.dat',action='read',status ='old')     
        !open(103, file='k_values.dat',action='read',status='old')
        !open(1, file='hi.txt',STATUS='NEW')     
        
        !open(200, file='rmsd.txt',action='write',STATUS='NEW')     
        open(200, file='int_ma_final.csv',action='write',STATUS='NEW')     
       open(201, file='int_ma_initial.csv',action='write',STATUS='NEW') 
        !open(3, file='matrix.txt',STATUS='NEW')     
        open(18, file='samples.txt',STATUS='NEW',action='write',position='append')

        !open(25, file='r_error.txt',STATUS='NEW',action='write',position='append')
       open(204, file='av_pop.csv',action='write',STATUS='NEW')
       open(205, file='final_pop.csv',action='write',STATUS='NEW')
       open(206, file='sample_input.dat',action='read',STATUS='old')


     ! close(12)
     do i = 1, tnob
        read(101,*) (ds(i,j), j=1, tnon)

        !write(*,*) k0(i)
        read(20,*) b_av(i)
        !read(206,*) b_fi(i)
     end do
        close(101)
        close(20)
     do i = 1, tnob

       read(102,*) (bs(i,j), j=1, tnob)
       read(14,*) (bs_ex(i,j), j=1, tnob)
     end do
       close(14)
        close(102)
           
     !do i=1,tnob
     !btav=btav+b_av(i)
     !end do

     do i = 1, tnob
  do j=1,tnob
       write(200,*) bs(i,j)
       write(201,*) bs_ex(i,j)
 end do
     end do

        close(200)
        close(201)



 !do i=1,tnob
    read(206,*) (b_fi(i),i=1,tnob)
 ! end do
 close(206)

do i=1,tnob
write(205,*) b_fi(i)
write(204,*) b_av(i)

end do


close(204)
close(205)












do f=1,tnob

write(6,*) f
!do i = 1, tnob

 !       write(6,*) (bs(i,j), j=1, tnob)

  !   end do




        open(12,file='sample_input.dat',action='read',status ='old')     
        !open(103, file='k_values.dat',action='read',status='old')
       do ss=1,sp      
       read(12,*) (b(i,ss), i=1,tnob)
        end do
       close(12)
             do ss=1,sp
                x(ss)=ss*dx
           
                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do           


!do m=1,10000

rr3=1.0


ts1loop:  do ts=1,tp
!write(6,*) ts
         ri1=int((ub-lb+1)*rand()+lb)
         ri2=int((ub-lb+1)*rand()+lb)
         ri4=int((ub-lb+1)*rand()+lb)
         ri5=int((ub-lb+1)*rand()+lb)
         ri3=int((ub1-lb1+1)*rand()+lb1)-2
        !mod(ts,ts1)
        if(mod(ts,21600) .eq. 0) then
     
               do ss=1,sp
                                    
                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0

               end do
        end if 
          


    do ss=1,sp

        bsum=0 
        do i=1,tnob
        bsum=bsum+b(i,ss)
       ! b_new(i,ss)=0
        end do

    
       do i=1,tnob
          pb(i,1) = delta*(c_p(ss)/(k_s+c_p(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          cb(i,2) = delta*(c_c(ss)/(k_s+c_c(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fb(i,3) = delta*(c_f(ss)/(k_s+c_f(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fib(i,4) =delta*(c_fi(ss)/(k_s+c_fi(ss)))*(1.0/(exp(-ds(i,4)+1.0)))

          alpha_p(i)=pb(i,1)/delta1
          alpha_c(i)=cb(i,2)/delta1
          alpha_f(i)=fb(i,3)/delta1
          alpha_fi(i)=fib(i,4)/delta1
         
          ph(i,1) =delta2*(c_p(ss)/(k_s+c_p(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          ch(i,2) =delta2*(c_c(ss)/(k_s+c_c(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fh(i,3) =delta2*(c_f(ss)/(k_s+c_f(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fih(i,4)=delta2*(c_fi(ss)/(k_s+c_fi(ss)))*(1.0/(exp(-ds(i,4)+1.0)))

          pt(i,1)  =   pb(i,1)+ph(i,1)
          ct(i,2)  =   cb(i,2)+ch(i,2)
          ft(i,3)  =   fb(i,3)+fh(i,3)
          fit(i,4) =   fib(i,4)+fih(i,4)

          asso=0

         
          asso1=0
          asso2=0
     
               do j=1,tnob
                    if (bs(i,j) < 0.0) then
                       asso1=asso1+bs(i,j)*b(j,ss)
                    end if
    
                    if (bs(i,j) > 0.0 ) then
                       asso2=asso2+bs(i,j)*b(j,ss)
                    end if         
               end do
          
          growth_p(i,1)=alpha_p(i)*dt*b(i,ss)*(1+(-sc*b(i,ss)+asso1+asso2)/k)
          growth_c(i,2)=alpha_c(i)*dt*b(i,ss)*(1+(-sc*b(i,ss)+asso1+asso2)/k)
          growth_f(i,3)=alpha_f(i)*dt*b(i,ss)*(1+(-sc*b(i,ss)+asso1+asso2)/k)
          growth_fi(i,4)=alpha_fi(i)*dt*b(i,ss)*(1+(-sc*b(i,ss)+asso1+asso2)/k)
         
          
           deathb(i) = death*dt*(b(i,ss))
           growth(i) = growth_p(i,1)+growth_c(i,2)+growth_f(i,3)+growth_fi(i,4)
        
           growthf(i)=growth(i)
           b_new(i,SS)=b(i,ss)+growthf(i)-deathb(i)
       
      end do


                            p_con=0.0
                            c_con=0.0
                            f_con=0.0
                            fi_con=0.0
                
                       do i=1,tnob
                            p_con=p_con+(b(i,ss))*pt(i,1)
                            c_con=c_con+(b(i,ss))*ct(i,2)
                            f_con=f_con+(b(i,ss))*ft(i,3)
                            fi_con=fi_con+(b(i,ss))*fit(i,4)
                             
 !                           b(i,ss)=b_new(i,ss)                            
!
                       end do
                
                            C_P_new(ss)=C_P(ss)-p_con
                            C_c_new(ss)=C_c(ss)-c_con
                            C_f_new(ss)=c_f(ss)-f_con
                            C_fi_new(ss)=c_fi(ss)-fi_con
                
               end do
                
                
                    do ss=1,sp   
                       C_P(ss)=c_P_new(ss)
                        C_C(ss)=c_C_new(ss)
                        C_F(ss)=c_F_new(ss)
                       C_Fi(ss)=c_Fi_new(ss)
                      end do
         do ss=1,sp
         do i=1,tnob
            if(i .eq. f) then
               b(i,ss)=0.0
            else
               b(i,ss)=b_new(i,ss)
            end if
 
         end do
         end do



         d=mod(ts,4320000)
         if (d .eq. 0) then
         !if (ts .eq. 2150000) then  
                do ss=1,sp

                    write(18,*)(b(i,ss), i=1,tnob)
                end do
         end if
end do ts1loop


end do
close(18)
end
 
