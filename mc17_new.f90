        program main 
 
           use constants
        implicit none 
       
        real,dimension(tnob,tnon)::ds,pb,cb,fb,fib,ph,ch,fh,fih,pt,ct,ft,fit
        real,dimension(tnob,tnon) ::growth_p,growth_c,growth_f,growth_fi
        real,dimension(tnob,tnon) ::growth_pn,growth_cn,growth_fn,growth_fin
        real, dimension(tnob,tnob) :: bs,ec,bs_new,bs1,corr_bu,corr_au,bs2
        real,dimension(tnob,sp)::b,b_new,br,br_new,br_newn,b_newn,b_com1,b_com2,b_new1
        integer :: i, j,ss,ts,xp1,xm1,un,m,um,q,mm
        real, dimension(sp) :: x, c_p,c_f,c_c,c_fi,rn,ln,c_p_new,c_f_new,c_c_new,c_fi_new,&
         & c_p_newn, c_f_newn,c_c_newn,c_fi_newn
        real ::asso,p_con,c_con,f_con,fi_con,a,asso1,asso2,asso3,asso4
        real,dimension(tnob)::alpha_p,alpha_c,alpha_f,alpha_fi,deathb,growth,growthf,k0,b_av
        real,dimension(tnob)::alpha_pn,alpha_cn,alpha_fn,alpha_fin,deathbn,growthn,growthfn,re,re1
        real :: bsum,t1,rmsd1,rmsd2,rmsd3,rmsd4
        real :: ri3,rmsd
        integer :: ri1,ri2,ri4,ri5
        real ::rr0, rrr0, rr1 , rr2, rr3,rr4,rr5,rr6,re2,re3,btav,ree
        integer :: cim, acim,d,tss

        cim=0
        btav=0
        acim=0
        un=7   ! Initial value is 7 is because 6 is used to write on screen 
        um=2
        cfl=c*(dt/dx)
        open(101, file='bd.dat',action='read',status='old')
        open(102, file='bs4.dat',action='read',status='old')
        open(20,file='av_bi.dat',action='read',status ='old')
       
        open(1, file='hi.txt',STATUS='NEW')     
        open(200, file='rmsd.txt',action='write',STATUS='NEW')     
        open(37, file='matrix.txt',STATUS='NEW')     
        open(18, file='pop.txt',STATUS='NEW',action='write',position='append')
        open(25, file='r_error.txt',STATUS='NEW',action='write',position='append')
        open(455, file='matrix_element.txt',status='new',action='write',position='append')
       
      do ss=1,sp
       do i=1,tnob
         b(i,ss)=(bt/tnob)
       end do
      end do

     do i = 1, tnob
        read(101,*) (ds(i,j), j=1, tnon)
        read(20,*) b_av(i)

     end do
        close(101)
        close(20)
     do i = 1, tnob

        read(102,*) (bs(i,j), j=1, tnob)

     end do
       
        close(102)

  do i=1 ,tnob
write(6,*) b_av(i)
 end do   
  

do i=1,tnob
     btav=btav+b_av(i)
     end do

              do ss=1,sp
                x(ss)=ss*dx
           
                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do           

!do m=1,1000
!write(6,*) m
rr3=1.0

open(19, file='pop_bu.txt',STATUS='NEW',action='write',position='append')


ts1loop:  do ts=1,tp

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

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
                              rr0=0.0
                              do i=1,tnob
                     
                                rr0=rr0+((b_new(i,ss)-b(i,ss)))**2
                     
                     
                              end do
                       
                                rrr0=sqrt(rr0)/float(tnob)
                             
                                  
                               if(rrr0 .lt. abs(0.001)) then        
                     
                
!                                     do ss=1,sp
                     write(19,*)(b(i,ss), i=1, tnob )
                                     !   end do
        
            
                                    exit ts1loop
    !                                tss=tp                     
                     
                               else 
                     
                              
                             
                                     do i=1,tnob
                     
                                          b(i,ss)=b_new(i,ss)
                                    end do
                      
                              end if   



                            p_con=0.0
                            c_con=0.0
                            f_con=0.0
                            fi_con=0.0
                
                       do i=1,tnob
                            p_con=p_con+(b(i,ss))*pt(i,1)
                            c_con=c_con+(b(i,ss))*ct(i,2)
                            f_con=f_con+(b(i,ss))*ft(i,3)
                            fi_con=fi_con+(b(i,ss))*fit(i,4)
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

          d=mod(ts,ts1)
         if (d .eq. 0) then

                do ss=1,sp

                    write(18,*)m,ts,(b(i,ss), i =1,tnob)
                end do
         end if


end do ts1loop



open(69, file='bs_bu.txt',STATUS='NEW',action='write',position='append')




   do i=1,tnob


  
    write(69,*) ( bs(i,j) , j=1,tnob)
   end do

!close(69)

close(69, status = 'delete')
!close(19)
close(19, status = 'delete')

rmsd1=0.0

 do i=1,tnob
  

      rmsd1=rmsd1+((b(i,ss)-b_av(i))/b_av(i))**2

  !    rmsd2=rmsd2+((b_newn(i,ss)-b_av(i))/b_av(i))**2


    end do

      rmsd1=sqrt(rmsd1/float(tnob))

  !    rmsd2=sqrt(rmsd2/float(tnob))

!call samples_bu()

!call p_corr_bu(corr_bu)


!call corr_rmsd_bu(rmsd3)


do m=1, tnos
 
write(6,*) m


     do i=1,tnob
        do j=1,tnob
           bs2(i,j)=bs(i,j)
        end do
    end do








         ri1=int((ub-lb+1)*rand()+lb)
         ri2=int((ub-lb+1)*rand()+lb)
         ri4=int((ub-lb+1)*rand()+lb)
         ri5=int((ub-lb+1)*rand()+lb)
         ri3=int((ub1-lb1+1)*rand()+lb1)-2

        
                           if (bs2(ri4,ri5) > 0) then
                           bs1(ri4,ri5)=bs2(ri4,ri5)+((rand()-0.5)*ctm)
                           !bs1(ri5,ri4)=bs(ri5,ri4)+((cor_t-0.5)*0.1)

                          ! bs1(ri4,ri5)=bs1(ri5,ri4)


!                             bs1(ri4,ri5)=bs(ri4,ri5)+ri3*r1

                           else if (bs2(ri4,ri5) < 0) then
                 !
                             bs1(ri4,ri5)=bs2(ri4,ri5)+((rand()-0.5)*ctm)
             
                !             bs1(ri4,ri5)=bs(ri4,ri5)+ri3*r1
                  !
                           end if
                             do i=1,tnob
                             
                                do j=1, tnob

                                 if ((i .eq. ri4) .and. (j .eq. ri5) ) then
                                      
                                     bs_new(i,j)=bs1(ri4,ri5)
                                                                          
        !                        else if ((i .eq. ri5) .and.  (j .eq. ri4)) then

         !                            bs_new(i,j)=bs1(ri4,ri5)
!
                                 else      ! this is not working
                                
                                     bs_new(i,j)=bs2(i,j)
                                           
                                 end if  
                             
                                end do
                          
                            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!Calculation with updated matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(file='pop.txt', action='read',position='append',status='old')

open(29, file='pop_au.txt',STATUS='NEW',action='write',position='append')

open(89, file='bs_au.txt',STATUS='NEW',action='write',position='append')

   do i=1,tnob
    write(89,*) ( bs_new(i,j) , j=1,tnob)
   end do




! do i=1,tnob
 !   write(6,*) ( bs_new(i,j) , j=1,tnob)
 !end do

close(89, status = 'delete')
!close(89)

 do ss=1,sp
!                x(ss)=ss*dx

                C_P_new(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F_new(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C_new(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi_new(ss)=fi0
              end do




   do ss=1,sp
        do i=1,tnob
         b_new1(i,ss)=b(i,ss)
        end do
    end do


ts2loop:           do ts=1,tp
!write(6,*) 'hi'

              if(mod(ts,21600) .eq. 0) then

               do ss=1,sp

                C_P_new(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F_new(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C_new(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi_new(ss)=fi0

               end do
            end if
       !do ss=,sp
       ! do i=1,tnob
        ! b_new(i,ss)=b(i,ss)
        !end do
       !end do

  do ss=1,sp
 
     do i=1,tnob
  
          pb(i,1) = delta*(c_p_new(ss)/(k_s+c_p_new(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          cb(i,2) = delta*(c_c_new(ss)/(k_s+c_c_new(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fb(i,3) = delta*(c_f_new(ss)/(k_s+c_f_new(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fib(i,4) =delta*(c_fi_new(ss)/(k_s+c_fi_new(ss)))*(1.0/(exp(-ds(i,4)+1.0)))

          alpha_p(i)=pb(i,1)/delta1
          alpha_c(i)=cb(i,2)/delta1
          alpha_f(i)=fb(i,3)/delta1
          alpha_fi(i)=fib(i,4)/delta1

          ph(i,1) =delta2*(c_p_new(ss)/(k_s+c_p_new(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          ch(i,2) =delta2*(c_c_new(ss)/(k_s+c_c_new(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fh(i,3) =delta2*(c_f_new(ss)/(k_s+c_f_new(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fih(i,4)=delta2*(c_fi_new(ss)/(k_s+c_fi_new(ss)))*(1.0/(exp(-ds(i,4)+1.0)))
!
          pt(i,1)  =   pb(i,1)+ph(i,1)
          ct(i,2)  =   cb(i,2)+ch(i,2)
          ft(i,3)  =   fb(i,3)+fh(i,3)
          fit(i,4) =   fib(i,4)+fih(i,4)

                       ! do i=1,tnob   
                              asso3=0
                              asso4=0

                              do j=1,tnob
                                 if (bs_new(i,j) < 0.0) then
                                    asso3=asso3+bs_new(i,j)*b_new1(j,ss)
                                 end if
                                 if (bs_new(i,j) > 0.0 ) then
                                    asso4=asso4+bs_new(i,j)*b_new1(j,ss)
                                 end if
                              end do
                        !  end do
          growth_pn(i,1)=alpha_p(i)*dt*b_new1(i,ss)*(1+(-sc*b_new1(i,ss)+asso3+asso4)/k)
          growth_cn(i,2)=alpha_c(i)*dt*b_new1(i,ss)*(1+(-sc*b_new1(i,ss)+asso3+asso4)/k)
          growth_fn(i,3)=alpha_f(i)*dt*b_new1(i,ss)*(1+(-sc*b_new1(i,ss)+asso3+asso4)/k)
          growth_fin(i,4)=alpha_fi(i)*dt*b_new1(i,ss)*(1+(-sc*b_new1(i,ss)+asso3+asso4)/k)


          deathbn(i) = death*dt*(b_new1(i,ss))
          growthn(i) = growth_pn(i,1)+growth_cn(i,2)+growth_fn(i,3)+growth_fin(i,4)

         growthfn(i)=growthn(i)
         b_newn(i,SS)=b_new1(i,ss)+growthfn(i)-deathbn(i)
  !write(6,*) b_newn(i,ss)     
       end do  

!!!!!!!!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!

                                    rr1=0.0
                                    do i=1,tnob

                                      rr1=rr1+((b_newn(i,ss)-b_new1(i,ss)))**2

                                    end do
  
                                      rr1=sqrt(rr1)/float(tnob)
                                   
                                        
                                    if(rr1 .lt. abs(0.001)) then        
                     
     write(29,*)(b_newn(i,ss), i=1,tnob)
                                                        
                                       exit ts2loop 

                                    else

                                       do i=1,tnob
                                          b_new1(i,ss)=b_newn(i,ss)
                                         end do

                                    end if  
        
                            p_con=0.0
                            c_con=0.0
                            f_con=0.0
                            fi_con=0.0
                
                       do i=1,tnob
                            p_con=p_con+(b_new1(i,ss))*pt(i,1)
                            c_con=c_con+(b_new1(i,ss))*ct(i,2)
                            f_con=f_con+(b_new1(i,ss))*ft(i,3)
                            fi_con=fi_con+(b_new1(i,ss))*fit(i,4)
                       end do
                
                            C_P_newn(ss)=C_P_new(ss)-p_con
                            C_c_newn(ss)=C_c_new(ss)-c_con
                            C_f_newn(ss)=c_f_new(ss)-f_con
                            C_fi_newn(ss)=c_fi_new(ss)-fi_con
                
                  end do
                
                
                    do ss=1,sp   
                        C_P_new(ss)=c_P_newn(ss)
                        C_C_new(ss)=c_C_newn(ss)
                        C_F_new(ss)=c_F_newn(ss)
                        C_Fi_new(ss)=c_Fi_newn(ss)
                     end do

         a =mod(ts,ts1)
         if (a .eq. 0) then
     
             do ss=1,sp

                    write(18,*) m,ts,(b_new1(i,ss),i=1,tnob)
        
             end do
     
         end if
    
   end do ts2loop

!close(29)
close(29, status = 'delete')




!call samples_au()

!call p_corr_au(corr_au)


!call corr_rmsd_au(rmsd4)





!!!!!!!!!!!!!!!!!!!!!!!!comparision of two RMSDs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !rmsd1=0
    rmsd2=0
   ! rr3=0
    
    do i=1,tnob
      !write(6,*) m, b(i,ss), b_newn(i,ss)
    
   !   rmsd1=rmsd1+((b(i,ss)-b_av(i))/b_av(i))**2
      
      rmsd2=rmsd2+((b_new1(i,ss)-b_av(i))/b_av(i))**2
      
    
    end do

      rmsd2=sqrt(rmsd2/float(tnob))

    cim=cim+1

 !write(6,*) m, rmsd1,rmsd2,cim,acim
 
 !           If (rmsd2 < rmsd1 .and.  rmsd4 < rmsd3) then
                 
          If (rmsd2 < rmsd1) then            
        !  If (rmsd4 < rmsd3) then            

                  acim=acim+1

          ! write(6,*) m, rmsd1, rmsd2,rmsd3,rmsd4, cim, acim
              write(6,*) m, rmsd1, rmsd2, cim, acim
          ! write(6,*) m, rmsd3, rmsd4, cim, acim
            

         ! write(200,*) m, rmsd1,rmsd2,rmsd3,rmsd4,cim,acim
           write(200,*) m, rmsd1,rmsd2,cim,acim
           !write(200,*) m, rmsd3,rmsd4,cim,acim
                  do i=1,tnob

                        do j =1,tnob
                           bs(i,j)=bs_new(i,j)
                        end do
                          ! b(i,ss)=b_newn(i,ss)
                  
                  end do

 
               do ss=1,sp

                    do i=1, tnob
   
                        b(i,ss)=b_newn(i,ss)
   
                    end do
   
                end do
             

               
               do i= 1, tnob
                    write(37,*) (bs(i,j), j=1, tnob)
               end do
 
           rmsd1=rmsd2
           rmsd3=rmsd4

          !  else 
     
           ! do ss=1,sp

            !   do i=1,tnob
             !     b(i,ss)=b_new(i,ss)
              ! end do   

          ! end do
!write(6,*) m
           
              end if 
                    
             
              !!!!!!!!!!!!!!!!!!!$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!


!write(45,*) m,b(1,2),b(1,4),b(2,1),b(2,3),b(3,2),b(3,4),b(4,1),b(4,3)
!write(45,*) m,bs(1,2),bs(1,4),bs(2,1),bs(2,3),bs(3,2),bs(3,4),bs(4,1),bs(4,3)

write(455,*)m,bs(1,1),bs(1,2),bs(1,3),bs(1,4),bs(2,2),bs(2,3),bs(2,4),bs(3,3),bs(3,4),bs(4,4),&
&bs(2,1),bs(3,1),bs(3,2),bs(4,1),bs(4,2),bs(4,3)


end do  ! number of steps
close(200)
close(445)
close(37)
end program main


subroutine samples_bu()
use constants

implicit none
real,dimension(tnob,tnon)::ds,pb,cb,fb,fib,ph,ch,fh,fih,pt,ct,ft,fit
        real,dimension(tnob,tnon) ::growth_p,growth_c,growth_f,growth_fi
        real,dimension(tnob,tnon) ::growth_pn,growth_cn,growth_fn,growth_fin
        real, dimension(tnob,tnob) :: bs,ec,bs_new,bs1
        real,dimension(tnob,sp):: b,b_bu,b_new,br,br_new,br_newn,b_newn,b_com1,b_com2
        integer :: i, j,ss,ts,xp1,xm1,un,m,um,q,mm
        real, dimension(sp) :: x,c_p,c_f,c_c,c_fi,rn,ln,c_p_new,c_f_new,c_c_new,c_fi_new,&
         & c_p_newn, c_f_newn,c_c_newn,c_fi_newn
        real ::asso,p_con,c_con,f_con,fi_con,a,asso1,asso2,asso3,asso4
        real,dimension(tnob)::alpha_p,alpha_c,alpha_f,alpha_fi,deathb,growth,growthf,k0,b_av
        real,dimension(tnob)::alpha_pn,alpha_cn,alpha_fn,alpha_fin,deathbn,growthn,growthfn,re,re1
        real :: bsum,t1,rmsd1,rmsd2
        real :: ri3,rmsd
        integer :: ri1,ri2,ri4,ri5
        real ::rr0, rrr0, rr1 , rr2, rr3,rr4,rr5,rr6,re2,re3,btav,ree
        integer :: cim, acim,d,tss,f
open(101, file='bd.dat',action='read',status='old')
open(21, file='pop_samples_bu.txt',STATUS='NEW',action='write',position='append')
open(69, file='bs_bu.txt',STATUS='old',action='read')

    do i = 1, tnob
       read(101,*) (ds(i,j), j=1, tnon)
     end do
        close(101)
      
      do i = 1, tnob
       read(69,*) (bs(i,j), j=1, tnob)
      end do
!     close(69,status ='delete')
     close(69)
!    do i = 1, tnob
!       write(6,*) (bs(i,j), j=1, tnob)
 !     end do



do f=1,tnob

     !write(6,*) f

            do ss=1,sp
                x(ss)=ss*dx

                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do

        open(12,file='pop_bu.txt',action='read',status ='old')
       
        do ss=1,sp
          read(12,*)(b(i,ss),i=1,tnob)
          !read(12,*)b(1,ss),b(2,ss),b(3,ss),b(4,ss),b(5,ss),b(6,ss),b(7,ss),&
           !&b(8,ss),b(9,ss),b(10,ss),b(11,ss),b(12,ss)
        end do

              !close(12)

  
             do ss=1,sp
                x(ss)=ss*dx

                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do

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
         !if (d .eq. 0) then
         if (ts .eq. 2160000) then  
                do ss=1,sp
                    write(21,*)(b(i,ss),i=1,tnob)
          !          write(21,*)b(1,ss),b(2,ss),b(3,ss),b(4,ss),b(5,ss),b(6,ss),b(7,ss),&
         !  &b(8,ss),b(9,ss),b(10,ss),b(11,ss),b(12,ss)
                    !write(21,*)b(1,ss),b(2,ss),b(3,ss),b(4,ss)
                end do
         end if
end do ts1loop

close(12)
end do

open(12,file='pop_bu.txt',action='read',status ='old')
close(12, status='delete')
close(21)
return
end subroutine samples_bu


subroutine p_corr_bu(corr_bu)
use constants
implicit none

         ! integer, parameter :: tnob=5
          real, dimension(tnob,tnob) :: r,rs,corr_bu
          integer :: i,j,q
          real,dimension(tnob) :: m,m_s,m_rs, sd ,m0
          real :: corr1,corr2
          open(55,file='pop_samples_bu.txt', action='read',status = 'old')
            open(33,file='pearson_bu.txt', action='write',status = 'new')

            !open(44,file='pearson1.txt', action='write',status = 'new')

          do i=1,tnob
             read(55,*) (r(i,j), j=1,tnob)
          end do
         ! close(55,status = 'delete')
           close(55)

          do i=1,tnob
             m(i)=0.0
             !sd(i)=0.0
             m_rs(i)=0.0
          end do

    !!!!!!!!!!!!!!!!mean!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i=1,tnob
            do j=1,tnob
               m(i)=m(i)+r(j,i)
            end do
           m(i)=m(i)/float(tnob)
          end do

      !     do i=1,tnob
!            write(6,*) m(i)
       !    end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!Standard deviation!!!!!!!!!!!!!!!!
            do i=1,tnob

               do j=1,tnob

                   rs(i,j)= r(i,j)**2

               end do

            end do

             do i=1,tnob
               do j=1,tnob
                m_rs(i)=m_rs(i)+rs(j,i)
            end do
               m_rs(i)=m_rs(i)/float(tnob)
          end do

           do i=1,tnob

             m_s(i)=m(i)**2

          end do

            do i=1,tnob
              sd(i)=sqrt(m_rs(i)-m_s(i))
            end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!Pearson correlation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
           do i=1,tnob
               do j=1,tnob
                  corr_bu(i,j)=0.0
               end do
           end do

 !          write(6,*) 'hi'
           corr1=0.0

           do i=1,tnob

            corr1=corr1+((r(i,1)-m(1))*(r(i,2)-m(2)))

          end do

           corr1=corr1/(sd(1)*sd(2)*tnob)

 !          write(6,*) corr1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          ! do k=1,tnob
               do i=1,tnob
                  do j=1, tnob
                     do q=1, tnob
                     corr_bu(i,j)=corr_bu(i,j)+((r(q,j)-m(j))*(r(q,i)-m(i)))
                     end do
                    ! write(6,*) corr2
                  end do
               end do
             ! end do 

               do i=1,tnob
                  do j=1,tnob
                    corr_bu(i,j)=corr_bu(i,j)/(sd(i)*sd(j)*tnob)
                  end do
               end do

           do i=1,tnob
             write(33,*) (corr_bu(i,j), j=1,tnob)
           end do

     !      do i=1,tnob
      !       do j=1, tnob
       !        write(44,*) corr_bu(i,j)
        !     end do
         !  end do

    close(33)
    !close(44)
    return
 end subroutine p_corr_bu




subroutine samples_au()
use constants

implicit none
 real,dimension(tnob,tnon)::ds,pb,cb,fb,fib,ph,ch,fh,fih,pt,ct,ft,fit
        real,dimension(tnob,tnon) ::growth_p,growth_c,growth_f,growth_fi
        real,dimension(tnob,tnon) ::growth_pn,growth_cn,growth_fn,growth_fin
        real, dimension(tnob,tnob) :: bs,ec,bs_new,bs1
        real,dimension(tnob,sp):: b,b_new,br,br_new,br_newn,b_newn,b_com1,b_com2
        integer :: i, j,ss,ts,xp1,xm1,un,m,um,q,mm
        real, dimension(sp) :: x,c_p,c_f,c_c,c_fi,rn,ln,c_p_new,c_f_new,c_c_new,c_fi_new,&
         & c_p_newn, c_f_newn,c_c_newn,c_fi_newn
        real ::asso,p_con,c_con,f_con,fi_con,a,asso1,asso2,asso3,asso4
        real,dimension(tnob)::alpha_p,alpha_c,alpha_f,alpha_fi,deathb,growth,growthf,k0,b_av
        real,dimension(tnob)::alpha_pn,alpha_cn,alpha_fn,alpha_fin,deathbn,growthn,growthfn,re,re1
        real :: bsum,t1,rmsd1,rmsd2
        real :: ri3,rmsd
        integer :: ri1,ri2,ri4,ri5
        real ::rr0, rrr0, rr1 , rr2, rr3,rr4,rr5,rr6,re2,re3,btav,ree
        integer :: cim, acim,d,tss,f
open(101, file='bd.dat',action='read',status='old')
open(21, file='pop_samples_au.txt',STATUS='NEW',action='write',position='append')

open(89, file='bs_au.txt',STATUS='old',action='read')
    do i = 1, tnob
       read(101,*) (ds(i,j), j=1, tnon)
     end do
        close(101)

     do i = 1, tnob
       read(89,*) (bs_new(i,j), j=1, tnob)

     end do
  !close(89)
     close(89, status ='delete')

  !  write(6,*) 'hi'

   !    do i = 1, tnob
    !   write(6,*) (bs_new(i,j), j=1, tnob)

     !end do




do f=1,tnob

    ! write(6,*) f

            do ss=1,sp
                x(ss)=ss*dx

                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do

        open(12,file='pop_au.txt',action='read',status ='old')
       do ss=1,sp
       read(12,*)(b_new(i,ss),i=1,tnob)
       !read(12,*)b_new(1,ss),b_new(2,ss),b_new(3,ss),b_new(4,ss),b_new(5,ss),b_new(6,ss),b_new(7,ss),&
        !   &b_new(8,ss),b_new(9,ss),b_new(10,ss),b_new(11,ss),b_new(12,ss)
        end do
!       close(12)
             do ss=1,sp
                x(ss)=ss*dx

                C_P(ss)=p0    !Concentration of Cnutrient in mili molar
                C_F(ss)=f0    !Concentration of nutrient Cin mili molar
                C_C(ss)=c0    !Concentration of nutrient Cin mili molar
                C_fi(ss)=fi0
              end do

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
                       asso1=asso1+bs_new(i,j)*b_new(j,ss)
                    end if

                    if (bs(i,j) > 0.0 ) then
                       asso2=asso2+bs_new(i,j)*b_new(j,ss)
                    end if
               end do

          growth_p(i,1)=alpha_p(i)*dt*b_new(i,ss)*(1+(-sc*b_new(i,ss)+asso1+asso2)/k)
          growth_c(i,2)=alpha_c(i)*dt*b_new(i,ss)*(1+(-sc*b_new(i,ss)+asso1+asso2)/k)
          growth_f(i,3)=alpha_f(i)*dt*b_new(i,ss)*(1+(-sc*b_new(i,ss)+asso1+asso2)/k)
          growth_fi(i,4)=alpha_fi(i)*dt*b_new(i,ss)*(1+(-sc*b_new(i,ss)+asso1+asso2)/k)


           deathb(i) = death*dt*(b_new(i,ss))
           growth(i) = growth_p(i,1)+growth_c(i,2)+growth_f(i,3)+growth_fi(i,4)

           growthf(i)=growth(i)
           b_newn(i,SS)=b_new(i,ss)+growthf(i)-deathb(i)

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
               b_new(i,ss)=0.0
            else
               b_new(i,ss)=b_newn(i,ss)
            end if

         end do
         end do
         !if (d .eq. 0) then
         if (ts .eq. 2160000) then  
                do ss=1,sp
           write(21,*)(b_new(i,ss),i=1,tnob)
           !write(21,*)b_new(1,ss),b_new(2,ss),b_new(3,ss),b_new(4,ss),b_new(5,ss),b_new(6,ss),b_new(7,ss),&
           !&b_new(8,ss),b_new(9,ss),b_new(10,ss),b_new(11,ss),b_new(12,ss)

                   ! write(21,*)b_new(1,ss),b_new(2,ss),b_new(3,ss),b_new(4,ss)
                end do
         end if
end do ts1loop

close(12)
end do

open(12,file='pop_au.txt',action='read',status ='old')
close(12, status='delete')
close(21)
return
end subroutine samples_au



subroutine p_corr_au(corr_au)
use constants
implicit none

          !integer, parameter :: tnob=5
          real, dimension(tnob,tnob) :: r,rs,corr_au
          integer :: i,j,q
          real,dimension(tnob) :: m,m_s,m_rs, sd ,m0
          real :: corr1,corr2
          open(55,file='pop_samples_au.txt', action='read',status = 'old')
            open(33,file='pearson_au.txt', action='write',status = 'new')

        !    open(44,file='pearson1_au.txt', action='write',status = 'new')

          do i=1,tnob
             read(55,*) (r(i,j), j=1,tnob)
          end do
          close(55,status = 'delete')


          do i=1,tnob
             m(i)=0.0
             !sd(i)=0.0
             m_rs(i)=0.0
          end do

    !!!!!!!!!!!!!!!!mean!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i=1,tnob
            do j=1,tnob
               m(i)=m(i)+r(j,i)
            end do
           m(i)=m(i)/float(tnob)
          end do

  !         do i=1,tnob
  !          write(6,*) m(i)
   !        end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!Standard deviation!!!!!!!!!!!!!!!!
            do i=1,tnob

               do j=1,tnob

                   rs(i,j)= r(i,j)**2

               end do

            end do

             do i=1,tnob
               do j=1,tnob
                m_rs(i)=m_rs(i)+rs(j,i)
            end do
               m_rs(i)=m_rs(i)/float(tnob)
          end do

           do i=1,tnob

             m_s(i)=m(i)**2

          end do

            do i=1,tnob
              sd(i)=sqrt(m_rs(i)-m_s(i))
            end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!Pearson correlation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
           do i=1,tnob
               do j=1,tnob
                  corr_au(i,j)=0.0
               end do
           end do

    !       write(6,*) 'hi'
           corr1=0.0

           do i=1,tnob

            corr1=corr1+((r(i,1)-m(1))*(r(i,2)-m(2)))

          end do

           corr1=corr1/(sd(1)*sd(2)*tnob)

     !      write(6,*) corr1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
          ! do k=1,tnob
               do i=1,tnob
                  do j=1, tnob
                     do q=1, tnob
                     corr_au(i,j)=corr_au(i,j)+((r(q,j)-m(j))*(r(q,i)-m(i)))
                     end do
                    ! write(6,*) corr2
                  end do
               end do
             ! end do 

               do i=1,tnob
                  do j=1,tnob
                    corr_au(i,j)=corr_au(i,j)/(sd(i)*sd(j)*tnob)
                  end do
               end do

           do i=1,tnob
             write(33,*) (corr_au(i,j), j=1,tnob)
           end do

     !      do i=1,tnob
      !       do j=1, tnob
       !        write(44,*) corr_au(i,j)
        !     end do
         !  end do

    close(33)
    !close(44)
    return
 end subroutine p_corr_au
 

  subroutine corr_rmsd_bu(rmsd3)
     use constants     
   implicit none
      
          integer :: i,j
        real, dimension(tnob,tnob) :: corr_ex,corr_au,corr_bu
        real :: rmsd3
      
       open(1, file= 'pearson_bu.txt', status = 'old')
         do i=1,tnob
           !do j=1,tnob
             read(1,*) (corr_bu(i,j) , j=1, tnob)
           !end do
         end do
       !close(1, status = 'delete')
 close(1)
     

       !open(2, file = 'pearson_au.txt' , status ='old')


         ! do i=1,tnob
            !do j=1,tnob
          !     read(2,*) (corr_au(i,j) , j=1,tnob)
            !end do
         !end do
       !close(2, status = 'delete')
     
       

        open(3, file = 'co_ma_ex.dat' , status ='old')


          do i=1,tnob
            !do j=1,tnob
               read(3,*) (corr_ex(i,j) , j=1,tnob)
            !end do
         end do
       close(3)
          rmsd3=0.0
     !     rmsd4=0.0


             do i=1,tnob
               do j=1, tnob
                 rmsd3=rmsd3+((corr_bu(i,j)-corr_ex(i,j))/corr_ex(i,j))**2
                 
            !     rmsd4=rmsd4+((corr_au(i,j)-corr_ex(i,j))/corr_ex(i,j))**2
                 
                end do
             end do 
              

                 rmsd3=sqrt(rmsd3)/float(tnob**2)
           !      rmsd4=sqrt(rmsd4)/float(tnob**2)
  return                            
  end subroutine corr_rmsd_bu


  subroutine corr_rmsd_au(rmsd4)
     use constants     
   implicit none
      
          integer :: i,j
        real, dimension(tnob,tnob) :: corr_ex,corr_au,corr_bu
        real :: rmsd4
      
      ! open(1, file= 'pearson_bu.txt', status = 'old')
       !  do i=1,tnob
           !do j=1,tnob
        !     read(1,*) (corr_bu(i,j) , j=1, tnob)
           !end do
        ! end do
       !close(1, status = 'delete')

     

       open(2, file = 'pearson_au.txt' , status ='old')


          do i=1,tnob
            !do j=1,tnob
               read(2,*) (corr_au(i,j) , j=1,tnob)
            !end do
         end do
       close(2, status = 'delete')
     
       

        open(3, file = 'co_ma_ex.dat' , status ='old')


          do i=1,tnob
            !do j=1,tnob
               read(3,*) (corr_ex(i,j) , j=1,tnob)
            !end do
         end do
       close(3)
   !       rmsd3=0.0
          rmsd4=0.0


             do i=1,tnob
               do j=1, tnob
         !        rmsd3=rmsd3+((corr_bu(i,j)-corr_ex(i,j))/corr_ex(i,j))**2
                 
                 rmsd4=rmsd4+((corr_au(i,j)-corr_ex(i,j))/corr_ex(i,j))**2
                 
                end do
             end do 
              

          !       rmsd3=sqrt(rmsd3)/float(tnob**2)
                 rmsd4=sqrt(rmsd4)/float(tnob**2)
  return                            
  end subroutine corr_rmsd_au

