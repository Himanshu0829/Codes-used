        program main
        implicit none
         integer , parameter :: tnob = 12   
         real, dimension (tnob,tnob) :: bs
         integer :: i,j
         real :: r1, r2, r3
       open(11, file= 'bs4.dat' , status = 'new', action = 'write')
       
              do i=1,tnob

          do j =1,tnob

             r1=rand()
             r2=rand()
             r3=rand() 
             if (j .gt. i ) then

                 bs(i,j)=(r1-0.5)

             else if (i .eq. j) then

                 bs(i,j) =1.0

             else if (j .lt. i) then

                 bs(i,j)=bs(j,i)
             end if
           end do

       end do





  write(6,*) 'hi'
 
      do i = 1, tnob
    
           write(11,*) (bs(i,j), j=1, tnob)

      end do

  close(11)


      do i = 1, tnob
          
           write(6,*) (bs(i,j), j=1, tnob)

      end do


   end 
