  program main
!use constants
implicit none

          integer, parameter :: tnob=12
          real, dimension(tnob,tnob) :: r,rs,corr_bu
          integer :: i,j,q
          real,dimension(tnob) :: m,m_s,m_rs, sd ,m0
          real :: corr1,corr2
          open(55,file='samples.txt', action='read',status = 'old')
            open(33,file='p_corr.txt', action='write',status = 'new')

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

           !do i=1,tnob
            ! write(33,*) (corr_bu(i,j), j=1,tnob)
           !end do

          do i=1,tnob
             do j=1, tnob
               write(33,*) corr_bu(i,j)
             end do
           end do

    close(33)
    !close(44)
   
 end 
