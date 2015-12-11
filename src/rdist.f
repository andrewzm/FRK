c fields, Tools for spatial data
c Copyright 2004-2013, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html

       subroutine rdist( nd,x1,n1,x2,n2, k)
       integer nd,n1,n2,ic   
       double precision  x1(n1,nd), x2(n2,nd), k(n1,n2)
       double precision xtemp, radfun
      
         do  j =1,n2
              xtemp= x2(j,1)
              do   i= 1, n1
c** accumulate squared differences
                 k(i,j)=  (x1(i,1)- xtemp)**2 
              enddo  
          enddo
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce memory swaps
       if( nd.ge.2) then  
       do  ic= 2, nd
          do  j =1,n2
              xtemp= x2(j,ic)
              do   i= 1, n1
c** accumulate squared differences
                 k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
              enddo  
          enddo
       enddo
       endif
c**** at this point k( i,j) is the squared distance between x1_i and x2_j
       do  j =1,n2
               do   i= 1, n1
                k(i,j)= sqrt( k(i,j))
               enddo
       enddo 
       return
       end



      subroutine rdist1( nd,x1,n1,k)
       integer nd,n1,ic   
       double precision  x1(n1,nd), k(n1,n1)
       double precision xtemp,  dtemp
      
         do  j =1,n1
              xtemp= x1(j,1)
              do   i= 1, j
c** accumulate squared differences
                 k(i,j)=  (x1(i,1)- xtemp)**2 
              enddo  
          enddo
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce memory swaps
       if( nd.ge.2) then  
       do  ic= 2, nd
          do  j =1,n1
              xtemp= x1(j,ic)
              do   i= 1, j
c** accumulate squared differences
                 k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
              enddo  
          enddo
       enddo
       endif
c**** at this point k( i,j) is the squared distance between x1_i and x2_j
c**** for the upper triangle of matrix
       do  j = 1,n1
               do   i= 1, j
                dtemp = sqrt( k(i,j))
                k(i,j)= dtemp
c                
c filling in lower part takes a substantial time and is omitted  
c This means the returned matrix k has indeterminant vlaues in the
c lower triangle.              
c                k(j,i)= dtemp
               enddo
       enddo 
       return
       end
