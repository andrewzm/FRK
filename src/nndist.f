C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine nndisg(pts,npts,dists,neighs)
c
c return nearest neighbour distances and the indices of the neighbours
c
      implicit real*8 (a-h,o-z)

      real*8 pts(2,npts)
      real*8 dists(npts)
      integer neighs(npts)

      do i=1,npts
        dsmin=1d40
        do j=1,npts
          if(i.ne.j)then
            dij=((pts(1,i)-pts(1,j))**2+(pts(2,i)-pts(2,j))**2)
            if(dij.lt.dsmin)then
              dsmin=dij
              nbr=j
            end if
          end if
        end do
        dists(i)=sqrt(dsmin)
        neighs(i)=nbr
      end do

      end

      subroutine n2dist(x1,y1,n1pts,x2,y2,n2pts,dists,neighs)
      implicit real*8 (a-h,o-z)

      dimension x1(n1pts),y1(n1pts),x2(n2pts),y2(n2pts)
      dimension dists(n2pts),neighs(n2pts)

      
c
c  return distances from the n2pts in x2,y2 to the nearest points
c in x1,y1
c

      do i=1,n2pts
        dsmin=1d40
	xp2=x2(i)
        yp2=y2(i)
        do j=1,n1pts
          dij=((x1(j)-xp2)**2+(y1(j)-yp2)**2)
          if(dij.lt.dsmin)then
            dsmin=dij
            nbr=j
          end if
        end do
        dists(i)=sqrt(dsmin)
        neighs(i)=nbr
      end do


      end


      subroutine nndisf(x1,y1,n1pts,x2,y2,n2pts,
     &                   dists)

      implicit real*8 (a-h,o-z)

      dimension x1(n1pts),y1(n1pts),x2(n2pts),y2(n2pts)
      dimension dists(n2pts)
c
c  return distances from the n2pts in x2,y2 to the nearest points
c in x1,y1
c

      do i=1,n2pts
        dsmin=1d40
	xp2=x2(i)
        yp2=y2(i)
        do j=1,n1pts
          dij=((x1(j)-xp2)**2+(y1(j)-yp2)**2)
          if(dij.lt.dsmin)then
            dsmin=dij
          end if
        end do
        dists(i)=sqrt(dsmin)
      end do

      end
