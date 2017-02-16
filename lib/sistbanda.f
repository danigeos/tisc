C ****************************************************************
C   SISTEMA DE n INCOGNITES I n EQUACIONS LINEALMENT INDEPENDENTS
C                MATRIU a EN BANDA
C                  SOLUCIO UNICA
C ****************************************************************

C                   a ˙ x = b
C  a: Matriu del sistema d'equacions, a(nfiles,nbanda)
C  nfiles: nß d'equacions = nß d'incognites x(nfiles) = nß de files d'a
C  nbanda: nß de columnes d'a   ( nbanda = Li+1+Ls )
C  b: terme independent, nfiles files
C  Li: banda inferior, termes per sota/esquerra de la diagonal
C  Ls: banda superior, termes per sobre/dreta de la diagonal
C      Els termes de la diagonal estan en la columna id

      SUBROUTINE sistbanda(a,b,nfiles,nbanda,Li,x)

      implicit double precision (a-h,o-z)
      dimension a(nfiles,nbanda),b(nfiles),x(nfiles)
      
      id=Li+1
      Ls=nbanda-(Li+1)

      do 5 j=1,nfiles-1
          if(a(j,id).eq.0.D0) then
             do 53 l=j+1,j+Li
                 jj=id+j-l
                 if(a(l,jj).ne.0.D0) then
                     do 8 k1=1,nbanda
                        k2=k1+j-l
                        ap=a(l,k2)
                        a(l,k2)=a(j,k1)
                        a(j,k1)=ap
8                    continue
                     bp=b(l)
                     b(l)=b(j)
                     b(j)=bp
                     goto 54
                 endif
53           continue
             print*,'no he trobat com substituir el terme',j
             print*,'  SISTEMA INDETERMINAT '
             stop 
           endif
54        continue

          do 7 i=j+1,j+Li
                IF(i.GT.nfiles) GOTO 5
                kj=id+j-i
                if(a(i,kj).eq.0.D0) goto 7
                fac=a(i,kj)/a(j,id)
                do 9 k1=kj,nbanda
                     k2=k1+i-j
                     if(k2.gt.nbanda) then
                          aaa=0.d0
                       else
                          aaa=a(j,k2)
                     endif
                     a(i,k1)=aaa*fac-a(i,k1)
9               continue
                b(i)=b(j)*fac-b(i)
7         continue
5     continue

c ******** Aãllar les incïgnites ********************

       if(a(nfiles,id).eq.0.D0) then
           print*,'L`ULTIM TERME DE LA DIAGONAL êS ZERO'
           stop
       endif
       x(nfiles)=(b(nfiles))/(a(nfiles,id))
       do 17 kk=1,nfiles-1
          k=nfiles-kk
          if(a(k,id).eq.0.D0) then
              print*,'EL TERME DE LA DIAGONAL DE LA FILA',k,' êS ZERO'
              stop
          endif
          c=0.d0
          do 19 l=id+1,nbanda
              ll=l+k-id
              c=c+a(k,l)*x(ll)
19        continue
         x(k)=(b(k)-c)/a(k,id)
17     continue

20    RETURN
      END
