include 'random_module.f90'
program WSbidimensionaleuc
  USE randomize

    implicit none

    integer, allocatable, dimension(:,:)     :: AA, dist, V, AA2
    real,allocatable,dimension(:)            :: CC
    integer,allocatable,dimension(:)         :: grado
    integer                                  :: i,j, i_p, j_p
    integer                                  :: N, L, m, k, num_args, mmm, idum, nn
    character (len=50), allocatable          :: args(:)        
    real (kind=4)                            :: ppi, pi(50), lmedio, Clust, lmedio_n, Clust_n, lmedio0, Clust0, lmedio_np, Clust_np


    !num_args = command_argument_count()
    !allocate(args(num_args))

    !do i=1,num_args
    !  call get_command_argument(i,args(i))
    !end do

    !read(args(:),*) k,L
    !open(99, file ='conexiones_4.dat')
    !open(86, file='WS.dat')
    !open(44, file='WS_conexiones.dat')
    open(23, file='cluster_30_3.dat', status='new')
    open(24, file='lmedio_30_3.dat', status='new')
    !write(*,*) "Ingresar el tamaño de la red"
    !read(*,*) L

    !write(*,*) "Ingresar el grado k de vecinos más cercanos"
    !read(*,*) k
    
    !Estos valores son solo las coordenadas de un sitio de la red para graficar
    !i_p = 4
    !j_p = 4

    idum=2873466      ! semilla generador de numeros aleatorios
    
    k = 3
    L = 30
    N = L**2

    mmm = mzranset(521288629,362436069,16163801,idum)
    
    call logspace_array(50,pi)
    
    write(23,*) '      pi               cl'
    write(24,*) '      pi             lmedio'

    do m=1,50
    !  lmedio_np=0
    !  Clust_np=0
    !  do nn=1,10
        allocate(AA(0:N-1,0:N-1))
        allocate(AA2(N,N))
        allocate(dist(N,N))
        allocate(V(N,N))
        allocate(grado(N))
        allocate(CC(N))
     
        ppi = pi(m)                !Setea la probabilidad p de reconexion
        call WS
         do i=1,N
          do j=1,N
           AA2(i,j)=AA(i-1,j-1)
          end do
         end do
        call minimum_distances
        call vecinos
        call clustering 
     
        Clust=0
        do i=1,N
          Clust=Clust+CC(i) 
        enddo
     
        lmedio = 0
        do i=1,N
          do j=1,N
            lmedio = lmedio+dist(i,j)
          enddo
        enddo
     
        lmedio = lmedio/(N*(N-1))
        Clust = Clust/N
        if (m==1) then
          lmedio0=lmedio
          Clust0=Clust
        end if 
        lmedio_n = lmedio/lmedio0
        Clust_n = Clust/Clust0
        !lmedio_np = lmedio_np+lmedio_n/10
        !Clust_np = Clust_np+Clust_n/10
        deallocate(AA)
        deallocate(AA2)
        deallocate(dist)
        deallocate(V)
        deallocate(grado)
        deallocate(CC)
      !end do 
      write(23,*) ppi, Clust_n
      write(24,*) ppi, lmedio_n
    end do

    close(23)
    close(24)
    

!------------------------------------------------------------------------------------------------------------------------
!ESTA ES OTRA FORMA DE HACER EL ALGORITMO DE W-S, LO DEJO POR LAS DUDAS
!       do z_i=0,N-1
!        do z_j=z_i+1,N-1 
!          ppi = 0.001
!          !zj = z_j
!          if (AA(z_i,z_j)==1) then 
!            r = rand()
!            if (r < ppi) then
!              AA(z_i,z_j)=0
!              AA(z_j,z_i)=0                                              
!        21    nuevositio=int(rand()*(N+1)) 
!              write(86,*) r, ppi, nuevositio
!                do while ((nuevositio/=z_j) .and. (nuevositio/=z_i))
!                    AA(z_i,nuevositio)=1
!                    AA(nuevositio,z_i)=1
!                    write(86,*) z_i, z_j, nuevositio
!                    write (86,*) ''
!                end do 
!              if ((nuevositio==z_j) .or. (nuevositio==z_i)) then 
!                go to 21
!              end if
!            end if
!          end if
!        end do
!       end do
!----------------------------------------------------------------------------------------------------------------------------

    contains 
!-----------------------------------------------------------------------------------------------------------------------------  
      subroutine print_matrix(b,n,m) !matriz b 
        integer :: n,m
        integer :: b(n,m) !n = # rows, m = # columns
        do i=1,n
           print ('(20i6.2)'), b(i,1:m)
         enddo
      end subroutine 
!----------------------------------------------------------------------------------------------------------------------------
      subroutine WS
        
        implicit none 
        integer          :: i_n, i_l, j_n, j_l, z_n, z_l, nuevositio, m_
        integer          :: d_i, d_i_per, d_j, d_j_per
        integer          :: z_i, z_j
        real(kind=4)     :: r

        AA=0                    !Inicializa la matriz de adyacencia con todas sus entradas nulas
  
        do z_n=0,N-1            !Corro para los elementos de fila de AA desde 0 a N-1 
         do z_l=z_n+1,N-1       !Para cada uno de los elementos de fila, corro para los elementos de columna  
          i_n = z_n/L           !Calculo las coordenadas de los nodos correspondientes a cada 
          i_l = z_l/L           !elemento de matriz   
          j_n = mod(z_n,L)
          j_l = mod(z_l,L)   
          do m_=1,k
              
              d_i = abs(i_n-i_l)                     !Calcula las distancias entre pares de nodos considerando condiciones
              d_j = abs(j_n-j_l)                     !periodicas de borde en las d_i(j)_per
              d_i_per = mod(abs(i_n-i_l),L)
              d_j_per = mod(abs(j_n-j_l),L)
          
              if (m_==1) then 
                  if (((d_i==1) .and. (d_j==0)) .or. ((d_i==0) .and. (d_j==1)) .or. ((d_i_per==L-1) .and. (d_j==0)) &
                  .or. ((d_i==0) .and. (d_j_per==L-1))) then 
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if 
              
              else if (m_==2) then 
                  if (((d_i==1) .and. (d_j==1)) .or. ((d_i_per==L-1) .and. (d_j==1)) .or. ((d_i==1) .and. (d_j_per==L-1)) & 
                  .or. ((d_i_per==L-1) .and. (d_j_per==L-1))) then 
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if 
              
              else if (m_==3) then 
                  if (((d_i==2) .and. (d_j==0)) .or. ((d_i==0) .and. (d_j==2)) .or. ((d_i_per==L-2) .and. (d_j==0)) &
                  .or. ((d_i==0) .and. (d_j_per==L-2))) then 
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if 
              
              else if (m_==4) then 
                  if (((d_i==1) .and. (d_j==2)) .or. ((d_i==2) .and. (d_j==1)) .or. ((d_i_per==L-1) .and. (d_j==2)) &
                  .or. ((d_i_per==L-2) .and. (d_j==1)) .or. ((d_i==1) .and. (d_j_per==L-2)) &
                  .or. ((d_i==2) .and. (d_j_per==L-1)) .or. ((d_i_per==L-1) .and. (d_j_per==L-2)) &
                  .or. ((d_i_per==L-2) .and. (d_j_per==L-1))) &
                  then
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if
              
              else if (m_==5) then 
                  if (((d_i==2) .and. (d_j==2)) .or. ((d_i_per==L-2) .and. (d_j==2)) .or. ((d_i==2) .and. (d_j_per==L-2)) & 
                  .or. ((d_i_per==L-2) .and. (d_j_per==L-2))) then 
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if 
              
              else if (m_==6) then 
                  if (((d_i==3) .and. (d_j==0)) .or. ((d_i==0) .and. (d_j==3)) .or. ((d_i_per==L-3) .and. (d_j==0)) &
                  .or. ((d_i==0) .and. (d_j_per==L-3))) then 
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if
              
              else if (m_==7) then 
                  if (((d_i==1) .and. (d_j==3)) .or. ((d_i==3) .and. (d_j==1)) .or. ((d_i_per==L-3) .and. (d_j==1)) &
                  .or. ((d_i_per==L-1) .and. (d_j==3)) .or. ((d_i==1) .and. (d_j_per==L-3)) &
                  .or. ((d_i==3) .and. (d_j_per==L-1)) .or. ((d_i_per==L-1) .and. (d_j_per==L-3)) &
                  .or. ((d_i_per==L-3) .and. (d_j_per==L-1))) &
                  then
                      AA(z_n,z_l) = 1               
                      AA(z_l,z_n) = 1
                  end if
              end if 
            end do 
          end do 
        end do
      
      !Algoritmo de Watts-Strogatz
        do z_i=0,N-1                                                  
          do z_j=z_i+1,N-1  
            r = random()                                                 !Sortea un numero random
            !write(86,'(ES16.8,2X,F8.5,2X,I5,3X,I5)') r, ppi, z_i, z_j  !Escribe el par de nodos donde estoy parada y el numero sorteado
            if (r < ppi) then                                          
              if (AA(z_i,z_j)==1) then                                 !Si los nodos estan conectados los desconecta     
                AA(z_i,z_j)=0
                AA(z_j,z_i)=0
                nuevositio=z_i                                              
                do while((nuevositio == z_i).or.(AA(z_i,nuevositio) == 1))  !Como nuevositio=z_i, siempre entra en este loop. La 
                   nuevositio=int(random()*N)                               !segunda condicion esta para que en el caso que debamos 
                   !if ((nuevositio==z_j).or.(nuevositio==N)) then           !volver a sortear nuevositio pues ya esta conectado a z_j,  
                   ! cycle                                                   !entonces vuelva a entrar en el do while pues AA(z_i,nuevositio)=1
                   !end if     
                end do                                                                     
                AA(z_i,nuevositio)=1                                        !Conectamos al nuevo sitio sorteado
                AA(nuevositio,z_i)=1
                !write(86,'(A,2X,I5)') "nodo reconectado, su nueva conexion es el nodo", nuevositio
                !write (86,*) ''
              end if
            end if 
          end do
        end do
         
      end subroutine
!--------------------------------------------------------------------------------------------------------------------
      subroutine logspace_array(n,x)
        implicit none
        
        integer, intent(in)           :: n ! number of elements
        real (kind=4), intent(out)    :: x(n)
        real, parameter               :: x0 = 1.0E-6
        real, parameter               :: x1 = 1.0
        integer                       :: i
        
        x(1) = 0.0
        x(2) = x0
        x(n) = x1
        
        do i = 3, n-1
          x(i) = 10.0**(log10(x(2)) + (log10(x(n)) - log10(x(2))) * (i-2) / (n-1))
        end do
        
      end subroutine logspace_array     
!--------------------------------------------------------------------------------------------------------------------
      subroutine linspace(from, to, array)
        
        implicit none
        real(kind=4), intent(in)  :: from, to
        real(kind=4), intent(out) :: array(:)
        real(kind=4)              :: range
        integer :: n, i
        n = size(array)
        range = to - from
    
        if (n == 0) return
    
        if (n == 1) then
            array(1) = from
            return
        end if
    
    
        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
      end subroutine
!--------------------------------------------------------------------------------------------------------------------
      subroutine minimum_distances ! algoritmo Floyd-Warshall
        
        implicit none
        integer        :: ii,jj,kk
        
        do ii=1,N
          do jj=1,N
            if (AA2(ii,jj) == 1) then
              dist(ii,jj)=AA2(ii,jj)
            else
              dist(ii,jj) = 100*N
            end if
          enddo
        enddo
        
        do ii=1,N
          dist(ii,ii)=0
        enddo
        
        do kk=1,N 
          do ii=1,N
            do jj=1,N
               dist(ii,jj) = min(dist(ii,jj), dist(ii,kk)+dist(kk,jj))
             enddo
          enddo
        enddo
      
      end subroutine minimum_distances
!----------------------------------------------------------------------------------------------------------------
      subroutine vecinos ! calcula lista de vecinos y grado de cada nodo
        
        implicit none
        integer        :: ii,jj,kk,aux
        
        do ii=1,N
          aux=0
          do jj=1,N
            aux=aux+AA2(ii,jj)
          enddo
          grado(ii)=aux
        enddo
        
        V=0
        do ii=1,N
          kk=1
          do jj=1,N
            if (AA2(ii,jj) == 1) then
              V(ii,kk)=jj
              kk=kk+1
            end if
          enddo
        enddo
      end subroutine vecinos 
!-----------------------------------------------------------------------------------------------------------------------
      subroutine clustering
        
        implicit none
        integer           :: ii,jj,kk,aux,TT,ll,rr
        
        do ii=1,N
          TT=0
          do ll=1,grado(ii)-1
            jj=V(ii,ll)
            do kk=ll+1,grado(ii)
              rr=V(ii,kk)
              if (AA2(jj,rr) /=0) TT=TT+1
            enddo
          enddo
          if ((grado(ii)==0).or.(grado(ii)==1)) then 
            CC(ii)=0
          else
            CC(ii)=2*TT/float((grado(ii)*(grado(ii)-1)))
          end if 
        enddo
      
      end subroutine clustering
!---------------------------------------------------------------------------------------------------------------------------
      subroutine GRAFICAR
     
       implicit none
       integer                        :: i_p, j_p, zp_ij, z_n, i_n, j_n
       integer                        :: dp_i, dp_i_per, dp_j, dp_j_per
       real (kind=4)                  :: d, d_1, d_2, d_3, d_4

       !ESCRITURA DE DATOS PARA GRAFICAR
       !Escribe las coordenadas de los nodos que estan conectados en un archivo de texto 
        do z_n = 0, N-1
          zp_ij = L*i_p+j_p
          if (AA(zp_ij,z_n)==1) then 
           i_n = z_n/L
           j_n = mod(z_n,L)
           dp_i = abs(i_p-i_n)
           dp_j = abs(j_p-j_n)
           dp_i_per = mod(abs(i_p-i_n),L)
           dp_j_per = mod(abs(j_p-j_n),L)
           d_1 = sqrt(real(dp_i**2+dp_j**2))
           d_2 = sqrt(real(dp_i_per**2+dp_j**2))
           d_3 = sqrt(real(dp_i**2+dp_j_per**2))
           d_4 = sqrt(real(dp_i_per**2+dp_j_per**2))
           d = min(d_1, d_2, d_3, d_4)
           write(99, *) i_p, j_p, 8
           write(99, *) i_n, j_n, d
           write(99, *) ""  
          end if
        end do
      end subroutine
!---------------------------------------------------------------------------------------------------------------
end program 