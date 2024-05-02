include 'random_module.f90'
program GHbid
    USE randomize

    implicit none
    integer                                        :: N, L, m, k
    integer                                        :: i,j,p,ii,jj,kk,idum,mmm
    integer                                        :: tmax,tran,Npuntos
    integer                                        :: num_clusters,max_cluster
    integer                                        :: max_cluster_size,maxlist
    real                                           :: ppi,inicial
    real                                           :: giant_m,smedio,smedio_m
    real                                           :: r,r1,r2,T,Tini,DT
    integer,allocatable,dimension(:,:)             :: AA,V, AA2
    integer, allocatable, dimension(:)             :: size_cluster
    real,allocatable,dimension(:,:)                :: WW
    integer (kind = 2),allocatable,dimension(:)    :: grado,x,xant
    integer,allocatable,dimension(:)               :: cluster,list

    open(86, file='smedio_10_0,2.dat')
    open(99, file='cluster_gigante_10_0,2.dat')

    idum=2873466      ! semilla generador de numeros aleatorios
    k=3
    L=10
    N=L**2
    ppi=0.2
    inicial=0.3
    r2 = 0.3
    r1 = 0.001
    tran = 1000
    Tini = 0.1
    tmax=10000
    DT = 0.001
    Npuntos = 500
    mmm = mzranset(521288629,362436069,16163801,idum)

    allocate(AA(0:N-1,0:N-1))
    allocate(AA2(N,N))
    allocate(WW(N,N))
    allocate(V(N,N))
    allocate(grado(N))
    allocate(x(N))
    allocate(xant(N))
    allocate(cluster(N))
    allocate(size_cluster(N))
    allocate(list(N))

    call WS
    do i=1,N
      do j=1,N
         AA2(i,j)=AA(i-1,j-1)
      end do
    end do
    call pesos
    call vecinos
    deallocate(AA)
    deallocate(AA2)

    do i=1,N                         !inicializa neuronas
      r = random()
      if (r < inicial) then
        x(i) = 1
      else
        x(i) = 0
      end if
    end do

    do k=1,tran                     !transitorio
      call dinamicaGH
    enddo

    write(99,*) '           t       cl'
    write(86,*) '           t     smedio'
    
    giant_m=0
    smedio_m=0

    T = Tini
    do kk=1,Npuntos
      do p=1,tmax                     
        call dinamicaGH
        call get_clusters
        giant_m=giant_m+max_cluster_size/float(N)
        call average_clusters_size
        smedio_m=smedio_m+smedio
      enddo
      giant_m=giant_m/float(tmax)
      smedio_m=smedio_m/float(tmax)  
      write(1,*) T, ' ',giant_m,' ',smedio_m
      T = T + DT
    enddo

    deallocate(WW)
    deallocate(V)
    deallocate(grado)
    deallocate(x)
    deallocate(cluster)
    deallocate(size_cluster)
    deallocate(list)
    deallocate(xant)

    contains
!-------------------------------------------------------------------------------------------------
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
            if (r < ppi) then                                          
              if (AA(z_i,z_j)==1) then                                 !Si los nodos estan conectados los desconecta     
                AA(z_i,z_j)=0
                AA(z_j,z_i)=0
                nuevositio=z_i                                              
                do while((nuevositio == z_i).or.(AA(z_i,nuevositio) == 1))  !Como nuevositio=z_i, siempre entra en este loop. La 
                   nuevositio=int(random()*N)                               !segunda condicion esta para que en el caso que debamos 
                end do                                                      !volver a sortear nuevositio pues ya esta conectado a z_j,               
                AA(z_i,nuevositio)=1                                        !entonces vuelva a entrar en el do while pues AA(z_i,nuevositio)=1
                AA(nuevositio,z_i)=1                                        !Conectamos al nuevo sitio sorteado
              end if
            end if 
          end do
        end do
         
      end subroutine
!----------------------------------------------------------------------------------------------------------
      function min(xx,yy)

        implicit none
        integer xx,yy,min
        
        if (xx <= yy) then
          min = xx
        else
          min=yy
        endif
        
      end function min
!------------------------------------------------------------------------------------------------------------------------
      function delta(xx,yy)

        implicit none
        integer (kind = 2) xx
        integer yy,delta
        
        if (xx == yy) then
          delta=1
        else
          delta=0
        endif
        
      end function delta
!------------------------------------------------------------------------------------------------------------------------
      function expran()

        implicit none
        real expran,r
        
        r = random()
        expran = -0.08*log(1-r)
        
      end function expran
!------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------
      subroutine pesos ! calcula matriz de pesos sinapticos excitatorios

        implicit none
        integer         :: ii,jj,kk
        
        WW = 0                     
        do ii=1,N
          do jj=ii+1,N
            if (AA2(ii,jj) == 1) then
              WW(ii,jj) = expran()
              WW(jj,ii) = WW(ii,jj)
            endif
          enddo
        enddo
      end subroutine pesos 
!-----------------------------------------------------------------------------------------------------------------------------
      subroutine dinamicaGH ! un paso de la dinamica de GH

        implicit none
        integer         :: ii,jj,kk
        real            :: sum
        
        xant = x                    
        
        do ii=1,N
          if (xant(ii) == 1) then
            x(ii) = 2
          elseif (xant(ii) == 2) then
            r = random()
            if (r < r2 ) x(ii) = 0
          else
            sum=0
            do kk=1,grado(ii)
              jj=V(ii,kk)
              sum = sum + WW(ii,jj)*delta(xant(jj),1)
            enddo
            r=random()
            if ((sum > T).or.(r < r1)) x(ii) = 1 
          endif
        enddo
      
      end subroutine dinamicaGH
!---------------------------------------------------------------------------------------------------------------------------
      subroutine paint(site,color)

        implicit none
        integer color,site
        ! la lista x tiene la informacion de la ocupacion de la red.
        ! la lista cluster es una particion de la red que identifica los clusters. 
        if((x(site)==1).and.(cluster(site)==0)) then
          cluster(site)=color
          maxlist=maxlist+1
          list(maxlist)=site 
        end if
        
      end subroutine paint
!-----------------------------------------------------------------------------------------------------------------------------
      subroutine genera_cluster(ii,color)
        ! generate a cluster starting from site ii 
        
        implicit none
        integer, intent(in)  :: color,ii
        integer              :: i,j,k,ll
        
        list=0
        cluster(ii)=color
        maxlist=1
        list(maxlist)=ii
        i=1
        !si el sitio vecino esta ocupado la subrutina 'paint'  pinta el sitio e incrementa maxlist
        !en 1 y guarda el nuevo sitio en list(maxlist) para explorar posteriormente sus vecinos.
        do while (i<=maxlist)
              ll = list(i)
              do j=1,grado(ll)
                 k = V(ll,j)
                 call paint(k,color)
              end do
          i=i+1
        end do
      end subroutine genera_cluster
!---------------------------------------------------------------------------------------------------------------------------------
      subroutine get_clusters

        implicit none
        integer color,ii,jj
        
        cluster=0
        color=0
        
        !se genera la particion que identifica los clusters
        do ii=1,N
          if ((x(ii)==1).and.(cluster(ii)==0)) then
            color=color+1
            call genera_cluster(ii,color)
          end if
        end do
        num_clusters=color
        
        size_cluster=0 
        do jj=1,num_clusters !extract clusters size
          do ii=1,N
            if (cluster(ii)==jj) then
              size_cluster(jj)=size_cluster(jj)+1
            end if  
          end do
        end do
        
        max_cluster_size=size_cluster(1) ! identify giant cluster
        max_cluster=1
        do ii=2,num_clusters
          if (size_cluster(ii) > max_cluster_size) then
            max_cluster_size=size_cluster(ii)
            max_cluster=ii
          end if
        end do
        
      end subroutine get_clusters
!--------------------------------------------------------------------------------------------------------------------------------
      subroutine average_clusters_size

        integer ii
        real    summ,sum2
        
        do ii=max_cluster,num_clusters-1      !extract giant cluster from list
          size_cluster(ii)=size_cluster(ii+1)
        end do
        num_clusters=num_clusters-1
        
        
        summ=0
        sum2=0
        do ii=1,num_clusters
          summ=summ+size_cluster(ii)/float(N)
          sum2=sum2+size_cluster(ii)**2/float(N)
        end do
        
        if (summ /= 0) then
          smedio=sum2/summ
        else
          smedio=0
        endif
      end subroutine average_clusters_size
!-----------------------------------------------------------------------------------------------------------------------------------
    end program 





    