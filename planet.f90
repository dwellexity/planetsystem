module var
  implicit none
  real, parameter :: G=1.0
  integer, parameter :: n=2
  real, parameter :: mEarth = 1.0
  real, parameter :: mSun = 332946.0
  real, parameter :: distPerih = 1.52e11
  real, dimension(3,n) :: xyz,vel,acc
  real, dimension(3,n) :: angmom
  real, parameter, dimension(n) :: m=[10000.0, 1.0]
  real, dimension(3) :: orig
  real, parameter :: tot_t=1000.0
  real, parameter :: dt=0.001
  real, parameter, dimension(n) :: rad=[1.0,1.0]
  real :: ke,pe
end module var
module initial
  use var, only: xyz,vel,acc,m,orig,angmom,G
  use cross_product, only:cross
  implicit none
  private
  public :: init
contains
  subroutine init()
    real :: temp
    xyz(:,1)=0.0
    xyz(:,2)=[10.0,0.0,0.0]
    orig=(m(1)*xyz(:,1)+m(2)*xyz(:,2))/(sum(m))
    vel(:,1)=[0,0,0]
    vel(:,2)=[0.0,30.0,0.0]
    angmom(:,1)=cross(xyz(:,1),m(1)*vel(:,1))
    angmom(:,2)=cross(xyz(:,2),m(2)*vel(:,2))
  end subroutine init
end module initial
module cross_product
  implicit none
contains
  function cross(a, b)
    real, dimension(3) :: cross
    real, dimension(3), intent(in) :: a, b
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
end module cross_product
module update
  use var, only: xyz,vel,acc,m,orig,angmom,G,n,dt
  implicit none
  private
  integer :: i,j
  real, dimension(3,n) :: tempacc
  real, dimension(3) :: dist
  public :: posUpd,velUpd,accUpd
contains
  subroutine posUpd()
    implicit none
    integer :: i,k
    do i=2,n
      do k=1,3
        xyz(k,i)=xyz(k,i) + vel(k,i)*dt + (acc(k,i)*(dt**2))/2
      enddo
    enddo
  end subroutine posUpd
  subroutine accUpd()
    implicit none
    integer :: i,j,k
    real :: r,temp
    do i=2,n
      tempacc(:,i)=acc(:,i)
      acc(:,i)=0.0
    enddo
    do i=1,n
      do j=i+1,n
        dist=xyz(:,i)-xyz(:,j)
        r=sqrt(sum(dist**2))
        temp=(G*m(i)*m(j))/(r**3)
        do k=1,3
          !acc(k,i)=acc(k,i)+ temp*dist(k)
          acc(k,j)=acc(k,j) + temp*dist(k)
        enddo
      enddo
    enddo
  end subroutine accUpd
  subroutine velUpd()
    implicit none
    integer :: i,k
    do i=2,n
      do k=1,3
        vel(k,i)=vel(k,i) + 0.5*(acc(k,i)+tempacc(k,i))*dt
      enddo
    enddo
  end subroutine velUpd
end module update
module energy
  use var, only: xyz,vel,m,ke,pe,n,G
  implicit none
  private
  real, dimension(3) ::dist
  real :: r
  integer :: i,j
  public :: kinetic,potential
contains
  subroutine kinetic()
    implicit none
    do i=1,n
      ke=0.5*m(i)*sum(vel(:,i)*vel(:,i))
    enddo
  end subroutine kinetic
  subroutine potential()
    real :: temp
    pe=0.0
    do i=1,n-1
      do j=i+1,n
        dist=xyz(:,i)-xyz(:,j)
        r=sqrt(sum(dist**2))
        pe= pe + (G*m(i)*m(j))/r       !r is relative distance.
      enddo
    enddo
  end subroutine potential
end module energy
program planet
  use var, only : xyz,vel,m,n,rad,ke,pe,tot_t,dt
  use energy, only : kinetic,potential
  use update, only : posUpd,velUpd,accUpd
  use initial, only: init
  implicit none
  integer :: i,t,iter
  character(len=30) :: fmt
  open(100,file="xyz.dat",status="replace")
  open(200,file="vel.dat",status="replace")
  open(300,file="acc.dat",status="replace")
  open(400,file="energy.dat",status="replace")
  open(500,file="params.dat",status="replace")
  call init()
  print*, "The Simulation is running."
  iter=int(tot_t/dt)
  do t=1,iter
    if(t==1) then
      call kinetic()
      call potential()
      write(400,*) t,ke,pe,ke+pe
    endif
    if(mod(t,100)==0) then
      call kinetic()
      call potential()
      write(400,*) t,ke,pe,ke+pe
    endif
    do i=1,n
      write(100,*) xyz(:,i)
    enddo
    call posUpd()
    call accUpd()
    call velUpd()
  enddo
  write(500,*) n
  write(500,*) rad
  call kinetic()
  print*, "The Kinetic Energy is ", ke
  call potential()
  print*, "The Potential Energy is ", pe
end program planet
