module var
  implicit none
  double precision, parameter :: mScale = 1.988544d30   ! Mass of Sun
  double precision, parameter :: lScale = 1.49597870700d11   ! 1 AU
  double precision, parameter :: tScale = 86400.d0    ! Mean Solar day
  double precision, parameter :: G=0.0002959122083d0
  double precision, parameter :: mEarth = 3.0024584d-6
  double precision, parameter :: mSun = 1.d0
  ! double precision, parameter :: mMoon = 0.07342d24/lScale
  integer, parameter :: n=2
  double precision, dimension(3,n) :: xyz,vel,acc
  double precision, dimension(3,n) :: angmom
  double precision, parameter, dimension(n) :: m=[mSun, mEarth]
  double precision, dimension(3) :: orig
  double precision, parameter :: tot_t=365.d0
  double precision, parameter :: dt=0.001d0
  double precision, parameter, dimension(n) :: rad=[0.01, 0.01]
  double precision :: ke,pe
end module var
module initial
  use var, only: xyz,vel,acc,m,orig,angmom,G,tScale,lScale,mEarth,mSun
  use cross_product, only:cross
  implicit none
  private
  public :: init
contains
  subroutine init()
    double precision :: temp
    !SUN
    xyz(:,1)=[(0.025d11)/lScale,0.d0,0.d0]
    vel(:,1)=[0.d0,(mEarth/mSun)*(30300.d0)*tScale/lScale,0.d0]
    acc(:,1)=[-G*m(2)/(1.496d11/lScale)**2,0.d0,0.d0]
    !EARTH
    xyz(:,2)=[(-1.471d11)/lScale,0.d0,0.d0]
    vel(:,2)=[0.d0,-(30300.d0)*tScale/lScale,0.d0]
    acc(:,2)=[G*m(1)/(1.496d11/lScale)**2,0.d0,0.d0]
    ! !MOON
    ! xyz(:,3)=xyz(:,2)-[(405400.0d3)/lScale,0.d0,0.d0]
    ! vel(:,3)=vel(:,2)-[0.d0,(0.964d3)*tScale/lScale,0.d0]
    ! acc(:,3)=[G*m(2)/(384748.0d3/lScale)**2,0.d0,0.d0]
    ! orig=(m(1)*xyz(:,1)+m(2)*xyz(:,2))/(sum(m))
    ! angmom(:,1)=cross(xyz(:,1),m(1)*vel(:,1))
    ! angmom(:,2)=cross(xyz(:,2),m(2)*vel(:,2))
  end subroutine init
end module initial
module cross_product
  implicit none
contains
  function cross(a, b)
    double precision, dimension(3) :: cross
    double precision, dimension(3), intent(in) :: a, b
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
  double precision, dimension(3,n) :: tempacc
  double precision, dimension(3) :: dist
  public :: posUpd,velUpd,accUpd
contains
  subroutine posUpd()
    implicit none
    integer :: i,k
    do i=1,n
      do k=1,3
        xyz(k,i)=xyz(k,i) + vel(k,i)*dt + (acc(k,i)*(dt**2))/2
      enddo
    enddo
  end subroutine posUpd
  subroutine accUpd()
    implicit none
    integer :: i,j,k
    double precision :: r,temp
    do i=1,n
      tempacc(:,i)=acc(:,i)
      acc(:,i)=0.0
    enddo
    do i=1,n
      do j=i+1,n
        dist=xyz(:,i)-xyz(:,j)
        r=sqrt(sum(dist**2))
        temp=(G*m(i)*m(j))/(r**3)
        do k=1,3
          acc(k,i)=acc(k,i) - temp*dist(k)/m(i)
          acc(k,j)=acc(k,j) + temp*dist(k)/m(j)
        enddo
      enddo
    enddo
  end subroutine accUpd
  subroutine velUpd()
    implicit none
    integer :: i,k
    do i=1,n
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
  double precision, dimension(3) ::dist
  double precision :: r
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
    double precision :: temp
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
  call execute_command_line("start /b python show.py")
end program planet
