module param_module
  implicit none
  save
  integer :: n
end module param_module

program main
  use mpi_lib_ours
  implicit none
  call prepare_mpi_lib
  call Prop
  call finalize_mpi_lib
end program main

subroutine  Prop

  use param_module, only : n
  implicit none
  !Define Variables.  Integer, then real*8, then complex*16

  integer                 :: st
  integer                 :: i,j
  integer                 :: nt                             ! number of time steps.  nx is n.
  real*8                  :: pi                             ! no need for three_sigma
  real*8                  :: dt
  real*8                  :: t_final
  real*8                  :: xmin, dx, box
  real*8                  :: dH, H_max, H_min, H_ave
  real*8,     allocatable :: coeff(:)
  real*8,     allocatable :: x(:)
  real*8,     allocatable :: kin(:,:)                 !real*8, not integer!!!
  real*8,     allocatable :: v(:,:)
  real*8,     allocatable :: e(:,:)                         ! E: unity
  real*8,     allocatable :: H(:,:), H_norm(:,:)
  real*8,     allocatable :: H_vc(:,:), H_evl(:)
  complex*16, parameter   :: ci=(0d0,1d0)
  complex*16, allocatable :: psi_0(:), psi_t(:)

  call read_param   
  call set_pi       
  call make_unity   
  call set_grid
  call make_kin     
  call make_v       
  call make_h       
  call open_evl
  call diag_h
  call make_coeff
  call make_psi0

contains
    
  subroutine read_param  
    implicit none
    open(7,file='param.inp',status='old',err=99) 
    rewind(7)
    call fetch_i(n     ,'Ngrid   ',7)    
    call fetch_r(dt    ,'dt      ',7)
    call fetch_i(nt    ,'nt      ',7)
    call fetch_r(xmin  ,'xmin    ',7)
    call fetch_r(dx    ,'dx      ',7)
    close(7)
    return
99  continue
    stop ' param.inp missing '
  end subroutine read_param

  subroutine make_unity
    implicit none
    allocate(e(n,n), stat=st); if(st/=0) stop ' e alloc problems '
    e=0d0
    do i=1,n
       e(i,i) = 1d0
    end do
  end subroutine make_unity

  subroutine set_grid
    implicit none
    allocate(x(n), stat=st); if(st/=0) stop 'x alloc problems '
    do i=1,n
      x(i) = xmin + (i-1)*dx
!      x(i) = 1.5d0 * cos( (2d0*i-1) * pi /2d0/n)
    end do
    box = n*dx
    write(6,*)' box size ',box
  end subroutine set_grid

  subroutine set_pi
    implicit none
    pi = dacos(-1d0)
  end subroutine set_pi

  subroutine make_kin
    implicit none
    allocate(kin(n,n), stat=st);    if(st/=0) stop ' kin alloc. problems '
    
    do i=1,n
      do j=1,n
        if(i==j)then
          kin(i,j) = 1d0/(dx ** 2)
        else if(abs(i-j)==1) then
          kin(i,j) = -0.5d0/(dx ** 2)
        else
          kin(i,j) = 0d0
        end if
      end do
    end do
  end subroutine make_kin

  subroutine make_v
    implicit none
    allocate(v(n,n), stat=st);    if(st/=0) stop ' v alloc. problems '
    v = 0d0
    do i=1, n
      v(i,i) = pot(x(i))
    end do
  end subroutine make_v

  subroutine make_h
    implicit none
    real*8  ::  diff
    allocate(H(n,n), stat=st);          if(st/=0) stop ' h alloc. problems '
    H = kin+v
    diff=maxval(abs(H-transpose(H)))
    if(diff>1d-8) then
      write(*,*)"H is not a Hermitian Matrix, Error"
      stop
    end if
  end subroutine make_h

  subroutine open_evl
    implicit none
    open(12,file='EigenValue.txt',status='replace')
  end subroutine open_evl
 
  subroutine diag_h
    use mat_module, only : mat_diag
    implicit none
    allocate(H_vc(n,n), H_evl(n), stat=st)
    if(st/=0) then
      stop 'H_diag problem'
    end if
    call mat_diag(H, H_vc,H_evl)
    write(12,*)  real(H_evl)
    call flush(12)
    close(12)
  end subroutine diag_h

  subroutine make_coeff
    implicit none
    allocate(coeff(n),stat=st); if(st/=0) stop 'coeff(n) alloc problem'
    do i=1,n
      coeff(i) = Psi(x(i))
    enddo
  end subroutine make_coeff

  subroutine make_psi0
    implicit none
    real*8    :: evl(n)
    allocate(psi_0(n),psi_t(n),stat=st); if(st/=0) stop 'psi_0 alloc problem'
    psi_0 = coeff
    evl = MATMUL(H,psi_0)
    do i=1,n
      write(*,*) 'evl/psi0:',evl(i)/psi_0(i)
    enddo
  end subroutine make_psi0

  subroutine open_plotfile  
    open(11,file='dens_t.txt',status='replace')
  end subroutine open_plotfile
  
  function pot(xv)
    implicit none
    real*8    :: xv, pot, kk
    kk = 4d0* pi**2
    pot = 0.5d0 * kk * xv**2
  end function pot
    
  
  function Psi(x)
    implicit none
    real*8      :: x,Psi
    
    Psi = 2d0**(0.25d0) * exp(-PI*x**2)
  end function Psi

end subroutine Prop
