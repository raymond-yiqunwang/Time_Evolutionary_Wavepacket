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

  integer                 :: st,it
  integer                 :: i,j
  integer                 :: nt                             ! number of time steps.  nx is n.
  real*8                  :: pi                             ! no need for three_sigma
  real*8                  :: dt
  real*8                  :: t_final
  real*8                  :: xmin, dx, box
  real*8                  :: dH, H_max, H_min, H_ave
  real*8,     allocatable :: x(:)
  real*8,     allocatable :: kin(:,:)                 !real*8, not integer!!!
  real*8,     allocatable :: v(:,:)
  real*8,     allocatable :: e(:,:)                         ! E: unity
  real*8,     allocatable :: H(:,:), H_norm(:,:)
  real*8,     allocatable :: H_vc(:,:), H_evl(:)
  complex*16, parameter   :: ci=(0d0,1d0)
  complex*16, allocatable :: u(:,:) 
  complex*16, allocatable :: Chebyshev(:,:)
  complex*16, allocatable :: prev1(:,:), prev2(:,:)
  complex*16, allocatable :: psi_0(:), psi_t(:)

  call read_param   
  call set_pi       
  call make_unity   
  call set_grid
  call set_times    
  call make_kin     
  call make_v       
  call make_h       
  call open_evl
  call diag_h
  call make_h_norm
  call make_u
  call make_psi0
  call normalize(psi_0,n,dx)
    do it = 1, 30 
    call prop_1step
    call normalize(psi_t,n,dx)
  enddo
  call open_plotfile
  call plot
  call dealloc_all

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

  subroutine fetch_i(num,str,file_id)
    implicit none
    integer          :: file_id, num
    integer          :: status
    character(len=8) :: str, string
    logical          :: success    
    status = 0
    success = .false.
    rewind(file_id)
    fetching : do    !but give it a name to make it easier to udnerstand what you "exit" from
      read(file_id,*,iostat=status) string
      if(status/=0) then
        rewind(file_id)
        exit fetching  
      end if
      if(TRIM(string)==TRIM(str))then  
        read(file_id,*)string
        read(string,*) num
        success = .true.
        exit fetching  
      end if
    end do fetching
    if(.not.success) then
       write(6,*)' did not succeed reading ',str,' therefore stopping '
       call flush(6)
       stop
    end if
  end subroutine fetch_i

  subroutine fetch_r(flt,str,file_id)   ! change in the same way fetch_i was changed.
    implicit none
    integer          :: file_id, status
	  real*8           :: flt
    character(len=8) :: str,string
    logical          :: success
    status = 0
    success = .false.
    rewind(file_id)
    fetching : do
      read(file_id,*,iostat=status) string
      if(status/=0) then
        rewind(file_id)
        exit fetching
      end if
      if(TRIM(string)==TRIM(str)) then
        read(file_id,*) string
        read(string,*) flt
        success = .true.
        exit fetching
      end if
    end do fetching
    if(.not. success) then
      write(6,*) 'did not succeed reading ',str,'therefore stopping'
      call flush(6)
      stop
    end if
  end subroutine fetch_r

  subroutine make_unity
    implicit none
    allocate(e(n,n), stat=st); if(st/=0) stop ' e alloc problems '
    e=0d0
    do i=1,n
       e(i,i) = 1d0
    end do
  end subroutine make_unity

  subroutine set_times
    implicit none
    t_final = dt*nt
    write(6,*)' t_final: ',t_final
  end subroutine set_times

  subroutine set_grid
    implicit none
    allocate(x(n), stat=st); if(st/=0) stop 'x alloc problems '
    do i=1,n
      x(i) = xmin + i*dx
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
 
  subroutine make_H_norm
    implicit none
    allocate(H_norm(n,n), stat=st); if(st/=0) stop ' h_norm alloc. problems '
    H_max = maxval(H_evl)
    H_min = minval(H_evl)
    dH = 0.5d0 * (H_max - H_min)
    H_ave = 0.5d0 * (H_max + H_min)
    H_norm = (H - H_ave * e)/dH
  end subroutine make_H_norm
 
  subroutine make_u
    implicit none
    allocate(u(n,n), Chebyshev(n,n), stat=st); if(st/=0) stop ' u alloc. problems '
    allocate(prev1(n,n), prev2(n,n), stat=st); if(st/=0) stop ' prev alloc. problems'
    u = 0d0; prev2 = e; prev1 = H_norm
    do i=0,n
      if(i == 0) then
        Chebyshev = 0.5d0*e
      elseif(i == 1) then
        Chebyshev = H_norm
      else
        call make_Chebyshev
      endif
      u = u + a_coeff(i,dH*dt) * Chebyshev
    enddo
    u = shift(H_ave, dt) * u
  end subroutine make_u

  subroutine make_Chebyshev
    implicit none
    Chebyshev = 2d0 * MATMUL(H_norm,prev1) - prev2
    prev2 = prev1
    prev1 = Chebyshev
  end subroutine make_Chebyshev

  subroutine make_psi0
    implicit none
    integer           :: i, check
    real*8            :: num
    character(len=8)  :: string
    open(13,file='random_vec',status='old')
    allocate(psi_0(n), psi_t(n), stat=st); if(st/=0) stop 'psi_0/psi_t alloc problems'
    i = 1
    fetching: do
      read(13,*,iostat=check) string
      read(string,*) num
      psi_0(i) = num
      if(check/=0) then
        rewind(13)
        exit fetching
      endif
      i = i + 1
    enddo fetching
    call flush(13)
    close(13)
    psi_t = psi_0
  end subroutine make_psi0

  subroutine prop_1step
    implicit none
    psi_t = MATMUL(u,psi_t)
  end subroutine prop_1step

  subroutine open_plotfile  
    open(11,file='dens_t.txt',status='replace')
  end subroutine open_plotfile

  subroutine plot
    implicit none
    write(11,*) abs(psi_t)**2
    call flush(11)
    close(11)
  end subroutine plot
  
  subroutine dealloc_all
    implicit none
    deallocate(e, kin, v, u, x)
    deallocate(H_evl, Chebyshev, prev1, prev2)
    deallocate(H, H_vc, H_norm)
    deallocate(psi_0, psi_t)
  end subroutine dealloc_all

  function pot(xv)
    implicit none
    real*8    :: xv, pot, kk
    kk = 4d0* pi**2
    pot = 0.5d0 * kk * xv**2
  end function pot
    
  subroutine normalize(a,n,dx)   ! lower case. use n,m etc. for integer, not b. and give dx
    implicit none
    integer    :: n
    real*8     :: dx  ! no need for L.  and also align the ::.
    complex*16 :: a(n)
    a(:) = a(:)/sqrt(dx*sum(abs(a)**2))    ! need to properly normalize with dx; this is more compact
  end subroutine normalize

  function shift(h, t)
    implicit none
    real*8      :: h, t
    complex*16  :: shift
    shift = exp(-ci*t*h)
  end function shift
    
  function a_coeff(i,rl)  
    implicit none
    integer       :: i
    real*8        :: rl
    complex*16    :: a_coeff
    a_coeff = 2d0 * (-ci)**i * BESSEL_JN(i,rl)
  end function a_coeff

end subroutine Prop
