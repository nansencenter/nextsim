!$id: $ 
module m_pseudo2D
implicit none

#ifdef AIX
   integer :: naux1,naux2,naux3,nn,s1
   double precision, dimension(:), allocatable,save :: aux1,aux2,aux3
#endif


contains
subroutine pseudo2D(Amat,nx,ny,lde,rh,n1,n2)
#if defined (QMPI)
   use qmpi, only: master, stop_mpi
#else
  implicit none

  logical, parameter :: master = .true.
  integer, parameter :: qmpi_num_proc = 1
  integer, parameter :: qmpi_proc_num = 0

contains

subroutine stop_mpi()
   stop
end subroutine stop_mpi

#endif
! This routine calculates the pseudo random fields using
! the procedure outlined in Evensen (1994) \cite{eve94a}.

#ifdef DEC
!   use mydxml
#endif
   use m_zeroin
#ifdef LINUX
   use m_cmplx2real
#endif
   implicit none
   integer, intent(in) :: nx,ny           ! horizontal dimensions
   integer, intent(in) :: lde             ! number of fields stored at the time
   real, intent(out)   :: Amat(nx,ny,lde) ! generated random fields
   real, intent(in)    :: rh              ! Horizontal decorrelation length
   integer, intent(in) :: n1,n2           ! horizontal dimensions in fft grid

#ifdef CRAY
   real, allocatable, dimension(:), save :: work 
   real, allocatable, dimension(:), save :: table 
#endif

#ifdef DEC
   integer status
   record /dxml_d_fft_structure_2d/ fft_struct
#endif

#ifdef SGI
   real, allocatable, dimension(:), save :: coeff
#endif

#if defined(IA32) && defined(FFTW)
   integer*8, save :: plan
   real*8 :: fftwy(n1,n2)
   include 'fftw3.f'
#endif


   real, save ::  rh_save=0.0  ! saving rh used in preprocessing.  Allow for new call to
                               ! zeroin if rh changes.

   real, save :: sigma,sigma2
   real, save :: c
   integer, save :: n1_save=0
   integer, save :: n2_save=0


   integer l,p,j,m,i
   real kappa2,lambda2,kappa,lambda
   real pi2,deltak,accum,factor
   real a1,b1,tol,fval

   real, allocatable    :: fampl(:,:,:)
   real, allocatable    :: phi(:,:)
   real, allocatable    :: y(:,:)   ! Physical field
   complex, allocatable :: x(:,:)   ! Fourier amplitudes
!#ifndef ICE2
#ifdef LINUX
   complex, allocatable :: speq(:)  ! Nyquist freqencies (LINUX)
#endif 

   real, parameter :: dx=1.0
   real, parameter :: pi=3.141592653589

   real, external :: func2D

   if (lde < 1)    stop 'pseudo2D: error lde < 1'
   if (rh <= 0.0)  stop 'pseudo2D: error, rh <= 0.0'
   if (n1 < nx)    stop 'pseudo2D: n1 < nx'
   if (n2 < ny)    stop 'pseudo2D: n2 < ny'

   allocate(fampl(0:n1/2,-n2/2:n2/2,2))
   allocate(phi(0:n1/2,-n2/2:n2/2))
#ifdef LINUX
   allocate(speq(n2))
   allocate(y(n1,n2))
   allocate(x(n1/2,n2))
#else
   allocate(y(0:n1+1,0:n2-1))
   allocate(x(0:n1/2,0:n2-1))
#endif

#ifndef LINUX
#ifndef CRAY
#ifndef AIX
#ifndef DEC
#ifndef SGI
#if !defined(IA32) || !defined(FFTW)
   print *,'ranfield is only running on the following machines:'
   print *,'   LINUX having Numerical Recipes'
   print *,'   LINUX having FFTW             '
   print *,'   CRAY'
   print *,'   AIX having essl'
   print *,'   DEC having dxml'
   print *,'   SGI'
   call stop_mpi()
#endif /*IA32 || FFTW*/
#endif /*SGI*/
#endif /*DEC*/
#endif /*AIX*/
#endif /*CRAY*/
#endif /*LINUX*/

   pi2=2.0*pi
   deltak=pi2**2/(float(n1*n2)*dx*dx)
   kappa=pi2/(float(n1)*dx)
   kappa2=kappa**2
   lambda=pi2/(float(n2)*dy)
   lambda2=lambda**2
   factor=1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization.
   if (rh /= rh_save .or. n1 /= n1_save .or. n2 /= n2_save) then
      rh_save=rh
      n1_save=n1
      n2_save=n2
#ifdef SGI
      if (mod(n1,2) /= 0) print *,'pseudo2D: mod(n1,2) must be zero, n1=',n1
      if (mod(n2,2) /= 0) print *,'pseudo2D: mod(n2,2) must be zero, n2=',n2
      if(allocated(coeff)) deallocate(coeff); allocate( coeff((n1+15) + 2*(n2+15)) )
      call dzfft2dui(n1,n2,coeff)
#endif /*SGI*/
#ifdef CRAY
      if(allocated(work)) deallocate(work)
      allocate(work(512*n1))

      if(allocated(table)) deallocate(table)
      allocate(table(100+2*(n1+n2)))

      call scfft2d(0,n1,n2,factor,x,n1/2+1,y,n1+2,table,work,0)
#endif /*CRAY*/
#ifdef AIX
      nn   = max(n1/2,n2)

      if (nn<=2048) then
         naux1= 42000 
      else
         naux1= ceiling(40000+1.64*n1+2.28*n2)
      end if

      if (n1 <= 4096 ) then
         naux2 = 20000
      else if (n1 > 4096 ) then 
         naux2 = ceiling(20000+1.14*n1)
      end if

      if ( n2 > 252) then 
         s1 = min(64, 1+n1/2)
         naux2 = naux2 + ceiling((2*n2+256)*(2.28+s1))
      end if

      naux3=1

      if(allocated(aux1)) deallocate(aux1); allocate(aux1(naux1))
      if(allocated(aux2)) deallocate(aux2); allocate(aux2(naux2))
      if(allocated(aux3)) deallocate(aux3); allocate(aux3(naux3))

      !print *,n1,n2
      !print *,factor
      !print *,naux1,naux2,naux3
      call dcrft2(1,x,n1/2+1,y,n1+2,n1,n2,-1,factor,&
               aux1,naux1,aux2,naux2,aux3,naux3)
#endif /*AIX*/
#ifdef DEC
      if (mod(n1,2) /= 0) then
         print *,'ranfield: n1 is not even. n1=',n1
      endif
      status=dfft_init_2d(n1,n2,fft_struct,.true.)
      if (status /= 0 ) print *,'status: dfft_init_2d',status
#endif /*DEC*/
#if defined(IA32) && defined(FFTW)
      if (master) then
         print *,'Using FFTW for fourier transform'
         print *,'Feel the power of the Fastest Fourier Transform in the West!'
      end if
#endif /*IA32 && FFTW*/

      rh_save=rh
      if (master) print '(a,2f6.2)','pseudo2D: Solving for sigma',rh,dx
      a1=0.1e-07
      b1=0.1e-06
      tol=0.1e-10
      if (master) print *,'pseudo2D: Go into  zeroin'
      call zeroin(func2D,sigma,a1,b1,tol,rh,dx,fval,n1,n2)
      if (master) print *,'pseudo2D: Leaving  zeroin'

      sigma2=sigma**2
      accum=0.0
      do p=-n2/2+1,n2/2
      do l=-n1/2+1,n1/2
#if defined(EXPCOV)
         accum=accum+1./((1.+(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)**3.)
#else
         accum=accum+exp(-2.0*(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)
#endif
      enddo
      enddo
      if (master) print*,'pseudo2D: Leving do loop'

      c=sqrt(1.0/(deltak*accum))

      if (master) print *,'pseudo2D: sigma  ',sigma
      if (master) print *,'pseudo2D: c=     ',c
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do j=1,lde
      ! Calculating the random wave phases
      call random_number(phi)
      phi=pi2*phi

      ! Calculating the wave amplitues
      do p=-n2/2,n2/2
      do l=0,n1/2 
#if defined(EXPCOV)
         fampl(l,p,1)=&
            ((1.+(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)**(-1.5))*&
            cos(phi(l,p))*sqrt(deltak)*c
         fampl(l,p,2)=&
            ((1.+(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)**(-1.5))*&
            sin(phi(l,p))*sqrt(deltak)*c
#else
         fampl(l,p,1)=&
            exp(-(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)*&
            cos(phi(l,p))*sqrt(deltak)*c
         fampl(l,p,2)=&
            exp(-(kappa2*float(l*l)+lambda2*float(p*p))/sigma2)*&
            sin(phi(l,p))*sqrt(deltak)*c
#endif
      enddo
      enddo
      fampl(0,0,2)=0.0

#ifdef LINUX
      do l=1,n1/2
        do p=1,n2/2+1
          x(l,p)=cmplx(fampl(l-1,p-1,1),fampl(l-1,p-1,2))
        enddo
        do p=n2/2+2,n2
          x(l,p)=cmplx(fampl(l-1,-n2+p-1,1),fampl(l-1,-n2+p-1,2))
        enddo
      enddo

      do p=1,n2/2+1
        speq(p)=cmplx(fampl(n1/2,p-1,1),fampl(n1/2,p-1,2))
      end do
      do p=n2/2+2,n2
        speq(p)=cmplx(fampl(n1/2,-n2+p-1,1),fampl(n1/2,-n2+p-1,2))
      end do

      call rlft3(x,speq,n1,n2,1,-1)
      call cmplx2real(x,y,(n1/2)*n2)

      do m=1,ny
       do i=1,nx
         Amat(i,m,j)=2.0*y(i,m)
       enddo
      enddo 

#endif /*LINUX*/

#ifndef LINUX
      do p=0,n2/2-1
         x(:,p)=cmplx(fampl(:,p,1),fampl(:,p,2))
      enddo

      do p=n2/2,n2-1
         x(:,p)=cmplx(fampl(:,-n2+p,1),fampl(:,-n2+p,2))
      enddo

#ifdef CRAY
      call csfft2d(-1,n1,n2,factor,x,n1/2+1,y,n1+2,table,work,0)
#endif
#ifdef SGI
      call zdfft2du(-1,n1,n2,x,n1+2,coeff)
      y=reshape(transfer(x,(/0.0 /) ),(/n1+2,n2/))
#endif

#ifdef AIX
      !print *,n1,n2
      !print *,factor
      !print *,naux1,naux2,naux3
      call dcrft2(0,x,n1/2+1,y,n1+2,n1,n2,-1,factor,&
               aux1,naux1,aux2,naux2,aux3,naux3)
      !print *,'ok'
#endif

#ifdef DEC
      status=dfft_apply_2d('C','R','B',x,y,n1+2,fft_struct,1,1)
      if (status /= 0 ) print *,'status: dfft_apply_2d',status
      y=y*float(n1*n2)
#endif
#if defined(IA32) && defined(FFTW)
      !print *,'IA32 fft ...',nx,ny,n1,n2
      call dfftw_plan_dft_c2r_2d(plan,n1,n2,x,fftwy,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
      !print *,'IA32 fft done...'
      y(0:n1-1 ,0:n2-1)=fftwy(1:n1,1:n2)
      y(n1:n1+1,0:n2-1)=fftwy(1:2 ,1:n2)
#endif

      do m=1,ny
      do i=1,nx
         Amat(i,m,j)=y(i-1,m-1)
      enddo
      enddo

!!! ifndef LINUX  !!!!!
#endif

   enddo


   deallocate(fampl, phi, y, x)
#ifdef LINUX
   deallocate(speq)
#endif

#ifdef DEC
   status=dfft_exit_2d(fft_struct)
   print *,'status: dfft_exit_2d',status
#endif

   ! Test
!   call tecfld('rand',Amat(:,:,1),nx,ny,Amat(:,:,lde))
!   call vario_regugrid('randomFFT',Amat(:,:,:3),nx,ny,3,rh)

end subroutine pseudo2D
end module m_pseudo2D
