program CrossValidationPCE
  use dimpce,only:id_proc
  implicit none  
  include 'mpif.h'

  integer,parameter::ndim=2
  integer::fct,initpts,ncyc,npts,stat,nseed,ierr
  real*8::DAT,cverrout
  real*8::Tsample(ndim,500)
  real*8::testpoint(ndim)
  real*8::CVE(500),MaxCVE
  real*8::MCVE
  integer::j,k,i,knot,nidx,nfunc
  double precision, allocatable, dimension(:,:) :: sample
  integer::nmodel ,order,nterms,nptstmp

  call MPI_START

  if (id_proc.eq.0) then
     print*,''
     print*,'=========================='
     print*,' CROSS VALIDATION PCE '
     print*,'=========================='
     print*,''
  end if

  DAT=77 !screen or file

  stat=0

  do nfunc=1,3

     if (nfunc.eq.1) fct=2
     if (nfunc.eq.2) fct=4
     if (nfunc.eq.3) fct=6

     if (id_proc.eq.0) then
        print*,''
        print*,'=========================='
        print*,'Function number ',fct
        print*,'=========================='
        print*,''
     end if

     !2=runge,4=exp,6=rosen

     if (id_proc.eq.0) then

        if (fct.eq.2) open(unit=37,file='CVPCEfct02dim2.his',form='formatted',status='replace')
        if (fct.eq.4) open(unit=37,file='CVPCEfct04dim2.his',form='formatted',status='replace')
        if (fct.eq.6) open(unit=37,file='CVPCEfct06dim2.his',form='formatted',status='replace')
     end if

     if (id_proc.eq.0) write(37,'(a)') 'npts Order  MeanCVE MaxCVE NMODELS'

     nmodel=0

     do order=2,10

        call combination(nDIM+order,nDIM,nterms)
        ! Get number of points based on stat,solver,oversamp ratios
        call getnpts(2,stat,ndim,nterms,2,nptstmp) 

        do npts=nptstmp,nptstmp

           allocate(sample(ndim,npts))

           if (id_proc.eq.0) then

              sample(:,1)=0.5
              call get_seed(nseed)
              call latin_random(ndim,npts-1,nseed,sample(:,2:npts))       

           end if

           call MPI_BCAST(sample,npts*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

           ! We now have a set of npts samples (first one is the center of domain)

           CVE(:)=0.0
           MaxCVE=0.0

           do k=1,npts !loop over all points to find the CVE and finally MCVE

              testpoint(:)=sample(:,k) !k-th point is the test point
              knot=k

              nidx=0

              do i=1,npts

                 ! Check for all points and see if they dont match the test point and store it in the training vector to be passed to kriging

                 if (i.ne.knot) then
                    nidx=nidx+1
                    Tsample(:,nidx)=sample(:,i)
                 end if

              end do

              ! Now we have the test and training points

              cverrout=0.0
              nmodel=nmodel+1

              call PCestimate(ndim,ndim,fct,DAT,order,npts-1,Tsample(:,1:npts-1),testpoint,cverrout)
              CVE(k)=cverrout
              if (CVE(k).gt.MaxCVE) MaxCVE=cve(k)
           end do

           MCVE=0.0
           do j=1,npts
              MCVE=MCVE+CVE(j)
           end do
           MCVE=MCVE/dble(npts)


           if (id_proc.eq.0) write(37,'(2i8,2e15.5,i8)'), npts, Order, MCVE, maxCVE,nmodel
           if (id_proc.eq.0) print*,'NPTS:',npts,'MCVE:',MCVE,maxCVE,nmodel


           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           deallocate(Sample)

        end do  ! loop npts

     end do !loop order

     if (id_proc.eq.0) close(37)

  end do  !  loop function

  call stop_all

end program CrossValidationPCE
