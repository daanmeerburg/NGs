
!this program aims to compute the S/N from the bispectra generated
!using alpha beta gamma delta from komatsu's code; PDM 2018
!Can be used to make forecast for SO, CMBS4, PICO and Planck Legacy 

!tested on CITA machine lobster/homard. Takes 14 seconds for lmax = 500, 6 minutes for lmax = 2000 and 190 minutes for lmax = 5000 for local
!longer for equilateral and orthogonal. So be wise...

program BS_SN

  !wignerJ lib
  use fwigxjpf
  implicit none

  !params
  integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359

  !opening file params 

  real(dl), allocatable:: rarray(:), deltar(:)
  integer :: dimr

  !to check 'heat map' of a,b,g,d integrals:
  logical :: test_rintegral  = .false. 

  character(120) :: temp_bisp, temp_cov, Covfile
  character(120) :: alphaTfile, betaTfile, alphaEfile, betaEfile, rfile
  character(120) :: gammaTfile, deltaTfile, gammaEfile, deltaEfile
  character(70) :: alphabetafile, gammadeltafile, header

  !alpha. beta files:

  real(dl), allocatable :: ar(:,:,:), br(:,:,:), gr(:,:,:), dr(:,:,:)
  integer :: nfields, minfields 
  character(6) :: form = '(I3)'


  !bispectra
  real(dl), allocatable :: BredTemp(:,:,:), Blens(:,:,:)
  real(dl) :: Btemp
  real(dl) :: abint
  integer :: local = 1, equil = 2, ortho = 3, folded = 4, flat = 5 !not yet implemented
  character(7) :: nshape
  integer :: shape


  !loop params
  integer :: L1, L2, L3
  integer :: i, p, j, k, q, r, s

  integer :: ell
  real(dl) :: rell !for cmbs4 noise file...
  character(6) :: clmax
  integer :: Fishlmax
  integer :: lmax, lmin, max_l, min_l
  integer :: FileUnit

  !limiting r-integration range
  real(dl) :: rmin, rmax
  integer :: irmin, irmax, ch1, ch2
  logical :: lim_r_samp = .true.

  !wigner 3j
  real(dl)  :: a3j(0:20000)
  real(dl)  :: a3j2(0:20000,3,2)
  !these are needed only when doing r-integral for speed up
  real(dl),allocatable  :: sa3j(:,:)
  integer :: l1l2index

  !version control
  character(3) :: vers, dat

  !comparing bispectrum covariance code to fisher
  logical :: comp_code = .false.
  real(dl) :: alpha, beta


  !cov data
  real(dl) :: l2CLTT, l2CLTE, l2CLEE
  real(dl), allocatable :: CL(:,:,:), CLtil(:,:,:)
  real(dl), allocatable :: invCov(:,:,:), detcov(:)
  real(dl), allocatable :: invCovCV(:,:,:), detcovCV(:)

  !experiments 
  integer :: exper, SO = 1, CMBS4 = 2, planck  = 4, PICO = 3
  character(6) :: cexp

  !noise
  real(dl), allocatable :: NL(:,:,:)
  real(dl) :: Noise(2), sigma2(2)
  character(120) :: NoiseFileT, NoiseFileEB, NoiseFileCMBS4, NoiseFilePlanck
  real(dl) :: CMB2COBEnorm = 7428350250000.d0
  real(dl) :: fsky  = 0.4
  real(dl) :: fskyCV = 1.
  character(40) :: DefT = 'SOV3_T_default1-4-2_noisecurves' 
  character(40) :: Defpol = 'SOV3_pol_default1-4-2_noisecurves'
  character(40) :: DefTS4 = 'S4_190604d_2LAT_T_default_noisecurves' !new curves, as of June 11, 2019
  character(40) :: DefpolS4 = 'S4_190604d_2LAT_pol_default_noisecurves'
  character(20) :: deproj, SENS, skyfrac

  !Fisher sums
  real(dl) :: TempCov, TempCovCV, Det
  real(dl) :: TotSum, TotSumCV

  !Fisher contributions
  real(dl) :: FishC(2,2,2,2,2,2), FishC_CV(2,2,2,2,2,2)

  !ISW -lensing
  real(dl) :: BLensISW, DetISWLens,DetLensCross
  real(dl) :: TotSumISWLens, TotSumCVISWLens
  real(dl) :: TotSumLensCross, TotSumCVLensCross
  real(dl) :: btlenstemp

  !lensing covariance
  logical :: want_lensingCV
  logical :: use_lensed_spectra

  real(dl) :: l2CLPP, l2CLPT, l2CLPE, l2CLBB
  real(dl), allocatable :: CLISW(:,:,:)
  character(120) :: ISWfile
  real(dl) :: TempCovISW,TempCovCVISW
  real(dl), allocatable :: invCovISW(:,:,:)
  real(dl), allocatable :: invCovCVISW(:,:,:)

  real(dl) :: LensNorm 

  real(dl) :: DetFish, DetFishCV
  real(dl) :: DefnlMar, DefnlMarCV

  !reconstruction noise params 
  character(120) :: TTreconstionNoiseFile
  logical :: use_recon
  real*8 val3j

  !SW approximation for fnl_local
  logical :: want_SW_template = .false.
  real(dl) :: norm


  !if you want to compare output of bispectrum covariance code
  !you get some extra prints

  comp_code =  .false. 
  
  !shape: local = 1, equilateral = 2, orthogonal  = 3, folded = 4, flat = 5
  shape = 2
  !experiment, SO = 1, CMBS4 = 2, PICO = 3, Planck = 4
  exper = 1

  !for SO and new S4:
  !deproj0, deproj1, deproj2, deproj3
  deproj = 'deproj0'
  !SENS0 (threshold), SENS1 (beseline), SENS2 (goal)
  sens = 'SENS2'
  !4000, 8000, 16000; 10-20-40%
  skyfrac = '16000'

  !lmax: 
  lmax = 5000 !==> File ell max (allocation)
  lmin = 2
  !you can change the Fisher-lmax here. Note that you cant change lmax since
  !it would ruiin reading the a,b,g, and d files
  Fishlmax = 5000 !==> Fisher lmax

  !version (change if you dont want to overwrite, all produced files will be labelled by this); up to 3 characters 
  vers = 'x1'


  !do you want to run the actual local shape or the SW approximation?
  !only works for local and when you run temperature only
  want_SW_template = .false.
  
  !do you want to include lensing covariance?
  !only works with temperature at the moment 
  want_lensingCV = .false.
  !iff you want to use the lensed spectra in covariance 
  use_lensed_spectra = .false.

  !for l < 51, do you want to use reconstruction noise?

  use_recon = .false.
  TTreconstionNoiseFile = 'noise/hatT2_S4_lmax5000.dat'

  !test r-integral instead of fisher?
  test_rintegral = .false.

  !once tested, you can then limit integration range
  !putting new integration range below:
  lim_r_samp = .true.
  
  if(test_rintegral) lim_r_samp = .false.

  FileUnit  = 10
  minfields  = 1 !2 for only E
  nfields  = 2 !2 for T and E, 1 for T only
  if (nfields .eq. 1) minfields  = 1

  if (nfields .eq. 1) then
     write(*,*) "temperature only <TTT>"
     dat = 'TT'
  elseif(minfields .eq. 2 .and. nfields .eq.2) then
     write(*,*) "polarization only <EEE>"
     dat = 'EE'
     want_lensingCV = .false.
  else
     write(*,*) "temperature and polarization <TTT>, <TTE>, <TEE>, <EEE>"
     dat = 'TE'
     want_lensingCV = .false.
  endif


  !currently using 2018 noise covariance directly from Planck Repo. Noise proabbly a little conservative (e.g. fsky not accounted for correctly)
  !lmax  = 1996, while for actual Planck analysis lmaxT = 2000 and lmaxE = 1500

  !if want to use 'blue book' Planck noise instead. Change code accordingly below. 
  Noise(1) = 1.07695E-017 !Dimensionless, for Planck T
  Noise(2) = 1.61543E-017 !Dimensionless, for Planck E
  !sigma2(1) = 9.765981E-007 !Planck beam
  Noise(2) = 1.07695E-017 !ACT
  sigma2(1) = 9.765981E-007/(3.333**2) !ACT beam
  if (exper .eq. 1) then
     cexp = 'SO'
     fsky  = 0.4 !should be consistent with noise choice above 
  elseif(exper .eq. 2) then
     cexp = 'S4'
     fsky  = 0.43
  elseif(exper .eq. 3) then
     cexp  = 'PICO'
     fsky  = 0.8
  elseif (exper .eq. 4) then
     cexp = 'planck'
     fsky = 0.4
     if(lmax .ge. 1995) Fishlmax = 3000 !cant be larger then this when using noise files/ can be larger if you use Blue Book
  else
     write(*,*) 'choose experiment 1 = SO, 2 = S4, 3 = PICO'
  endif

  !tested these below. Orthogonal becomes better constrained for narrower range of r. 
  if (shape .eq. 1) then
     nshape = 'local'
     !tested using  test_rintegral = .true.
     !will make some plots/currently tested for lmax  = 2000.
     !gives 1% accuracy 
     rmin = 13805.1
     rmax = 14067.6
  elseif(shape .eq. 2) then
     nshape = 'equil'
     rmin = 13728.1
     rmax = 14137.6
  elseif(shape .eq. 3) then
     nshape = 'ortho'
     rmin = 13675.6
     rmax = 14172.6
  elseif (shape .eq. 4) then
     nshape  = 'folded'
     !most conservative, but needs to be checked
     rmin = 13675.6
     rmax = 14172.6
  elseif (shape .eq. 5) then
     nshape  = 'flat'
     !most conservative, but needs to be checked
     rmin = 13675.6
     rmax = 14172.6
  endif

  !if(lmax .ge. 1000) form = '(I4)'
  !write(clmax, form)  Fishlmax

  !convert integer to character string 
  if(lmax .ge. 1000) form = '(I4)'
  write(clmax, form)  lmax

  !locate alpha,beta,gamma,delta files
  alphabetafile = 'alphabeta/l_r_alpha_beta_new_Lmax'//trim(clmax)//'.txt' 
  gammadeltafile = 'alphabeta/l_r_gamma_delta_new_Lmax'//trim(clmax)//'.txt' 

  !neutrino mass = 0
  if(use_lensed_spectra) then     
     !Covfile ='CAMB/Lmax5500_lensedCls.dat'
     !Covfile ='CAMB/cosmo2017_10K_acc3_lensedCls.dat'
     Covfile = '/mnt/raid-cita/meerburg/Bispectrum_Covariance/SOspectra/SOspectra_lensedCls.dat'
  else
     Covfile ='CAMB/Lmax'//trim(clmax)//'_nu0_scalCls.dat'
  endif

  !ISW-lensing:
  !ISWfile = 'CAMB/Lmax'//trim(clmax)//'_nu0_lenspotentialCls.dat'
  ISWfile = '/mnt/raid-cita/meerburg/Bispectrum_Covariance/SOspectra/SOspectra_lenspotentialCls.dat'
  !Covfile ='CAMB/Lmax5000_nu0_Lmax5000_nu0_lensedCls.dat'
  open(unit=FileUnit,file = alphabetafile, status='old')
  open(unit=FileUnit+1,file = gammadeltafile, status='old')

  !reading in header (nice illustration of the limitations of Fortran)
  !blank line
  write(*,*)
  do i=1,3
     !alpha/beta file, show header:
     read(FileUnit,'(1A70)') header
     if(i .eq. 1 .or. i .eq. 2) write(*,*) trim(header)
  enddo
  write(*,*)
  !uncomment if you want to make code non-interactive. Need to specify lmax and dimr though
  !write(*,'(A60)',advance='no') '# of radial points (enter the number shown above: eg, 603):'
  !read(*,*) dimr
  !write(*,'(A68)',advance='no') 'maximum multipole (enter the number shown above or less: eg, 5000):'
  !read(*,*) lmax

  !if you know what it is an are always using the same file (like me...)
  dimr = 603

  do i=1,3
     !gamma/delta file, skip header:
     read(FileUnit+1,'(1A70)') header
  enddo

  allocate(rarray(dimr), deltar(dimr))

  write(*,'(A6,X,I4)') 'lmax:', lmax
  write(*,'(A6,X,I4)') 'lmin:', lmin
  write(*,'(A6,X,F4.2)') 'fsky:', fsky
  write(*,'(A13,X,I4)') 'Fisher-lmax:', Fishlmax

  allocate(ar(dimr,2:lmax,2))
  allocate(br(dimr,2:lmax,2))
  allocate(gr(dimr,2:lmax,2))
  allocate(dr(dimr,2:lmax,2))

  !read
  write(*,*) 'reading alpha/beta/gamma/delta files ...'

  ch1 = 1
  ch2 = 1
  do l1 = 2, lmax
     do i  = 1, dimr
        read(FileUnit,*) ell, rarray(i), deltar(i), ar(i,l1,1), br(i,l1,1) , ar(i,l1,2), br(i,l1,2)
        read(FileUnit+1,*) ell, rarray(i), deltar(i), gr(i,l1,1), dr(i,l1,1) , gr(i,l1,2), dr(i,l1,2)
        !polarization alpha,beta,gamma,delta need factor ((l-2)!/(l+2)!^)1/2
        br(i,l1,2) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*br(i,l1,2)
        ar(i,l1,2) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*ar(i,l1,2)
        gr(i,l1,2) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*gr(i,l1,2)
        dr(i,l1,2) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*dr(i,l1,2)
        !write(*,*) ell, rarray(i), deltar(i), ar(i,l1,1), br(i,l1,1) , ar(i,l1,2), br(i,l1,2)



     enddo
  enddo
  !if not limiting integration range:
  irmin  = dimr
  irmax = 1
  !if limiting, change irmin and irmax
  if(lim_r_samp .and. .not. test_rintegral) then
     write(*,*) 'limiting integration in r. 99% accuracy'
     do i  = 1, dimr
        !locate rmin and rmax
        if (rarray(i) .le. rmin .and. ch1 .eq. 1) then
           irmin = i
           ch1 = 0
        endif
        if (rarray(i) .le. rmax .and. ch2 .eq. 1) then
           irmax = i
           ch2 = 0
        endif
        !write(*,*) 'using only', irmin - irmax, 'instead of', dimr, 'points' 
     enddo
     write(*,*) 'using only', irmin - irmax, 'instead of', dimr, 'points'
  endif


  !reverse indexing when doing r-check 
!!$  allocate(ar(lmin:lmax,2,dimr))
!!$  allocate(br(lmin:lmax,2,dimr))
!!$  allocate(gr(lmin:lmax,2,dimr))
!!$  allocate(dr(lmin:lmax,2,dimr))
!!$  
!!$  do l1 = lmin, lmax
!!$     do i  = 1, dimr
!!$        read(FileUnit,*) ell, rarray(i), deltar(i), ar(l1,1,i), br(l1,1,i) , ar(l1,2,i), br(l1,2,i)
!!$        read(FileUnit+1,*) ell, rarray(i), deltar(i), gr(l1,1,i), dr(l1,1,i) , gr(l1,2,i), dr(l1,2,i)
!!$        !polarization alpha,beta,gamma,delta need factor ((l-2)!/(l+2)!^)1/2
!!$        br(l1,2,i) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*br(l1,2,i)
!!$        ar(l1,2,i) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*ar(l1,2,i)
!!$        gr(l1,2,i) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*gr(l1,2,i)
!!$        dr(l1,2,i) = sqrt((l1-1d0)*l1*(l1+1d0)*(l1+2d0))*dr(l1,2,i)
!!$        !write(*,*) ell, rarray(i), deltar(i), ar(i,l1,1), br(i,l1,1) , ar(i,l1,2), br(i,l1,2)
!!$     enddo
!!$  enddo

  close(FileUnit) !alpha & beta
  close(FileUnit + 1) !gamma & delta

  allocate(CL(2:lmax,2,2))
  open(unit=FileUnit + 6,file = Covfile, status='old')
  write(*,*) 'reading in power spectra for covariance ...'
  do i = 1, lmax
     if (i .eq. 1) then
        read(FileUnit + 6,*)
        l1 = 2
     else

        read(FileUnit + 6,*) ell, l2ClTT, l2CLEE, l2CLTE
        !convert l(l+1)C_l/2pi [muK^2] to dimensionless C_l
        Cl(l1,1,1) = 2.*pi*l2ClTT/ell/(ell+1.)/CMB2COBEnorm
        Cl(l1,2,2) = 2.*pi*l2ClEE/ell/(ell+1.)/CMB2COBEnorm
        Cl(l1,2,1) = 2.*pi*l2ClTE/ell/(ell+1.)/CMB2COBEnorm
        Cl(l1,1,2) = 2.*pi*l2ClTE/ell/(ell+1.)/CMB2COBEnorm
        l1 = l1 + 1
     endif
  enddo
  close(FileUnit + 6)

  !high ell noise
  if(exper .eq. SO) then 
     !NoiseFileT = 'noise/'//trim(defT)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_TT_yy.txt'
     !NoiseFileEB = 'noise/'//trim(defpol)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_EE_BB.txt'
     NoiseFileT = 'noise/'//trim(defT)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_TT_yy_39GHzfix_ext15yrs.txt'
     NoiseFileEB = 'noise/'//trim(defpol)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_EE_BB_ext15yrs.txt'

     !use for ell < 40 in T en ell < 10 in E. For consistency only
     NoiseFileCMBS4 = 'noise/CMBS4noise.txt'


     open(unit=FileUnit + 7,file = NoiseFileT, status='old')
     open(unit=FileUnit + 8,file = NoiseFileEB, status='old')
     open(unit=FileUnit + 9,file = NoiseFileCMBS4, status='old')

     allocate(NL(2:lmax,2,2))
     allocate(CLtil(2:lmax,2,2))
     write(*,*) 'reading noise SO ...'
     do l1 = 2, lmax
        if (l1 .le. 39) then
           !Planck noise
           !NL(l1,1,1) = Noise(1)*exp((l1+1)*(l1)*sigma2(1))
           read(FileUnit + 9,*) rell, NL(l1,1,1)
           NL(l1,2,1) = 0.d0
           NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
           !adding noise to Cl's
           CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
           CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1)
        else
           !from the ground
           read(FileUnit + 7,*) ell, NL(l1,1,1)
           NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
           NL(l1,2,1) = 0.d0

           !adding noise to Cl's
           CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
           CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1)
        endif
     enddo
     !for SO E en T files have differnt lmin...
     do l1 = 2, lmax
        if (l1 .le. 9) then
           !Planck noise
           !NL(l1,2,2) = Noise(2)*exp((l1+1)*(l1)*sigma2(2))
           read(FileUnit + 9,*) rell, NL(l1,1,1), NL(l1,2,2)
           NL(l1,1,2) = 0.d0
           NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm
           !adding noise to Cl's
           CLtil(l1,2,2) = CL(l1,2,2) + NL(l1,2,2)
           CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2)
        else
           !from the ground
           read(FileUnit + 8,*) ell, NL(l1,2,2)
           NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm
           NL(l1,1,2) = 0.d0
           !adding noise to Cl's
           CLtil(l1,2,2) = CL(l1,2,2) + NL(l1,2,2)
           CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2)
        endif
     enddo
     close(FileUnit + 7)
     close(FileUnit + 8)
     close(FileUnit + 9)
  elseif(exper .eq. CMBS4) then

     NoiseFileT = 'noise/'//trim(defTS4)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_TT_yy.txt'
     NoiseFileEB = 'noise/'//trim(defpolS4)//'_'//trim(deproj)//'_'//trim(SENS)//'_mask_'//trim(skyfrac)//'_ell_EE_BB.txt'

     NoiseFileCMBS4 = 'noise/CMBS4noise_v2.txt'

     open(unit=FileUnit + 7,file = NoiseFileT, status='old')
     open(unit=FileUnit + 8,file = NoiseFileEB, status='old')
     open(unit=FileUnit + 9,file = NoiseFileCMBS4, status='old')

     allocate(NL(2:lmax,2,2))
     allocate(CLtil(2:lmax,2,2))
     write(*,*) 'reading noise S4 ...'
!!$     do l1 = lmin, lmax
!!$        read(FileUnit + 7,*) rell, NL(l1,1,1), NL(l1,2,2)
!!$        NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm 
!!$        NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm 
!!$        NL(l1,2,1) = 0.d0
!!$        NL(l1,1,2) = 0.d0
!!$
!!$        !adding noise to Cl's
!!$        CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
!!$        CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1) 
!!$        CLtil(l1,2,2) = CL(l1,2,2) + NL(l1,2,2) 
!!$        CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2) 
!!$        !write(*,*) ell, NL(l1,1,1), NL(l1,2,2)
!!$     enddo
     do l1 = 2, lmax
        if (l1 .le. 39) then
           !Planck noise
           read(FileUnit + 9,*) rell, NL(l1,1,1)
           NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
           NL(l1,2,1) = 0.d0
           !adding noise to Cl's
           CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
           CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1)
        else
           !from the ground
           read(FileUnit + 7,*) ell, NL(l1,1,1)
           NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
           NL(l1,2,1) = 0.d0

           !adding noise to Cl's
           CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
           CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1)
        endif
     enddo
     !for S4 E en T files have differnt lmin...
     do l1 = 2, lmax
        if (l1 .le. 9) then
           !Planck noise; using Joel Meyers old file 
           read(FileUnit + 9,*) rell, NL(l1,1,1), NL(l1,2,2)
           NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm
           NL(l1,1,2) = 0.d0
           !adding noise to Cl's
           CLtil(l1,2,2) = CL(l1,2,2) + NL(l1,2,2)
           CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2)
        else
           !from the ground
           read(FileUnit + 8,*) ell, NL(l1,2,2)
           !after discussion with S4 noise tiger Colin Hill, using Planck noise up to ell = 30 in E:
           if (l1 .le. 29) then
              read(FileUnit + 9,*) rell, NL(l1,1,1), NL(l1,2,2)
           endif
           NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm
           NL(l1,1,2) = 0.d0
           !adding noise to Cl's
           CLtil(l1,2,2) = CL(l1,2,2) + NL(l1,2,2)
           CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2)
        endif
     enddo

     close(FileUnit + 7)
     close(FileUnit + 8)
     close(FileUnit + 9)
  elseif(exper .eq. planck) then
     allocate(NL(2:lmax,2,2))
     allocate(CLtil(2:lmax,2,2))

     !generated using the 
     NoiseFilePlanck = 'noise/PlanckRealNoise.tsv'

     open(unit=FileUnit + 7,file = NoiseFilePlanck, status='old')
     write(*,*) 'getting noise Planck ...'
     do l1 = 2, lmax !min(lmax,1996) !change if you want to use the file instead of the Blue Book values
        !read(FileUnit + 7,*) ell, NL(l1,1,1), NL(l1,2,2)
        !NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
        NL(l1,1,1) = Noise(1)*exp((l1+1)*(l1)*sigma2(1))
        !NL(l1,2,2) = NL(l1,2,2)/CMB2COBEnorm
        NL(l1,2,2) = Noise(2)*exp((l1+1)*(l1)*sigma2(1))

        NL(l1,2,1) = 0.d0
        NL(l1,1,2) = 0.d0 
        !write(*,*) ell, NL(l1,1,1), NL(l1,2,2)
        !adding noise to Cl's (allready includes Cosmic Variance)
        CLtil(l1,1,1) = NL(l1,1,1) + CL(l1,1,1)
        CLtil(l1,2,1) = CL(l1,2,1) + NL(l1,2,1)
        CLtil(l1,2,2) = NL(l1,2,2) + CL(l1,2,2)
        CLtil(l1,1,2) = CL(l1,1,2) + NL(l1,1,2)
     enddo
     !if(lmax .ge. 1996) CLtil(1996:lmax,1:2,1:2) = 10E10 !large number

     close(FileUnit + 7)
  endif
  !replace low ell with reconstruction noise:

  if (use_recon) then
     write(*,*) 'Beware: using reconstruction noise at ell < 50 ...'
     open(unit=FileUnit + 7,file = TTreconstionNoiseFile, status='old')
     do l1 = 2, 50
        read(FileUnit + 7,*) ell, NL(l1,1,1)
        NL(l1,1,1) = NL(l1,1,1)/CMB2COBEnorm
        CLtil(l1,1,1) = CL(l1,1,1) + NL(l1,1,1)
     enddo
     close(FileUnit + 7)
  endif

  !
  allocate(invcov(2:lmax,2,2))
  allocate(invcovCV(2:lmax,2,2))

  allocate(detCov(2:lmax))
  allocate(detCovCV(2:lmax))

  do l1  = 2, lmax
     if (nfields .eq. 2 .and. minfields .eq. 1) then 
        !determinant of covariance 
        detCov(l1) = (CLtil(l1,1,1)*CLtil(l1,2,2)-CLtil(l1,2,1)**2)
        detCovCV(l1) = (CL(l1,1,1)*CL(l1,2,2)-CL(l1,2,1)**2) 
        !inverse covariance 
        invcov(l1,1,1) = CLtil(l1,2,2)/detCov(l1)
        invcov(l1,1,2) = -CLtil(l1,1,2)/detCov(l1)
        invcov(l1,2,1) = -CLtil(l1,2,1)/detCov(l1)
        invcov(l1,2,2) = CLtil(l1,1,1)/detCov(l1)
        !inverse covariance in CV limit 
        invcovCV(l1,1,1) = CL(l1,2,2)/detCovCV(l1)
        invcovCV(l1,1,2) = -CL(l1,1,2)/detCovCV(l1)
        invcovCV(l1,2,1) = -CL(l1,2,1)/detCovCV(l1)
        invcovCV(l1,2,2) = CL(l1,1,1)/detCovCV(l1)
     elseif(minfields .eq. 2 .and. nfields .eq. 2) then 
        !when only using E-mode Polarization
        invcov(l1,2,2) =  1./CLtil(l1,2,2)
        invcovCV(l1,2,2)  = 1./CL(l1,2,2)
     else
        !when using only temperature 
        invcov(l1,1,1) =  1./CLtil(l1,1,1)
        invcovCV(l1,1,1)  = 1./CL(l1,1,1)
     endif

  enddo

  !deallocating things I no longer need:
  deallocate(NL)
  !deallocate(CL)

  deallocate(detcovCV)
  deallocate(detcov)

  !compute ISW bispectrum:
  allocate(CLISW(2:lmax,5,5))
  open(unit=28,file=ISWfile, status="old")
  do l1 = 1, lmax
     if (l1 .eq. 1) then 
        read(28,*)
     else
        !11 = TT
        !21 = TP
        read(28,*) ell, l2ClTT, l2CLEE, l2CLBB, l2CLTE, l2CLPP, l2CLPT, l2CLPE
        !convert l(l+1)C_l/2pi [muK^2] to dimensionless C_l
        CLISW(l1,3,3) = 2.*pi*l2ClBB/ell/(ell+1.) /CMB2COBEnorm

        CLISW(l1,1,1) = CL(l1,1,1)
        CLISW(l1,2,2) = CL(l1,2,2)

        CLISW(l1,1,2) = CL(l1,1,2)
        CLISW(l1,2,1) = CL(l1,2,1)

        !see https://camb.info/readme.html
        CLISW(l1,4,1) = 2.*pi*l2ClPT/(real(ell,dl)*(real(ell,dl)+1.))**(3.d0/2.d0) /CMB2COBEnorm**(1./2)
        CLISW(l1,1,4) = 2.*pi*l2ClPT/(real(ell,dl)*(real(ell,dl)+1.))**(3.d0/2.d0) /CMB2COBEnorm**(1./2)
        CLISW(l1,4,4) = 2.*pi*l2CLPP/(ell*(ell+1.))**2
        !including reionziation-lensing 
        CLISW(l1,1,5) = 2.*pi*l2ClPE/(real(ell,dl)*(real(ell,dl)+1.))**(3.d0/2.d0) /CMB2COBEnorm**(1./2)
        CLISW(l1,5,1) = 2.*pi*l2ClPE/(real(ell,dl)*(real(ell,dl)+1.))**(3.d0/2.d0) /CMB2COBEnorm**(1./2)
        !write(*,*) ell, ClISW(l1,1:2,1)
     endif
  enddo
  close(28)


  !print lensing Bispectra in squeezed limit:
  !see arXiv: 1101.2234v2; looks good. 
  LensNorm = 1d6*CMB2COBEnorm**(-3./2)
  allocate(BLens(2,2,2))
  open(unit=FileUnit + 10,file = 'sigmafnl/LensingBispectrum_'//trim(dat)//'.txt', status='replace')
  write(FileUnit + 10,'(A23)') '#l^4b_4ll+4/10^6 muK^3'
  if(nfields .eq. 1) then
     write(FileUnit + 10,'(2A5)') '#ell','TTT'
  elseif(minfields .eq. 1 .and. nfields .eq. 2) then 
     !note that TTE and TET as well as ETE and EET are almost identical (except the l2 and l3 are interchanged)
     write(FileUnit + 10,'(9A5)') '#ell','TTT','TTE','TET','TEE','ETT','ETE','EET','EEE'
  else
     write(FileUnit + 10,'(2A5)') '#ell','EEE'
  endif
  do l2 =  max(lmin,4), 2000,2
     l1 = 4
     l3 = l2+4
     a3j2(:,:,1)=0.5d0
     call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
     if (nfields > 1) then 
        call GetThreeJs(a3j2(max(2,abs(l2-l1)),1,2),l1,l2,2,0)   
        call GetThreeJs(a3j2(max(2,abs(l2-l1)),2,2),l1,l2,0,2)   
        call GetThreeJs(a3j2(max(0,abs(l2-l1)),3,2),l1,l2,2,-2)

        a3j2(l3,:,2) = a3j2(l3,:,2)/a3j(l3)*0.5d0           

     endif
     do p  = minfields,nfields !T,E (8 terms only)
        do j = minfields,nfields !T,E
           do k = minfields,nfields !T,E
              BLens(k,j,p) = fc1(l1,l2,l3) *  &
                   (a3j2(l3,1,p)*ClISW(l2,j+3,1)*ClISW(l3,p,k) + &
                   a3j2(l3,3,p)*ClISW(l3,k+3,1)*ClISW(l2,p,j)) + &       
                   fc1(l2,l3,l1) * &
                   (a3j2(l3,3,j)*ClISW(l3,k+3,1)*ClISW(l1,j,p) + &
                   a3j2(l3,2,j)*ClISW(l1,p+3,1)*ClISW(l3,j,k)) + &       
                   fc1(l3,l1,l2) * &
                   (a3j2(l3,2,k)*ClISW(l1,p+3,1)*ClISW(l2,k,j) + &
                   a3j2(l3,1,k)*ClISW(l2,j+3,1)*ClISW(l1,k,p))
           enddo
        enddo
     enddo

     !storing squeezed bispectra
     write(FileUnit + 10,'(9E17.6)') real(l2,dl), 2.*real(l2,dl)**4*BLens(1:nfields,1:nfields,1:nfields)/LensNorm

  enddo
  close(FileUnit + 10)

  !lensing 4pt function, relevant for covariance:
  !using Eq. 76 in https://arxiv.org/pdf/astro-ph/0105117.pdf



  !now set Fisher lmax:
  if (Fishlmax .ne. lmax) then 
     lmax = Fishlmax
     if(lmax .lt. 1000) form = '(I3)'
     write(clmax, form)  lmax
  endif

  !allocating temp reduced bispectrum 
  allocate(BredTemp(2,2,2))

  !reduced bispectrum 
  BredTemp = 0.d0
  lmin  = 2

  write(*,*) 'FISHER lmin:', lmin
  !below we use:
  !* symmetry: l1=<l2=<l3
  !* even parity: l1 + l2 + l3 = even
  !* triangle condition: abs(l2-l1) =< l3 =< l1+l2

  if (.not. test_rintegral) then
     write(*,*) 'computing reduced',' ',nshape,'bispectrum ...'
     write(*,*) 'fisher error computed in loop ...'
     !filename for Fisher result:
     open(unit=FileUnit + 9,file = 'sigmafnl/Fisher_fnl_'//trim(nshape)//'_lmax_'//trim(clmax)//'_'//trim(dat)//'_'//trim(cexp)//'_'//trim(vers)//'.txt', status='replace')
     !putting sums to zero 
     TotSum = 0.d0
     TotSumCV = 0.d0
     TotSumISWLens = 0.d0
     TotSumCVISWLens = 0.d0
     TotSumLensCross = 0.d0
     TotSumCVLensCross = 0.d0
     FishC = 0.d0
     FishC_CV = 0.d0

     !parallel loop l1:
     !lmin = 4
     !call fwig_table_init(2*500,9)
     !call fwig_temp_init(2*500)
     !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
     !$OMP PRIVATE(l1,l2,l3, min_l,max_l,BredTemp, Blens), &
     !$OMP PRIVATE(p,j,k,q,r,s,i,Det,TempCov,TempCovCV,abint,Btemp,a3j,a3j2,val3j) &
     !$OMP PRIVATE(BLensISW,DetISWLens,DetLensCross) &
     !$OMP REDUCTION(+:TotSum,TotSumCV,FishC,FishC_CV,TotSumISWLens,TotSumCVISWLens,TotSumLensCross,TotSumCVLensCross) 

     do l1 = lmin, lmax

        do l2 =  max(lmin,l1), lmax
           min_l = max(abs(l1-l2),l2)
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if
           max_l = min(lmax,l1+l2)
           a3j2(:,:,1)=0.5d0
           call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
           if (nfields > 1) then 
              call GetThreeJs(a3j2(max(2,abs(l2-l1)),1,2),l1,l2,2,0)   
              call GetThreeJs(a3j2(max(2,abs(l2-l1)),2,2),l1,l2,0,2)   
              call GetThreeJs(a3j2(max(0,abs(l2-l1)),3,2),l1,l2,2,-2)


              do l3=min_l,max_l ,2
                 a3j2(l3,:,2) = a3j2(l3,:,2)/a3j(l3)*0.5d0           
              end do
           endif
           do l3=min_l,max_l, 2 !sum has to be even

              do p  = minfields,nfields !T,E (8 terms only)
                 do j = minfields,nfields !T,E
                    do k = minfields,nfields !T,E
                       Btemp = 0.d0

                       do i = irmax, irmin !rmini, rmaxi !
                          if (shape .eq. local) then 
                             abint = br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)
                          elseif (shape .eq. equil) then
                             abint = -3.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  -6.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k)) + &
                                  +3.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  +3.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  +3.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. ortho) then
                             abint = -9.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  -24.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k)) + &
                                  +9.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  +9.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  +9.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. folded) then
                             abint = 3.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  9.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k))  &
                                  -3.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  -3.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  -3.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. flat) then
                             abint = dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k) !1/(k1*k2*k3)**2
                          endif

                          !1D heat map of integral to determine r values to include
                          !cum_rint(i) = cum_rint(i) + rarray(i)**2*deltar(i)*abint

                          Btemp = Btemp + rarray(i)**2*deltar(i)*abint

                       enddo !r loop

                       !val3j = fwig3jj(2* l1 , 2* l2 , 2* l3 , &
                       !     2*(0), 2*(0) , 2*(0))

                       !means you have to recompute bispectrum elements every time, but requires low memory 
                       BredTemp(k,j,p)  = prefactor(l1,l2,l3)*a3j(l3)*Btemp
                       !write(*,*)  BredTemp(k,j,p)
                       

                       
                       !BredTemp(k,j,p)  = prefactor(l1,l2,l3)*val3j*Btemp
                       !if you write in this loop, it severaly slows down the inner loop.
                       !write(*,*) BredTemp(k,j,p)

                       !parity even bispectra from lensing-ISW and reioniation-ISW
                       !see e.g. https://arxiv.org/pdf/1302.5799.pdf Eq 2.6
                       !also Rio-phi for E
                       !p = field1
                       !j  = field2
                       !k = fierld3
                       BLens(k,j,p) = fc1(l1,l2,l3) *  &
                            (a3j2(l3,1,p)*ClISW(l2,j+3,1)*ClISW(l3,p,k) + &
                            a3j2(l3,3,p)*ClISW(l3,k+3,1)*ClISW(l2,p,j)) + &       
                            fc1(l2,l3,l1) * &
                            (a3j2(l3,3,j)*ClISW(l3,k+3,1)*ClISW(l1,j,p) + &
                            a3j2(l3,2,j)*ClISW(l1,p+3,1)*ClISW(l3,j,k)) + &       
                            fc1(l3,l1,l2) * &
                            (a3j2(l3,2,k)*ClISW(l1,p+3,1)*ClISW(l2,k,j) + &
                            a3j2(l3,1,k)*ClISW(l2,j+3,1)*ClISW(l1,k,p))



                       BLens(k,j,p)  = a3j(l3)*BLens(k,j,p)*prefactor(l1,l2,l3)
                       !if you want semi analytical only for local Temp and cross correlation with true shape:
                       if(want_SW_template .and. (nfields .eq. 1) .and. (shape .eq. 1)) then
                          BLens(k,j,p) = -prefactor(l1,l2,l3)*a3j(l3)*floc(l1,l2,l3)                     
                       endif


                    enddo !T, E loop
                 enddo  !T, E loop
              enddo  !T, E loop

              !now compute Fisher elements by looping over T/E
              do p  = minfields,nfields !T,E
                 do j = minfields,nfields !T,E
                    do k = minfields,nfields !T,E
                       do q  = minfields,nfields !T,E
                          do r = minfields,nfields !T,E
                             do s = minfields,nfields !T,E

                                !auto primordial 
                                Det = BredTemp(k,j,p)*BredTemp(s,r,q)
                                TempCov = invcov(l1,p,q)*invcov(l2,j,r)*invcov(l3,k,s)
                                TempCovCV = invcovCV(l1,p,q)*invcovCV(l2,j,r)*invcovCV(l3,k,s)

                                TotSum = TotSum + Det*TempCov/tr(l1,l2,l3)
                                TotSumCV = TotSumCV + Det*TempCovCV/tr(l1,l2,l3)

                                !auto lensing 
                                DetISWLens = BLens(k,j,p)*BLens(s,r,q) 

                                TotSumISWLens = TotSumISWLens + DetISWLens*TempCov/tr(l1,l2,l3)
                                TotSumCVISWLens = TotSumCVISWLens + DetISWLens*TempCovCV/tr(l1,l2,l3)

                                !ISW-lensing x primordial and ISW-reinization x primordial (any shape)
                                DetLensCross = BLens(k,j,p)*BredTemp(s,r,q)

                                TotSumLensCross = TotSumLensCross + DetLensCross*TempCov/tr(l1,l2,l3)
                                TotSumCVLensCross = TotSumCVLensCross + DetLensCross*TempCovCV/tr(l1,l2,l3)                             

                                !Seperating Fisher contributions; not completely correct. Need to use Wick theorem for off-diagnoals. To do.                               
                                FishC(s,r,q,k,j,p) = FishC(s,r,q,k,j,p) + Det/tr(l1,l2,l3)/Cltil(l1,p,q)/Cltil(l2,j,r)/Cltil(l3,k,s)  
                                FishC_CV(s,r,q,k,j,p) = FishC_CV(s,r,q,k,j,p) + Det/tr(l1,l2,l3)/Cl(l1,p,q)/Cl(l2,j,r)/Cl(l3,k,s)


                             enddo !T, E loop
                          enddo  !T, E loop
                       enddo  !T, E loop!

                    enddo !T, E loop
                 enddo  !T, E loop
              enddo  !T, E loop
!!$
           enddo !l3 loop
        enddo !l2 loop

     enddo !L1 loop
     !$OMP END PARAllEl DO
     !call fwig_temp_free();
     !call fwig_table_free();
     write(*,*)
     write(*,*) 'fisher error weighted by spectrum:'
     write(*,'(A4,X,A11,X,A11)') '<x>', 'experiment', 'cosmic var'
     if (minfields .eq. 1) then
        write(*,'(A4,X,F11.3,X,F11.3)') 'TTT', sqrt(1./FishC(1,1,1,1,1,1)/fsky), sqrt(1./FishC_CV(1,1,1,1,1,1))
     endif

     if (nfields .gt. 1 .and. minfields .eq. 1) then 
        write(*,'(A4,X,F11.3,X,F11.3)') 'TTE', sqrt(1./FishC(1,1,2,1,1,2)/fsky), sqrt(1./FishC_CV(1,1,2,1,1,2))
        write(*,'(A4,X,F11.3,X,F11.3)') 'TEE', sqrt(1./FishC(1,2,2,1,2,2)/fsky), sqrt(1./FishC_CV(1,2,2,1,2,2))
        write(*,'(A4,X,F11.3,X,F11.3)') 'EEE', sqrt(1./FishC(2,2,2,2,2,2)/fsky), sqrt(1./FishC_CV(2,2,2,2,2,2))
        write(*,'(A4,X,F11.3,X,F11.3)') 'ETT', sqrt(1./FishC(2,1,1,2,1,1)/fsky), sqrt(1./FishC_CV(2,1,1,2,1,1))
        write(*,'(A4,X,F11.3,X,F11.3)') 'EET', sqrt(1./FishC(2,2,1,2,2,1)/fsky), sqrt(1./FishC_CV(2,2,1,2,1,1))
     endif
     if (nfields .gt. 1 .and. minfields .eq. 2) then
        write(*,'(A4,X,F11.3,X,F11.3)') 'EEE', sqrt(1./FishC(2,2,2,2,2,2)/fsky), sqrt(1./FishC_CV(2,2,2,2,2,2))
     endif
 
     
     
     if (want_SW_template) then
        write(*,'(A12,X,I4,X,A31,X,F11.3,X,A9,X,F11.3,X,A1)') 'For lmax = ', lmax, 'the error on fnl (SW approx) = ', sqrt(1./TotSumISWLens/fsky), '(CV-limit', sqrt(1./TotSumCVISWLens/fskyCV), ')'
        write(*,'(A25,X,F11.3,X,A9,X,F11.3,X,A1)') 'correlation coefficient:', TotSumLensCross/TotSum**(1./2)/TotSumISWLens**(1./2), &
          '(CV-limit', TotSumCVLensCross/TotSumCV**(1./2)/TotSumCVISWLens**(1./2), ')'
     else
        
        !marginalization over lensing contributions
        DetFish = TotSumISWLens*TotSum -TotSumLensCross**2
        DetFishCV = TotSumCVISWLens*TotSumCV -TotSumCVLensCross**2

        DefnlMar = TotSumISWLens/DetFish
        DefnlMarCV = TotSumCVISWLens/DetFishCV

        if(comp_code) then 
        alpha = TotSumCVISWLens/DetFishCV !C/det
        beta = -TotSumCVLensCross/DetFishCV !A/det 

        write(*,*) 'lensing-ISW-fnl_local correlation coefficient:', TotSumCVLensCross/TotSumCV**(1./2)/TotSumCVISWLens**(1./2)
        write(*,*) 'Fisher local error:', 1/TotSumCV**(1./2)
        write(*,*) 'Fisher ISW-lensing error:', 1/TotSumCVISWLens**(1./2)
        write(*,*) 'alpha', alpha, 'beta', beta
        endif  
        
        write(*,'(A12,X,I4,X,A19,X,F11.3,X,A9,X,F11.3,X,A1)') 'For lmax = ', lmax, 'ISW-lensing bias = ', &
             TotSumLensCross/TotSum, '(CV-limit', TotSumCVLensCross/TotSumCV, ')'
        write(*,'(A12,X,I4,X,A31,X,F11.3,X,A9,X,F11.3,X,A1)') 'For lmax = ', lmax, 'the error on fnl-ISW-lensing = ', &
             sqrt(1./TotSumISWLens/fsky), '(CV-limit', sqrt(1./TotSumCVISWLens/fskyCV), ')'
        write(*,'(A25,X,F11.3,X,A9,X,F11.3,X,A1)') 'correlation coefficient:', TotSumLensCross/TotSum**(1./2)/TotSumISWLens**(1./2), &
          '(CV-limit', TotSumCVLensCross/TotSumCV**(1./2)/TotSumCVISWLens**(1./2), ')'
     endif
     
     write(*,'(A12,X,I4,X,A19,X,F11.3,X,A9,X,F11.3,X,A1)') 'For lmax = ', lmax, 'the error on fnl = ', sqrt(1./TotSum/fsky), '(CV-limit', sqrt(1./TotSumCV/fskyCV), ')'
     write(*,'(A28,X,F11.3,X,A9,X,F11.3,X,A1)') 'Lensing Marginalized Error:', sqrt(DefnlMar/fsky), '(CV-limit', sqrt(DefnlMarCV), ')'

     write(FileUnit + 9,'(A6,X,A21,F4.2,X,A25,X,A12,X,A15)') '#lmax', 'delta fnl/mar fsky = ', fsky, 'delta_fnl/mar CV fsky = 1', 'lensing bias', 'lensing bias CV'
     write(FileUnit + 9,'(I4,8E17.8)') lmax, sqrt(1./TotSum/fsky),sqrt(DefnlMar/fsky),sqrt(1./TotSumCV), sqrt(DefnlMarCV), TotSumISWLens/TotSum, TotSumCVISWLens/TotSumCV, &
          TotSumLensCross/TotSum, TotSumCVLensCross/TotSumCV
     !ix = ix + 1
     !enddo !lmax loop 
     close(FileUnit + 9)
  endif !fisher error done


  !for tesing r-integration/should need this only once per shape per cosomology per instrument 
  if (test_rintegral) then 
     write(*,*) 'getting S/N as function of r ...'
     !signal per rbin (to determine how many r bins to use for each shape)

     !if you want heat map:
     open(unit=FileUnit + 13,file = 'HeatMap/rint_'//trim(nshape)//'_lmax_'//trim(clmax)//'_'//trim(dat)//'_'//trim(vers)//'.txt', status='replace')

     allocate(sa3j(0:2*lmax+500,lmax**2/2))
     l1l2index = 1
     do l1 = lmin, lmax

        do l2 =  max(lmin,l1), lmax
           min_l = max(abs(l1-l2),l2)
           if (mod(l1+l2+min_l,2)/=0) then
              min_l = min_l+1 !l3 should only lead to parity even numbers
           end if
           max_l = min(lmax,l1+l2)
           call GetThreeJs(sa3j(abs(l2-l1),l1l2index),l1,l2,0,0)
           l1l2index = l1l2index + 1
           !write(*,*) l1l2index
        enddo
     enddo

     !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
     !$OMP PRIVATE(l1,l2,l3, min_l,max_l,BredTemp, Blens), &
     !$OMP PRIVATE(p,j,k,q,r,s,i,Det,TempCov,TempCovCV,abint,Btemp,a3j,a3j2) &
     !$OMP PRIVATE(TotSum,TotSumCV,l1l2index) 

     do i = 1, dimr-4
        TotSum = 0.d0
        TotSumCV = 0.d0
        TotSumISWLens = 0.d0
        TotSumCVISWLens = 0.d0
        TotSumLensCross = 0.d0
        TotSumCVLensCross = 0.d0
        l1l2index = 1
        do l1 = lmin, lmax

           do l2 =  max(lmin,l1), lmax
              min_l = max(abs(l1-l2),l2)
              if (mod(l1+l2+min_l,2)/=0) then
                 min_l = min_l+1 !l3 should only lead to parity even numbers
              end if
              max_l = min(lmax,l1+l2)

              do l3=min_l,max_l, 2 !sum has to be even

                 do p  = minfields,nfields !T,E (8 terms only)
                    do j = minfields,nfields !T,E
                       do k = minfields,nfields !T,E
                          Btemp = 0.d0

                          !do i = 1, dimr-4 !rmini, rmaxi !
                          if (shape .eq. local) then 
                             abint = br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)
                          elseif (shape .eq. equil) then
                             abint = -3.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  -6.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k)) + &
                                  +3.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  +3.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  +3.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. ortho) then
                             abint = -9.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  -24.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k)) + &
                                  +9.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  +9.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  +9.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. folded) then
                             abint = 3.*( br(i,l1,p)*br(i,l2,j)*ar(i,l3,k) + &
                                  br(i,l3,k)*br(i,l1,p)*ar(i,l2,j) + &
                                  br(i,l2,j)*br(i,l3,k)*ar(i,l1,p)) + &
                                  9.*(dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k))  &
                                  -3.*(gr(i,l1,p)*dr(i,l2,j) + gr(i,l2,j)*dr(i,l1,p))*br(i,l3,k) + &
                                  -3.*(dr(i,l1,p)*br(i,l2,j) + dr(i,l2,j)*br(i,l1,p))*gr(i,l3,k) + &
                                  -3.*(gr(i,l1,p)*br(i,l2,j) + gr(i,l2,j)*br(i,l1,p))*dr(i,l3,k) 
                          elseif (shape .eq. flat) then
                             abint = dr(i,l1,p)*dr(i,l2,j)*dr(i,l3,k) !1/(k1*k2*k3)**2
                          endif


                          Btemp = Btemp + rarray(i)**2*deltar(i)*abint

                          !means you have to recompute bispectrum elements every time, but requires low memory 
                          BredTemp(k,j,p)  = prefactor(l1,l2,l3)*sa3j(l3,l1l2index)*Btemp
                          !if you write in this loop, it severaly slows down the inner loop.
                          !write(*,*) BredTemp(k,j,p)



                       enddo !T, E loop
                    enddo  !T, E loop
                 enddo  !T, E loop

                 !now compute Fisher elements by looping over T/E
                 do p  = minfields,nfields !T,E
                    do j = minfields,nfields !T,E
                       do k = minfields,nfields !T,E
                          do q  = minfields,nfields !T,E
                             do r = minfields,nfields !T,E
                                do s = minfields,nfields !T,E

                                   !appears that taking sqrt leads to better results
                                   !not sure about full T,E bispectrum. But certinaly for seperate case. 
                                   Det = abs(BredTemp(k,j,p)*BredTemp(s,r,q))
                                   TempCov = abs(invcov(l1,p,q)*invcov(l2,j,r)*invcov(l3,k,s))
                                   TempCovCV = abs(invcovCV(l1,p,q)*invcovCV(l2,j,r)*invcovCV(l3,k,s))

                                   TotSum = TotSum + sqrt(Det*TempCov/tr(l1,l2,l3))
                                   TotSumCV = TotSumCV + sqrt(Det*TempCovCV/tr(l1,l2,l3))

                                enddo !T, E loop
                             enddo  !T, E loop
                          enddo  !T, E loop!

                       enddo !T, E loop
                    enddo  !T, E loop
                 enddo  !T, E loop

              enddo !l3 loop
              l1l2index = l1l2index + 1
           enddo !l2 loop
        enddo !L1 loop
        write(*,'(4E18.7)') rarray(i), deltar(i), sqrt(1./TotSum/fsky), sqrt(1./TotSumCV/fskyCV)
        write(FileUnit + 13,'(4E18.7)') rarray(i), deltar(i), sqrt(1./TotSum/fsky), sqrt(1./TotSumCV/fskyCV)
     enddo !r loop
     !$OMP END PARAllEl DO

     close(FileUnit + 13)
     deallocate(sa3j)
  endif




  !deallocating remaining arrays:

  deallocate(CL)
  deallocate(CLISW)
  deallocate(ar)
  deallocate(br)
  deallocate(gr)
  deallocate(dr)
  deallocate(rarray,deltar)
  deallocate(CLtil)
  deallocate(BredTemp, BLens)
  deallocate(invCov)
  deallocate(invCovCV)


  !end of prgram

  !subroutines and functions needed:
contains


  real(dl) function floc(l1,l2,l3)
    !SW approximation 
    integer :: l1, l2, l3
    real(dl) :: amp
    real(dl) :: As = 2.1056d-9
    !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
    amp = (2.d0/27./pi**2)*As
    floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
         1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
    floc = floc*amp*2.E-7 !2.E-7  is introduced to get roughly same amplitude at l_max = 500 to full fnl_local
  end function floc
  
  real function tr(l1,l2,l3)
    integer, intent(in) :: l1,l2,l3
    if ((l1.eq.l2).and.(l2.eq.l3)) then
       tr  = 6.d0
    elseif ((l1.eq.l2).or.(l2.eq.l3).or.(l3.eq.l1)) then 
       tr = 2d0
    else
       tr = 1.d0
    endif

  end function tr

  real function prefactor(l1,l2,l3)
    integer, intent(in) :: l1,l2,l3

    prefactor = 2.0*sqrt((1./4.)*((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.))/pi)
  end function prefactor

  real function fc1(l1,l2,l3)
    integer, intent(in) :: l1,l2,l3

    fc1 = (-1.d0*l1*(l1+1.d0)+1.d0*l2*(l2+1.d0)+1.d0*l3*(l3+1.d0))/2.d0

  end function fc1

  real function parity(n)
    integer, intent(in) :: n

    if (MOD(n,2)==0) then
       parity = 1.0
    else
       parity = 0.0
    endif
  end function parity

  real function fac(l)
    integer, intent(in) :: l
    fac = 1.d0*l*(l+1.)/(2.*pi)
  end function fac

  real function Noise_l(sigmab,w0,ell)
    real(dl) :: sigmab
    real(dl) :: w0
    integer :: ell
    real(dl) :: arcmin
    real(dL) :: fac
    arcmin = 0.000290888
    fac = arcmin / 2.d0*sqrt(2.d0*log(2.d0))

    Noise_l = w0**2*arcmin**2*exp(ell**2*fac**2*sigmab**2)

  end function Noise_l

  subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input
    !parameters.
    implicit none
    integer, parameter :: dl = KIND(1.d0)
    integer, intent(in) :: l2in,l3in, m2in,m3in
    real(dl), dimension(*) :: thrcof
    INTEGER, PARAMETER :: i8 = selected_int_kind(18)
    integer(i8) :: l2,l3,m2,m3
    integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2

    real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1,sumuni
    real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
    integer i,ier, index, nlim, sign2
    integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
    real(dl), parameter :: zero = 0._dl, one = 1._dl
    real(dl), parameter ::  tiny = 1.0d-30, srtiny=1.0d-15, huge = 1.d30,srhuge = 1.d15

    ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

    ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
    !                to l1max = l2+l3
    ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from l1min and l1max
    ! respectively.
    !
    ! lmatch is the l1-value at which forward and backward recursion are
    ! matched.
    !
    ! ndim is the length of the array thrcof
    !
    ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
    ! ier = -2 if possible 3j's exceed ndim
    ! ier >= 0 otherwise

    l2=l2in
    l3=l3in
    m2=m2in
    m3=m3in
    newfac = 0
    lmatch = 0
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    ier = 0

    if (l2 < abs(m2) .or. l3 < m3) then
       ier = -1
       ! call MpiStop('error ier = -1')
       print*, 'error ier = -1'
       stop
       return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3),abs(m1))
    l1max = l2+l3

    if (l1min >= l1max) then
       if (l1min/=l1max) then
          ier = -1

          !call MpiStop('error ier = -1')
          print*, 'error ier = -1' 
          stop
          return
       end if

       ! reached if l1 can take only one value, i.e.l1min=l1max
       thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
       return

    end if

    nfin = l1max-l1min+1

    ! starting forward recursion from l1min taking nstep1 steps
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1 + 1)*tiny

    lstep = 1

30  lstep = lstep+1
    l1 = l1+1

    oldfac = newfac
    a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
    a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
    newfac = sqrt(a2*real(a1,dl))
    if (l1 == 1) then
       !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
       c1 = -(2*l1-1)*l1*(m3-m2)/newfac
    else

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
       denom = (l1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*l1-1)*dv/denom

    end if

    if (lstep<= 2) then

       ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       thrcof(2) = x
       sum1 = sum1+tiny*(2*l1+1)*c1*c1
       if(lstep==nfin) then
          sumuni=sum1
          go to 230
       end if
       goto 30

    end if

    c2 = -l1*oldfac/denom

    ! recursion to the next 3j-coeff x  
    x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
    thrcof(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (2*l1+1)*x*x
    if (lstep/=nfin) then

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(x) >= srhuge) then

          ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
          ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
          ! HAS TO BE RESCALED TO PREVENT OVERFLOW

          ier = ier+1
          do i = 1, lstep
             if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
             thrcof(i) = thrcof(i)/srhuge
          end do

          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge

       end if

       ! as long as abs(c1) is decreasing, the recursion proceeds towards
       ! increasing
       ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
       ! detected, the recursion direction is reversed.

       if (c1old > abs(c1)) goto 30

    end if !lstep/=nfin

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion

    lmatch = l1-1
    x1 = x
    x2 = thrcof(lstep-1)
    x3 = thrcof(lstep-2)
    nstep2 = nfin-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from l1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! l1 = lmatch-1, lmatch, lmatch+1

    nfinp1 = nfin+1
    nfinp2 = nfin+2
    nfinp3 = nfin+3
    l1 = l1max
    thrcof(nfin) = srtiny
    sum2 = tiny*(2*l1+1)

    l1 = l1+2
    lstep=1

    do
       lstep = lstep + 1
       l1= l1-1

       oldfac = newfac
       a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
       a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
       newfac = sqrt(a1*real(a2,dl))

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

       denom = l1*newfac
       c1 = -(2*l1-1)*dv/denom
       if (lstep <= 2) then

          ! if l2=l2max+1, the third term in the recursion vanishes

          y = srtiny*c1
          thrcof(nfin-1) = y
          sumbac = sum2
          sum2 = sum2 + tiny*(2*l1-3)*c1*c1

          cycle

       end if

       c2 = -(l1-1)*oldfac/denom

       ! recursion to the next 3j-coeff y
       y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

       if (lstep==nstep2) exit

       thrcof(nfinp1-lstep) = y
       sumbac = sum2
       sum2 = sum2+(2*l1-3)*y*y

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(y) >= srhuge) then

          ! reached if 3j-coeff larger than srhuge so that the recursion series
          ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent
          ! overflow

          ier=ier+1
          do i = 1, lstep
             index=nfin-i+1
             if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
             thrcof(index) = thrcof(index)/srhuge
          end do

          sum2=sum2/huge
          sumbac=sumbac/huge

       end if

    end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

    y3 = y
    y2 = thrcof(nfinp2-lstep)
    y1 = thrcof(nfinp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal
    ! error

    ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
    nlim = nfin-nstep2+1

    if (abs(ratio) >= 1) then

       thrcof(1:nlim) = ratio*thrcof(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

    else

       nlim = nlim+1
       ratio = 1/ratio
       do n = nlim, nfin
          thrcof(n) = ratio*thrcof(n)
       end do
       sumuni = sumfor + ratio*ratio*sumbac

    end if
    ! normalise 3j-coeffs

230 cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

    sign1 = sign(one,thrcof(nfin))
    sign2 = (-1)**(abs(l2+m2-l3+m3))
    if (sign1*sign2 <= 0) then
       cnorm = -cnorm
    end if
    if (abs(cnorm) >= one) then
       thrcof(1:nfin) = cnorm*thrcof(1:nfin)
       return
    end if

    thresh = tiny/abs(cnorm)

    do n = 1, nfin
       if (abs(thrcof(n)) < thresh) thrcof(n) = zero
       thrcof(n) = cnorm*thrcof(n)
    end do
    return 

  end subroutine GetThreeJs

end program BS_SN
