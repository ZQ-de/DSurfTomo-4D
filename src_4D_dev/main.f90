! CODE FOR SURFACE WAVE TOMOGRAPHY USING DISPERSION MEASUREMENTS
! VERSION:
!      1.0 
! AUTHOR:
!      HONGJIAN FANG. fanghj@mail.ustc.edu.cn
! PURPOSE:
!      DIRECTLY INVERT SURFACE WAVE DISPERSION MEASUREMENTS FOR 3-D
! STUCTURE WITHOUT THE INTERMEDIATE STEP OF CONSTUCTION THE PHASE
! OR GROUP VELOCITY MAPS.
! REFERENCE:
! Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D.
! (2015). Direct inversion of surface wave dispersion for
! three-dimensional shallow crustal structure based on ray tracing:
! methodology and application. Geophysical Journal International,
! 201(3), 1251-1263.
! HISTORY:
!       2015/01/31 START TO REORGONIZE THE MESSY CODE   
!

program SurfTomo
        use lsmrModule, only:lsmr
        use  lsmrblasInterface, only : dnrm2
        use omp_lib
        implicit none

        ! VARIABLE DEFINE

        character inputfile*80
        character logfile*100
        character outmodel*100
        character outsyn*100
        logical ex
        character dummy*40
        character datafile*80
        ! ZHAOQIAN
        character datafile_cur*80
        character true_model_file_t1*80
        character true_model_file_t2*80

        integer nx,ny,nz
        real goxd,gozd
        real dvxd,dvzd
        integer nsrc,nrc
        real weight,weight0
        real damp
        real minthk
        integer kmax,kmaxRc,kmaxRg,kmaxLc,kmaxLg
        real*8,dimension(:),allocatable:: tRc,tRg,tLc,tLg
        real,dimension(:),allocatable:: depz
        ! ZHAOQIAN
        real,dimension(:),allocatable:: depz_true_t1,depz_true_t2
        integer itn
        integer nout
        integer localSize
        real mean,std_devs,balances,balanceb
        integer msurf
        ! ZHAOQIAN
        real,dimension(:),allocatable:: obst,obst_cur,obstdiff,dsyn,dsyn_ref,&
                cbst,wt,dtres,dist,datweight
        real,dimension(:),allocatable:: ref_src_lat_deg,ref_src_lon_deg,&
                ref_rec_lat_deg,ref_rec_lon_deg
        integer,dimension(:),allocatable:: ref_period_id,ref_wavetype_id,&
                ref_veltype_id
        
        real,dimension(:),allocatable:: pvall,depRp,pvRp
        real sta1_lat,sta1_lon,sta2_lat,sta2_lon
        real sta1_lat_deg,sta1_lon_deg,sta2_lat_deg,sta2_lon_deg
        real dist1
        integer dall
        integer istep
        real,parameter :: pi=3.1415926535898
        integer checkstat
        integer ii,jj,kk
        real, dimension (:,:), allocatable :: scxf,sczf
        real, dimension (:,:,:), allocatable :: rcxf,rczf
        integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
        integer,dimension(:),allocatable::nsrc1
        integer,dimension(:,:),allocatable::periods
        real,dimension(:),allocatable::rw
        integer,dimension(:),allocatable::iw,col,nrow
        real,dimension(:),allocatable::dv,norm,dvsub,dvstd,dvall
!        real,dimension(:),allocatable::dvall
        !ZHAOQIAN
        real,dimension(:,:,:),allocatable::vsf,vsf_init
        real,dimension(:,:,:),allocatable::vsftrue_t1,vsftrue_t2,&
                dvs_true
        character strf
        integer veltp,wavetp
        real velvalue
        integer knum,knumo,err
        integer istep1,istep2
        integer period
        ! ZHAOQIAN
        integer knumi,srcnum,count1,count_cur
        integer HorizonType,VerticalType
        character line*200
        integer iter,maxiter
        integer maxnar
        real acond
        real anorm
        real arnorm
        real rnorm
        real xnorm
        character str1
        real atol,btol
        real conlim
        integer istop
        integer itnlim
        integer lenrw,leniw
        integer nar,nar_tmp,nars
        integer count3,nvz,nvx
        integer m,maxvp,n
        integer i,j,k
        real Minvel,MaxVel
        real spfra
        real noiselevel
        integer ifsyn
        real averdws
        real maxnorm
        real threshold0,q25,q75


        ! OPEN FILES FIRST TO OUTPUT THE PROCESS
        !nout=36
        !open(nout,file='lsmr.txt')

        ! OUTPUT PROGRAM INFOMATION            
        write(*,*)
        write(*,*) '                             DSurfTomo-4D (v1.0, based on DSurfTomo v1.4)'
        !write(*,*) 'PLEASE contact Hongjain Fang &
        !        (fanghj@mail.ustc.edu.cn) if you find any bug'
        write(*,*) 'For bug report, PLEASE contact Hongjain Fang &
                (fanghj1990@gmail.com)'
        write(*,*)

        ! READ INPUT FILE
        if (iargc() < 1) then
                write(*,*) 'input file [DSurfTomo.in(default)]:'
                read(*,'(a)') inputfile
                if (len_trim(inputfile) <=1 ) then
                        inputfile = 'DSurfTomo.in'
                else
                        inputfile = inputfile(1:len_trim(inputfile))
                endif
        else
                call getarg(1,inputfile)
        endif
        inquire(file = inputfile, exist = ex)
        if (.not. ex) stop 'unable to open the inputfile'

        open(10,file=inputfile,status='old')
        read(10,'(a30)')dummy
        read(10,'(a30)')dummy
        read(10,'(a30)')dummy
        read(10,*) datafile
        read(10,*) datafile_cur
        read(10,*) nx,ny,nz
        read(10,*) goxd,gozd
        read(10,*) dvxd,dvzd
        read(10,*) nsrc
        read(10,*) weight0,damp
        read(10,*) minthk
        read(10,*) Minvel,Maxvel
        read(10,*) maxiter
        read(10,*) spfra
        read(10,*) kmaxRc
        write(*,*)  'model origin:latitude,longitue'
        write(*,'(2f10.5)') goxd,gozd
        write(*,*) 'grid spacing:latitude,longitue'
        write(*,'(2f10.5)') dvxd,dvzd
        write(*,*) 'model dimension:nx,ny,nz'
        write(*,'(3i5)') nx,ny,nz
        write(logfile,'(a,a)')trim(inputfile),'.log' 
        open(66,file=logfile)
        write(66,*)
        write(66,*)'                         S U R F  T O M O'
        write(66,*)'PLEASE contact Hongjain Fang &
                (fanghj@mail.ustc.edu.cn) if you find any bug'
        write(66,*)
        write(66,*) 'model origin:latitude,longitue'
        write(66,'(2f10.5)') goxd,gozd
        write(66,*) 'grid spacing:latitude,longitue'
        write(66,'(2f10.5)') dvxd,dvzd
        write(66,*) 'model dimension:nx,ny,nz'
        write(66,'(3i5)') nx,ny,nz
        if(kmaxRc.gt.0)then
                allocate(tRc(kmaxRc),&
                        stat=checkstat)
                if (checkstat > 0) stop 'error allocating RP'
                read(10,*)(tRc(i),i=1,kmaxRc)
                write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
                write(*,'(50f7.2)')(tRc(i),i=1,kmaxRc)
                write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
                write(66,'(50f7.2)')(tRc(i),i=1,kmaxRc)
        endif
        read(10,*)kmaxRg
        if(kmaxRg.gt.0)then
                allocate(tRg(max(kmaxRg,1)), stat=checkstat)
                if (checkstat > 0) stop 'error allocating RP'
                read(10,*)(tRg(i),i=1,kmaxRg)
                write(*,*)'Rayleigh wave group velocity used,periods:(s)'
                write(*,'(50f7.2)')(tRg(i),i=1,kmaxRg)
                write(66,*)'Rayleigh wave group velocity used,periods:(s)'
                write(66,'(50f7.2)')(tRg(i),i=1,kmaxRg)
        endif
        read(10,*)kmaxLc
        if(kmaxLc.gt.0)then
                allocate(tLc(max(kmaxLc,1)), stat=checkstat)
                if (checkstat > 0) stop 'error allocating RP'
                read(10,*)(tLc(i),i=1,kmaxLc)
                write(*,*)'Love wave phase velocity used,periods:(s)'
                write(*,'(50f7.2)')(tLc(i),i=1,kmaxLc)
                write(66,*)'Love wave phase velocity used,periods:(s)'
                write(66,'(50f7.2)')(tLc(i),i=1,kmaxLc)
        endif
        read(10,*)kmaxLg
        if(kmaxLg.gt.0)then
                allocate(tLg(max(kmaxLg,1)), stat=checkstat)
                if (checkstat > 0) stop 'error allocating RP'
                read(10,*)(tLg(i),i=1,kmaxLg)
                write(*,*)'Love wave group velocity used,periods:(s)'
                write(*,'(50f7.2)')(tLg(i),i=1,kmaxLg)
                write(66,*)'Love wave group velocity used,periods:(s)'
                write(66,'(50f7.2)')(tLg(i),i=1,kmaxLg)
        endif
        read(10,*)ifsyn
        ! ZHAOQIAN
        if(ifsyn.eq.1) then
                read(10,*) true_model_file_t1
                read(10,*) true_model_file_t2
        endif
        read(10,*)noiselevel
        read(10,*) threshold0

        close(10)
        nrc=nsrc
        kmax=kmaxRc+kmaxRg+kmaxLc+kmaxLg+100

        ! READ MEASUREMENTS            
        open(unit=87,file=datafile,status='old')
        allocate(scxf(nsrc,kmax),sczf(nsrc,kmax),&
                rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate'
        endif
        allocate(periods(nsrc,kmax),wavetype(nsrc,kmax),&
                nrc1(nsrc,kmax),nsrc1(kmax),&
                igrt(nsrc,kmax),stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate'
        endif
        ! ZHAOQIAN
        allocate(obst(nrc*nsrc*kmax),obst_cur(nrc*nsrc*kmax),&
                obstdiff(nrc*nsrc*kmax),dist(nrc*nsrc*kmax),&
                stat=checkstat)
        allocate(ref_src_lat_deg(nrc*nsrc*kmax),&
                ref_src_lon_deg(nrc*nsrc*kmax),&
                ref_rec_lat_deg(nrc*nsrc*kmax),&
                ref_rec_lon_deg(nrc*nsrc*kmax),&
                ref_period_id(nrc*nsrc*kmax),&
                ref_wavetype_id(nrc*nsrc*kmax),&
                ref_veltype_id(nrc*nsrc*kmax),&
                stat=checkstat)
        if(checkstat > 0) stop 'error allocating epoch-pair key arrays'
        allocate(pvall(nrc*nsrc*kmax),depRp(nrc*nsrc*kmax),&
                pvRp(nrc*nsrc*kmax),stat=checkstat)
        IF(checkstat > 0)THEN
                write(6,*)'error with allocate'
        ENDIF
        istep=0
        istep2=0
        dall=0
        knumo=12345
        knum=0
        istep1=0
        do 
        read(87,'(a)',iostat=err) line
        if(err.eq.0) then
                if(line(1:1).eq.'#') then
                        read(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
                        sta1_lat_deg=sta1_lat
                        sta1_lon_deg=sta1_lon
                        if(wavetp.eq.2.and.veltp.eq.0) knum=period
                        if(wavetp.eq.2.and.veltp.eq.1) knum=period
                        if(wavetp.eq.1.and.veltp.eq.0) knum=kmaxRg+period
                        if(wavetp.eq.1.and.veltp.eq.1) knum=kmaxLc+kmaxRg+&
                                period
                        if(knum.ne.knumo) then
                                istep=0
                                istep2=istep2+1
                        endif
                        istep=istep+1
                        istep1=0
                        sta1_lat=(90.0-sta1_lat)*pi/180.0
                        sta1_lon=sta1_lon*pi/180.0
                        scxf(istep,knum)=sta1_lat
                        sczf(istep,knum)=sta1_lon
                        periods(istep,knum)=period
                        wavetype(istep,knum)=wavetp
                        igrt(istep,knum)=veltp
                        nsrc1(knum)=istep
                        knumo=knum
                else
                        read(line,*) sta2_lat,sta2_lon,velvalue
                        sta2_lat_deg=sta2_lat
                        sta2_lon_deg=sta2_lon
                        istep1=istep1+1
                        dall=dall+1

                        ref_src_lat_deg(dall)=sta1_lat_deg
                        ref_src_lon_deg(dall)=sta1_lon_deg
                        ref_rec_lat_deg(dall)=sta2_lat_deg
                        ref_rec_lon_deg(dall)=sta2_lon_deg
                        ref_period_id(dall)=period
                        ref_wavetype_id(dall)=wavetp
                        ref_veltype_id(dall)=veltp

                        sta2_lat=(90.0-sta2_lat)*pi/180.0
                        sta2_lon=sta2_lon*pi/180.0
                        rcxf(istep1,istep,knum)=sta2_lat
                        rczf(istep1,istep,knum)=sta2_lon
                        call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
                        dist(dall)=dist1
                        obst(dall)=dist1/velvalue
                        pvall(dall)=velvalue
                        nrc1(istep,knum)=istep1
                endif
        else
                exit
        endif
        enddo
        close(87)
        
        ! ZHAOQIAN
                ! READ CURRENT-PERIOD DATA FOR 4D PAIRWISE INVERSION
        ! Reference-period data have already been stored in obst(:)
        ! Current-period data file is specified in the input parameter file

        inquire(file=datafile_cur,exist=ex)
        if(.not.ex) then
                write(*,*) 'Current-period data file does not exist: ',trim(datafile_cur)
                stop
        endif

        open(89,file=datafile_cur,status='old')

        count_cur = 0

        do
        read(89,'(a)',iostat=err) line

        if(err.eq.0) then
                if(line(1:1).eq.'#') then
                        read(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
                        sta1_lat_deg=sta1_lat
                        sta1_lon_deg=sta1_lon
                        cycle
                else
                        read(line,*) sta2_lat,sta2_lon,velvalue
                        count_cur = count_cur + 1

                        if(count_cur.gt.dall) then
                                write(*,*) 'Error: current data has more measurements than reference data'
                                write(*,*) 'count_cur = ',count_cur,' dall = ',dall
                                stop
                        endif

                        if(abs(sta1_lat_deg-ref_src_lat_deg(count_cur)).gt.1.0e-5 .or. &
                           abs(sta1_lon_deg-ref_src_lon_deg(count_cur)).gt.1.0e-5 .or. &
                           abs(sta2_lat-ref_rec_lat_deg(count_cur)).gt.1.0e-5 .or. &
                           abs(sta2_lon-ref_rec_lon_deg(count_cur)).gt.1.0e-5 .or. &
                           period.ne.ref_period_id(count_cur) .or. &
                           wavetp.ne.ref_wavetype_id(count_cur) .or. &
                           veltp.ne.ref_veltype_id(count_cur)) then

                                write(*,*) 'Error: epoch observation key mismatch'
                                write(*,*) 'measurement index = ',count_cur
                                write(*,*) 'epoch-1 source = ',&
                                        ref_src_lat_deg(count_cur),&
                                        ref_src_lon_deg(count_cur)
                                write(*,*) 'epoch-2 source = ',&
                                        sta1_lat_deg,sta1_lon_deg
                                write(*,*) 'epoch-1 receiver = ',&
                                        ref_rec_lat_deg(count_cur),&
                                        ref_rec_lon_deg(count_cur)
                                write(*,*) 'epoch-2 receiver = ',&
                                        sta2_lat,sta2_lon
                                write(*,*) 'epoch-1 period/wave/velocity type = ',&
                                        ref_period_id(count_cur),&
                                        ref_wavetype_id(count_cur),&
                                        ref_veltype_id(count_cur)
                                write(*,*) 'epoch-2 period/wave/velocity type = ',&
                                        period,wavetp,veltp
                                error stop 1
                        endif

                        if(velvalue.le.0.0) then
                                write(*,*) 'Error: non-positive velocity in current data'
                                write(*,*) 'line index = ',count_cur,' velocity = ',velvalue
                                stop
                        endif

                        obst_cur(count_cur) = dist(count_cur) / velvalue
                endif
        else
                exit
        endif
        enddo

        close(89)

        if(count_cur.ne.dall) then
                write(*,*) 'Error: current and reference data have different number of measurements'
                write(*,*) 'current count = ',count_cur,' reference count = ',dall
                stop
        endif

        do i = 1,dall
        obstdiff(i) = obst_cur(i) - obst(i)
        enddo

        mean = sum(obstdiff(1:dall))/dall
        std_devs = sqrt(sum(obstdiff(1:dall)**2)/dall - mean**2)

        ! ZHAOQIAN
        if(ifsyn.eq.0) then
                write(*,*) '4D pairwise data loaded successfully'
                write(*,*) 'Reference data file: ',trim(datafile)
                write(*,*) 'Current data file: ',trim(datafile_cur)
                write(*,'(a,f10.4,a,f10.4,a,f10.4,a)') &
                        'DeltaTobs mean/std/rms: ',&
                        mean*1000.0,' ms ',std_devs*1000.0,' ms ',&
                        sqrt(sum(obstdiff(1:dall)**2)/dall)*1000.0,' ms'
        endif



        allocate(depz(nz), stat=checkstat)
        ! ZHAOQIAN
        allocate(depz_true_t1(nz),depz_true_t2(nz),stat=checkstat)
        if(checkstat.gt.0) stop 'error allocating true-model depth arrays'
        maxnar = spfra*dall*nx*ny*nz!sparsity fraction
        if (maxnar<0) print*, 'number overflow, decrease your sparsefrac'
        maxvp = (nx-2)*(ny-2)*(nz-1)
        allocate(dv(maxvp),dvsub(maxvp),dvstd(maxvp),dvall(maxvp), stat=checkstat)
!        allocate(dvall(maxvp*nrealizations),stats=checkstat)
        allocate(norm(maxvp), stat=checkstat)
        ! ZHAOQIAN
        allocate(vsf(nx,ny,nz),vsf_init(nx,ny,nz), stat=checkstat)
        allocate(vsftrue_t1(nx,ny,nz),vsftrue_t2(nx,ny,nz),&
                dvs_true(nx,ny,nz),stat=checkstat)

        allocate(rw(maxnar), stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate: real rw'
        endif
        allocate(iw(2*maxnar+1), stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate: integer iw'
        endif
        allocate(col(maxnar),nrow(dall), stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate:  integer iw'
        endif
        ! ZHAOQIAN
        allocate(cbst(dall+maxvp),dsyn(dall),dsyn_ref(dall),&
                datweight(dall),wt(dall+maxvp),dtres(dall+maxvp),&
                stat=checkstat)

        ! MEASUREMENTS STATISTICS AND READ INITIAL MODEL           
        write(*,'(a,i7)') ' Number of all measurements',dall

        open(10,file='MOD',status='old')
        read(10,*) (depz(i),i=1,nz)
        do k = 1,nz
        do j = 1,ny
        read(10,*)(vsf(i,j,k),i=1,nx)
        enddo
        enddo
        close(10)
        
        ! ZHAOQIAN,Save reference model for 4D pairwise inversion
        vsf_init = vsf
        
        write(*,*) 'grid points in depth direction:(km)'
        write(*,'(50f7.2)') depz

        ! ZHAOQIAN, TWO-EPOCH 4D SYNTHETIC TEST
        if (ifsyn == 1) then
                write(*,*) 'Two-epoch synthetic test begins'

                ! Check whether both true-model files exist
                inquire(file=true_model_file_t1,exist=ex)
                if(.not.ex) then
                        write(*,*) 'Error: epoch-1 true model does not exist: ',&
                                trim(true_model_file_t1)
                        error stop 1
                endif

                inquire(file=true_model_file_t2,exist=ex)
                if(.not.ex) then
                        write(*,*) 'Error: epoch-2 true model does not exist: ',&
                                trim(true_model_file_t2)
                        error stop 1
                endif

                ! Read epoch-1 true model, including the depth line
                open(11,file=true_model_file_t1,status='old')
                read(11,*) (depz_true_t1(i),i=1,nz)
                do k = 1,nz
                do j = 1,ny
                        read(11,*) (vsftrue_t1(i,j,k),i=1,nx)
                enddo
                enddo
                close(11)

                ! Read epoch-2 true model, including the depth line
                open(12,file=true_model_file_t2,status='old')
                read(12,*) (depz_true_t2(i),i=1,nz)
                do k = 1,nz
                do j = 1,ny
                        read(12,*) (vsftrue_t2(i,j,k),i=1,nx)
                enddo
                enddo
                close(12)

                ! Verify that both true models use the same depth grid as MOD
                if(any(abs(depz_true_t1-depz).gt.1.0e-5)) then
                        write(*,*) 'Error: epoch-1 true-model depths do not match MOD'
                        error stop 1
                endif

                if(any(abs(depz_true_t2-depz).gt.1.0e-5)) then
                        write(*,*) 'Error: epoch-2 true-model depths do not match MOD'
                        error stop 1
                endif

                if(minval(vsftrue_t1).le.0.0 .or. &
                   minval(vsftrue_t2).le.0.0) then
                        write(*,*) 'Error: true models contain non-positive Vs'
                        error stop 1
                endif

                ! True time-lapse shear-wave velocity change
                dvs_true = vsftrue_t2 - vsftrue_t1

                ! ZHAOQIAN, check whether T1 true model equals MOD
                write(*,'(a,f12.8)') &
                        'T1 true/MOD maximum Vs difference: ',&
                        maxval(abs(vsftrue_t1-vsf_init))

                ! Generate epoch-1 synthetic travel times
                write(*,*) 'Generating epoch-1 synthetic travel times...'
                call synthetic(nx,ny,nz,maxvp,vsftrue_t1,obst,&
                        goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                        tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                        scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
                        nsrc,nrc,noiselevel)

                ! Generate epoch-2 synthetic travel times independently
                write(*,*) 'Generating epoch-2 synthetic travel times...'
                call synthetic(nx,ny,nz,maxvp,vsftrue_t2,obst_cur,&
                        goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                        tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                        scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
                        nsrc,nrc,noiselevel)

                ! Recompute the differential observations after both forward runs
                do i = 1,dall
                        obstdiff(i) = obst_cur(i) - obst(i)
                enddo

                mean = sum(obstdiff(1:dall))/dall
                std_devs = sqrt(sum(obstdiff(1:dall)**2)/dall - mean**2)

                write(*,*) 'Two-epoch synthetic data generated successfully'
                write(*,*) 'Epoch-1 true model: ',trim(true_model_file_t1)
                write(*,*) 'Epoch-2 true model: ',trim(true_model_file_t2)
                write(*,'(a,2f10.5)') 'True delta Vs min/max: ',&
                        minval(dvs_true),maxval(dvs_true)
                write(*,'(a,f10.3,a)') 'Synthetic differential travel-time RMS: ',&
                        sqrt(sum(obstdiff(1:dall)**2)/dall)*1000.0,' ms'
        endif


        ! ============================================================
        ! ZHAOQIAN 4D PAIRWISE:
        ! Compute reference theoretical travel times using MOD_T1.
        ! At this point, vsf = vsf_init = MOD_T1.
        ! The result dsyn_ref is Tcalc_ref.
        ! ============================================================
        iw = 0
        rw = 0.0
        col = 0

        write(*,*) 'computing reference theoretical travel times...'
        call CalSurfG(nx,ny,nz,maxvp,vsf_init,iw,rw,col,dsyn,&
                goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
                nsrc,nrc,nar)

        dsyn_ref = dsyn

        ! ZHAOQIAN, validate the epoch-1 synthetic forward calculation
        if(ifsyn.eq.1) then
                write(*,'(a,f12.6,a)') &
                        'T1 synthetic/reference forward RMS mismatch: ',&
                        sqrt(sum((obst(1:dall)-dsyn_ref(1:dall))**2)/dall)*&
                        1000.0,' ms'
        endif

        ! ZHAOQIAN, check the differential forward-operator mismatch
        if(ifsyn.eq.1) then
                iw = 0
                rw = 0.0
                col = 0
                nar = 0

                call CalSurfG(nx,ny,nz,maxvp,vsftrue_t2,iw,rw,col,dsyn,&
                        goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                        tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                        scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
                        nsrc,nrc,nar)

                write(*,'(a,f12.6,a)') &
                        'Differential forward RMS mismatch: ',&
                        sqrt(sum(((obst_cur(1:dall)-obst(1:dall))-&
                        (dsyn(1:dall)-dsyn_ref(1:dall)))**2)/dall)*&
                        1000.0,' ms'
        endif

        ! Clear work arrays before normal inversion iterations
        iw = 0
        rw = 0.0
        col = 0
        nar = 0


        ! ITERATE UNTILL CONVERGE
        do iter = 1,maxiter
        iw = 0
        rw = 0.0
        col = 0

        ! COMPUTE SENSITIVITY MATRIX
        write(*,*) 'computing sensitivity matrix...'
        call CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
                goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
                nsrc,nrc,nar)

        ! ZHAOQIAN 4D PAIRWISE DIFFERENTIAL RESIDUAL
        ! obstdiff(i) = Tobs_current(i) - Tobs_reference(i)
        ! dsyn(i)     = Tcalc_current(i)
        ! dsyn_ref(i) = Tcalc_reference(i)
        ! residual    = dTobs - dTcalc

        do i = 1,dall
        cbst(i) = obstdiff(i) - (dsyn(i) - dsyn_ref(i))
        enddo

        call getpercentile(dall,cbst,q25,q75)
        datweight = 1.0
        
       
       ! ===== ZhaoQian 4D pseudo double-difference modification =====
       ! Disable residual outlier rejection.
       ! Keep all valid data residuals and unit data weights.
       ! Original code:
       ! do i = 1,dall
       ! if (cbst(i)<q25*threshold0 .or. cbst(i)>q75*threshold0) then
       !    datweight(i) = 0.0
       !     cbst(i) = 0
       ! endif
       ! enddo

       ! do i = 1,dall
       ! datweight(i) = 0.01+1.0/(1+0.05*exp(cbst(i)**2*threshold0))
       ! cbst(i) = cbst(i)*datweight(i)
       ! enddo

        do i = 1,nar
        rw(i) = rw(i)*datweight(iw(1+i))
        enddo

       ! norm=0
       ! do i = 1,dall
       ! datweight(i) = 0.01+1.0/(1+0.05*exp(cbst(i)**2*threshold0))
       ! cbst(i) = cbst(i)*datweight(i)
       ! enddo

       ! do i = 1,nar
       ! rw(i) = rw(i)*datweight(iw(1+i))
       ! enddo

        norm=0
        do i=1,nar
        norm(col(i))=norm(col(i))+abs(rw(i))
        enddo
        
        averdws=0
        maxnorm=0
        do i=1,maxvp
        averdws =  averdws+norm(i)
        if(norm(i)>maxnorm) maxnorm=norm(i)
        enddo
        averdws=averdws/maxvp
        write(66,*)'Maximum and Average DWS values:',maxnorm,averdws

        ! WRITE OUT RESIDUAL FOR THE FIRST AND LAST ITERATION
        if(iter.eq.1) then
                open(88,file='residualFirst.dat')
                do i=1,dall
                write(88,*) dist(i),dsyn(i),obst(i), &
                        dsyn(i)*datweight(i),obst(i)*datweight(i),datweight(i)
                enddo
                close(88)
        endif
        if(iter.eq.maxiter) then
                open(88,file='residualLast.dat')
                do i=1,dall
                write(88,*) dist(i),dsyn(i),obst(i), &
                        dsyn(i)*datweight(i),obst(i)*datweight(i),datweight(i)
                enddo
                close(88)
        endif


        weight=weight0
        nar_tmp=nar
        nars=0

        count3=0
        nvz=ny-2
        nvx=nx-2
        do k=1,nz-1
        do j=1,nvz
        do i=1,nvx
        if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
                count3=count3+1
                col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
                rw(nar+1)=2.0*weight
                iw(1+nar+1)=dall+count3
                cbst(dall+count3)=0
                nar=nar+1
        else
                count3=count3+1
                col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
                rw(nar+1)=6.0*weight
                iw(1+nar+1)=dall+count3
                rw(nar+2)=-1.0*weight
                iw(1+nar+2)=dall+count3
                col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1
                rw(nar+3)=-1.0*weight
                iw(1+nar+3)=dall+count3
                col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1
                rw(nar+4)=-1.0*weight
                iw(1+nar+4)=dall+count3
                col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i
                rw(nar+5)=-1.0*weight
                iw(1+nar+5)=dall+count3
                col(nar+5)=(k-1)*nvz*nvx+j*nvx+i
                rw(nar+6)=-1.0*weight
                iw(1+nar+6)=dall+count3
                col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i
                rw(nar+7)=-1.0*weight
                iw(1+nar+7)=dall+count3
                col(nar+7)=k*nvz*nvx+(j-1)*nvx+i
                cbst(dall+count3)=0
                nar=nar+7
        endif
        enddo
        enddo
        enddo
        m = dall + count3
        n = maxvp

        iw(1)=nar
        do i=1,nar
        iw(1+nar+i)=col(i)
        enddo
        if (nar > maxnar) stop 'increase sparsity fraction(spfra)'

        ! CALLING IRLS TO SOLVE THE PROBLEM

        leniw = 2*nar+1
        lenrw = nar
        dv = 0
        atol = 1e-6
        btol = 1e-6
        !atol = 1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
        !btol = 1e-3/(dvxd*nx*111.19/3.0)!1e-3
        conlim = 100
        itnlim = 400
        istop = 0
        anorm = 0.0
        acond = 0.0
        arnorm = 0.0
        xnorm = 0.0
        localSize = 10

        call LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
                atol, btol, conlim, itnlim, localSize, nout,&
                dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
        
        mean = sum(cbst(1:dall))/dall
        std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
        write(*,'(i2,a)')iter,'th iteration...'
!        write(*,'(a,f7.3)')'weight is:',weight
        write(*,'(a,f8.1,a,f8.2,a,f8.3)')' mean,std_devs and rms of &
                residual after weighting: ',mean*1000,'ms ',1000*std_devs,'ms ',&
                dnrm2(dall,cbst,1)/sqrt(real(dall))

        !do i =1,dall
        !cbst(i)=cbst(i)/datweight(i)
        !enddo

        !mean = sum(cbst(1:dall))/dall
        !std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
        !write(*,'(a,f8.1,a,f8.2,a,f8.3)')' residual before weighting: ',mean*1000,'ms ',1000*std_devs,'ms ',&
        !        dnrm2(dall,cbst,1)/sqrt(real(dall))
        write(66,'(i2,a)')iter,'th iteration...'
!        write(66,'(a,f7.3)')'weight is:',weight
        write(66,'(a,f8.1,a,f8.2,a,f8.3)')'mean,std_devs and rms of &
                residual: ',mean*1000,'ms ',1000*std_devs,'ms ',&
                dnrm2(dall,cbst,1)/sqrt(real(dall))

        write(*,'(a,2f7.4)')' min and max velocity variation ',&
                minval(dv),maxval(dv)
        write(66,'(a,2f7.4)')'min and max velocity variation ',&
                minval(dv),maxval(dv)

        do k=1,nz-1
        do j=1,ny-2
        do i=1,nx-2
        if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i).ge.0.500) then
                dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)=0.500
        endif
        if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i).le.-0.500) then
                dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)=-0.500
        endif
        vsf(i+1,j+1,k)=vsf(i+1,j+1,k)+dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
        if(vsf(i+1,j+1,k).lt.Minvel) vsf(i+1,j+1,k)=Minvel
        if(vsf(i+1,j+1,k).gt.Maxvel) vsf(i+1,j+1,k)=Maxvel
        enddo
        enddo
        enddo

        write(outmodel,'(a,a,i3.3)') trim(inputfile),'Measure.dat.iter',iter
        open(64,file=outmodel)
        do k=1,nz-1
        do j=1,ny-2
        do i=1,nx-2
        write(64,'(5f10.5)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
        enddo
        enddo
        enddo
        close(64)


        enddo !end iteration

        ! OUTPUT THE VELOCITY MODEL

        write(*,*)'Program finishes successfully'
        write(66,*)'Program finishes successfully'

        if(ifsyn == 1) then
                open(65,file='Vs_model.real')
                write(outsyn,'(a,a)') trim(inputfile),'Syn.dat'
                open(63,file=outsyn)
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2
                write(65,'(5f10.5)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsftrue_t1(i+1,j+1,k)
                write(63,'(5f10.5)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
                enddo
                enddo
                enddo
                close(65)
                close(63)

                ! ZHAOQIAN, output true, recovered and error delta Vs
                open(67,file='DSurfTomo.trueDeltaVs.dat')
                open(68,file='DSurfTomo.deltaVs.dat')
                open(69,file='DSurfTomo.deltaVsError.dat')
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2
                        write(67,'(5f10.5)') gozd+(j-1)*dvzd,&
                                goxd-(i-1)*dvxd,depz(k),&
                                dvs_true(i+1,j+1,k)
                        write(68,'(5f10.5)') gozd+(j-1)*dvzd,&
                                goxd-(i-1)*dvxd,depz(k),&
                                vsf(i+1,j+1,k)-vsf_init(i+1,j+1,k)
                        write(69,'(5f10.5)') gozd+(j-1)*dvzd,&
                                goxd-(i-1)*dvxd,depz(k),&
                                vsf(i+1,j+1,k)-vsf_init(i+1,j+1,k)-&
                                dvs_true(i+1,j+1,k)
                enddo
                enddo
                enddo
                close(67)
                close(68)
                close(69)

                write(*,'(a,2f10.5)') 'Recovered delta Vs min/max: ',&
                        minval(vsf(2:nx-1,2:ny-1,1:nz-1)-&
                        vsf_init(2:nx-1,2:ny-1,1:nz-1)),&
                        maxval(vsf(2:nx-1,2:ny-1,1:nz-1)-&
                        vsf_init(2:nx-1,2:ny-1,1:nz-1))
                write(*,'(a,f10.6)') 'Delta Vs RMSE: ',&
                        sqrt(sum(((vsf(2:nx-1,2:ny-1,1:nz-1)-&
                        vsf_init(2:nx-1,2:ny-1,1:nz-1))-&
                        dvs_true(2:nx-1,2:ny-1,1:nz-1))**2)/&
                        real((nx-2)*(ny-2)*(nz-1)))

                write(*,*)'Output True velocity model &
                        to Vs_model.real'
                write(*,*)'Output inverted shear velocity model &
                        to ',outsyn
                write(66,*)'Output True velocity model &
                        to Vs_model.real'
                write(66,*)'Output inverted shear velocity model &
                        to ',outsyn
        else
                write(outmodel,'(a,a)') trim(inputfile),'Measure.dat'
                open(64,file=outmodel)
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2
                write(64,'(5f10.5)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
                enddo
                enddo
                enddo
                close(64)

                ! ZHAOQIAN 4D PAIRWISE:
                ! Output delta Vs model relative to the reference model MOD_T1.
                ! DeltaVs = Vs_current - Vs_reference
                open(67,file='DSurfTomo.deltaVs.dat')
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2
                write(67,'(5f10.5)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),&
                        vsf(i+1,j+1,k)-vsf_init(i+1,j+1,k)
                enddo
                enddo
                enddo
                close(67)

                write(*,*)'Output inverted shear velocity model &
                        to ',outmodel
                write(66,*)'Output inverted shear velocity model &
                        to ',outmodel

                write(*,*)'Output 4D delta Vs model to DSurfTomo.deltaVs.dat'
                write(66,*)'Output 4D delta Vs model to DSurfTomo.deltaVs.dat'
             
        endif
        
        !close(40)
        !close(nout) !close lsmr.txt
        close(66) !close surf_tomo.log

        ! ZHAOQIAN
        deallocate(obst,obst_cur,obstdiff)
        deallocate(dsyn,dsyn_ref)
        deallocate(dist)
        deallocate(depz)
        ! ZHAOQIAN
        deallocate(depz_true_t1,depz_true_t2)
        deallocate(scxf,sczf)
        deallocate(rcxf,rczf)
        deallocate(wavetype,igrt,nrc1)
        deallocate(nsrc1,periods)
        deallocate(rw)
        deallocate(iw,col,nrow)
        ! ZHAOQIAN
        deallocate(cbst,wt,dtres,datweight)
        deallocate(dv,dvsub,dvstd,dvall)
        deallocate(norm)
        ! ZHAOQIAN
        deallocate(vsf,vsf_init)
        deallocate(vsftrue_t1,vsftrue_t2,dvs_true)
        if(kmaxRc.gt.0) then
                deallocate(tRc)
        endif
        if(kmaxRg.gt.0) then
                deallocate(tRg)
        endif
        if(kmaxLc.gt.0) then
                deallocate(tLc)
        endif
        if(kmaxLg.gt.0) then
                deallocate(tLg)
        endif

end program            

