! Electrospinning: axisymmetric electrostatics + steady slender Newtonian jet.
! Fortran 2008, double precision, SI units. See README_RU.md for model limits.
! Reduced research model with an optional transient air-plasma / electron-energy block.
module spinning_base
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: pi=acos(-1.0_dp), eps0=8.8541878128e-12_dp, gas_R=8.314462618_dp
    character(1024) :: failure_status_path=''
contains
    subroutine require(ok,msg)
        logical, intent(in) :: ok
        character(*), intent(in) :: msg
        integer :: status_unit,ios
        if (.not.ok) then
            if(len_trim(failure_status_path)>0) then
                open(newunit=status_unit,file=trim(failure_status_path),status='replace',iostat=ios)
                if(ios==0) then
                    write(status_unit,'(A)') 'FAILED: '//trim(msg)
                    write(status_unit,'(A)') 'Any final CSV from an earlier run with this prefix is stale.'
                    close(status_unit)
                end if
            end if
            write(error_unit,'(A)') 'ERROR: '//trim(msg)
            error stop 1
        end if
    end subroutine

    real(dp) function interp(x,y,t) result(v)
        real(dp), intent(in) :: x(:),y(:),t
        integer :: lo,hi,k
        lo=1; hi=size(x)
        if(t<=x(1)) then
            v=y(1); return
        end if
        if(t>=x(hi)) then
            v=y(hi); return
        end if
        do while(hi-lo>1)
            k=(hi+lo)/2
            if(x(k)>t) then
                hi=k
            else
                lo=k
            end if
        end do
        v=y(lo)+(y(hi)-y(lo))*(t-x(lo))/(x(hi)-x(lo))
    end function

    real(dp) function slope3(x0,x1,x2,y0,y1,y2) result(v)
        ! Derivative at the middle of three nonuniform points.
        real(dp), intent(in) :: x0,x1,x2,y0,y1,y2
        real(dp) :: h0,h1
        h0=x1-x0; h1=x2-x1
        v=-h1*y0/(h0*(h0+h1))+(h1-h0)*y1/(h0*h1)+h0*y2/(h1*(h0+h1))
    end function

    subroutine stretched_grid(x,length,stretch)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in) :: length,stretch
        integer :: k,n
        n=size(x)
        do k=1,n
            if(abs(stretch)<1e-12_dp) then
                x(k)=length*real(k-1,dp)/real(n-1,dp)
            else
                x(k)=length*(exp(stretch*real(k-1,dp)/real(n-1,dp))-1)/(exp(stretch)-1)
            end if
        end do
        x(1)=0; x(n)=length
    end subroutine

    subroutine tridiagonal(lo,di,up,b,x)
        real(dp), intent(in) :: lo(:),di(:),up(:),b(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: piv(size(b)),rhs(size(b)),f
        integer :: k,n
        n=size(b); piv(1)=di(1); rhs(1)=b(1)
        call require(piv(1)>0,'Non-positive tridiagonal pivot.')
        do k=2,n
            f=lo(k)/piv(k-1)
            piv(k)=di(k)-f*up(k-1); rhs(k)=b(k)-f*rhs(k-1)
            call require(piv(k)>0 .and. ieee_is_finite(piv(k)),'Singular jet system.')
        end do
        x(n)=rhs(n)/piv(n)
        do k=n-1,1,-1
            x(k)=(rhs(k)-up(k)*x(k+1))/piv(k)
        end do
    end subroutine
end module

module spinning_config
    use spinning_base
    implicit none
    ! Original physical names are retained where meaningful. All quantities SI.
    real(dp) :: V_voltage=3500, Z_max=0.10_dp, R_needle=2.5e-4_dp
    real(dp) :: Gamma=0.055_dp, Viscosity=0.2_dp, Q_pump=1.38888888888889e-10_dp
    real(dp) :: Epsilon_r=40, Conductivity=2e-6_dp, Rho_fluid=-1
    real(dp) :: T_air=25, Humidity=40, Cd_drag=0.01_dp, h_conv=100
    real(dp) :: R_max=0.15_dp, jet_radius_in=3e-6_dp, polymer_fraction=0.12_dp
    real(dp) :: rho_polymer=1200, rho_solvent=997, cp_polymer=1700, cp_solvent=4180
    real(dp) :: thermal_conductivity_fluid=0.6_dp
    real(dp) :: latent_heat=2.44e6_dp, molar_mass=0.01801528_dp, k_mass=0.01_dp
    real(dp) :: psat_ref=3169, T_psat_ref=298.15_dp, pressure_air=101325
    real(dp) :: viscosity_exponent=3, activation_energy=20000
    real(dp) :: conductivity_exponent=0, Jet_current_A=-1, current_coefficient=1
    real(dp) :: g_z=9.81_dp, air_velocity=0, rho_air=1.184_dp
    real(dp) :: stretch_r=11, stretch_z=7, stretch_jet=6
    real(dp) :: poisson_tol=1e-10_dp, jet_tol=2e-7_dp, outer_tol=2e-5_dp
    real(dp) :: relaxation=0.25_dp, charge_relaxation=0.35_dp
    integer :: Nr=96, Nz=160, Nj=401, max_pcg=5000, max_jet=1200, max_outer=100
    logical :: space_charge=.true., evaporation=.true., joule_heat=.true., viscous_heat=.true.
    logical :: polarization_force=.true.
    character(16) :: side_boundary='grounded'
    character(32) :: solvent_name='water'
    character(512) :: output_prefix=''
    real(dp) :: mp_in,ms_in,u_in,total_current
    namelist /spinning/ V_voltage,Z_max,R_needle,Gamma,Viscosity,Q_pump,Epsilon_r,Conductivity, &
        Rho_fluid,T_air,Humidity,Cd_drag,h_conv,R_max,jet_radius_in,polymer_fraction,rho_polymer, &
        rho_solvent,cp_polymer,cp_solvent,thermal_conductivity_fluid,latent_heat,molar_mass,k_mass,psat_ref,T_psat_ref, &
        pressure_air,viscosity_exponent,activation_energy,conductivity_exponent,Jet_current_A, &
        current_coefficient,g_z,air_velocity,rho_air,stretch_r,stretch_z,stretch_jet,poisson_tol, &
        jet_tol,outer_tol,relaxation,charge_relaxation,Nr,Nz,Nj,max_pcg,max_jet,max_outer, &
        space_charge,evaporation,joule_heat,viscous_heat,polarization_force,side_boundary,solvent_name,output_prefix
contains
    subroutine read_configuration(path)
        character(*), intent(in) :: path
        integer :: u,ios,k,n
        character(2048) :: line,msg
        real(dp) :: old(14)
        open(newunit=u,file=path,status='old',action='read',iostat=ios,iomsg=msg)
        call require(ios==0,'Cannot read configuration: '//trim(path)//' '//trim(msg))
        line=''
        do
            read(u,'(A)',iostat=ios) line
            call require(ios==0,'Empty configuration.')
            line=adjustl(line)
            if(len_trim(line)==0) cycle
            if(line(1:1)=='!' .or. line(1:1)=='#') cycle
            exit
        end do
        rewind(u)
        if(line(1:1)=='&') then
            read(u,nml=spinning,iostat=ios,iomsg=msg)
            call require(ios==0,'Invalid namelist: '//trim(msg))
        else
            ! Legacy: exactly the original 11 required + 3 optional numbers.
            old=0; old(12)=0; old(13)=Cd_drag; old(14)=h_conv; n=0
            do
                read(u,'(A)',iostat=ios) line
                if(ios<0) exit
                call require(ios==0,'Legacy config read error.')
                k=index(line,'!'); if(k>0) line=line(:k-1)
                k=index(line,'#'); if(k>0) line=line(:k-1)
                if(len_trim(line)==0) cycle
                n=n+1
                call require(n<=14,'Legacy config has more than 14 values; use namelist.')
                read(line,*,iostat=ios) old(n)
                call require(ios==0,'Invalid number in legacy config.')
            end do
            call require(n>=11,'Legacy config needs at least 11 values.')
            V_voltage=old(1); Z_max=old(2); R_needle=old(3); Gamma=old(4)
            Viscosity=old(5); Q_pump=old(6); Epsilon_r=old(7); Conductivity=old(8)
            Rho_fluid=old(9); T_air=old(10); Humidity=old(11)
            Cd_drag=old(13); h_conv=old(14); jet_radius_in=R_needle
            write(*,'(A)') 'WARNING: legacy inlet radius = needle radius; Taylor cone is unresolved.'
            write(*,'(A)') 'WARNING: legacy alpha_whip is ignored; it violates mass conservation.'
            write(*,'(A)') 'WARNING: solvent and polymer data use defaults; supply measured material properties.'
        end if
        close(u)
    end subroutine

    subroutine set_preset(name)
        character(*), intent(in) :: name
        select case(trim(name))
        case('quick')
            Nr=64; Nz=96; Nj=201
        case('standard')
            Nr=96; Nz=160; Nj=401
        case('fine')
            Nr=160; Nz=256; Nj=801
        case('heavy')
            Nr=256; Nz=400; Nj=1601
        case default
            call require(.false.,'Preset must be quick, standard, fine or heavy.')
        end select
    end subroutine

    subroutine validate_config()
        real(dp) :: vals(46),rho0
        vals=[V_voltage,Z_max,R_needle,Gamma,Viscosity,Q_pump,Epsilon_r,Conductivity,Rho_fluid, &
            T_air,Humidity,Cd_drag,h_conv,R_max,jet_radius_in,polymer_fraction,rho_polymer, &
            rho_solvent,cp_polymer,cp_solvent,latent_heat,molar_mass,k_mass,psat_ref,T_psat_ref, &
            pressure_air,viscosity_exponent,activation_energy,conductivity_exponent,Jet_current_A, &
            current_coefficient,g_z,air_velocity,rho_air,stretch_r,stretch_z,stretch_jet,poisson_tol, &
            jet_tol,outer_tol,relaxation,charge_relaxation,real(Nr,dp),real(Nz,dp),real(Nj,dp), &
            thermal_conductivity_fluid]
        call require(all(ieee_is_finite(vals)),'All input numbers must be finite.')
        call require(Z_max>0 .and. R_max>0,'Domain lengths must be positive.')
        call require(R_needle>0 .and. R_needle<R_max,'Need 0 < R_needle < R_max.')
        call require(jet_radius_in>0 .and. jet_radius_in<=R_needle,'Need 0 < inlet jet radius <= needle radius.')
        call require(Q_pump>0 .and. Viscosity>0 .and. Gamma>=0,'Invalid flow, viscosity or surface tension.')
        call require(Epsilon_r>=1 .and. Conductivity>=0,'Invalid permittivity or conductivity.')
        call require(polymer_fraction>0 .and. polymer_fraction<1,'Polymer mass fraction must be in (0,1).')
        call require(rho_polymer>0 .and. rho_solvent>0,'Densities must be positive.')
        call require(Rho_fluid==-1 .or. Rho_fluid>0,'Rho_fluid must be -1 or positive.')
        call require(cp_polymer>0 .and. cp_solvent>0,'Heat capacities must be positive.')
        call require(thermal_conductivity_fluid>0,'Thermal conductivity must be positive.')
        call require(latent_heat>0 .and. molar_mass>0 .and. psat_ref>0,'Invalid solvent properties.')
        call require(T_psat_ref>=200 .and. T_psat_ref<=420,'T_psat_ref is in kelvins, range 200..420.')
        call require(T_air+273.15_dp>200 .and. T_air+273.15_dp<420,'Ambient temperature outside model range.')
        call require(Humidity>=0 .and. Humidity<=100,'Humidity/saturation ratio must be 0..100 percent.')
        call require(k_mass>=0 .and. h_conv>=0 .and. Cd_drag>=0 .and. rho_air>=0,'Negative transfer coefficient.')
        call require(pressure_air>0,'Air pressure must be positive.')
        call require(viscosity_exponent>=0 .and. activation_energy>=0 .and. conductivity_exponent>=0, &
            'Material exponents and activation energy must be nonnegative.')
        call require(Jet_current_A==-1 .or. Jet_current_A>=0,'Jet_current_A must be -1 (estimate) or >=0.')
        call require(current_coefficient>=0,'Current coefficient must be nonnegative.')
        call require(.not.(V_voltage==0 .and. Jet_current_A>0),'Specified jet current needs nonzero voltage.')
        call require(Nr>=8 .and. Nz>=8 .and. Nj>=8,'Each mesh needs at least 8 points/cells.')
        call require(Nr<=4096 .and. Nz<=4096 .and. Nj<=200000,'Mesh exceeds supported memory guard.')
        call require(min(stretch_r,stretch_z,stretch_jet)>=0 .and. &
            max(stretch_r,stretch_z,stretch_jet)<=16,'Grid stretch must lie in [0,16].')
        call require(poisson_tol>=1e-14_dp .and. poisson_tol<0.01_dp,'Invalid Poisson tolerance.')
        call require(jet_tol>=1e-12_dp .and. jet_tol<0.01_dp,'Invalid jet tolerance.')
        call require(outer_tol>=1e-10_dp .and. outer_tol<0.01_dp,'Invalid outer tolerance.')
        call require(relaxation>0 .and. relaxation<=1 .and. charge_relaxation>0 .and. &
            charge_relaxation<=1,'Relaxations must lie in (0,1].')
        call require(min(max_pcg,max_jet,max_outer)>0,'Iteration limits must be positive.')
        call require(trim(side_boundary)=='grounded' .or. trim(side_boundary)=='insulating', &
            'side_boundary must be grounded or insulating.')
        if(Rho_fluid>0) then
            rho_polymer=Rho_fluid; rho_solvent=Rho_fluid
        end if
        rho0=1/(polymer_fraction/rho_polymer+(1-polymer_fraction)/rho_solvent)
        mp_in=rho0*Q_pump*polymer_fraction; ms_in=rho0*Q_pump*(1-polymer_fraction)
        u_in=Q_pump/(pi*jet_radius_in**2)
        total_current=Jet_current_A
        if(Jet_current_A<0) then
            total_current=current_coefficient*sqrt(Gamma*Conductivity*Q_pump)
            write(*,'(A)') 'WARNING: current uses an uncalibrated cone-jet scaling; measured Jet_current_A is preferred.'
        end if
        if(V_voltage==0) then
            total_current=0
        else
            total_current=sign(total_current,V_voltage)
        end if
        call require(all(ieee_is_finite([mp_in,ms_in,u_in,total_current])),'Input scaling overflow.')
    end subroutine
end module

module axisymmetric_field
    use spinning_base
    implicit none
    type :: field_t
        integer :: nr,nz,last_iterations=0
        real(dp) :: radius,length,voltage,last_residual=1
        logical :: insulating,full_top
        real(dp), allocatable :: rf(:),zf(:),r(:),z(:),area(:),vol(:,:)
        real(dp), allocatable :: gr(:,:),gz(:,:),diag(:,:),invpiv(:,:),bc(:,:)
        real(dp), allocatable :: phi(:,:),rho(:,:),er(:,:),ez(:,:)
    end type
contains
    subroutine init_field(f,nr,nz,radius,length,sr,sz,voltage,electrode_radius,insulating,full_top)
        type(field_t), intent(out) :: f
        integer, intent(in) :: nr,nz
        real(dp), intent(in) :: radius,length,sr,sz,voltage,electrode_radius
        logical, intent(in) :: insulating,full_top
        real(dp) :: conduct,top_area,pivot,cross_term
        integer :: i,j
        f%nr=nr; f%nz=nz; f%radius=radius; f%length=length; f%voltage=voltage
        f%insulating=insulating; f%full_top=full_top
        allocate(f%rf(nr+1),f%zf(nz+1),f%r(nr),f%z(nz),f%area(nr),f%vol(nr,nz))
        allocate(f%gr(nr-1,nz),f%gz(nr,nz-1),f%diag(nr,nz),f%invpiv(nr,nz),f%bc(nr,nz))
        allocate(f%phi(nr,nz),f%rho(nr,nz),f%er(nr,nz),f%ez(nr,nz))
        call stretched_grid(f%rf,radius,sr); call stretched_grid(f%zf,length,sz)
        f%r=(f%rf(:nr)+f%rf(2:))/2; f%z=(f%zf(:nz)+f%zf(2:))/2
        f%area=pi*(f%rf(2:)**2-f%rf(:nr)**2)
        f%diag=0; f%bc=0; f%rho=0; f%phi=0
        do j=1,nz
            f%vol(:,j)=f%area*(f%zf(j+1)-f%zf(j))
            do i=1,nr-1
                conduct=2*pi*f%rf(i+1)*(f%zf(j+1)-f%zf(j))/(f%r(i+1)-f%r(i))
                f%gr(i,j)=conduct
                f%diag(i,j)=f%diag(i,j)+conduct; f%diag(i+1,j)=f%diag(i+1,j)+conduct
            end do
            if(.not.insulating) f%diag(nr,j)=f%diag(nr,j)+ &
                2*pi*radius*(f%zf(j+1)-f%zf(j))/(radius-f%r(nr))
        end do
        do j=1,nz-1
            f%gz(:,j)=f%area/(f%z(j+1)-f%z(j))
            f%diag(:,j)=f%diag(:,j)+f%gz(:,j)
            f%diag(:,j+1)=f%diag(:,j+1)+f%gz(:,j)
        end do
        do i=1,nr
            top_area=max(0.0_dp,pi*(min(f%rf(i+1),electrode_radius)**2-f%rf(i)**2))
            if(full_top) top_area=f%area(i)
            conduct=top_area/f%z(1)
            f%diag(i,1)=f%diag(i,1)+conduct; f%bc(i,1)=conduct*voltage
            f%diag(i,nz)=f%diag(i,nz)+f%area(i)/(length-f%z(nz))
        end do
        ! Modified incomplete Cholesky (tau=0.97) for the SPD integrated FV matrix.
        ! Dropped cross terms are compensated on the diagonal; both grid directions enter.
        do j=1,nz
            do i=1,nr
                pivot=f%diag(i,j)
                if(i>1) then
                    cross_term=0
                    if(j<nz) cross_term=f%gz(i-1,j)
                    pivot=pivot-f%gr(i-1,j)*(f%gr(i-1,j)+0.97_dp*cross_term)*f%invpiv(i-1,j)
                end if
                if(j>1) then
                    cross_term=0
                    if(i<nr) cross_term=f%gr(i,j-1)
                    pivot=pivot-f%gz(i,j-1)*(f%gz(i,j-1)+0.97_dp*cross_term)*f%invpiv(i,j-1)
                end if
                call require(pivot>0,'Incomplete Cholesky pivot failure; reduce mesh stretch.')
                f%invpiv(i,j)=1/pivot
            end do
        end do
    end subroutine

    subroutine apply_matrix(f,x,y)
        type(field_t), intent(in) :: f
        real(dp), intent(in) :: x(:,:)
        real(dp), intent(out) :: y(:,:)
        integer :: j
        y=f%diag*x
        y(:f%nr-1,:)=y(:f%nr-1,:)-f%gr*x(2:,:)
        y(2:,:)=y(2:,:)-f%gr*x(:f%nr-1,:)
        do j=1,f%nz-1
            y(:,j)=y(:,j)-f%gz(:,j)*x(:,j+1)
            y(:,j+1)=y(:,j+1)-f%gz(:,j)*x(:,j)
        end do
    end subroutine

    subroutine precondition(f,b,x)
        type(field_t), intent(in) :: f
        real(dp), intent(in) :: b(:,:)
        real(dp), intent(out) :: x(:,:)
        integer :: i,j
        x(1,1)=b(1,1)
        do i=2,f%nr
            x(i,1)=b(i,1)+f%gr(i-1,1)*f%invpiv(i-1,1)*x(i-1,1)
        end do
        do j=2,f%nz
            x(1,j)=b(1,j)+f%gz(1,j-1)*f%invpiv(1,j-1)*x(1,j-1)
            do i=2,f%nr
                x(i,j)=b(i,j)+f%gr(i-1,j)*f%invpiv(i-1,j)*x(i-1,j) &
                    +f%gz(i,j-1)*f%invpiv(i,j-1)*x(i,j-1)
            end do
        end do
        do j=f%nz,1,-1
            do i=f%nr,1,-1
                if(i<f%nr) x(i,j)=x(i,j)+f%gr(i,j)*x(i+1,j)
                if(j<f%nz) x(i,j)=x(i,j)+f%gz(i,j)*x(i,j+1)
                x(i,j)=x(i,j)*f%invpiv(i,j)
            end do
        end do
    end subroutine

    subroutine solve_field(f,tol,maxit)
        type(field_t), intent(inout) :: f
        real(dp), intent(in) :: tol
        integer, intent(in) :: maxit
        real(dp) :: b(f%nr,f%nz),r(f%nr,f%nz),p(f%nr,f%nz),q(f%nr,f%nz),w(f%nr,f%nz)
        real(dp) :: bn,rr,rrnew,denom,alpha,res
        integer :: k
        b=f%bc+f%rho*f%vol/eps0
        bn=sqrt(sum(b*b))
        if(bn<tiny(1.0_dp)) then
            f%phi=0; f%last_residual=0; f%last_iterations=0; return
        end if
        call apply_matrix(f,f%phi,q); r=b-q
        res=sqrt(sum(r*r))/bn
        if(res<=tol) then
            f%last_residual=res; f%last_iterations=0; return
        end if
        call precondition(f,r,w); p=w; rr=sum(r*w)
        do k=1,maxit
            call apply_matrix(f,p,q); denom=sum(p*q)
            call require(denom>0 .and. ieee_is_finite(denom),'Poisson PCG lost positive definiteness.')
            alpha=rr/denom; f%phi=f%phi+alpha*p; r=r-alpha*q
            res=sqrt(sum(r*r))/bn
            if(res<=tol) then
                ! Check the actual residual, not only the recursively updated one.
                call apply_matrix(f,f%phi,q); r=b-q; res=sqrt(sum(r*r))/bn
                if(res<=5*tol) exit
                call precondition(f,r,w); p=w; rr=sum(r*w); cycle
            end if
            call precondition(f,r,w); rrnew=sum(r*w)
            p=w+(rrnew/rr)*p; rr=rrnew
        end do
        f%last_iterations=min(k,maxit); f%last_residual=res
        call require(res<=5*tol,'Poisson did not converge. Increase max_pcg or check the mesh.')
        call require(all(ieee_is_finite(f%phi)),'Non-finite electrostatic potential.')
    end subroutine

    subroutine field_gradients(f,electrode_radius)
        type(field_t), intent(inout) :: f
        real(dp), intent(in) :: electrode_radius
        real(dp) :: top
        integer :: i,j,n,m
        n=f%nr; m=f%nz
        do j=1,m
            f%er(1,j)=-2*f%r(1)*(f%phi(2,j)-f%phi(1,j))/(f%r(2)**2-f%r(1)**2)
            do i=2,n-1
                f%er(i,j)=-slope3(f%r(i-1),f%r(i),f%r(i+1), &
                    f%phi(i-1,j),f%phi(i,j),f%phi(i+1,j))
            end do
            if(f%insulating) then
                f%er(n,j)=0
            else
                f%er(n,j)=-slope3(f%r(n-1),f%r(n),f%radius,f%phi(n-1,j),f%phi(n,j),0.0_dp)
            end if
        end do
        do i=1,n
            if(f%r(i)<electrode_radius .or. f%full_top) then
                top=f%voltage
                f%ez(i,1)=-slope3(0.0_dp,f%z(1),f%z(2),top,f%phi(i,1),f%phi(i,2))
            else
                f%ez(i,1)=-slope3(-f%z(1),f%z(1),f%z(2), &
                    f%phi(i,1),f%phi(i,1),f%phi(i,2))
            end if
            do j=2,m-1
                f%ez(i,j)=-slope3(f%z(j-1),f%z(j),f%z(j+1), &
                    f%phi(i,j-1),f%phi(i,j),f%phi(i,j+1))
            end do
            f%ez(i,m)=-slope3(f%z(m-1),f%z(m),f%length,f%phi(i,m-1),f%phi(i,m),0.0_dp)
        end do
    end subroutine

    subroutine sample_jet_field(f,z,a,e)
        type(field_t), intent(in) :: f
        real(dp), intent(in) :: z(:),a(:)
        real(dp), intent(out) :: e(:)
        real(dp) :: ezline(f%nz),h1,h2,top_e,bottom_e,xp(f%nz+2),yp(f%nz+2)
        integer :: k,j
        xp(1)=0; xp(2:f%nz+1)=f%z; xp(f%nz+2)=f%length
        do k=1,size(z)
            do j=1,f%nz
                ezline(j)=interp(f%r,f%ez(:,j),a(k))
            end do
            ! Second-order one-sided surface derivatives; entry lies under electrode.
            h1=f%z(1); h2=f%z(2)-f%z(1)
            top_e=(2*h1+h2)/(h1*(h1+h2))*f%voltage &
                -(h1+h2)/(h1*h2)*interp(f%r,f%phi(:,1),a(k)) &
                +h1/(h2*(h1+h2))*interp(f%r,f%phi(:,2),a(k))
            h1=f%length-f%z(f%nz); h2=f%z(f%nz)-f%z(f%nz-1)
            bottom_e=(h1+h2)/(h1*h2)*interp(f%r,f%phi(:,f%nz),a(k)) &
                -h1/(h2*(h1+h2))*interp(f%r,f%phi(:,f%nz-1),a(k))
            yp(1)=top_e; yp(2:f%nz+1)=ezline; yp(f%nz+2)=bottom_e
            e(k)=interp(xp,yp,z(k))
        end do
    end subroutine

    subroutine deposit_charge(f,z,a,lambda,error)
        ! Cross-section averaged source. Exact annular volumes conserve line charge.
        type(field_t), intent(inout) :: f
        real(dp), intent(in) :: z(:),a(:),lambda(:)
        real(dp), intent(out) :: error
        real(dp) :: aa,ll,fraction,check
        integer :: i,j
        error=0; f%rho=0
        do j=1,f%nz
            aa=interp(z,a,f%z(j)); ll=interp(z,lambda,f%z(j))
            call require(aa>0 .and. aa<f%radius,'Jet left the radial domain.')
            do i=1,f%nr
                fraction=max(0.0_dp,min(f%rf(i+1),aa)**2-min(f%rf(i),aa)**2)/aa**2
                f%rho(i,j)=ll*fraction/f%area(i)
            end do
            check=sum(f%rho(:,j)*f%area)
            error=max(error,abs(check-ll)/max(abs(ll),1e-30_dp))
        end do
    end subroutine

    subroutine field_self_tests()
        type(field_t) :: f
        real(dp) :: exact,err,prev,ratio,dummy
        real(dp) :: zz(2),aa(2),ll(2)
        integer :: i,j,k,n
        call init_field(f,20,32,0.1_dp,0.2_dp,2.0_dp,3.0_dp,1000.0_dp,0.1_dp,.true.,.true.)
        call solve_field(f,1e-12_dp,10000)
        err=0
        do j=1,f%nz
            err=max(err,maxval(abs(f%phi(:,j)-1000*(1-f%z(j)/f%length)))/1000)
        end do
        write(*,'(A,ES12.4)') 'TEST planar capacitor relative error: ',err
        call require(err<1e-8_dp,'Planar capacitor test failed.')
        prev=0
        do k=1,3
            n=12*2**(k-1)
            call init_field(f,n,n,0.1_dp,0.2_dp,0.0_dp,0.0_dp,0.0_dp,0.1_dp,.false.,.true.)
            do j=1,n
                do i=1,n
                    f%rho(i,j)=eps0*((4/f%radius**2)*(f%z(j)/f%length)*(1-f%z(j)/f%length) &
                        +2/f%length**2*(1-(f%r(i)/f%radius)**2))
                end do
            end do
            call solve_field(f,1e-12_dp,10000)
            err=0
            do j=1,n
                do i=1,n
                    exact=(1-(f%r(i)/f%radius)**2)*(f%z(j)/f%length)*(1-f%z(j)/f%length)
                    err=max(err,abs(f%phi(i,j)-exact))
                end do
            end do
            write(*,'(A,I4,A,ES12.4)') 'TEST manufactured Poisson n=',n,' max error: ',err
            if(k>1) then
                ratio=prev/err
                call require(ratio>3.0_dp,'Poisson grid convergence is below expected order.')
            end if
            prev=err
        end do
        zz=[0.0_dp,f%length]; aa=[1e-7_dp,2e-3_dp]; ll=[-1e-9_dp,3e-9_dp]
        call deposit_charge(f,zz,aa,ll,dummy)
        write(*,'(A,ES12.4)') 'TEST annular charge conservation: ',dummy
        call require(dummy<1e-12_dp,'Charge deposition test failed.')
        write(*,'(A)') 'ALL FIELD SELF-TESTS PASSED'
    end subroutine
end module

module slender_jet
    use spinning_config
    implicit none
    type :: jet_t
        integer :: n,last_iterations=0
        real(dp) :: momentum_residual=1,energy_residual=1,mass_residual=1
        real(dp), allocatable :: z(:),u(:),a(:),ms(:),temp(:),eta(:),sigma(:),area(:),e(:)
        real(dp), allocatable :: lambda(:),evap_loss(:),heat_in(:)
    end type
contains
    real(dp) function saturation_pressure(t) result(p)
        real(dp), intent(in) :: t
        p=psat_ref*exp(latent_heat*molar_mass/gas_R*(1/T_psat_ref-1/t))
    end function

    subroutine local_properties(ms,t,u,a,eta,sigma)
        real(dp), intent(in) :: ms,t,u
        real(dp), intent(out) :: a,eta,sigma
        real(dp) :: w,sv,sv0,volume
        volume=mp_in/rho_polymer+ms/rho_solvent
        a=sqrt(volume/(pi*u)); w=mp_in/(mp_in+ms)
        eta=Viscosity*(w/polymer_fraction)**viscosity_exponent &
            *exp(activation_energy/gas_R*(1/t-1/(T_air+273.15_dp)))
        sv=(ms/rho_solvent)/volume
        sv0=(ms_in/rho_solvent)/Q_pump
        if(conductivity_exponent==0) then
            sigma=Conductivity
        else
            sigma=Conductivity*(sv/sv0)**conductivity_exponent
        end if
    end subroutine

    real(dp) function evaporation_flux(ms,t) result(flux)
        real(dp), intent(in) :: ms,t
        real(dp) :: solvent_volume_fraction,ambient_conc
        flux=0
        if(.not.evaporation .or. k_mass==0) return
        solvent_volume_fraction=(ms/rho_solvent)/(mp_in/rho_polymer+ms/rho_solvent)
        ambient_conc=Humidity/100*saturation_pressure(T_air+273.15_dp)/(T_air+273.15_dp)
        ! Ideal volume-fraction activity; evaporation only (no condensation model).
        flux=max(0.0_dp,k_mass*molar_mass/gas_R*(solvent_volume_fraction*saturation_pressure(t)/t-ambient_conc))
    end function

    subroutine implicit_mass(t,u,ds,ms_old,ms,a,eta,sigma)
        real(dp), intent(in) :: t,u,ds,ms_old
        real(dp), intent(out) :: ms,a,eta,sigma
        real(dp) :: lo,hi,mid,res,volume,flux,df,deriv,trial
        integer :: k
        if(.not.evaporation .or. k_mass==0 .or. ms_old==0) then
            ms=ms_old; call local_properties(ms,t,u,a,eta,sigma); return
        end if
        lo=0; hi=ms_old; mid=ms_old
        do k=1,70
            volume=mp_in/rho_polymer+mid/rho_solvent
            a=sqrt(volume/(pi*u)); flux=evaporation_flux(mid,t)
            res=mid-ms_old+ds*2*pi*a*flux
            if(abs(res)<2e-14_dp*ms_old) exit
            if(res>0) then
                hi=mid
            else
                lo=mid
            end if
            deriv=1
            if(flux>0) then
                df=k_mass*molar_mass/gas_R*saturation_pressure(t)/t &
                    *(mp_in/rho_polymer)/(rho_solvent*volume**2)
                deriv=1+ds*2*pi*a*(flux/(2*rho_solvent*volume)+df)
            end if
            trial=mid-res/deriv
            if(trial<=lo .or. trial>=hi) trial=(lo+hi)/2
            mid=trial
        end do
        call require(k<=70,'Implicit solvent balance did not converge.')
        ms=mid
        call local_properties(ms,t,u,a,eta,sigma)
    end subroutine

    subroutine thermal_residual(t,u,du,ds,ef,ms_old,t_old,res,ms,a,eta,sigma,heat)
        real(dp), intent(in) :: t,u,du,ds,ef,ms_old,t_old
        real(dp), intent(out) :: res,ms,a,eta,sigma,heat
        real(dp) :: power,area
        call implicit_mass(t,u,ds,ms_old,ms,a,eta,sigma)
        area=pi*a*a
        power=2*pi*a*h_conv*(T_air+273.15_dp-t)
        if(joule_heat) power=power+sigma*area*ef**2
        if(viscous_heat) power=power+3*eta*area*du**2
        heat=power*ds
        res=(mp_in*cp_polymer+ms_old*cp_solvent)*(t-t_old)-heat+(ms_old-ms)*latent_heat
    end subroutine

    subroutine init_jet(jet)
        type(jet_t), intent(out) :: jet
        integer :: n
        n=Nj; jet%n=n
        allocate(jet%z(n),jet%u(n),jet%a(n),jet%ms(n),jet%temp(n),jet%eta(n),jet%sigma(n), &
            jet%area(n),jet%e(n),jet%lambda(n),jet%evap_loss(n),jet%heat_in(n))
        call stretched_grid(jet%z,Z_max,stretch_jet)
        ! Ballistic gravity provides a positive initial guess for very slow feed speeds.
        ! It is only an initial guess; the converged solution satisfies the viscous system.
        jet%u=sqrt(u_in**2+2*max(g_z,0.0_dp)*jet%z)
        jet%a=sqrt(Q_pump/(pi*jet%u)); jet%e=0; jet%lambda=0
        jet%ms=ms_in; jet%temp=T_air+273.15_dp; jet%eta=Viscosity; jet%sigma=Conductivity
        jet%area=pi*jet%a**2; jet%evap_loss=0; jet%heat_in=0
    end subroutine

    subroutine thermal_sweep(jet)
        type(jet_t), intent(inout) :: jet
        real(dp) :: lo,hi,t,fl,fh,res,ds,du,ms,a,eta,sigma,heat,scale,merr,eerr,rplus,deriv,trial
        real(dp) :: heat_capacity
        integer :: i,k
        jet%ms(1)=ms_in; jet%temp(1)=T_air+273.15_dp
        call local_properties(ms_in,jet%temp(1),jet%u(1),jet%a(1),jet%eta(1),jet%sigma(1))
        jet%evap_loss=0; jet%heat_in=0; jet%mass_residual=0; jet%energy_residual=0
        do i=2,jet%n
            ds=jet%z(i)-jet%z(i-1); du=(jet%u(i)-jet%u(i-1))/ds
            lo=200; hi=420
            call thermal_residual(lo,jet%u(i),du,ds,jet%e(i),jet%ms(i-1),jet%temp(i-1), &
                fl,ms,a,eta,sigma,heat)
            call thermal_residual(hi,jet%u(i),du,ds,jet%e(i),jet%ms(i-1),jet%temp(i-1), &
                fh,ms,a,eta,sigma,heat)
            call require(fl<=0 .and. fh>=0,'Temperature leaves 200..420 K; heat/material model invalid.')
            t=min(hi,max(lo,jet%temp(i)))
            heat_capacity=mp_in*cp_polymer+jet%ms(i-1)*cp_solvent
            do k=1,70
                call thermal_residual(t,jet%u(i),du,ds,jet%e(i),jet%ms(i-1),jet%temp(i-1), &
                    res,ms,a,eta,sigma,heat)
                if(abs(res)<1e-11_dp*heat_capacity) exit
                if(res>0) then
                    hi=t
                else
                    lo=t
                end if
                call thermal_residual(t+0.001_dp,jet%u(i),du,ds,jet%e(i),jet%ms(i-1),jet%temp(i-1), &
                    rplus,ms,a,eta,sigma,heat)
                deriv=(rplus-res)/0.001_dp
                trial=t-res/max(deriv,heat_capacity)
                if(trial<=lo .or. trial>=hi) trial=(lo+hi)/2
                t=trial
            end do
            call require(k<=70,'Implicit energy balance did not converge.')
            call thermal_residual(t,jet%u(i),du,ds,jet%e(i),jet%ms(i-1),jet%temp(i-1), &
                res,ms,a,eta,sigma,heat)
            jet%ms(i)=ms; jet%a(i)=a; jet%temp(i)=t; jet%eta(i)=eta; jet%sigma(i)=sigma
            jet%evap_loss(i)=jet%ms(i-1)-ms; jet%heat_in(i)=heat
            scale=(mp_in*cp_polymer+ms_in*cp_solvent)*(T_air+273.15_dp)
            eerr=abs(res)/scale
            merr=abs(ms-jet%ms(i-1)+ds*2*pi*a*evaporation_flux(ms,t))/(mp_in+ms_in)
            jet%energy_residual=max(jet%energy_residual,eerr)
            jet%mass_residual=max(jet%mass_residual,merr)
        end do
        jet%area=pi*jet%a**2
        jet%lambda=(total_current-jet%sigma*jet%area*jet%e)/jet%u
        call require(all(ieee_is_finite(jet%temp)) .and. all(ieee_is_finite(jet%eta)), &
            'Non-finite jet temperature or viscosity.')
    end subroutine

    subroutine momentum_system(jet,lo,di,up,b)
        type(jet_t), intent(in) :: jet
        real(dp), intent(out) :: lo(:),di(:),up(:),b(:)
        real(dp) :: cw,ce,dzcv,adv,drag,force,ha,hb,east_e,west_e
        integer :: i,n
        n=jet%n; lo=0; up=0; di=0; b=0
        di(1)=1; b(1)=u_in
        do i=2,n-1
            dzcv=(jet%z(i+1)-jet%z(i-1))/2
            ha=jet%eta(i-1)*jet%area(i-1); hb=jet%eta(i)*jet%area(i)
            cw=3*2*ha*hb/(ha+hb)/(jet%z(i)-jet%z(i-1))
            ha=hb; hb=jet%eta(i+1)*jet%area(i+1)
            ce=3*2*ha*hb/(ha+hb)/(jet%z(i+1)-jet%z(i))
            adv=(mp_in+jet%ms(i))*dzcv/(jet%z(i)-jet%z(i-1))
            drag=pi*rho_air*Cd_drag*jet%a(i)*abs(jet%u(i)-air_velocity)
            force=jet%lambda(i)*jet%e(i)+((mp_in+jet%ms(i))/jet%u(i)-rho_air*jet%area(i))*g_z
            if(polarization_force) then
                ! Leading slender-jet dielectric Maxwell term: A*(eps-eps0)*E*dE/dz.
                ! The field remains an air-domain approximation; see model assumptions.
                east_e=(jet%e(i)+jet%e(i+1))/2; west_e=(jet%e(i-1)+jet%e(i))/2
                force=force+0.5_dp*eps0*(Epsilon_r-1)*jet%area(i)*(east_e**2-west_e**2)/dzcv
            end if
            lo(i)=-cw-adv; di(i)=cw+ce+adv+drag*dzcv; up(i)=-ce
            b(i)=(force+drag*air_velocity)*dzcv+pi*Gamma*(jet%a(i+1)-jet%a(i-1))/2
        end do
        ! Free outflow: vanishing Newtonian extensional stress (du/dz=0).
        lo(n)=-1; di(n)=1; b(n)=0
    end subroutine

    subroutine solve_jet(jet)
        type(jet_t), intent(inout) :: jet
        real(dp) :: lo(jet%n),di(jet%n),up(jet%n),b(jet%n),trial(jet%n),previous(jet%n)
        real(dp) :: resvec(jet%n),scale,change,step,res
        integer :: k,i,n
        n=jet%n; res=huge(1.0_dp)
        do k=1,max_jet
            call thermal_sweep(jet)
            call momentum_system(jet,lo,di,up,b)
            call tridiagonal(lo,di,up,b,trial)
            previous=jet%u; step=relaxation
            do i=2,n
                if(trial(i)<=0) step=min(step,0.8_dp*jet%u(i)/(jet%u(i)-trial(i)))
            end do
            call require(step>1e-10_dp,'No positive steady jet; check current, field and inlet conditions.')
            jet%u=previous+step*(trial-previous); jet%u(1)=u_in
            call require(all(jet%u>0) .and. all(ieee_is_finite(jet%u)),'Invalid jet velocity.')
            change=maxval(abs(jet%u-previous)/max(abs(jet%u),u_in))
            if(change<jet_tol) then
                call thermal_sweep(jet)
                call momentum_system(jet,lo,di,up,b)
                resvec=di*jet%u-b
                resvec(2:)=resvec(2:)+lo(2:)*jet%u(:n-1)
                resvec(:n-1)=resvec(:n-1)+up(:n-1)*jet%u(2:)
                scale=max(maxval(abs(b(2:n-1))),maxval(abs(di(2:n-1)*jet%u(2:n-1))),1e-30_dp)
                res=maxval(abs(resvec(2:n-1)))/scale
                if(res<10*jet_tol) exit
            end if
        end do
        jet%last_iterations=min(k,max_jet)
        call require(k<=max_jet,'Steady jet did not converge; inspect material data or increase max_jet.')
        do i=1,n
            call require(saturation_pressure(jet%temp(i))<pressure_air, &
                'Boiling condition reached in the converged jet; this model excludes boiling.')
        end do
        jet%momentum_residual=res
    end subroutine
end module

module spinning_output
    use axisymmetric_field
    use slender_jet
    implicit none
contains
    integer function output_file(name) result(u)
        character(*), intent(in) :: name
        integer :: ios
        character(1024) :: msg
        open(newunit=u,file=trim(output_prefix)//name,status='replace',action='write',iostat=ios,iomsg=msg)
        call require(ios==0,'Cannot write '//trim(output_prefix)//name//': '//trim(msg))
    end function

    subroutine write_results(f,jet,outer_iterations,coupling_error,deposition_error)
        type(field_t), intent(in) :: f
        type(jet_t), intent(in) :: jet
        integer, intent(in) :: outer_iterations
        real(dp), intent(in) :: coupling_error,deposition_error
        integer :: u,i,j,n
        real(dp) :: delta,ek,etip,current_err,poly_err,flow_err,w,relaxation_time,flight_time,remaining
        real(dp) :: dry_diameter,energy_global,biot,pscale
        logical :: corona_possible
        n=jet%n
        current_err=maxval(abs(jet%lambda*jet%u+jet%sigma*jet%area*jet%e-total_current)) &
            /max(abs(total_current),1e-30_dp)
        poly_err=maxval(abs((mp_in/(mp_in+jet%ms))*(mp_in+jet%ms)-mp_in))/mp_in
        flow_err=maxval(abs(jet%area*jet%u-(mp_in/rho_polymer+jet%ms/rho_solvent)))/Q_pump
        energy_global=0
        do i=2,n
            energy_global=energy_global+(mp_in*cp_polymer+jet%ms(i-1)*cp_solvent) &
                *(jet%temp(i)-jet%temp(i-1))+jet%evap_loss(i)*latent_heat-jet%heat_in(i)
        end do
        pscale=max(abs(V_voltage*total_current),latent_heat*ms_in,1e-30_dp)
        energy_global=abs(energy_global)/pscale
        ! Peek is a cylindrical-air screening estimate, not a discharge solution.
        delta=(pressure_air/101325.0_dp)*298.15_dp/(T_air+273.15_dp)
        ek=3.1e6_dp*delta*(1+0.0308_dp/sqrt(delta*R_needle))
        etip=abs(jet%e(1)); corona_possible=etip>ek
        flight_time=0
        do i=2,n
            flight_time=flight_time+(jet%z(i)-jet%z(i-1))*0.5_dp*(1/jet%u(i)+1/jet%u(i-1))
        end do
        dry_diameter=2*sqrt(mp_in/(rho_polymer*pi*jet%u(n)))
        remaining=jet%ms(n)/ms_in
        u=output_file('jet_profile.csv')
        write(u,'(A)') 'z_m,r_m,diameter_m,V_m_per_s,T_C,eta_Pa_s,polymer_mass_fraction,'// &
            'solvent_mass_flow_kg_s,area_m2,lambda_C_m,q_s_C_m2,E_axial_V_m,'// &
            'I_conv_A,I_cond_A,I_total_A,conductivity_S_m,charge_relaxation_s,'// &
            'evaporated_segment_kg_s,segment_heat_W'
        do i=1,n
            w=mp_in/(mp_in+jet%ms(i)); relaxation_time=huge(1.0_dp)
            if(jet%sigma(i)>0) relaxation_time=eps0*Epsilon_r/jet%sigma(i)
            write(u,'(ES24.16E3,18(",",ES24.16E3))') jet%z(i),jet%a(i),2*jet%a(i),jet%u(i), &
                jet%temp(i)-273.15_dp,jet%eta(i),w,jet%ms(i),jet%area(i),jet%lambda(i), &
                jet%lambda(i)/(2*pi*jet%a(i)),jet%e(i),jet%lambda(i)*jet%u(i), &
                jet%sigma(i)*jet%area(i)*jet%e(i),total_current,jet%sigma(i),relaxation_time, &
                jet%evap_loss(i),jet%heat_in(i)
        end do
        close(u)
        u=output_file('potential_rz.csv')
        write(u,'(A)') 'r_m,z_m,Phi_V,E_r_V_m,E_z_V_m,rho_C_m3,cell_volume_m3'
        do j=1,f%nz
            do i=1,f%nr
                write(u,'(ES24.16E3,6(",",ES24.16E3))') f%r(i),f%z(j),f%phi(i,j), &
                    f%er(i,j),f%ez(i,j),f%rho(i,j),f%vol(i,j)
            end do
        end do
        close(u)
        u=output_file('corona_diagnostic.csv')
        write(u,'(A)') 'voltage_V,entry_E_abs_V_m,Peek_estimate_V_m,relative_air_density,threshold_exceeded'
        write(u,'(ES24.16E3,3(",",ES24.16E3),",",I1)') V_voltage,etip,ek,delta,merge(1,0,corona_possible)
        close(u)
        u=output_file('diagnostics.csv')
        write(u,'(A)') 'metric,value'
        call metric('poisson_relative_residual',f%last_residual)
        call metric('coupling_relative_error',coupling_error)
        call metric('momentum_relative_residual',jet%momentum_residual)
        call metric('mass_segment_relative_residual',jet%mass_residual)
        call metric('energy_segment_relative_residual',jet%energy_residual)
        call metric('thermal_global_relative_residual',energy_global)
        call metric('polymer_relative_error',poly_err)
        call metric('volume_flow_relative_error',flow_err)
        call metric('current_relative_error',current_err)
        call metric('charge_deposition_relative_error',deposition_error)
        call metric('jet_diameter_end_m',2*jet%a(n))
        call metric('dry_diameter_same_speed_m',dry_diameter)
        call metric('velocity_end_m_s',jet%u(n))
        call metric('temperature_end_C',jet%temp(n)-273.15_dp)
        call metric('solvent_remaining_fraction',remaining)
        call metric('flight_time_s',flight_time)
        call metric('jet_current_A',total_current)
        call metric('jet_electric_power_W',abs(V_voltage*total_current))
        call metric('radial_first_cell_m',f%rf(2))
        call metric('minimum_jet_radius_m',minval(jet%a))
        call metric('minimum_radius_over_first_cell',minval(jet%a)/f%rf(2))
        call metric('entry_E_abs_V_m',etip)
        call metric('Peek_estimate_V_m',ek)
        call metric('outer_iterations',real(outer_iterations,dp))
        ! Thermal diffusion is not solved radially. A large Biot number invalidates lumping.
        biot=maxval(h_conv*jet%a/thermal_conductivity_fluid)
        call metric('Biot_screen',biot)
        close(u)
        u=output_file('run_summary.txt')
        write(u,'(A)') 'CONVERGED: steady axisymmetric-field / slender Newtonian-jet model'
        write(u,'(A)') 'Geometry: emitting disk at z=0, grounded collector at z=L.'
        write(u,'(A)') 'Lateral boundary: '//trim(side_boundary)
        write(u,'(A)') 'Liquid stage excludes Taylor cone, 3-D whipping, viscoelasticity and radial solvent diffusion.'
        write(u,'(A)') 'Optional gas-transient results are written separately to plasma_summary.txt.'
        write(u,'(A,I0,A,I0,A,I0)') 'Mesh: ',Nr,' x ',Nz,'; jet nodes: ',Nj
        write(u,'(A,ES16.8)') 'Final jet diameter [m]: ',2*jet%a(n)
        write(u,'(A,ES16.8)') 'Dry diameter at same terminal velocity [m]: ',dry_diameter
        write(u,'(A,ES16.8)') 'Terminal velocity [m/s]: ',jet%u(n)
        write(u,'(A,ES16.8)') 'Terminal temperature [C]: ',jet%temp(n)-273.15_dp
        write(u,'(A,ES16.8)') 'Fraction of initial solvent remaining: ',remaining
        write(u,'(A,ES16.8)') 'Specified/estimated jet current [A]: ',total_current
        write(u,'(A)') 'Dry diameter assumes later solvent removal at fixed speed; it is not a whipping prediction.'
        if(Jet_current_A<0) write(u,'(A)') 'WARNING: current is an uncalibrated scaling estimate.'
        if(minval(jet%a)<3*f%rf(2)) write(u,'(A)') &
            'WARNING: jet cross-section is poorly resolved in field mesh; refine radial grid.'
        if(maxval(mp_in/(mp_in+jet%ms))>0.65_dp) write(u,'(A)') &
            'WARNING: concentrated jet; Newtonian viscosity/ideal activity need experimental validation.'
        if(any(jet%lambda*sign(1.0_dp,V_voltage)<0) .and. V_voltage/=0) write(u,'(A)') &
            'WARNING: local surface charge reverses sign; check current, conductivity and inlet radius.'
        if(corona_possible) write(u,'(A)') &
            'WARNING: approximate corona threshold exceeded; no gas-discharge feedback is modeled.'
        if(biot>0.1_dp) write(u,'(A)') 'WARNING: radial temperature gradients may be significant.'
        close(u)
        u=output_file('effective_config.nml')
        write(u,nml=spinning)
        close(u)
        u=output_file('run_status.txt')
        write(u,'(A)') 'NUMERICALLY_CONVERGED'
        write(u,'(A)') 'Read run_summary.txt for physical applicability warnings.'
        close(u)
        write(*,'(A)') 'CONVERGED. See run_summary.txt and diagnostics.csv (with the chosen output prefix).'
        write(*,'(A,ES14.6,A)') 'Jet diameter at collector = ',2*jet%a(n)*1e6_dp,' um'
        write(*,'(A,ES14.6,A)') 'Velocity at collector     = ',jet%u(n),' m/s'
        write(*,'(A,ES14.6,A)') 'Temperature at collector  = ',jet%temp(n)-273.15_dp,' C'
        write(*,'(A,ES14.6)') 'Initial solvent remaining = ',remaining
        if(corona_possible) write(*,'(A)') 'WARNING: corona screening threshold exceeded; discharge physics is outside this model.'
        if(minval(jet%a)<3*f%rf(2)) write(*,'(A)') 'WARNING: refine radial mesh to resolve the jet charge distribution.'
    contains
        subroutine metric(name,value)
            character(*), intent(in) :: name
            real(dp), intent(in) :: value
            write(u,'(A,",",ES24.16E3)') name,value
        end subroutine
    end subroutine
end module

! BEGIN AIR PLASMA MODULE
! Transient local-mean-energy air model on a nested near-axis patch.
! Swarm transport data: BOLSIG+/Phelps, artificial dry air 80% N2 + 20% O2.
! Electron energy losses are reconstructed from the swarm equilibrium balance.
! This is not a complete humid-air corona/streamer chemistry; see PLASMA_RU.md.
module air_plasma_solver
    use spinning_output
    implicit none
    real(dp), parameter :: elementary_charge=1.602176634e-19_dp, boltzmann=1.380649e-23_dp
    real(dp), parameter :: electron_mass=9.1093837139e-31_dp
    integer :: p_nr=48,p_nz=96,p_max_steps=100000,p_log_every=20
    real(dp) :: p_radius=2.0e-4_dp,p_length=2.0e-3_dp,p_stretch_r=2,p_stretch_z=2
    real(dp) :: p_end_time=2.0e-10_dp,p_dt_max=5.0e-13_dp
    real(dp) :: p_background_ne=1e8_dp,p_seed_ne=1e17_dp,p_seed_z=1.5e-4_dp
    real(dp) :: p_seed_r_width=3e-5_dp,p_seed_z_width=5e-5_dp
    real(dp) :: p_positive_mobility=2e-4_dp,p_negative_mobility=2e-4_dp
    real(dp) :: p_beta_ei=2e-13_dp,p_beta_ii=1e-13_dp,p_energy_transport_factor=5.0_dp/3
    character(512) :: p_transport_file='air_transport.csv'
    namelist /air_plasma/ p_nr,p_nz,p_max_steps,p_log_every,p_radius,p_length,p_stretch_r,p_stretch_z, &
        p_end_time,p_dt_max,p_background_ne,p_seed_ne,p_seed_z,p_seed_r_width,p_seed_z_width, &
        p_positive_mobility,p_negative_mobility,p_beta_ei,p_beta_ii,p_energy_transport_factor,p_transport_file

    type :: swarm_t
        integer :: n
        real(dp) :: N_gas,thermal_energy
        real(dp),allocatable :: en(:),energy(:),muN(:),DN(:),ki(:),ka2(:),ka3(:),loss(:)
    end type
    type :: plasma_t
        type(field_t) :: f
        logical,allocatable :: gas(:,:)
        real(dp),allocatable :: ne(:,:),positive(:,:),negative(:,:),energy(:,:)
        real(dp),allocatable :: mu(:,:),diff(:,:),ion(:,:),attach(:,:),jet_rho(:,:),base_phi(:,:)
        real(dp),allocatable :: outer_phi(:),bottom_phi(:)
    end type
contains
    subroutine read_plasma_config(path)
        character(*),intent(in) :: path
        integer :: u,ios
        character(1024) :: msg
        real(dp) :: values(21)
        if(len_trim(path)>0) then
            open(newunit=u,file=trim(path),status='old',action='read',iostat=ios,iomsg=msg)
            call require(ios==0,'Cannot read plasma config: '//trim(path)//' '//trim(msg))
            read(u,nml=air_plasma,iostat=ios,iomsg=msg)
            close(u)
            call require(ios==0,'Invalid &air_plasma namelist: '//trim(msg))
        end if
        values=[p_radius,p_length,p_stretch_r,p_stretch_z,p_end_time,p_dt_max,p_background_ne, &
            p_seed_ne,p_seed_z,p_seed_r_width,p_seed_z_width,p_positive_mobility,p_negative_mobility, &
            p_beta_ei,p_beta_ii,p_energy_transport_factor,real(p_nr,dp),real(p_nz,dp), &
            real(p_max_steps,dp),real(p_log_every,dp),pressure_air]
        call require(all(ieee_is_finite(values)),'Plasma inputs must be finite.')
        call require(p_nr>=8.and.p_nz>=8.and.p_nr<=2048.and.p_nz<=4096,'Invalid plasma mesh.')
        call require(p_radius>jet_radius_in.and.p_radius<0.98_dp*R_needle, &
            'Near-axis plasma patch requires jet_radius_in < p_radius < 0.98 R_needle (excludes disk edge).')
        call require(p_length>0.and.p_length<Z_max,'Plasma patch must end before the collector.')
        call require(min(p_stretch_r,p_stretch_z)>=0.and.max(p_stretch_r,p_stretch_z)<=8, &
            'Plasma grid stretches must lie in [0,8].')
        call require(p_end_time>0.and.p_dt_max>0.and.p_max_steps>0.and.p_log_every>0, &
            'Plasma time, step and output limits must be positive.')
        call require(min(p_background_ne,p_seed_ne)>=0,'Negative plasma seed density.')
        call require(p_seed_r_width>0.and.p_seed_z_width>0.and.p_seed_z>=0.and.p_seed_z<=p_length, &
            'Invalid plasma seed position/width.')
        call require(min(p_positive_mobility,p_negative_mobility,p_beta_ei,p_beta_ii)>=0, &
            'Ion mobility and recombination coefficients must be nonnegative.')
        call require(p_energy_transport_factor>=1.and.p_energy_transport_factor<=3, &
            'Invalid electron energy transport factor.')
        call require(abs(pressure_air/101325.0_dp-1)<0.01_dp, &
            'Bundled air swarm closure is restricted to approximately one atmosphere.')
        call require(T_air+273.15_dp>=280.and.T_air+273.15_dp<=320, &
            'Bundled dry-air transport data approximate room-temperature gas (280..320 K).')
    end subroutine

    subroutine read_swarm(s,path,number_density,gas_temperature)
        type(swarm_t),intent(out) :: s
        character(*),intent(in) :: path
        real(dp),intent(in) :: number_density,gas_temperature
        integer :: u,ios,n,k
        character(2048) :: line
        real(dp) :: row(7),growth,elastic
        open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
        call require(ios==0,'Missing air_transport.csv; use the complete archive or set p_transport_file.')
        read(u,'(A)',iostat=ios) line
        call require(ios==0.and.index(line,'EN_Td')>0,'Invalid swarm CSV header.')
        n=0
        do
            read(u,'(A)',iostat=ios) line
            if(ios<0) exit
            call require(ios==0,'Swarm table read failure.')
            if(len_trim(line)>0) n=n+1
        end do
        call require(n>=3,'Swarm table has too few rows.')
        rewind(u); read(u,'(A)') line
        s%n=n+2; s%N_gas=number_density; s%thermal_energy=1.5_dp*boltzmann*gas_temperature/elementary_charge
        allocate(s%en(n+2),s%energy(n+2),s%muN(n+2),s%DN(n+2),s%ki(n+2), &
            s%ka2(n+2),s%ka3(n+2),s%loss(n+2))
        do k=3,n+2
            read(u,*,iostat=ios) row
            call require(ios==0.and.all(ieee_is_finite(row)),'Invalid swarm coefficient row.')
            call require(all(row>=0).and.minval(row(1:4))>0,'Invalid swarm coefficient sign.')
            s%en(k)=row(1); s%energy(k)=row(2); s%muN(k)=row(3); s%DN(k)=row(4)
            s%ki(k)=row(5); s%ka2(k)=row(6); s%ka3(k)=row(7)
        end do
        close(u)
        call require(all(s%en(4:)>s%en(3:s%n-1)).and.all(s%energy(4:)>s%energy(3:s%n-1)), &
            'Swarm field and energy coordinates must increase strictly.')
        call require(s%thermal_energy<s%energy(3),'Swarm low-energy join is invalid.')
        ! Add a thermal low-field limit. The original published data begin at 1 Td.
        s%en(1:2)=0; s%energy(1)=0; s%energy(2)=s%thermal_energy
        s%muN(1:2)=s%muN(3); s%DN(1)=0
        s%DN(2)=s%muN(3)*boltzmann*gas_temperature/elementary_charge
        s%ki(1:2)=0; s%ka2(1:2)=s%ka2(3); s%ka3(1:2)=s%ka3(3)
        do k=3,s%n
            growth=s%ki(k)-s%ka2(k)-number_density*s%ka3(k)
            ! Equilibrium of mean energy, including electron production/loss:
            ! L(eps) = muN*(E/N)^2 - eps*(ki-ka2-N*ka3).
            s%loss(k)=s%muN(k)*(s%en(k)*1e-21_dp)**2-s%energy(k)*growth
        end do
        elastic=s%muN(3)*(s%en(3)*1e-21_dp)**2/(s%energy(3)-s%thermal_energy)
        s%loss(1)=-elastic*s%thermal_energy
        s%loss(2)=s%thermal_energy*(s%ka2(2)+number_density*s%ka3(2))
        call require(all(s%loss(2:)>s%loss(:s%n-1)), &
            'Reconstructed energy loss is not monotone; this table needs another energy solver.')
    end subroutine

    subroutine lookup_swarm(s,energy,mu,diff,ion,attach)
        type(swarm_t),intent(in) :: s
        real(dp),intent(in) :: energy
        real(dp),intent(out) :: mu,diff,ion,attach
        real(dp) :: w,e
        integer :: lo,hi,k
        e=max(0.0_dp,min(energy,s%energy(s%n)))
        lo=1; hi=s%n
        do while(hi-lo>1)
            k=(hi+lo)/2
            if(s%energy(k)>e) then
                hi=k
            else
                lo=k
            end if
        end do
        w=(e-s%energy(lo))/(s%energy(hi)-s%energy(lo))
        mu=((1-w)*s%muN(lo)+w*s%muN(hi))/s%N_gas
        diff=((1-w)*s%DN(lo)+w*s%DN(hi))/s%N_gas
        ion=s%N_gas*((1-w)*s%ki(lo)+w*s%ki(hi))
        attach=s%N_gas*((1-w)*s%ka2(lo)+w*s%ka2(hi))+ &
            s%N_gas**2*((1-w)*s%ka3(lo)+w*s%ka3(hi))
    end subroutine

    real(dp) function swarm_equilibrium(s,en) result(energy)
        type(swarm_t),intent(in) :: s
        real(dp),intent(in) :: en
        if(en<=s%en(3)) then
            energy=s%thermal_energy+(s%energy(3)-s%thermal_energy)*max(en,0.0_dp)/s%en(3)
        else
            energy=interp(s%en(3:),s%energy(3:),en)
        end if
    end function

    real(dp) function bernoulli(x) result(v)
        real(dp),intent(in) :: x
        if(abs(x)<1e-4_dp) then
            v=1-x/2+x*x/12-x**4/720
        else if(x>50) then
            v=x*exp(-x)
        else if(x< -50) then
            v=-x
        else
            v=x/(exp(x)-1)
        end if
    end function

    subroutine sg_face(metric,delta_phi,charge_sign,mu_left,mu_right,d_left,d_right,c_left,c_right)
        real(dp),intent(in) :: metric,delta_phi,charge_sign,mu_left,mu_right,d_left,d_right
        real(dp),intent(out) :: c_left,c_right
        real(dp) :: d,pe,v_metric
        d=(d_left+d_right)/2
        if(d>tiny(1.0_dp)) then
            pe=-charge_sign*(mu_left+mu_right)/2*delta_phi/d
            ! Gradient of D*n includes thermal diffusion for variable electron D.
            c_left=metric*d_left*bernoulli(-pe); c_right=metric*d_right*bernoulli(pe)
        else
            v_metric=-metric*charge_sign*(mu_left+mu_right)/2*delta_phi
            c_left=max(v_metric,0.0_dp); c_right=max(-v_metric,0.0_dp)
        end if
    end subroutine

    real(dp) function reference_phi(f,r,z) result(v)
        type(field_t),intent(in) :: f
        real(dp),intent(in) :: r,z
        real(dp) :: line(f%nz)
        integer :: j
        do j=1,f%nz
            line(j)=interp(f%r,f%phi(:,j),r)
        end do
        if(z<f%z(1)) then
            v=f%voltage+(line(1)-f%voltage)*z/f%z(1)
        else
            v=interp(f%z,line,z)
        end if
    end function

    subroutine init_plasma(p,s,reference,jet)
        type(plasma_t),intent(out) :: p
        type(swarm_t),intent(in) :: s
        type(field_t),intent(in) :: reference
        type(jet_t),intent(in) :: jet
        integer :: i,j,n,m
        real(dp) :: dummy,a,seed,en
        n=p_nr; m=p_nz
        call init_field(p%f,n,m,p_radius,p_length,p_stretch_r,p_stretch_z,V_voltage,p_radius,.false.,.true.)
        allocate(p%gas(n,m),p%ne(n,m),p%positive(n,m),p%negative(n,m),p%energy(n,m), &
            p%mu(n,m),p%diff(n,m),p%ion(n,m),p%attach(n,m),p%jet_rho(n,m),p%base_phi(n,m), &
            p%outer_phi(m),p%bottom_phi(n))
        do j=1,m
            p%outer_phi(j)=reference_phi(reference,p_radius,p%f%z(j))
            p%f%bc(n,j)=p%f%bc(n,j)+2*pi*p_radius*(p%f%zf(j+1)-p%f%zf(j)) &
                /(p_radius-p%f%r(n))*p%outer_phi(j)
            a=interp(jet%z,jet%a,p%f%z(j))
            do i=1,n
                p%f%phi(i,j)=reference_phi(reference,p%f%r(i),p%f%z(j))
                p%gas(i,j)=p%f%r(i)>a
                seed=p_seed_ne*exp(-(p%f%r(i)/p_seed_r_width)**2 &
                    -((p%f%z(j)-p_seed_z)/p_seed_z_width)**2)
                p%ne(i,j)=merge(p_background_ne+seed,0.0_dp,p%gas(i,j))
            end do
        end do
        do i=1,n
            p%bottom_phi(i)=reference_phi(reference,p%f%r(i),p_length)
            p%f%bc(i,m)=p%f%bc(i,m)+p%f%area(i)/(p_length-p%f%z(m))*p%bottom_phi(i)
        end do
        call require(count(p%gas)>0,'Plasma patch contains no gas cells.')
        call deposit_charge(p%f,jet%z,jet%a,jet%lambda,dummy)
        if(.not.space_charge) p%f%rho=0
        p%jet_rho=p%f%rho; p%positive=p%ne; p%negative=0
        call solve_field(p%f,poisson_tol,max_pcg)
        call patch_gradients(p)
        p%base_phi=p%f%phi
        do j=1,m
            do i=1,n
                en=hypot(p%f%er(i,j),p%f%ez(i,j))/s%N_gas/1e-21_dp
                call require(en<=s%en(s%n),'Local E/N exceeds swarm table range; change geometry/data, not a coefficient cap.')
                p%energy(i,j)=p%ne(i,j)*swarm_equilibrium(s,en)
            end do
        end do
    end subroutine

    subroutine patch_gradients(p)
        type(plasma_t),intent(inout) :: p
        integer :: i,j,n,m
        n=p%f%nr; m=p%f%nz
        call field_gradients(p%f,p_radius)
        do j=1,m
            p%f%er(n,j)=-slope3(p%f%r(n-1),p%f%r(n),p_radius, &
                p%f%phi(n-1,j),p%f%phi(n,j),p%outer_phi(j))
        end do
        do i=1,n
            p%f%ez(i,m)=-slope3(p%f%z(m-1),p%f%z(m),p_length, &
                p%f%phi(i,m-1),p%f%phi(i,m),p%bottom_phi(i))
        end do
    end subroutine

    subroutine update_plasma_field(p,s)
        type(plasma_t),intent(inout) :: p
        type(swarm_t),intent(in) :: s
        real(dp) :: en
        p%f%rho=p%jet_rho+elementary_charge*(p%positive-p%ne-p%negative)
        call solve_field(p%f,poisson_tol,max_pcg); call patch_gradients(p)
        en=maxval(hypot(p%f%er,p%f%ez),mask=p%gas)/s%N_gas/1e-21_dp
        call require(en<=s%en(s%n),'E/N exceeds supplied swarm data during plasma evolution.')
        call require(maxval(p%ne)<1e-3_dp*s%N_gas,'Weak-ionization assumption has failed.')
    end subroutine

    subroutine plasma_coefficients(p,s)
        type(plasma_t),intent(inout) :: p
        type(swarm_t),intent(in) :: s
        real(dp) :: energy
        integer :: i,j
        do j=1,p%f%nz
            do i=1,p%f%nr
                energy=s%thermal_energy
                if(p%ne(i,j)>1e-100_dp) energy=p%energy(i,j)/p%ne(i,j)
                call require(energy>=0.and.energy<=s%energy(s%n)*(1+1e-10_dp), &
                    'Electron mean energy outside transport data; reduce step or supply a wider table.')
                call lookup_swarm(s,energy,p%mu(i,j),p%diff(i,j),p%ion(i,j),p%attach(i,j))
            end do
        end do
    end subroutine

    subroutine transport(f,gas,outer_phi,bottom_phi,x,mu,diff,charge_sign,dt,absorbing,wall_loss,balance_error,wall_speed)
        type(field_t),intent(in) :: f
        logical,intent(in) :: gas(:,:)
        real(dp),intent(in) :: outer_phi(:),bottom_phi(:)
        real(dp),intent(inout) :: x(:,:)
        real(dp),intent(in) :: mu(:,:),diff(:,:),charge_sign,dt
        logical,intent(in) :: absorbing
        real(dp),intent(out) :: wall_loss,balance_error
        real(dp),intent(in),optional :: wall_speed(:,:)
        real(dp) :: lo(max(f%nr,f%nz)),di(max(f%nr,f%nz)),up(max(f%nr,f%nz))
        real(dp) :: rhs(max(f%nr,f%nz)),answer(max(f%nr,f%nz))
        real(dp) :: cl(f%nr,f%nz),cr(f%nr,f%nz),left,right,metric,before,after,loss,wall_flux
        integer :: i,j,n,m
        n=f%nr; m=f%nz
        before=sum(x*f%vol); wall_loss=0
        ! Backward Euler in each direction, Scharfetter-Gummel flux, positive M matrices.
        do j=1,m
            lo(:n)=0; up(:n)=0; di(:n)=1; rhs(:n)=x(:,j)
            do i=1,n-1
                call sg_face(f%gr(i,j),f%phi(i+1,j)-f%phi(i,j),charge_sign, &
                    mu(i,j),mu(i+1,j),diff(i,j),diff(i+1,j),cl(i,j),cr(i,j))
                if(present(wall_speed)) then
                    if(gas(i,j).and..not.gas(i+1,j)) then
                        wall_flux=2*pi*f%rf(i+1)*(f%zf(j+1)-f%zf(j))*wall_speed(i,j)
                        cl(i,j)=cl(i,j)*wall_flux/max(cr(i,j)+wall_flux,tiny(1.0_dp)); cr(i,j)=0
                    else if(.not.gas(i,j).and.gas(i+1,j)) then
                        wall_flux=2*pi*f%rf(i+1)*(f%zf(j+1)-f%zf(j))*wall_speed(i+1,j)
                        cr(i,j)=cr(i,j)*wall_flux/max(cl(i,j)+wall_flux,tiny(1.0_dp)); cl(i,j)=0
                    end if
                end if
                left=dt*cl(i,j); right=dt*cr(i,j)
                di(i)=di(i)+left/f%vol(i,j); up(i)=up(i)-right/f%vol(i,j)
                di(i+1)=di(i+1)+right/f%vol(i+1,j); lo(i+1)=lo(i+1)-left/f%vol(i+1,j)
            end do
            cl(n,j)=0; cr(n,j)=0
            if(absorbing) then
                metric=2*pi*p_radius*(f%zf(j+1)-f%zf(j))/(p_radius-f%r(n))
                call sg_face(metric,outer_phi(j)-f%phi(n,j),charge_sign, &
                    mu(n,j),mu(n,j),diff(n,j),diff(n,j),cl(n,j),cr(n,j))
                if(present(wall_speed)) then
                    wall_flux=2*pi*p_radius*(f%zf(j+1)-f%zf(j))*wall_speed(n,j)
                    cl(n,j)=cl(n,j)*wall_flux/max(cr(n,j)+wall_flux,tiny(1.0_dp)); cr(n,j)=0
                end if
                di(n)=di(n)+dt*cl(n,j)/f%vol(n,j)
            end if
            do i=1,n
                if(.not.gas(i,j)) then
                    lo(i)=0; up(i)=0; di(i)=1; rhs(i)=0
                end if
            end do
            call tridiagonal(lo(:n),di(:n),up(:n),rhs(:n),answer(:n)); x(:,j)=answer(:n)
            loss=cl(n,j)*x(n,j)
            do i=1,n-1
                if(gas(i,j).neqv.gas(i+1,j)) then
                    loss=loss+abs(cl(i,j)*x(i,j)-cr(i,j)*x(i+1,j))
                end if
            end do
            wall_loss=wall_loss+dt*loss
        end do
        do i=1,n
            lo(:m)=0; up(:m)=0; di(:m)=1; rhs(:m)=x(i,:)
            do j=1,m-1
                call sg_face(f%gz(i,j),f%phi(i,j+1)-f%phi(i,j),charge_sign, &
                    mu(i,j),mu(i,j+1),diff(i,j),diff(i,j+1),cl(i,j),cr(i,j))
                if(present(wall_speed)) then
                    if(gas(i,j).and..not.gas(i,j+1)) then
                        wall_flux=f%area(i)*wall_speed(i,j)
                        cl(i,j)=cl(i,j)*wall_flux/max(cr(i,j)+wall_flux,tiny(1.0_dp)); cr(i,j)=0
                    else if(.not.gas(i,j).and.gas(i,j+1)) then
                        wall_flux=f%area(i)*wall_speed(i,j+1)
                        cr(i,j)=cr(i,j)*wall_flux/max(cl(i,j)+wall_flux,tiny(1.0_dp)); cl(i,j)=0
                    end if
                end if
                left=dt*cl(i,j); right=dt*cr(i,j)
                di(j)=di(j)+left/f%vol(i,j); up(j)=up(j)-right/f%vol(i,j)
                di(j+1)=di(j+1)+right/f%vol(i,j+1); lo(j+1)=lo(j+1)-left/f%vol(i,j+1)
            end do
            left=0; right=0
            if(absorbing) then
                metric=f%area(i)/f%z(1)
                call sg_face(metric,f%phi(i,1)-V_voltage,charge_sign, &
                    mu(i,1),mu(i,1),diff(i,1),diff(i,1),left,right)
                if(present(wall_speed)) then
                    wall_flux=f%area(i)*wall_speed(i,1)
                    right=right*wall_flux/max(left+wall_flux,tiny(1.0_dp)); left=0
                end if
                di(1)=di(1)+dt*right/f%vol(i,1)
            end if
            cl(i,m)=0; cr(i,m)=0
            if(absorbing) then
                metric=f%area(i)/(p_length-f%z(m))
                call sg_face(metric,bottom_phi(i)-f%phi(i,m),charge_sign, &
                    mu(i,m),mu(i,m),diff(i,m),diff(i,m),cl(i,m),cr(i,m))
                if(present(wall_speed)) then
                    wall_flux=f%area(i)*wall_speed(i,m)
                    cl(i,m)=cl(i,m)*wall_flux/max(cr(i,m)+wall_flux,tiny(1.0_dp)); cr(i,m)=0
                end if
                di(m)=di(m)+dt*cl(i,m)/f%vol(i,m)
            end if
            do j=1,m
                if(.not.gas(i,j)) then
                    lo(j)=0; up(j)=0; di(j)=1; rhs(j)=0
                end if
            end do
            call tridiagonal(lo(:m),di(:m),up(:m),rhs(:m),answer(:m)); x(i,:)=answer(:m)
            loss=right*x(i,1)+cl(i,m)*x(i,m)
            do j=1,m-1
                if(gas(i,j).neqv.gas(i,j+1)) then
                    loss=loss+abs(cl(i,j)*x(i,j)-cr(i,j)*x(i,j+1))
                end if
            end do
            wall_loss=wall_loss+dt*loss
        end do
        call require(all(x>=0).and.all(ieee_is_finite(x)),'Invalid density in implicit plasma transport.')
        after=sum(x*f%vol)
        balance_error=abs(after-before+wall_loss)/max(before,after,wall_loss,1e-100_dp)
        call require(balance_error<1e-8_dp,'Plasma finite-volume particle flux balance failed.')
    end subroutine

    subroutine joule_power(p,power,wall_speed)
        type(plasma_t),intent(in) :: p
        real(dp),intent(out) :: power(:,:)
        real(dp),intent(in) :: wall_speed(:,:)
        real(dp) :: left,right,flux,work,metric,delta,wall_flux
        integer :: i,j,n,m
        n=p%f%nr; m=p%f%nz; power=0
        do j=1,m
            do i=1,n-1
                delta=p%f%phi(i+1,j)-p%f%phi(i,j)
                call sg_face(p%f%gr(i,j),delta,-1.0_dp,p%mu(i,j),p%mu(i+1,j), &
                    p%diff(i,j),p%diff(i+1,j),left,right)
                if(p%gas(i,j).and..not.p%gas(i+1,j)) then
                    wall_flux=2*pi*p%f%rf(i+1)*(p%f%zf(j+1)-p%f%zf(j))*wall_speed(i,j)
                    left=left*wall_flux/max(right+wall_flux,tiny(1.0_dp)); right=0
                else if(.not.p%gas(i,j).and.p%gas(i+1,j)) then
                    wall_flux=2*pi*p%f%rf(i+1)*(p%f%zf(j+1)-p%f%zf(j))*wall_speed(i+1,j)
                    right=right*wall_flux/max(left+wall_flux,tiny(1.0_dp)); left=0
                end if
                flux=left*p%ne(i,j)-right*p%ne(i+1,j); work=flux*delta
                if(p%gas(i,j).and.p%gas(i+1,j)) then
                    power(i,j)=power(i,j)+work/2; power(i+1,j)=power(i+1,j)+work/2
                else if(p%gas(i,j)) then
                    power(i,j)=power(i,j)+work
                else if(p%gas(i+1,j)) then
                    power(i+1,j)=power(i+1,j)+work
                end if
            end do
            delta=p%outer_phi(j)-p%f%phi(n,j)
            metric=2*pi*p_radius*(p%f%zf(j+1)-p%f%zf(j))/(p_radius-p%f%r(n))
            call sg_face(metric,delta,-1.0_dp,p%mu(n,j),p%mu(n,j),p%diff(n,j),p%diff(n,j),left,right)
            wall_flux=2*pi*p_radius*(p%f%zf(j+1)-p%f%zf(j))*wall_speed(n,j)
            left=left*wall_flux/max(right+wall_flux,tiny(1.0_dp))
            if(p%gas(n,j)) power(n,j)=power(n,j)+left*p%ne(n,j)*delta
        end do
        do i=1,n
            do j=1,m-1
                delta=p%f%phi(i,j+1)-p%f%phi(i,j)
                call sg_face(p%f%gz(i,j),delta,-1.0_dp,p%mu(i,j),p%mu(i,j+1), &
                    p%diff(i,j),p%diff(i,j+1),left,right)
                if(p%gas(i,j).and..not.p%gas(i,j+1)) then
                    wall_flux=p%f%area(i)*wall_speed(i,j)
                    left=left*wall_flux/max(right+wall_flux,tiny(1.0_dp)); right=0
                else if(.not.p%gas(i,j).and.p%gas(i,j+1)) then
                    wall_flux=p%f%area(i)*wall_speed(i,j+1)
                    right=right*wall_flux/max(left+wall_flux,tiny(1.0_dp)); left=0
                end if
                work=(left*p%ne(i,j)-right*p%ne(i,j+1))*delta
                if(p%gas(i,j).and.p%gas(i,j+1)) then
                    power(i,j)=power(i,j)+work/2; power(i,j+1)=power(i,j+1)+work/2
                else if(p%gas(i,j)) then
                    power(i,j)=power(i,j)+work
                else if(p%gas(i,j+1)) then
                    power(i,j+1)=power(i,j+1)+work
                end if
            end do
            delta=p%f%phi(i,1)-V_voltage; metric=p%f%area(i)/p%f%z(1)
            call sg_face(metric,delta,-1.0_dp,p%mu(i,1),p%mu(i,1),p%diff(i,1),p%diff(i,1),left,right)
            wall_flux=p%f%area(i)*wall_speed(i,1)
            right=right*wall_flux/max(left+wall_flux,tiny(1.0_dp))
            if(p%gas(i,1)) power(i,1)=power(i,1)-right*p%ne(i,1)*delta
            delta=p%bottom_phi(i)-p%f%phi(i,m); metric=p%f%area(i)/(p_length-p%f%z(m))
            call sg_face(metric,delta,-1.0_dp,p%mu(i,m),p%mu(i,m),p%diff(i,m),p%diff(i,m),left,right)
            wall_flux=p%f%area(i)*wall_speed(i,m)
            left=left*wall_flux/max(right+wall_flux,tiny(1.0_dp))
            if(p%gas(i,m)) power(i,m)=power(i,m)+left*p%ne(i,m)*delta
        end do
        power=power/p%f%vol
    end subroutine

    real(dp) function pair_loss(a,b,beta_dt) result(loss)
        real(dp),intent(in) :: a,b,beta_dt
        real(dp) :: discriminant
        discriminant=1+2*beta_dt*(a+b)+(beta_dt*(a-b))**2
        loss=2*beta_dt*a*b/(1+beta_dt*(a+b)+sqrt(discriminant))
    end function

    subroutine chemical_step(ne,positive,negative,ion,attach,beta_ei,beta_ii,dt,energy_factor)
        real(dp),intent(inout) :: ne,positive,negative
        real(dp),intent(in) :: ion,attach,beta_ei,beta_ii,dt
        real(dp),intent(out) :: energy_factor
        real(dp) :: loss,pre_recomb
        ne=ne/(1-dt*(ion-attach))
        positive=positive+dt*ion*ne; negative=negative+dt*attach*ne
        pre_recomb=ne
        loss=pair_loss(ne,positive,dt*beta_ei); ne=ne-loss; positive=positive-loss
        loss=pair_loss(negative,positive,dt*beta_ii); negative=negative-loss; positive=positive-loss
        energy_factor=1
        if(pre_recomb>0) energy_factor=ne/pre_recomb
    end subroutine

    subroutine energy_root(s,ne,rhs,dt,energy,loss,ok)
        type(swarm_t),intent(in) :: s
        real(dp),intent(in) :: ne,rhs,dt
        real(dp),intent(out) :: energy,loss
        logical,intent(out) :: ok
        real(dp) :: target,low,high,w
        integer :: lo,hi,k
        ok=.true.; energy=0; loss=0
        if(ne<=0) then
            ok=abs(rhs)<1e-100_dp; return
        end if
        target=rhs/ne
        lo=1; hi=s%n
        low=s%energy(lo)+dt*s%N_gas*s%loss(lo)
        high=s%energy(hi)+dt*s%N_gas*s%loss(hi)
        if(target<low.or.target>high) then
            ok=.false.; return
        end if
        do while(hi-lo>1)
            k=(lo+hi)/2
            if(s%energy(k)+dt*s%N_gas*s%loss(k)>target) then
                hi=k
            else
                lo=k
            end if
        end do
        low=s%energy(lo)+dt*s%N_gas*s%loss(lo)
        high=s%energy(hi)+dt*s%N_gas*s%loss(hi)
        w=(target-low)/(high-low)
        energy=ne*((1-w)*s%energy(lo)+w*s%energy(hi))
        loss=dt*ne*s%N_gas*((1-w)*s%loss(lo)+w*s%loss(hi))
    end subroutine

    subroutine plasma_step(p,s,dt,wall_charge,wall_energy,collision_energy,recomb_energy,joule_energy,ok,max_flux_error)
        type(plasma_t),intent(inout) :: p
        type(swarm_t),intent(in) :: s
        real(dp),intent(in) :: dt
        real(dp),intent(out) :: wall_charge,wall_energy,collision_energy,recomb_energy,joule_energy,max_flux_error
        logical,intent(out) :: ok
        real(dp) :: muion(p%f%nr,p%f%nz),dion(p%f%nr,p%f%nz),power(p%f%nr,p%f%nz)
        real(dp) :: wall_speed(p%f%nr,p%f%nz)
        real(dp) :: le,lp,lm,err,fac,loss,rhs,new_energy,gas_temperature
        integer :: i,j
        logical :: cell_ok
        gas_temperature=T_air+273.15_dp
        wall_speed=sqrt(elementary_charge*(2.0_dp/3)* &
            max(p%energy/max(p%ne,tiny(1.0_dp)),0.0_dp)/(2*pi*electron_mass))
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,p%ne,p%mu,p%diff,-1.0_dp,dt,.true.,le,err,wall_speed); max_flux_error=err
        muion=p_positive_mobility; dion=p_positive_mobility*boltzmann*gas_temperature/elementary_charge
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,p%positive,muion,dion,1.0_dp,dt,.true.,lp,err); max_flux_error=max(max_flux_error,err)
        muion=p_negative_mobility; dion=p_negative_mobility*boltzmann*gas_temperature/elementary_charge
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,p%negative,muion,dion,-1.0_dp,dt,.true.,lm,err); max_flux_error=max(max_flux_error,err)
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,p%energy,p_energy_transport_factor*p%mu,p_energy_transport_factor*p%diff, &
            -1.0_dp,dt,.true.,wall_energy,err,(4.0_dp/3)*wall_speed)
        max_flux_error=max(max_flux_error,err)
        wall_charge=lp-le-lm; recomb_energy=0; collision_energy=0; ok=.true.
        do j=1,p%f%nz
            do i=1,p%f%nr
                if(.not.p%gas(i,j)) cycle
                call chemical_step(p%ne(i,j),p%positive(i,j),p%negative(i,j),p%ion(i,j),p%attach(i,j), &
                    p_beta_ei,p_beta_ii,dt,fac)
                recomb_energy=recomb_energy+p%energy(i,j)*(1-fac)*p%f%vol(i,j)
                p%energy(i,j)=p%energy(i,j)*fac
            end do
        end do
        call joule_power(p,power,wall_speed)
        joule_energy=dt*sum(power*p%f%vol)
        do j=1,p%f%nz
            do i=1,p%f%nr
                if(.not.p%gas(i,j)) cycle
                rhs=p%energy(i,j)+dt*power(i,j)
                call energy_root(s,p%ne(i,j),rhs,dt,new_energy,loss,cell_ok)
                if(.not.cell_ok) then
                    ok=.false.; return
                end if
                p%energy(i,j)=new_energy
                collision_energy=collision_energy+loss*p%f%vol(i,j)
            end do
        end do
    end subroutine

    subroutine plasma_map(p,s,t,name)
        type(plasma_t),intent(in) :: p
        type(swarm_t),intent(in) :: s
        real(dp),intent(in) :: t
        character(*),intent(in) :: name
        integer :: u,i,j
        real(dp) :: en,mean_energy,te,conductivity
        u=output_file(name)
        write(u,'(A)') 'time_s,r_m,z_m,ne_m3,n_positive_m3,n_negative_m3,energy_eV_m3,'// &
            'mean_energy_eV,Te_eV,Phi_V,E_r_V_m,E_z_V_m,EN_Td,gas_cell,Te_defined,conductivity_S_m'
        do j=1,p%f%nz
            do i=1,p%f%nr
                en=hypot(p%f%er(i,j),p%f%ez(i,j))/s%N_gas/1e-21_dp
                mean_energy=0
                if(p%ne(i,j)>0) mean_energy=p%energy(i,j)/p%ne(i,j)
                te=2*mean_energy/3
                conductivity=elementary_charge*(p%mu(i,j)*p%ne(i,j)+ &
                    p_positive_mobility*p%positive(i,j)+p_negative_mobility*p%negative(i,j))
                write(u,'(ES24.16E3,12(",",ES24.16E3),2(",",I1),",",ES24.16E3)') &
                    t,p%f%r(i),p%f%z(j),p%ne(i,j),p%positive(i,j),p%negative(i,j),p%energy(i,j), &
                    mean_energy,te,p%f%phi(i,j),p%f%er(i,j),p%f%ez(i,j),en, &
                    merge(1,0,p%gas(i,j)),merge(1,0,p%ne(i,j)>0),conductivity
            end do
        end do
        close(u)
    end subroutine

    subroutine run_air_plasma(reference,jet,path)
        type(field_t),intent(in) :: reference
        type(jet_t),intent(in) :: jet
        character(*),intent(in) :: path
        type(plasma_t) :: p
        type(swarm_t) :: s
        integer :: step,retry,u,history,n,m,pcg_sum,retries_total,history_count,k
        real(dp) :: t,dt,number_density,q0,w0,q,w,wallq,wallw,collw,recw,joulew
        real(dp) :: total_wallq,total_wallw,total_collw,total_recw,total_joulew,max_flux_error,flux_error
        real(dp) :: qerr,werr,rate_max,charge_rate,energy_scale,charge_scale,frozen_limit,flight,t0,t1
        real(dp),allocatable :: ne0(:,:),np0(:,:),nm0(:,:),energy0(:,:),mean_energy(:,:)
        real(dp),allocatable :: history_values(:,:)
        logical :: ok
        call read_plasma_config(path)
        flight=sum((jet%z(2:)-jet%z(:jet%n-1))/jet%u(2:))
        frozen_limit=0.01_dp*flight
        if(Conductivity>0) frozen_limit=min(frozen_limit,0.01_dp*eps0*Epsilon_r/Conductivity)
        call require(p_end_time<frozen_limit,'Plasma duration is too long for the frozen-jet approximation.')
        number_density=pressure_air/(boltzmann*(T_air+273.15_dp))
        call read_swarm(s,trim(p_transport_file),number_density,T_air+273.15_dp)
        call init_plasma(p,s,reference,jet)
        call plasma_coefficients(p,s)
        call plasma_map(p,s,0.0_dp,'plasma_initial.csv')
        n=p%f%nr; m=p%f%nz
        allocate(ne0(n,m),np0(n,m),nm0(n,m),energy0(n,m),mean_energy(n,m))
        allocate(history_values(11,4+p_max_steps/p_log_every)); history_count=0
        q0=sum((p%positive-p%ne-p%negative)*p%f%vol); w0=sum(p%energy*p%f%vol)
        total_wallq=0; total_wallw=0; total_collw=0; total_recw=0; total_joulew=0
        max_flux_error=0; retries_total=0; t=0; pcg_sum=0
        call cpu_time(t0)
        write(*,'(A)') 'AIR PLASMA: transient local mean energy, 80% N2/20% O2, one atmosphere; frozen jet.'
        write(*,'(A,I0,A,I0,A,ES12.4)') 'Plasma patch: ',n,' x ',m,'; end time [s]: ',p_end_time
        do step=1,p_max_steps
            call plasma_coefficients(p,s)
            rate_max=maxval(p%ion+p%attach+p_beta_ei*p%positive,mask=p%gas)
            charge_rate=elementary_charge/eps0*maxval(p%mu*p%ne+ &
                p_positive_mobility*p%positive+p_negative_mobility*p%negative,mask=p%gas)
            dt=min(p_dt_max,p_end_time-t,0.15_dp/max(rate_max,1e-100_dp), &
                0.15_dp/max(charge_rate,1e-100_dp))
            ne0=p%ne; np0=p%positive; nm0=p%negative; energy0=p%energy
            do retry=1,25
                call plasma_step(p,s,dt,wallq,wallw,collw,recw,joulew,ok,flux_error)
                if(ok) exit
                p%ne=ne0; p%positive=np0; p%negative=nm0; p%energy=energy0
                dt=dt/2; retries_total=retries_total+1
            end do
            call require(ok,'Electron energy step cannot stay within the available swarm table; inspect plasma setup.')
            call require(dt>1e-20_dp,'Plasma step collapsed; this reduced model is outside its numerical range.')
            t=t+dt
            total_wallq=total_wallq+wallq; total_wallw=total_wallw+wallw
            total_collw=total_collw+collw; total_recw=total_recw+recw; total_joulew=total_joulew+joulew
            max_flux_error=max(max_flux_error,flux_error)
            call update_plasma_field(p,s); pcg_sum=pcg_sum+p%f%last_iterations
            q=sum((p%positive-p%ne-p%negative)*p%f%vol); w=sum(p%energy*p%f%vol)
            charge_scale=max(sum((p%positive+p%ne+p%negative)*p%f%vol),1e-100_dp)
            energy_scale=max(w0,w,abs(total_wallw),abs(total_collw),abs(total_joulew),1e-100_dp)
            qerr=abs(q+total_wallq-q0)/charge_scale
            werr=abs(w-w0+total_wallw+total_collw+total_recw-total_joulew)/energy_scale
            call require(qerr<1e-7_dp.and.werr<1e-7_dp,'Plasma global charge/energy balance failed.')
            if(mod(step,p_log_every)==0.or.step==1.or.t>=p_end_time*(1-1e-13_dp)) then
                mean_energy=0
                where(p%ne>0) mean_energy=p%energy/max(p%ne,tiny(1.0_dp))
                history_count=history_count+1
                history_values(:,history_count)=[t,dt,maxval(p%ne), &
                    2*maxval(mean_energy,mask=p%gas)/3,sum(p%ne*p%f%vol), &
                    sum(p%positive*p%f%vol),sum(p%negative*p%f%vol),q*elementary_charge,qerr,werr, &
                    maxval(abs(p%f%phi-p%base_phi))/max(abs(V_voltage),1.0_dp)]
            end if
            if(t>=p_end_time*(1-1e-13_dp)) exit
        end do
        call require(step<=p_max_steps,'Plasma end time was not reached; increase p_max_steps or inspect step restrictions.')
        history=output_file('plasma_history.csv')
        write(history,'(A)') 'time_s,dt_s,ne_peak_m3,Te_peak_eV,electrons,positive_ions,negative_ions,'// &
            'gas_charge_C,charge_balance_relative,energy_balance_relative,phi_feedback_relative'
        do k=1,history_count
            write(history,'(ES24.16E3,10(",",ES24.16E3))') history_values(:,k)
        end do
        close(history)
        call plasma_coefficients(p,s)
        call plasma_map(p,s,t,'plasma_final.csv')
        u=output_file('effective_plasma_config.nml'); write(u,nml=air_plasma); close(u)
        u=output_file('plasma_summary.txt')
        write(u,'(A)') 'TRANSIENT_COMPLETED: electron/ion transport + dynamic electron energy + self-consistent gas Poisson.'
        write(u,'(A)') 'Gas approximation: 80% N2 + 20% O2, room temperature, one atmosphere; no water/solvent vapor chemistry.'
        write(u,'(A)') 'Short-time nested near-axis patch; prescribed outer potential; liquid velocity and charge are frozen.'
        write(u,'(A)') 'This is not a claim of steady corona or spontaneous inception; initial seed/background are prescribed.'
        write(u,'(A)') 'No photoionization, detachment, secondary emission, detailed ion conversion, or full external circuit.'
        write(u,'(A)') 'Mean energy loss closure is reconstructed from published BOLSIG+/Phelps swarm balance.'
        write(u,'(A)') 'Electron energy flux uses the 5/3 Maxwellian factor by default; the EEDF is not solved in space/time.'
        write(u,'(A)') 'Electron absorbing walls use thermal Robin flux; energy per lost particle is 4/3 of local mean energy.'
        write(u,'(A,ES24.16E3)') 'end_time_s=',t
        write(u,'(A,I0)') 'time_steps=',step
        write(u,'(A,I0)') 'energy_step_retries=',retries_total
        write(u,'(A,ES24.16E3)') 'electron_peak_m3=',maxval(p%ne)
        write(u,'(A,ES24.16E3)') 'effective_Te_peak_eV=',2*maxval(mean_energy,mask=p%gas)/3
        write(u,'(A,ES24.16E3)') 'charge_balance_relative=',qerr
        write(u,'(A,ES24.16E3)') 'energy_balance_relative=',werr
        write(u,'(A,ES24.16E3)') 'max_transport_flux_balance_relative=',max_flux_error
        write(u,'(A,ES24.16E3)') 'poisson_relative_residual=',p%f%last_residual
        write(u,'(A,ES24.16E3)') 'potential_feedback_relative=', &
            maxval(abs(p%f%phi-p%base_phi))/max(abs(V_voltage),1.0_dp)
        write(u,'(A,ES24.16E3)') 'maximum_EN_Td=', &
            maxval(hypot(p%f%er,p%f%ez),mask=p%gas)/number_density/1e-21_dp
        if(maxval(hypot(p%f%er,p%f%ez),mask=p%gas)/number_density/1e-21_dp>500) &
            write(u,'(A)') 'WARNING: E/N exceeds 500 Td locally; assess the drift-diffusion/mean-energy closure carefully.'
        close(u)
        call cpu_time(t1)
        write(*,'(A,I0,A,ES12.4,A,F10.3)') 'Plasma complete: steps=',step, &
            '; peak ne [m^-3]=',maxval(p%ne),'; CPU seconds=',t1-t0
        write(*,'(A,ES12.4,A,ES12.4)') 'Plasma charge error=',qerr,'; energy error=',werr
        write(*,'(A)') 'See plasma_summary.txt for physical limits; Te is an effective mean-energy temperature.'
    end subroutine

    subroutine plasma_self_tests()
        type(swarm_t) :: s
        type(plasma_t) :: p
        real(dp) :: N,mu,d,ion,attach,en,eq,ne,np,nm,factor,loss,rhs,energy,cl,cr,pe,err,wall
        real(dp) :: qbefore,qafter,x(12,16),mm(12,16),dd(12,16),before,after
        integer :: i,j,k
        logical :: ok
        N=101325/(boltzmann*298.15_dp)
        call read_swarm(s,trim(p_transport_file),N,298.15_dp)
        do k=-100,100
            pe=real(k,dp)
            call require(abs(bernoulli(-pe)-bernoulli(pe)-pe)<1e-12_dp,'Bernoulli identity failed.')
        end do
        call sg_face(1.0_dp,0.04_dp,-1.0_dp,1.0_dp,1.0_dp,0.1_dp,0.1_dp,cl,cr)
        call require(abs(cl-cr*exp(0.4_dp))<1e-14_dp,'Scharfetter-Gummel Boltzmann equilibrium failed.')
        write(*,'(A)') 'TEST Bernoulli stability and SG Boltzmann equilibrium: passed'
        do k=1,3
            en=real(k*100,dp); eq=swarm_equilibrium(s,en)
            ! Use an actual knot so the reconstructed loss and coefficients balance exactly.
            i=minloc(abs(s%en-en),dim=1); en=s%en(i); eq=s%energy(i)
            call lookup_swarm(s,eq,mu,d,ion,attach)
            ne=1e12_dp; np=ne; nm=0; energy=ne*eq
            qbefore=np-ne-nm
            call chemical_step(ne,np,nm,ion,attach,0.0_dp,0.0_dp,1e-13_dp,factor)
            rhs=energy+1e-13_dp*ne*mu*(en*1e-21_dp*N)**2
            call energy_root(s,ne,rhs,1e-13_dp,energy,loss,ok)
            call require(ok.and.abs(energy/ne-eq)<1e-11_dp,'Homogeneous electron energy equilibrium failed.')
            qafter=np-ne-nm
            call require(abs(qafter-qbefore)<1e-3_dp,'Ionization/attachment charge conservation failed.')
        end do
        ne=1e12_dp; np=ne; nm=0
        call chemical_step(ne,np,nm,0.0_dp,1e8_dp,0.0_dp,0.0_dp,1e-9_dp,factor)
        call require(abs(ne-1e12_dp/1.1_dp)<1e-3_dp.and.abs(np-ne-nm)<1e-3_dp,'Attachment transfer test failed.')
        write(*,'(A)') 'TEST homogeneous energy closure and charge-conserving reactions: passed'
        call init_field(p%f,12,16,2e-4_dp,2e-3_dp,1.0_dp,2.0_dp,0.0_dp,2e-4_dp,.true.,.true.)
        allocate(p%gas(12,16),p%outer_phi(16),p%bottom_phi(12))
        p%gas=.true.; p%outer_phi=0; p%bottom_phi=0; p%f%phi=0
        do j=1,16
            do i=1,12
                x(i,j)=1+exp(-((p%f%r(i)-1e-4_dp)/3e-5_dp)**2-((p%f%z(j)-1e-3_dp)/1e-4_dp)**2)
            end do
        end do
        mm=0; dd=0.1_dp; before=sum(x*p%f%vol)
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,x,mm,dd,-1.0_dp,1e-10_dp,.false.,wall,err)
        after=sum(x*p%f%vol)
        call require(abs(after-before)/before<1e-12_dp.and.wall==0,'Closed FV transport conservation failed.')
        p%gas(1:2,:)=.false.; where(.not.p%gas) x=0
        call transport(p%f,p%gas,p%outer_phi,p%bottom_phi,x,mm,dd,-1.0_dp,1e-10_dp,.true.,wall,err)
        call require(wall>0.and.err<1e-12_dp,'Absorbing wall flux balance failed.')
        write(*,'(A)') 'TEST nonuniform cylindrical FV transport, positivity and wall accounting: passed'
        write(*,'(A)') 'ALL PLASMA SELF-TESTS PASSED'
    end subroutine
end module air_plasma_solver
! END AIR PLASMA MODULE

program electrospinning_rz
    use spinning_output
    use air_plasma_solver
    implicit none
    type(field_t) :: field
    type(jet_t) :: jet
    character(512) :: config_path,arg,preset,cli_prefix,plasma_config_path
    integer :: k,nargs,it,logunit,pcg_total
    real(dp) :: err,errq,erre,erru,deperr,t0,t1,lambda_scale
    real(dp), allocatable :: source_lambda(:),source_radius(:),old_e(:),old_u(:),candidate_e(:)
    logical :: exists,selftest,help,plasma_on,plasma_test
    config_path=''; preset=''; cli_prefix=''; plasma_config_path=''
    selftest=.false.; help=.false.; plasma_on=.false.; plasma_test=.false.
    nargs=command_argument_count(); k=1
    do while(k<=nargs)
        call get_command_argument(k,arg)
        select case(trim(arg))
        case('--config','--preset','--output','--plasma-config')
            call require(k<nargs,'Missing value after '//trim(arg))
            k=k+1
            select case(trim(arg))
            case('--config'); call get_command_argument(k,config_path)
            case('--preset'); call get_command_argument(k,preset)
            case('--output'); call get_command_argument(k,cli_prefix)
            case('--plasma-config'); call get_command_argument(k,plasma_config_path); plasma_on=.true.
            end select
        case('--self-test'); selftest=.true.
        case('--plasma'); plasma_on=.true.
        case('--plasma-test'); plasma_test=.true.
        case('--help','-h'); help=.true.
        case default
            call require(.false.,'Unknown argument: '//trim(arg))
        end select
        k=k+1
    end do
    if(help) then
        write(*,'(A)') 'Usage: electrospinning [--config file] [--preset quick|standard|fine|heavy] [--output prefix]'
        write(*,'(A)') '       electrospinning --self-test'
        write(*,'(A)') '       electrospinning --config config_spinning.nml --plasma-config plasma_air.nml'
        write(*,'(A)') '       electrospinning --plasma-test (requires air_transport.csv)'
        write(*,'(A)') 'Config: &spinning namelist in SI units, or original 11..14-line format.'
        write(*,'(A)') 'Without --config: config_spinning.nml, then config_spinning.txt, otherwise illustrative defaults.'
        stop
    end if
    if(selftest) then
        call field_self_tests(); stop
    end if
    if(plasma_test) then
        call plasma_self_tests(); stop
    end if
    if(len_trim(config_path)==0) then
        inquire(file='config_spinning.nml',exist=exists)
        if(exists) then
            config_path='config_spinning.nml'
        else
            inquire(file='config_spinning.txt',exist=exists)
            if(exists) config_path='config_spinning.txt'
        end if
    end if
    if(len_trim(config_path)>0) then
        call read_configuration(trim(config_path))
    else
        write(*,'(A)') 'WARNING: using illustrative defaults, not measured installation data.'
    end if
    if(len_trim(preset)>0) call set_preset(trim(preset))
    if(len_trim(cli_prefix)>0) output_prefix=cli_prefix
    call validate_config()
    failure_status_path=trim(output_prefix)//'run_status.txt'
    logunit=output_file('run_status.txt')
    write(logunit,'(A)') 'RUNNING: final CSV files are valid only after NUMERICALLY_CONVERGED.'
    close(logunit)
    call cpu_time(t0)
    write(*,'(A,I0,A,I0,A,I0)') 'Field cells: ',Nr,' x ',Nz,'; jet nodes: ',Nj
    write(*,'(A,ES12.4,A,ES12.4)') 'Inlet speed [m/s]: ',u_in,'; signed current [A]: ',total_current
    call init_field(field,Nr,Nz,R_max,Z_max,stretch_r,stretch_z,V_voltage,R_needle, &
        trim(side_boundary)=='insulating',.false.)
    call init_jet(jet)
    allocate(source_lambda(Nj),source_radius(Nj),old_e(Nj),old_u(Nj),candidate_e(Nj))
    source_lambda=0; source_radius=jet%a; old_e=0; old_u=jet%u; deperr=0; pcg_total=0
    logunit=output_file('convergence.csv')
    write(logunit,'(A)') 'outer_iteration,field_change,velocity_change,source_change,poisson_residual,pcg_iterations,jet_iterations'
    do it=1,max_outer
        if(space_charge) then
            call deposit_charge(field,jet%z,source_radius,source_lambda,deperr)
        else
            field%rho=0
        end if
        call solve_field(field,poisson_tol,max_pcg)
        pcg_total=pcg_total+field%last_iterations
        call field_gradients(field,R_needle)
        call sample_jet_field(field,jet%z,jet%a,candidate_e)
        jet%e=candidate_e
        call solve_jet(jet)
        erre=maxval(abs(jet%e-old_e))/max(maxval(abs(jet%e)),abs(V_voltage)/Z_max,1e-20_dp)
        erru=maxval(abs(jet%u-old_u)/max(abs(jet%u),u_in))
        lambda_scale=max(maxval(abs(jet%lambda)),abs(total_current)/u_in,1e-20_dp)
        errq=maxval(abs(jet%lambda-source_lambda))/lambda_scale
        if(.not.space_charge) errq=0
        err=max(erre,erru,errq)
        write(logunit,'(I0,4(",",ES16.8E3),2(",",I0))') it,erre,erru,errq,field%last_residual, &
            field%last_iterations,jet%last_iterations
        flush(logunit)
        write(*,'(A,I3,A,ES10.3,A,I5,A,I5)') 'Outer ',it,' error=',err, &
            ' PCG=',field%last_iterations,' jet=',jet%last_iterations
        if(it>=2 .and. err<outer_tol) exit
        source_lambda=(1-charge_relaxation)*source_lambda+charge_relaxation*jet%lambda
        source_radius=(1-charge_relaxation)*source_radius+charge_relaxation*jet%a
        old_e=jet%e; old_u=jet%u
    end do
    close(logunit)
    call require(it<=max_outer,'Field/jet coupling did not converge; no final physical result was exported.')
    call write_results(field,jet,it,err,deperr)
    if(plasma_on) then
        logunit=output_file('run_status.txt')
        write(logunit,'(A)') 'RUNNING_PLASMA: liquid stage complete; plasma transient has not finished.'
        close(logunit)
        call run_air_plasma(field,jet,trim(plasma_config_path))
        logunit=output_file('run_status.txt')
        write(logunit,'(A)') 'NUMERICALLY_CONVERGED'
        write(logunit,'(A)') 'PLASMA_TRANSIENT_COMPLETED: read plasma_summary.txt for scope and assumptions.'
        close(logunit)
    end if
    call cpu_time(t1)
    write(*,'(A,I0,A,F12.3)') 'Total PCG iterations: ',pcg_total,'; CPU seconds: ',t1-t0
end program
