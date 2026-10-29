! Transient Lagrangian, lumped-mass Maxwell filament. All quantities are SI.
! See docs/PHYSICS_RU.md. This is a reduced model, not a resolved fluid simulation.
module jet3d_model
 use iso_fortran_env, only: dp=>real64, int64
 use ieee_arithmetic, only: ieee_is_finite
 implicit none
 real(dp),parameter :: pi=acos(-1.0_dp),eps0=8.8541878128e-12_dp,ke=1/(4*pi*eps0),rgas=8.314462618_dp
 real(dp) :: voltage=10000, gap=.1_dp, flow_rate=1.388888888888889e-10_dp, jet_radius=3e-6_dp
 real(dp) :: density_polymer=1200,density_solvent=997,polymer_fraction=.12_dp,viscosity=.2_dp
 real(dp) :: relaxation_time=.01_dp,surface_tension=.055_dp,current=2e-8_dp,conductivity=2e-6_dp
 real(dp) :: relative_permittivity=40,temperature=298.15_dp,ambient_temperature=298.15_dp,humidity=.4_dp
 real(dp) :: mass_transfer=.01_dp,heat_transfer=100,latent_heat=2.44e6_dp,molar_mass=.01801528_dp
 real(dp) :: psat_ref=3169,psat_temperature=298.15_dp,cp_polymer=1700,cp_solvent=4180
 real(dp) :: viscosity_exponent=3,activation_energy=20000,air_density=1.184_dp,air_viscosity=1.85e-5_dp
 real(dp) :: drag_coefficient=.01_dp,air_velocity(3)=0,gravity(3)=[0._dp,0._dp,9.81_dp]
 real(dp) :: duration=.03_dp,dt_max=2e-6_dp,dt_min=1e-12_dp,injection_interval=2e-4_dp
 real(dp) :: output_interval=2e-4_dp,perturbation_amplitude=1e-6_dp,perturbation_frequency=250
 real(dp) :: collector_velocity(2)=0,coulomb_softening=3e-6_dp
 integer :: max_beads=500,max_steps=2000000
 logical :: evaporation=.true.,coulomb=.true.,image_force=.false.
 ! Legacy defaults remain available for comparison. Desktop presets opt in.
 character(32) :: perturbation_mode='periodic',constitutive_model='maxwell',drag_model='effective'
 integer :: random_seed=42
 real(dp) :: noise_correlation_time=4e-4_dp,finite_extensibility=100._dp,noise_sample_interval=1e-5_dp
 logical :: charge_transport=.false.
 logical :: surface_contact=.false.
 character(32) :: transport_model='constant',activity_model='volume_fraction'
 real(dp) :: gas_boundary_radius=.001_dp,vapor_diffusivity=2.5e-5_dp
 real(dp) :: air_thermal_conductivity=.0262_dp,air_heat_capacity=1005,air_molar_mass=.02897_dp
 real(dp) :: ambient_pressure=101325,flory_huggins_chi=.5_dp
 logical :: liquid_diffusion=.false.
 real(dp) :: liquid_diffusivity=1e-9_dp,diffusion_concentration_exponent=3
 character(1024) :: field_file=''
 namelist /jet3d/ voltage,gap,flow_rate,jet_radius,density_polymer,density_solvent,polymer_fraction, &
 viscosity,relaxation_time,surface_tension,current,conductivity,relative_permittivity,temperature, &
 ambient_temperature,humidity,mass_transfer,heat_transfer,latent_heat,molar_mass,psat_ref, &
 psat_temperature,cp_polymer,cp_solvent,viscosity_exponent,activation_energy,air_density,air_viscosity, &
 drag_coefficient,air_velocity,gravity,duration,dt_max,dt_min,injection_interval,output_interval, &
 max_beads,max_steps,perturbation_amplitude,perturbation_frequency,collector_velocity,coulomb_softening, &
 evaporation,coulomb,image_force,field_file,perturbation_mode,random_seed,noise_correlation_time, &
 constitutive_model,finite_extensibility,drag_model,charge_transport,surface_contact, &
 transport_model,activity_model,gas_boundary_radius,vapor_diffusivity,air_thermal_conductivity, &
 air_heat_capacity,air_molar_mass,ambient_pressure,flory_huggins_chi,liquid_diffusion, &
 liquid_diffusivity,diffusion_concentration_exponent,noise_sample_interval
 real(dp),allocatable :: pos(:,:),vel(:,:),force(:,:),mp(:),ms(:),q(:),temp(:),stress(:)
 real(dp),allocatable :: ell(:),radius(:),eta(:),oldell(:),evrate(:),heatrate(:),land_time(:)
 real(dp),allocatable :: fr(:),fz(:),fer(:,:),fez(:,:)
 real(dp),allocatable :: fphi(:,:),azz(:),arr(:),qrate(:),pmat(:,:),gcond(:),drag_damping(:)
 real(dp),allocatable :: travel(:),birth(:)
 real(dp),allocatable :: km_air(:),h_air(:),re_air(:),phi_surface(:),transfer_valid(:),flux_air(:)
 real(dp),allocatable :: noise_history(:,:)
 real(dp) :: charge_step=huge(1._dp),ground_current=0,noise_xy(2)=0
 integer(int64) :: rng_state=42
 logical,allocatable :: landed(:)
 integer :: n=0,nland=0,steps=0,nr=0,nz=0,ut,ud,ub
 real(dp) :: t=0,dt=0,mp0,ms0,m0,q0,u0,l0,evaporated=0,delivered=0,next_output=0
contains
 subroutine require(ok,msg)
  logical,intent(in)::ok
  character(*),intent(in)::msg
  integer::u
  if(ok)return
  open(newunit=u,file='run_status.txt',status='replace');write(u,'(a)')'FAILED: '//msg;close(u)
  write(*,'(a)')'FAILED: '//msg
  error stop 1
 end subroutine
 subroutine load_config(path)
  character(*),intent(in)::path
  integer::u,ios
  character(1024)::msg
  real(dp)::vals(49)
  open(newunit=u,file=trim(path),status='old',iostat=ios,iomsg=msg)
  call require(ios==0,'Cannot open namelist: '//trim(msg))
  read(u,nml=jet3d,iostat=ios,iomsg=msg);close(u)
  call require(ios==0,'Invalid &jet3d namelist: '//trim(msg))
  vals=[voltage,gap,flow_rate,jet_radius,density_polymer,density_solvent,polymer_fraction,viscosity, &
   relaxation_time,surface_tension,current,conductivity,relative_permittivity,temperature,ambient_temperature, &
   humidity,mass_transfer,heat_transfer,latent_heat,molar_mass,psat_ref,psat_temperature,cp_polymer,cp_solvent, &
   viscosity_exponent,activation_energy,air_density,air_viscosity,drag_coefficient,air_velocity,gravity,duration, &
   dt_max,dt_min,injection_interval,output_interval,perturbation_amplitude,perturbation_frequency, &
   collector_velocity,coulomb_softening,real(max_beads,dp),real(max_steps,dp),0._dp,0._dp]
  call require(all(ieee_is_finite(vals)),'Nonfinite input parameter')
  call require(gap>0.and.flow_rate>0.and.jet_radius>0,'gap, flow_rate and jet_radius must be positive')
  call require(min(density_polymer,density_solvent,viscosity,relaxation_time)>0,'Invalid density/viscosity/relaxation time')
  call require(polymer_fraction>0.and.polymer_fraction<1,'polymer_fraction must be a mass fraction between 0 and 1')
  call require(temperature>150.and.ambient_temperature>150.and.psat_temperature>150,'Temperatures must be in kelvin (>150)')
  call require(humidity>=0.and.humidity<=1,'humidity must be a fraction in [0,1]')
  call require(min(surface_tension,conductivity,mass_transfer,heat_transfer,drag_coefficient,air_density)>=0,'Negative transport coefficient')
  call require(min(relative_permittivity,cp_polymer,cp_solvent,molar_mass,psat_ref,latent_heat)>0,'Invalid thermal/electrical property')
  call require(duration>0.and.dt_min>0.and.dt_max>=dt_min.and.output_interval>0,'Invalid time controls')
  call require(injection_interval>0.and.injection_interval<duration,'injection_interval must be less than duration')
  call require(coulomb_softening>0.and.perturbation_amplitude>=0,'Invalid softening or perturbation')
  call require(max_beads>=2.and.max_steps>=1,'Invalid allocation or step limit')
  call require(int(duration/injection_interval+1e-8_dp)<=max_beads,'max_beads too small; increase to ceil(duration/injection_interval)')
  call require(voltage*current>=0,'voltage and current signs must agree (both polarities supported)')
  call require(trim(perturbation_mode)=='periodic'.or.trim(perturbation_mode)=='colored_noise','Invalid perturbation_mode')
  call require(trim(constitutive_model)=='maxwell'.or.trim(constitutive_model)=='fene_p','Invalid constitutive_model')
  call require(trim(drag_model)=='effective'.or.trim(drag_model)=='cylinder','Invalid drag_model')
  call require(ieee_is_finite(noise_correlation_time).and.noise_correlation_time>0,'Invalid noise correlation time')
  call require(ieee_is_finite(noise_sample_interval).and.noise_sample_interval>0,'Invalid noise sample interval')
  call require(duration/noise_sample_interval<1000000,'Noise history too large; increase noise sample interval')
  call require(ieee_is_finite(finite_extensibility).and.finite_extensibility>3,'FENE L squared must exceed 3')
  call require(air_viscosity>0,'Air dynamic viscosity must be positive')
  call require(trim(transport_model)=='constant'.or.trim(transport_model)=='gas_film'.or. &
   trim(transport_model)=='crossflow','Invalid transport_model')
  call require(trim(activity_model)=='volume_fraction'.or.trim(activity_model)=='flory_huggins','Invalid activity_model')
  call require(all(ieee_is_finite([gas_boundary_radius,vapor_diffusivity,air_thermal_conductivity, &
   air_heat_capacity,air_molar_mass,ambient_pressure,flory_huggins_chi,liquid_diffusivity, &
   diffusion_concentration_exponent])),'Nonfinite transport property')
  call require(min(gas_boundary_radius,vapor_diffusivity,air_thermal_conductivity,air_heat_capacity, &
   air_molar_mass,ambient_pressure,liquid_diffusivity)>0,'Invalid gas/liquid transport property')
  call require(flory_huggins_chi>=0.and.flory_huggins_chi<=.5_dp,'Single-phase Flory-Huggins closure requires 0<=chi<=0.5')
  call require(diffusion_concentration_exponent>=0,'Diffusion exponent must be >=0')
  if(trim(transport_model)=='gas_film')call require(gas_boundary_radius>jet_radius,'Gas boundary must exceed jet radius')
  if(trim(transport_model)/='constant')call require(air_density>0,'Gas transfer requires positive air density')
  call require(random_seed>0,'Seed must be a positive integer')
  rng_state=modulo(int(random_seed,int64)-1,2147483646_int64)+1
  ! Additive phase volumes define inlet solution density from MASS fraction.
  m0=flow_rate*injection_interval/(polymer_fraction/density_polymer+(1-polymer_fraction)/density_solvent)
  mp0=polymer_fraction*m0;ms0=(1-polymer_fraction)*m0
  q0=current*injection_interval;u0=flow_rate/(pi*jet_radius**2);l0=u0*injection_interval
  call require(l0<gap*.2_dp,'Injection segment too long for gap; reduce injection_interval')
  call require(perturbation_amplitude<.1_dp*l0,'Inlet perturbation must be small versus injection segment length')
  allocate(pos(3,max_beads),vel(3,max_beads),force(3,max_beads),mp(max_beads),ms(max_beads),q(max_beads), &
   temp(max_beads),stress(max_beads),ell(max_beads),radius(max_beads),eta(max_beads),oldell(max_beads), &
   evrate(max_beads),heatrate(max_beads),land_time(max_beads),landed(max_beads))
  pos=0;vel=0;force=0;mp=0;ms=0;q=0;temp=temperature;stress=0;landed=.false.;land_time=0
  allocate(azz(max_beads),arr(max_beads),qrate(max_beads),pmat(max_beads,max_beads), &
   gcond(max_beads),drag_damping(max_beads))
  azz=1;arr=1;qrate=0;pmat=0;gcond=0;drag_damping=0
  allocate(travel(max_beads),birth(max_beads));travel=0;birth=0
  allocate(km_air(max_beads),h_air(max_beads),re_air(max_beads),phi_surface(max_beads),transfer_valid(max_beads),flux_air(max_beads))
  km_air=0;h_air=0;re_air=0;phi_surface=0;transfer_valid=1;flux_air=0
  if(trim(perturbation_mode)=='colored_noise')call prepare_noise()
  if(len_trim(field_file)>0)call load_field()
  open(newunit=u,file='effective_jet3d.nml',status='replace');write(u,nml=jet3d);close(u)
 end subroutine
 subroutine prepare_noise()
  integer::j,last
  real(dp)::decay
  last=ceiling(duration/noise_sample_interval)+1
  allocate(noise_history(2,0:last));noise_history(:,0)=0
  decay=exp(-noise_sample_interval/noise_correlation_time)
  do j=1,last
   noise_history(:,j)=decay*noise_history(:,j-1)+sqrt(1-decay**2)*[normal_random(),normal_random()]/sqrt(2._dp)
  enddo
 end subroutine
 function inlet_noise(time)result(xy)
  real(dp),intent(in)::time
  real(dp)::xy(2),fraction
  integer::j
  j=min(int(time/noise_sample_interval),size(noise_history,2)-2)
  fraction=(time-j*noise_sample_interval)/noise_sample_interval
  xy=(1-fraction)*noise_history(:,j)+fraction*noise_history(:,j+1)
 end function
 subroutine load_field()
  integer::u,ios,count,i,j
  real(dp)::row(7),zfirst
  character(4096)::line
  open(newunit=u,file=trim(field_file),status='old',iostat=ios)
  call require(ios==0,'Cannot open external field CSV')
  read(u,'(a)',iostat=ios)line
  call require(index(line,'E_r_V_m')>0.and.index(line,'E_z_V_m')>0,'External field CSV header incompatible')
  count=0;nr=0;zfirst=huge(1._dp)
  do
   read(u,*,iostat=ios)row
   if(ios<0)exit
   call require(ios==0,'Malformed field CSV row')
   call require(all(ieee_is_finite(row)),'Nonfinite external field')
   count=count+1
   if(count==1)zfirst=row(2)
   if(row(2)==zfirst)nr=nr+1
  enddo
  call require(nr>=2.and.mod(count,nr)==0,'External field must be tensor-product CSV, r varying fastest')
  nz=count/nr;call require(nz>=2,'External field needs >=2 axial points')
  allocate(fr(nr),fz(nz),fer(nr,nz),fez(nr,nz),fphi(nr,nz))
  rewind(u);read(u,'(a)')line
  do j=1,nz
   do i=1,nr
    read(u,*)row
    if(j==1)fr(i)=row(1)
    if(i==1)fz(j)=row(2)
    call require(abs(row(1)-fr(i))<1e-12_dp.and.abs(row(2)-fz(j))<1e-12_dp,'Inconsistent field tensor axes')
    fer(i,j)=row(4);fez(i,j)=row(5)
    fphi(i,j)=row(3)
   enddo
  enddo
  close(u)
  call require(all(fr(2:)>fr(:nr-1)).and.all(fz(2:)>fz(:nz-1)),'Field axes must increase')
  call require(fr(1)>=0.and.minval(fz)<=.01_dp*gap.and.maxval(fz)>=.9_dp*gap,'Field must span nozzle to collector')
 end subroutine
 integer function lower(a,x)result(i)
  real(dp),intent(in)::a(:),x
  integer::lo,hi,mid
  lo=1;hi=size(a)
  do while(hi-lo>1)
   mid=(lo+hi)/2
   if(a(mid)<=x)then
    lo=mid
   else
    hi=mid
   endif
  enddo
  i=min(lo,size(a)-1)
 end function
 function potential(x)result(v)
  real(dp),intent(in)::x(3)
  real(dp)::v,r,z,ar,az
  integer::i,j
  if(nr==0)then
   v=voltage*(1-x(3)/gap);return
  endif
  r=max(fr(1),sqrt(sum(x(1:2)**2)));z=max(fz(1),min(fz(nz),x(3)))
  i=lower(fr,r);j=lower(fz,z)
  ar=(r-fr(i))/(fr(i+1)-fr(i));az=(z-fz(j))/(fz(j+1)-fz(j))
  v=(1-az)*((1-ar)*fphi(i,j)+ar*fphi(i+1,j))+az*((1-ar)*fphi(i,j+1)+ar*fphi(i+1,j+1))
 end function
 function normal_random()result(v)
  real(dp)::v,u1,u2
  ! Park-Miller with 64-bit arithmetic; portable seed, no compiler RNG dependence.
  rng_state=modulo(16807_int64*rng_state,2147483647_int64);u1=real(rng_state,dp)/2147483647._dp
  rng_state=modulo(16807_int64*rng_state,2147483647_int64);u2=real(rng_state,dp)/2147483647._dp
  v=sqrt(-2*log(u1))*cos(2*pi*u2)
 end function
 subroutine fene_update(a,b,strain,h)
  real(dp),intent(inout)::a,b
  real(dp),intent(in)::strain,h
  real(dp)::aa,bb,tr,rate,total,disc,den
  ! Exact uniaxial affine deformation, implicit Peterlin relaxation.
  ! Positive root evaluation avoids subtracting nearly equal large numbers.
  aa=a*exp(2*strain*h);bb=b*exp(-strain*h);rate=h/relaxation_time
  total=aa+2*bb+3*rate
  den=finite_extensibility+total+rate*(finite_extensibility-3)
  disc=max(0._dp,den**2-4*total*finite_extensibility)
  tr=2*total*finite_extensibility/(den+sqrt(disc))
  call require(tr<finite_extensibility.and.tr>0,'FENE conformation trace invalid; reduce timestep')
  den=1+rate*(finite_extensibility-3)/(finite_extensibility-tr)
  a=(aa+rate)/den;b=(bb+rate)/den
 end subroutine
 function electric(x)result(e)
  real(dp),intent(in)::x(3)
  real(dp)::e(3),r,ar,az,er,ez,z,rr
  integer::i,j
  if(nr==0)then
   e=[0._dp,0._dp,voltage/gap];return
  endif
  r=sqrt(sum(x(1:2)**2));z=max(fz(1),min(fz(nz),x(3)))
  call require(r<=fr(nr),'Jet left external field radial domain; enlarge field grid/chamber')
  rr=max(fr(1),r);i=lower(fr,rr);j=lower(fz,z)
  ar=(rr-fr(i))/(fr(i+1)-fr(i));az=(z-fz(j))/(fz(j+1)-fz(j))
  er=(1-az)*((1-ar)*fer(i,j)+ar*fer(i+1,j))+az*((1-ar)*fer(i,j+1)+ar*fer(i+1,j+1))
  ez=(1-az)*((1-ar)*fez(i,j)+ar*fez(i+1,j))+az*((1-ar)*fez(i,j+1)+ar*fez(i+1,j+1))
  if(r<fr(1).and.fr(1)>0)er=er*r/fr(1)
  e=[0._dp,0._dp,ez]
  if(r>0)e(1:2)=er*x(1:2)/r
 end function
 subroutine geometry()
  real(dp)::up(3),w,expo
  integer::i
  do i=1,n
   ! Once both endpoints have landed, store the last material segment geometry.
   ! A zero planar chord is not a zero filament volume/length (unresolved slack).
   if(landed(i))then
    if(i==n)cycle
    if(landed(i+1))cycle
   endif
   up=0
   if(i<n)up=pos(:,i+1)
   ell(i)=sqrt(sum((pos(:,i)-up)**2))
   call require(ell(i)>1e-10_dp*l0,'Collapsed segment; refine timestep/injection discretization')
   radius(i)=sqrt((mp(i)/density_polymer+ms(i)/density_solvent)/(pi*ell(i)))
   w=mp(i)/(mp(i)+ms(i));expo=activation_energy/rgas*(1/temp(i)-1/temperature)
   call require(abs(expo)<500,'Constitutive viscosity outside numerical range')
   eta(i)=viscosity*(w/polymer_fraction)**viscosity_exponent*exp(expo)
   call require(ieee_is_finite(eta(i)).and.eta(i)>0,'Invalid constitutive viscosity')
  enddo
 end subroutine
 subroutine forces()
  real(dp)::up(3),d(3),f(3),tension,r2,vrel(3),speed,area,vol,activity,surf
  real(dp)::tangent(3),vnormal(3),vparallel(3),vn,re,cd,normal_coef,parallel_coef,factor,dissipation,flux
  integer::i,j
  force=0;evrate=0;heatrate=0;drag_damping=0
  do i=1,n
   if(landed(i))cycle
   force(:,i)=q(i)*electric(pos(:,i))+(mp(i)+ms(i))*gravity
   vrel=vel(:,i)-air_velocity;speed=sqrt(sum(vrel**2));area=2*radius(i)*ell(i)
   ! Effective quadratic cylinder drag (Cd requires calibration for actual airflow).
   if(trim(drag_model)=='effective')then
    drag_damping(i)=.5_dp*air_density*drag_coefficient*area*speed
    force(:,i)=force(:,i)-drag_damping(i)*vrel
   else
    up=0;if(i<n)up=pos(:,i+1)
    tangent=(pos(:,i)-up)/ell(i);vparallel=dot_product(vrel,tangent)*tangent
    vnormal=vrel-vparallel;vn=sqrt(sum(vnormal*vnormal))
    re=2*radius(i)*air_density*vn/air_viscosity
    cd=1.2_dp+10/max(re,1e-20_dp)**(2._dp/3)
    normal_coef=.5_dp*air_density*cd*area*vn
    parallel_coef=2*pi*air_viscosity*ell(i)/max(log(2*ell(i)/radius(i))-.5_dp,1._dp)
    force(:,i)=force(:,i)-normal_coef*vnormal-parallel_coef*vparallel
    drag_damping(i)=max(normal_coef,parallel_coef)
   endif
  enddo
  ! The newest cell is an open inlet reservoir. Incoming momentum enters at u_in.
  if(n>0)then
   if(.not.landed(n))force(:,n)=force(:,n)+m0/injection_interval*([0._dp,0._dp,u0]-vel(:,n))
  endif
  ! Maxwell/capillary segment tension, applied equally and oppositely.
  do i=1,n
   if(landed(i))then
    if(i==n)cycle
    if(landed(i+1))cycle
   endif
   up=0
   if(i<n)up=pos(:,i+1)
   d=(pos(:,i)-up)/ell(i)
   tension=pi*radius(i)**2*stress(i)+pi*surface_tension*radius(i)
   f=-tension*d
   if(.not.landed(i))force(:,i)=force(:,i)+f
   if(i<n)then
    if(.not.landed(i+1))force(:,i+1)=force(:,i+1)-f
   endif
  enddo
  if(coulomb)then
   do i=1,n-1
    if(landed(i))cycle
    do j=i+1,n
     if(landed(j))cycle
     d=pos(:,i)-pos(:,j);r2=sum(d*d)+coulomb_softening**2
     f=ke*q(i)*q(j)*d/(r2*sqrt(r2))
     force(:,i)=force(:,i)+f;force(:,j)=force(:,j)-f
    enddo
   enddo
  endif
  if(image_force)then
   ! Images in conducting plane z=gap, opposite charge. Self-image included.
   do i=1,n
    if(landed(i))cycle
    do j=1,n
     if(landed(j))cycle
     up=pos(:,j);up(3)=2*gap-up(3);d=pos(:,i)-up
     r2=sum(d*d)+coulomb_softening**2
     force(:,i)=force(:,i)-ke*q(i)*q(j)*d/(r2*sqrt(r2))
    enddo
   enddo
  endif
  do i=1,n
   if(landed(i))cycle
   vol=mp(i)/density_polymer+ms(i)/density_solvent
   surf=2*pi*radius(i)*ell(i)
   activity=(ms(i)/density_solvent)/vol
   up=0;if(i<n)up=pos(:,i+1)
   tangent=(pos(:,i)-up)/ell(i);vrel=vel(:,i)-air_velocity
   vnormal=vrel-dot_product(vrel,tangent)*tangent;vn=sqrt(sum(vnormal*vnormal))
   call evaporation_flux(radius(i),temp(i),activity,vn,flux,km_air(i),h_air(i),re_air(i), &
    phi_surface(i),transfer_valid(i))
   flux_air(i)=merge(flux,0._dp,evaporation)
   if(evaporation.and.ms(i)>0)evrate(i)=flux*surf
   heatrate(i)=h_air(i)*surf*(ambient_temperature-temp(i))-latent_heat*evrate(i) &
    +stress(i)**2/(3*eta(i))*vol
   if(trim(constitutive_model)=='fene_p')then
    factor=(finite_extensibility-3)/(finite_extensibility-azz(i)-2*arr(i))
    dissipation=eta(i)/(2*relaxation_time**2)*max(0._dp, &
     factor**2*(azz(i)+2*arr(i))-6*factor+1/azz(i)+2/arr(i))
    heatrate(i)=h_air(i)*surf*(ambient_temperature-temp(i))-latent_heat*evrate(i)+dissipation*vol
   endif
   if(i==n)heatrate(i)=heatrate(i)+(mp0*cp_polymer+ms0*cp_solvent)/injection_interval*(temperature-temp(i))
  enddo
  if(charge_transport)call conduction()
  call require(all(ieee_is_finite(force(:,:n))),'Nonfinite force; model diverged')
 end subroutine
 function cb_transfer(re,p)result(nu)
  real(dp),intent(in)::re,p
  real(dp)::nu
  nu=.3_dp+.62_dp*sqrt(re)*p**(1._dp/3)/(1+(.4_dp/p)**(2._dp/3))**.25_dp &
   *(1+(re/282000._dp)**.625_dp)**.8_dp
 end function
 function solvent_activity(phi)result(activity)
  real(dp),intent(in)::phi
  real(dp)::activity
  activity=0
  if(phi<=0)return
  if(phi>=1)then
   activity=1;return
  endif
  activity=phi
  if(trim(activity_model)=='flory_huggins')activity=exp(log(phi)+(1-phi)+flory_huggins_chi*(1-phi)**2)
 end function
 function vapor_fraction(pressure)result(y)
  real(dp),intent(in)::pressure
  real(dp)::x,y
  call require(pressure>=0.and.pressure<ambient_pressure,'Vapor pressure outside gas-film model: check temperature/pressure')
  x=pressure/ambient_pressure;y=molar_mass*x/(molar_mass*x+air_molar_mass*(1-x))
 end function
 function gas_evaporation(phi,ps,ya,km)result(flux)
  real(dp),intent(in)::phi,ps,ya,km
  real(dp)::flux,ys
  ys=vapor_fraction(solvent_activity(phi)*ps)
  flux=max(0._dp,air_density*km*log((1-ya)/(1-ys)))
 end function
 subroutine evaporation_flux(a,tk,phi,vn,flux,km,hc,re,phiw,valid)
  real(dp),intent(in)::a,tk,phi,vn
  real(dp),intent(out)::flux,km,hc,re,phiw,valid
  real(dp)::resistance,pr,sc,ps,pa,ya,phi0,diffusivity,conductance,lo,hi,mid
  integer::iteration
  km=mass_transfer;hc=heat_transfer;re=2*a*air_density*vn/air_viscosity;valid=1;phiw=phi
  ps=psat_ref*exp(latent_heat*molar_mass/rgas*(1/psat_temperature-1/tk))
  pa=psat_ref*exp(latent_heat*molar_mass/rgas*(1/psat_temperature-1/ambient_temperature))
  if(trim(transport_model)=='constant')then
   flux=max(0._dp,km*molar_mass/rgas*(solvent_activity(phi)*ps/tk-humidity*pa/ambient_temperature))
   return
  endif
  if(trim(transport_model)=='gas_film')then
   call require(gas_boundary_radius>a,'Jet reached gas diffusion boundary; gas-film geometry invalid')
   resistance=a*log(gas_boundary_radius/a)
   km=vapor_diffusivity/resistance;hc=air_thermal_conductivity/resistance
  else
   sc=air_viscosity/(air_density*vapor_diffusivity)
   pr=air_viscosity*air_heat_capacity/air_thermal_conductivity
   km=cb_transfer(re,sc)*vapor_diffusivity/(2*a);hc=cb_transfer(re,pr)*air_thermal_conductivity/(2*a)
   if(re*min(pr,sc)<.2_dp)valid=0
  endif
  ya=vapor_fraction(humidity*pa);flux=gas_evaporation(phi,ps,ya,km)
  if(liquid_diffusion.and.flux>0.and.phi>0)then
   phi0=((1-polymer_fraction)/density_solvent)/(polymer_fraction/density_polymer+(1-polymer_fraction)/density_solvent)
   diffusivity=liquid_diffusivity*(phi/phi0)**diffusion_concentration_exponent
   conductance=4*density_solvent*diffusivity/a;lo=0;hi=phi
   do iteration=1,45
    mid=(lo+hi)/2
    if(gas_evaporation(mid,ps,ya,km)>conductance*(phi-mid))then
     hi=mid
    else
     lo=mid
    endif
   enddo
   phiw=(lo+hi)/2;flux=gas_evaporation(phiw,ps,ya,km)
  endif
 end subroutine
 subroutine conduction()
  ! dq/dt=-L(Pq+phi_ext). L is the nearest-neighbour conductance graph;
  ! P uses the same positive-definite regularized kernel on AND off diagonal.
  ! Effective bead self-capacitance is softening/ke. Changing only the diagonal
  ! to a line-segment capacitance breaks positivity for endpoint lumped charges.
  ! Grounded deposited nodes have fixed potential 0.
  ! This is a reduced slender-conductor closure, not full leaky-dielectric CFD.
  real(dp)::phi(max_beads),d(3),im(3),flux,g,degree,bound,rowval,joule
  integer::i,j
  pmat=0;qrate=0;gcond=0;ground_current=0;charge_step=huge(1._dp)
  do i=1,n
   phi(i)=0
   if(landed(i))cycle
   phi(i)=potential(pos(:,i))
   do j=1,n
    if(landed(j))cycle
    d=pos(:,i)-pos(:,j)
    pmat(i,j)=ke/sqrt(sum(d*d)+coulomb_softening**2)
    if(image_force)then
     im=pos(:,j);im(3)=2*gap-im(3);d=pos(:,i)-im
     pmat(i,j)=pmat(i,j)-ke/sqrt(sum(d*d)+coulomb_softening**2)
    endif
    phi(i)=phi(i)+pmat(i,j)*q(j)
   enddo
  enddo
  do i=1,n-1
   if(landed(i).and.landed(i+1))cycle
   g=conductivity*pi*radius(i)**2/ell(i);gcond(i)=g
   flux=g*(phi(i+1)-phi(i))
   if(.not.landed(i))then
    qrate(i)=qrate(i)+flux
   else
    ground_current=ground_current+flux
   endif
   if(.not.landed(i+1))then
    qrate(i+1)=qrate(i+1)-flux
   else
    ground_current=ground_current-flux
   endif
   joule=g*(phi(i+1)-phi(i))**2
   if(.not.landed(i))heatrate(i)=heatrate(i)+.5_dp*joule
   if(.not.landed(i+1))heatrate(i+1)=heatrate(i+1)+.5_dp*joule
  enddo
  ! Infinity-norm bound for explicit charge relaxation; no charge/force clipping.
  do i=1,n
   if(landed(i))cycle
   degree=gcond(i);if(i>1)degree=degree+gcond(i-1)
   bound=0
   do j=1,n
    if(landed(j))cycle
    rowval=degree*pmat(i,j)
    if(i>1)rowval=rowval-gcond(i-1)*pmat(i-1,j)
    if(i<n)rowval=rowval-gcond(i)*pmat(i+1,j)
    bound=bound+abs(rowval)
   enddo
   if(bound>0)charge_step=min(charge_step,.2_dp/bound)
  enddo
 end subroutine
 function stable_step()result(h)
  real(dp)::h,scale,acc,speed,cap,elastic,kernel
  integer::i,j
  h=min(dt_max,charge_step)
  do i=1,n
   if(landed(i))cycle
   scale=ell(i)
   if(i>1)scale=min(scale,ell(i-1))
   acc=sqrt(sum(force(:,i)**2))/(mp(i)+ms(i))
   speed=sqrt(sum(vel(:,i)**2))
   h=min(h,.15_dp*scale/max(speed,1e-30_dp),.15_dp*sqrt(scale/max(acc,1e-30_dp)))
   elastic=3*eta(i)/relaxation_time*pi*radius(i)**2/ell(i)
   h=min(h,.15_dp*sqrt((mp(i)+ms(i))/max(elastic,1e-30_dp)))
   if(trim(constitutive_model)=='fene_p')then
    elastic=elastic*(finite_extensibility-3)/max(finite_extensibility-azz(i)-2*arr(i),1e-12_dp)
    h=min(h,.1_dp*sqrt((mp(i)+ms(i))/max(elastic,1e-30_dp)))
   endif
   if(drag_damping(i)>0)h=min(h,.15_dp*(mp(i)+ms(i))/drag_damping(i))
   if(evrate(i)>0)h=min(h,.05_dp*ms(i)/evrate(i))
   cap=mp(i)*cp_polymer+ms(i)*cp_solvent
   if(h_air(i)>0)h=min(h,.15_dp*cap/(h_air(i)*2*pi*radius(i)*ell(i)))
   if(evrate(i)>0)h=min(h,.15_dp*cap/(latent_heat*evrate(i)*latent_heat*molar_mass/(rgas*temp(i)**2)))
   if(abs(heatrate(i))>0)h=min(h,2._dp*cap/abs(heatrate(i)))
   if(image_force)h=min(h,.15_dp*max(gap-pos(3,i),coulomb_softening)/max(speed,1e-30_dp), &
    .15_dp*sqrt(max(gap-pos(3,i),coulomb_softening)/max(acc,1e-30_dp)))
  enddo
  ! Resolve regularized pair-force curvature, including close non-neighbours.
  if(coulomb)then
   do i=1,n-1
    if(landed(i))cycle
    do j=i+1,n
     if(landed(j))cycle
     kernel=ke*abs(q(i)*q(j))/(sum((pos(:,i)-pos(:,j))**2)+coulomb_softening**2)**1.5_dp
     if(kernel>0)h=min(h,.15_dp*sqrt(min(mp(i)+ms(i),mp(j)+ms(j))/kernel))
    enddo
   enddo
  endif
  call require(h>=dt_min,'Required timestep fell below dt_min; refine model/parameters, do not clip forces')
 end function
 subroutine inject()
  real(dp)::phase,fraction,newms,newq,oldmass,newmass,oldvelocity(3),oldtemperature
  call require(n<max_beads,'Bead allocation exhausted')
  ! The growing inlet cell is split conservatively into an old material cell
  ! and a new inlet reservoir. Solvent is well mixed within the reservoir.
  if(n>0)then
   call require(pos(3,n)>l0,'Inlet segment cannot accept next parcel; jet stalled/reversed near nozzle')
   oldmass=mp(n)+ms(n);oldvelocity=vel(:,n);oldtemperature=temp(n)
   fraction=mp0/mp(n)
   call require(fraction>0.and.fraction<1,'Invalid inlet reservoir split')
   newms=ms(n)*fraction;newq=q(n)*fraction;newmass=mp0+newms
   mp(n)=mp(n)-mp0;ms(n)=ms(n)-newms;q(n)=q(n)-newq
   vel(:,n)=(oldmass*oldvelocity-newmass*[0._dp,0._dp,u0])/(oldmass-newmass)
  else
   newms=ms0;newq=q0;oldtemperature=temperature
  endif
  n=n+1;phase=2*pi*perturbation_frequency*t
  if(trim(perturbation_mode)=='colored_noise')then
   ! Fixed physical noise clock: same path under BOTH time and injection refinement.
   noise_xy=inlet_noise(t)
   pos(:,n)=[perturbation_amplitude*noise_xy(1),perturbation_amplitude*noise_xy(2),l0]
  else
   pos(:,n)=[perturbation_amplitude*cos(phase),perturbation_amplitude*sin(phase),l0]
  endif
  vel(:,n)=[0._dp,0._dp,u0];mp(n)=mp0;ms(n)=newms;q(n)=newq
  temp(n)=oldtemperature;stress(n)=0;azz(n)=1;arr(n)=1;landed(n)=.false.
  travel(n)=sqrt(sum(pos(:,n)**2));birth(n)=t
 end subroutine
 subroutine advance(h)
  real(dp),intent(in)::h
  real(dp)::before(3),fraction,hit(3),hitv(3),cross_time,dryr,dm,capacity,strain,decay,up(3),w,contact_z
  integer::i
  oldell(:n)=ell(:n)
  if(charge_transport)then
   do i=1,n
    if(.not.landed(i))q(i)=q(i)+h*qrate(i)
   enddo
   delivered=delivered+h*ground_current
  endif
  do i=1,n
   if(landed(i))then
    pos(1:2,i)=pos(1:2,i)+collector_velocity*h
    cycle
   endif
   before=pos(:,i)
   contact_z=gap
   if(surface_contact)contact_z=gap-radius(i)
   vel(:,i)=vel(:,i)+.5_dp*h*force(:,i)/(mp(i)+ms(i))
   pos(:,i)=pos(:,i)+h*vel(:,i)
   call require(pos(3,i)>0,'Parcel crossed emitter plane; no steady emission for these parameters')
   dm=evrate(i)*h
   call require(dm<=ms(i)*(1+1e-12_dp),'Evaporation timestep consumed unavailable solvent')
   capacity=mp(i)*cp_polymer+ms(i)*cp_solvent
   temp(i)=temp(i)+h*heatrate(i)/capacity
   ms(i)=ms(i)-dm;evaporated=evaporated+dm
   if(i==n)then
    mp(i)=mp(i)+mp0*h/injection_interval
    ms(i)=ms(i)+ms0*h/injection_interval
    q(i)=q(i)+q0*h/injection_interval
   endif
   call require(temp(i)>150.and.temp(i)<1000,'Parcel temperature outside material model [150,1000] K')
   if(pos(3,i)>=contact_z)then
    call require(i<n,'Inlet reservoir reached collector before splitting; reduce injection_interval')
    fraction=max(0._dp,min(1._dp,(contact_z-before(3))/max(pos(3,i)-before(3),1e-30_dp)))
    cross_time=t+fraction*h
    travel(i)=travel(i)+fraction*sqrt(sum((pos(:,i)-before)**2))
    hit=before+fraction*(pos(:,i)-before);hit(3)=gap;hitv=vel(:,i)
    ! The deposited point subsequently follows the translating collector.
    pos(:,i)=hit;pos(1:2,i)=hit(1:2)+collector_velocity*(t+h-cross_time)
    landed(i)=.true.;nland=nland+1;land_time(i)=cross_time;delivered=delivered+q(i)
    up=0;if(i<n)up=pos(:,i+1)
    w=sqrt(sum((hit-up)**2));dryr=sqrt(mp(i)/(density_polymer*pi*max(w,1e-30_dp)))
    write(ud,'(es24.16e3,",",i0,10(",",es24.16e3))')cross_time,i,hit(1),hit(2), &
     hit(1)-collector_velocity(1)*cross_time,hit(2)-collector_velocity(2)*cross_time, &
     sqrt((mp(i)/density_polymer+ms(i)/density_solvent)/(pi*max(w,1e-30_dp))),dryr, &
     mp(i),ms(i),q(i),sqrt(sum(hitv**2))
    vel(:,i)=[collector_velocity(1),collector_velocity(2),0._dp]
   else
    travel(i)=travel(i)+sqrt(sum((pos(:,i)-before)**2))
   endif
  enddo
  call geometry()
  decay=exp(-h/relaxation_time)
  do i=1,n
   if(landed(i))cycle
   strain=log(ell(i)/oldell(i))/h
   ! Open inlet contains incoming material: subtract inlet material velocity.
   ! Using geometric growth against a fixed origin would create spurious strain.
   if(i==n)strain=dot_product(vel(:,i)-[0._dp,0._dp,u0],pos(:,i)/ell(i))/ell(i)
   if(trim(constitutive_model)=='maxwell')then
    stress(i)=stress(i)*decay+3*eta(i)*strain*(1-decay)
   else
    call fene_update(azz(i),arr(i),strain,h)
    if(i==n)then
     w=m0*h/injection_interval/(mp(i)+ms(i))
     azz(i)=(1-w)*azz(i)+w;arr(i)=(1-w)*arr(i)+w
    endif
    stress(i)=eta(i)/relaxation_time*(finite_extensibility-3) &
     /(finite_extensibility-azz(i)-2*arr(i))*(azz(i)-arr(i))
   endif
  enddo
  call forces()
  do i=1,n
   if(.not.landed(i))vel(:,i)=vel(:,i)+.5_dp*h*force(:,i)/(mp(i)+ms(i))
  enddo
  call require(all(ieee_is_finite(vel(:,:n))).and.all(ieee_is_finite(stress(:n))),'Nonfinite state')
 end subroutine
 subroutine outputs()
  integer::i
  real(dp)::retained,pmass,qactive,err,p_err,c_err,maxv,maxs,injected_count
  do i=1,n
   write(ut,'(es24.16e3,",",i0,14(",",es24.16e3),",",i0,8(",",es24.16e3))')t,i,pos(:,i),vel(:,i),radius(i), &
    mp(i)+ms(i),mp(i),ms(i),q(i),temp(i),stress(i),ell(i),merge(1,0,landed(i)), &
    travel(i),merge(land_time(i),t,landed(i))-birth(i),km_air(i),h_air(i),re_air(i),phi_surface(i),transfer_valid(i),flux_air(i)
  enddo
  retained=sum(mp(:n)+ms(:n));pmass=sum(mp(:n));qactive=sum(q(:n),mask=.not.landed(:n))
  injected_count=0
  if(n>0)injected_count=t/injection_interval
  err=abs(retained+evaporated-injected_count*m0)/max(injected_count*m0,1e-30_dp)
  p_err=abs(pmass-injected_count*mp0)/max(injected_count*mp0,1e-30_dp)
  c_err=abs(qactive+delivered-injected_count*q0)/max(abs(injected_count*q0),1e-30_dp)
  call require(max(err,p_err,c_err)<1e-8_dp,'Mass/polymer/charge balance failed')
  maxv=0;maxs=0
  if(n>0)then
   maxv=maxval(sqrt(sum(vel(:,:n)**2,dim=1)));maxs=maxval(abs(stress(:n)))
  endif
  write(ub,'(es24.16e3,2(",",i0),9(",",es24.16e3))')t,n,nland,injected_count*m0,retained,evaporated,err,p_err,c_err,dt,maxv,maxs
  flush(ut);flush(ud);flush(ub)
 end subroutine
 subroutine run()
  integer::u
  real(dp)::event_time,last_output
  open(newunit=ut,file='trajectory.csv',status='replace')
  write(ut,'(a)')'time_s,bead_id,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,radius_m,mass_kg,polymer_mass_kg,solvent_mass_kg,charge_C,temperature_K,stress_Pa,segment_length_m,deposited,path_length_m,resolved_age_s,mass_transfer_m_s,heat_transfer_W_m2_K,normal_Re,surface_solvent_volume_fraction,transfer_correlation_in_domain,evaporation_flux_kg_m2_s'
  open(newunit=ud,file='deposition.csv',status='replace')
  write(ud,'(a)')'time_s,bead_id,x_m,y_m,collector_x_m,collector_y_m,radius_m,dry_radius_m,polymer_mass_kg,solvent_mass_kg,charge_C,speed_m_s'
  open(newunit=ub,file='diagnostics.csv',status='replace')
  write(ub,'(a)')'time_s,beads,deposited,injected_mass_kg,retained_mass_kg,evaporated_mass_kg,mass_relative_error,polymer_relative_error,charge_relative_error,dt_s,max_speed_m_s,max_stress_Pa'
  call outputs();last_output=t;next_output=output_interval
  do while(t<duration*(1-1e-13_dp))
   event_time=min(duration,(n+1)*injection_interval)
   if(n==0)then
    t=min(event_time,next_output)
    if(abs(t-event_time)<1e-12_dp*duration)then
     call inject();call geometry();call forces()
    endif
   else
    call geometry();call forces();dt=stable_step()
    dt=min(dt,event_time-t,next_output-t)
    if(dt>1e-15_dp*duration)then
     steps=steps+1;call require(steps<=max_steps,'max_steps exhausted before requested duration')
     call advance(dt);t=t+dt
    else
     t=min(event_time,next_output)
    endif
    if(abs(t-event_time)<1e-12_dp*duration.and.t<duration*(1-1e-13_dp))then
     call inject();call geometry();call forces()
    endif
   endif
   if(t>=next_output-1e-12_dp*duration)then
    call outputs();last_output=t;next_output=next_output+output_interval
   endif
  enddo
  if(abs(last_output-t)>1e-12_dp*duration)call outputs()
  close(ut);close(ud);close(ub)
  open(newunit=u,file='summary.txt',status='replace')
  write(u,'(a)')'COMPLETED: reduced transient 3-D viscoelastic bead-chain model; not experimentally validated.'
  write(u,'(a,a)')'Constitutive law: ',trim(constitutive_model)
  write(u,'(a,a)')'Inlet perturbation: ',trim(perturbation_mode)
  write(u,'(a,a)')'Aerodynamic law: ',trim(drag_model)
  write(u,'(a,a)')'Heat/mass transfer model: ',trim(transport_model)
  write(u,'(a,a)')'Solvent activity model: ',trim(activity_model)
  write(u,'(a,i0)')'Steps: ',steps
  write(u,'(a,i0)')'Injected parcels: ',n
  write(u,'(a,i0)')'Deposited parcels: ',nland
  write(u,'(a,es16.8)')'Duration [s]: ',t
  write(u,'(a,es16.8)')'Inlet velocity [m/s]: ',u0
  write(u,'(a,es16.8)')'Initial segment length [m]: ',l0
  write(u,'(a,es16.8)')'Evaporated solvent mass [kg]: ',evaporated
  write(u,'(a,es16.8)')'Delivered charge [C]: ',delivered
  if(conductivity>0)write(u,'(a,es16.8)')'Charge relaxation screening time [s]: ',eps0*relative_permittivity/conductivity
  if(charge_transport)then
   write(u,'(a)')'Prescribed inlet TOTAL current; conservative longitudinal ohmic redistribution, grounded collector.'
   write(u,'(a)')'Regularized bead self-capacitance softening/ke; check injection mesh and softening sensitivity.'
  else
   write(u,'(a)')'Charge is prescribed CONVECTIVE current; conductivity/permittivity are diagnostic.'
  endif
  write(u,'(a)')'Fixed external field; no 3-D fluid/Poisson feedback. Softened pair charge forces, no liquid cone breakup/contact mechanics.'
  write(u,'(a)')'Activity and transfer laws are listed above; deposited parcels stop evaporating and drain their charge.'
  if(trim(transport_model)=='gas_film')write(u,'(a)')'Quasi-steady finite cylindrical gas film: outer radius is an uncertain physical boundary assumption.'
  if(trim(drag_model)=='effective')write(u,'(a)')'Effective quadratic drag coefficient requires calibration.'
  if(trim(drag_model)=='cylinder')write(u,'(a)')'Slender cylinder axial drag and empirical normal Cd(Re); not a resolved airflow solution.'
  write(u,'(a)')'Newest cell grows by continuous material/charge inflow and splits conservatively; inlet solvent is well mixed.'
  write(u,'(a)')'First parcel enters after injection_interval at z=u_in*injection_interval: unresolved inlet flight is excluded.'
  if(nr==0)write(u,'(a)')'WARNING: uniform electric field selected; no needle enhancement is represented.'
  if(nland==0)write(u,'(a)')'WARNING: no collector crossings during simulation; increase duration if trajectory remains admissible.'
  close(u)
  open(newunit=u,file='run_status.txt',status='replace');write(u,'(a)')'COMPLETED';close(u)
  write(*,'(a,i0,a,i0,a,i0)')'Completed 3-D integration: ',steps,' steps, ',n,' parcels, deposited ',nland
 end subroutine
 subroutine self_test()
  real(dp)::d(3),f(3),opposite(3),exact,actual,relaxed,volume0,a,length0,tension0,eps
  d=[2._dp,0._dp,0._dp];eps=1e-6_dp
  f=ke*1e-18_dp*d/(sum(d*d)+eps**2)**1.5_dp;opposite=-f
  call require(maxval(abs(f+opposite))<1e-30_dp,'Coulomb pair antisymmetry')
  exact=ke*1e-18_dp/4
  call require(abs(f(1)/exact-1)<1e-11_dp,'Coulomb far-field limit')
  d=[0._dp,0._dp,-.2_dp]
  actual=-ke*1e-18_dp*d(3)/sum(d*d)**1.5_dp;exact=ke*1e-18_dp/(4*.1_dp**2)
  call require(abs(actual/exact-1)<1e-12_dp,'Ground-plane self-image limit')
  relaxed=123._dp*exp(-.02_dp/.01_dp)
  call require(abs(relaxed-16.64623983810336_dp)<1e-12_dp,'Maxwell analytical relaxation')
  volume0=1e-12_dp;length0=.001_dp;a=sqrt(volume0/(pi*length0));tension0=pi*.055_dp*a
  exact=(2*pi*.055_dp*sqrt(volume0/(pi*(length0+1e-8_dp)))*(length0+1e-8_dp) &
   -2*pi*.055_dp*sqrt(volume0/(pi*(length0-1e-8_dp)))*(length0-1e-8_dp))/(2e-8_dp)
  call require(abs(tension0/exact-1)<1e-8_dp,'Capillary tension as surface-energy derivative')
  write(*,'(a)')'SELF-TEST PASSED: Coulomb pair/limit, plane self-image, Maxwell relaxation, capillary energy derivative.'
 end subroutine
end module
program electrospinning_3d
 use jet3d_model
 implicit none
 character(1024)::arg
 call get_command_argument(1,arg)
 if(trim(arg)=='--self-test')then
  call self_test()
 else
  call require(len_trim(arg)>0,'Usage: jet3d path/to/config.nml (cwd receives outputs)')
  call load_config(trim(arg));call run()
 endif
end program
