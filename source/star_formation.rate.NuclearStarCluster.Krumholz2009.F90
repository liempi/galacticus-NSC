!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !!{
  Implementation of the \cite{krumholz_star_2009} star formation rate law for galactic NSCs.
  !!}
  use :: Abundances_Structure, only : abundances

  !![
  <starFormationRateNSC name="starFormationRateNSCKrumholz2009">
   <description>
    A star formation rate surface density class implementing the model of \citep{krumholz_star_2009}:
    \begin{equation}
     \dot{\Sigma}_\star(R) = \nu_\mathrm{SF} f_\mathrm{H_2}(R)\Sigma_\mathrm{HI, NSC}(R) \left\{ \begin{array}{ll}
     (\Sigma_\mathrm{HI}/\Sigma_0)^{-1/3}, &amp; \hbox{ if } \Sigma_\mathrm{HI}/\Sigma_0 \le 1 \\
     (\Sigma_\mathrm{HI}/\Sigma_0)^{1/3}, &amp; \hbox{ if } \Sigma_\mathrm{HI}/\Sigma_0 &gt; 1 \end{array} \right. ,
    \end{equation}
    where $\nu_\mathrm{SF}=${\normalfont \ttfamily [frequencyStarFormation]} is a frequency and $\Sigma_0=85 M_\odot
    \hbox{pc}^{-2}$. The molecular fraction is given by
    \begin{equation}
     f_\mathrm{H_2} = 1 - \left( 1 + \left[ { 3 s \over 4 (1+\delta)} \right]^{-5} \right)^{-1/5},
    \end{equation}
    where
    \begin{equation}
     \delta = 0.0712 \left[ 0.1 s^{-1} + 0.675 \right]^{-2.8},
    \end{equation}
    and
    \begin{equation}
     s = {\ln(1+0.6\chi+0.01\chi^2) \over 0.04 \Sigma_\mathrm{comp,0} Z^\prime},
    \end{equation}
    with
    \begin{equation}
     \chi = 0.77 \left[ 1 + 3.1 Z^{\prime 0.365} \right],
    \end{equation}
    and $\Sigma_\mathrm{comp,0}=c \Sigma_\mathrm{HI}/M_\odot \hbox{pc}^{-2}$ where $c=${\normalfont \ttfamily
    [clumpingFactorMolecularComplex]} is a density enhancement factor relating the surface density of molecular complexes to
    the gas density on larger scales. Alternatively, if {\normalfont \ttfamily [molecularFractionFast]} is set to true, the
    molecular fraction will be computed using the faster (but less accurate at low molecular fraction) formula
    \begin{equation}
     f_\mathrm{H_2} = 1 - { 3s/4 \over (1 + s/4)}.
    \end{equation}
   </description>
  </starFormationRateNSC>
  !!]
  type, extends(starFormationRateNSCClass) :: starFormationRateNSCKrumholz2009
     !!{
     Implementation of the \cite{krumholz_star_2009} star formation rate surface density law for galactic NSCs.
     !!}
     private
     double precision                                          ::  metallicityRelativeToSolar, frequencyStarFormation, &
          &                                                        s                         , chi                            
     contains
     procedure :: rate                  => krumholz2009Rate
  end type starFormationRateNSCKrumholz2009

  interface starFormationRateNSCKrumholz2009
     !!{
     Constructors for the {\normalfont \ttfamily krumholz2009} star formation surface density rate in NSCs class.
     !!}
     module procedure krumholz2009ConstructorParameters
     module procedure krumholz2009ConstructorInternal
  end interface starFormationRateNSCKrumholz2009
    
contains

  function krumholz2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily krumholz2009} star formation surface density rate in NSCs class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateNSCKrumholz2009)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: frequencyStarFormation         

    !![
    <inputParameter>
      <name>frequencyStarFormation</name>
      <defaultSource>\citep{krumholz_star_2009}</defaultSource>
      <defaultValue>2.36d0</defaultValue>
      <description>The star formation frequency (in units of Gyr) in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateNSCKrumholz2009(frequencyStarFormation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function krumholz2009ConstructorParameters

  function krumholz2009ConstructorInternal(frequencyStarFormation) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily krumholz2009} star formation surface density rate from NSCs class.
    !!}
    implicit none
    type            (starFormationRateNSCKrumholz2009)                 :: self
    double precision                                  , intent(in   )  :: frequencyStarFormation
    !![
    <constructorAssign variables="frequencyStarFormation"/>
    !!]
    return
  end function krumholz2009ConstructorInternal

  double precision function krumholz2009Rate(self,node)
    !!{
    Returns the star formation rate (in $M_\odot$ Gyr$^{-1}$) for star formation
    in the galactic NSC of {\normalfont \ttfamily node}. The NSC is assumed to obey the
    \cite{krumholz_star_2009} star formation rule.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentNSC                , treeNode
    use :: Abundances_Structure, only : metallicityTypeLinearByMassSolar
    implicit none
    class           (starFormationRateNSCKrumholz2009), intent(inout), target :: self
    type            (treeNode                        ), intent(inout)         :: node
    class           (nodeComponentNSC                ), pointer               :: NSC
    type            (abundances                      ), save                  :: abundancesFuel
    double precision                                                          :: molecularFractiona, radiusNSC, &
         &                                                                       massGas           , t_SF                           
    double precision                                  , parameter             :: Sigma_th = 85.0d0
    double precision                                                          :: Sigma_res, Sigma_1
                        
    !$omp threadprivate(abundancesFuel)


    NSC       => node%NSC    ()
    massGas   =  NSC %massGas()
    radiusNSC =  NSC %radius ()*10d6 !pc

    if     (                                             &
         &   massGas                            <= 0.0d0 &
         &  .or.                                         &
         &   radiusNSC                          <= 0.0d0 &
         & ) then
       ! It is not, so return zero rate.
       krumholz2009Rate=0.0d0
       return
    else 
       Sigma_res = SurfaceDensityGas(radiusNSC,massGas)
       Sigma_1   = Sigma_res/1.0d0

       ! Find the hydrogen fraction in the NSC gas of the fuel supply.
       abundancesFuel=NSC%abundancesGas()
       
       call abundancesFuel%massToMassFraction(massGas)

       ! Get the metallicity in Solar units, and related quantities.
       self%metallicityRelativeToSolar=abundancesFuel%metallicity(metallicityTypeLinearByMassSolar)
       if (self%metallicityRelativeToSolar /= 0.0d0) then
           self%chi                        = 0.77d0*(1.0d0+3.1d0*self%metallicityRelativeToSolar**0.365d0)
           self%s                          = log(1.0d0+0.6d0*self%chi)/(0.04d0*self%metallicityRelativeToSolar*Sigma_1)
       else
            krumholz2009Rate=0.0d0
       end if 
 
       if (Sigma_res > Sigma_th ) then 
           t_SF = (1.0d0/self%frequencyStarFormation)*(Sigma_res/Sigma_th)**(-0.33d0)  !t_sf^-1
       else
           t_SF = (1.0d0/self%frequencyStarFormation)*(Sigma_res/Sigma_th)**0.34d0 !t_sf^-1
       end if 
      ! Compute the molecular fraction.
       molecularFractiona  = MolecularFraction(self%s)
       !Compute the star formation rate density.
       krumholz2009Rate= +massGas          &
            &            *molecularFractiona*t_SF
    end if
    return
  end function krumholz2009Rate
  
  double precision function SurfaceDensityGas(radius,massGas)
    !!{
    Compute surface density of the NSC.
    !!}
    use :: Numerical_Constants_Math , only : Pi
    implicit none
    double precision                , intent(in   ) :: radius, massGas
    ! Get gas surface density in units of M☉/Mpc²
    SurfaceDensityGas = massGas / (2.0d0*Pi*radius**2.0)
    return
  end function SurfaceDensityGas

  double precision function MolecularFraction(s)
    !!{
    Slow (but more accurate at low molecular fraction) fitting function from \cite{krumholz_star_2009} for the molecular
    hydrogen fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: s
    double precision, parameter     :: sTiny        =1.000000d-06
    double precision, parameter     :: sHuge        =1.000000d+10
    double precision, parameter     :: deltaInfinity=0.214008d+00 ! The value of δ for s → ∞.
    double precision, parameter     :: sMaximum     =10.0d0 
    double precision                :: delta

    if      (s <  sTiny   ) then
       ! Series expansion for very small s.
       MolecularFraction=1.0d0-0.75d0*s
    else if (s >= sHuge   ) then
       ! Truncate to zero for extremely large s.
       MolecularFraction=0.0d0
    else if (s >= sMaximum) then
       ! Simplified form for very large s.
       MolecularFraction=1.0d0/(0.75d0/(1.0d0+deltaInfinity))**5/5.0d0/s**5
    else
       ! Full expression.
       delta            =0.0712d0/((0.1d0/s+0.675d0)**2.8d0)
       MolecularFraction=1.0d0-1.0d0/((1.0d0+(((1.0d0+delta)/0.75d0/s)**5))**0.2d0)
    end if

    if (MolecularFraction > 0.02) then
        return
    else 
       MolecularFraction = 0.02
    end if 
    return
  end function MolecularFraction 