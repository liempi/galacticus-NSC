!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implementation of the \cite{antonini_coevolution_2015} star formation rate law for galactic \glspl{nsc}.
  !!}
  use :: Cosmology_Functions           , only : cosmologyFunctionsClass
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass

  !![
  <nuclearStarClusterGrowthRates name="nuclearStarClusterGrowthRatesAntonini2015Enhanced">
   <description>
    A gas inflow rate implementing the model of \citep{antonini_coevolution_2015} for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{gas}^\mathrm{NSC} = A_\mathrm{res}(1+z)^{\alpha}\dot{M}_\star^\mathrm{spheroid},
    \end{equation}    
    where $A_\mathrm{res}=${\normalfont \ttfamily [efficiency]} is a free parameter that controls the amount of gas that inflows 
    into the center, $\alpha$={\normalfont \ttfamily [exponent]} is the exponent of the $(1+z)$ term that enhance the gas inflow at early redshifts,
    and $\dot{M}_\star^\mathrm{spheroid}$ is the star formation rate of the spheroid component.
   </description>
  </nuclearStarClusterGrowthRates>
  !!]
  type, extends(nuclearStarClusterGrowthRatesClass) :: nuclearStarClusterGrowthRatesAntonini2015Enhanced
     !!{
     Implementation of the \cite{antonini_coevolution_2015} gas inflow rate with an enhance term for galactic \glspl{nsc}.
     !!}
     private
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_         => null()
     class           (starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     double precision                                           :: efficiency                           , exponent
     contains
     final     ::          antonini2015EnhancedDestructor
     procedure :: rate  => antonini2015EnhancedRate
  end type nuclearStarClusterGrowthRatesAntonini2015Enhanced

  interface nuclearStarClusterGrowthRatesAntonini2015Enhanced
     !!{
     Constructors for the \refClass{nuclearStarClusterGrowthRatesAntonini2015Enhanced} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure antonini2015EnhancedConstructorParameters
     module procedure antonini2015EnhancedConstructorInternal
  end interface nuclearStarClusterGrowthRatesAntonini2015Enhanced
    
contains

  function antonini2015EnhancedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015Enhanced} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nuclearStarClusterGrowthRatesAntonini2015Enhanced)                :: self
    type            (inputParameters                                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                          ), pointer       :: cosmologyFunctions_
    class           (starFormationRateSpheroidsClass                  ), pointer       :: starFormationRateSpheroids_
    double precision                                                                   :: efficiency                 , exponent

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultSource>\citep{antonini_coevolution_2015}</defaultSource>
      <defaultValue>1.0d-2</defaultValue>
      <description>Parameter controlling the rate of the gas inflow onto the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <defaultValue>2.0d0</defaultValue>
      <description>Value of the exponent in the power law that enhance the gas inflow onto the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"         name="cosmologyFunctions_"         source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=nuclearStarClusterGrowthRatesAntonini2015Enhanced(efficiency,exponent,cosmologyFunctions_,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"        />
    <objectDestructor name="cosmologyFunctions_"        />
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function antonini2015EnhancedConstructorParameters

  function antonini2015EnhancedConstructorInternal(efficiency,exponent,cosmologyFunctions_,starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015Enhanced} gas inflow rate from NSCs class.
    !!}
    implicit none
    type            (nuclearStarClusterGrowthRatesAntonini2015Enhanced)                        :: self
    class           (cosmologyFunctionsClass                          ), intent(in   ), target :: cosmologyFunctions_
    class           (starFormationRateSpheroidsClass                  ), intent(in   ), target :: starFormationRateSpheroids_
    double precision                                                   , intent(in   )         :: efficiency
    double precision                                                   , intent(in   )         :: exponent

    !![
    <constructorAssign variables="efficiency, exponent, *cosmologyFunctions_, *starFormationRateSpheroids_"/>
    !!]
    return
  end function antonini2015EnhancedConstructorInternal

  subroutine antonini2015EnhancedDestructor(self)
    !!{
    Destructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015Enhanced} class
    !!}
    implicit none
    type(nuclearStarClusterGrowthRatesAntonini2015Enhanced), intent(inout) :: self
    !![
    <objectDestructor name="self%cosmologyFunctions_"        />
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]
    return
  end subroutine antonini2015EnhancedDestructor

  double precision function antonini2015EnhancedRate(self,node,time) result(rate)
    !!{
    Returns the gas inflow rate (in $M_\odot$ Gyr$^{-1}$) onto the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The \gls{nsc} is assumed to obey the \cite{antonini_coevolution_2015} gas inflow rate model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentSpheroid
    implicit none
    class           (nuclearStarClusterGrowthRatesAntonini2015Enhanced), intent(inout), target  :: self
    type            (treeNode                                         ), intent(inout)          :: node
    double precision                                                   , intent(in   )          :: time    
    class           (nodeComponentSpheroid                            ),                pointer :: spheroid
    double precision                                                                            :: rateStarFormationSpheroid, redshift, &
         &                                                                                         expansionFactor

    ! Get the spheroid component.
    spheroid                  => node                            %spheroid(    )
    ! Get the star formation rate of the spheroid component.
    rateStarFormationSpheroid =  self%starFormationRateSpheroids_%rate    (node)

    expansionFactor =+self%cosmologyFunctions_%expansionFactor            (time           )
    redshift        =+self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)  
    ! Find the rate of mass accretion onto the nuclear star cluster.
    if (rateStarFormationSpheroid <= 0.0d0) then
      rate    =+0.0d0
    else
      ! Gas accretion rate model from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
       rate   =+self%efficiency                &
            &  *(1+redshift)**self%exponent    &
            &  *rateStarFormationSpheroid
    end if 
    return
  end function antonini2015EnhancedRate

  
