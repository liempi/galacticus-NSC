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

!!{
Contains custom functions for the standard dark core component.
!!}

double precision function Node_Component_Dark_Core_Standard_Velocity(self)
   !!{
   Return the circular velocity of the standard dark Core.
   !!}
   use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
   implicit none
   class(nodeComponentDarkCoreStandard), intent(inout) :: self
   double precision                                    :: radiusDarkCore, massDarkCore

   radiusDarkCore= self%radius()
   massDarkCore  = self%massStellar()

   !Units corrected, G in galacticus is in units of Mpc M⊙⁻¹km²s⁻².
   !The velocity is returned in km/s units
   if (massDarkCore > 0.0d0 .and. radiusDarkCore >0.0d0) then
      Node_Component_Dark_Core_Standard_Velocity = sqrt(gravitationalConstant_internal*massDarkCore/radiusDarkCore)
   else 
      Node_Component_Dark_Core_Standard_Velocity = 0.0
   end if 
   return
end function Node_Component_Dark_Core_Standard_Velocity

double precision function Node_Component_Dark_Core_Standard_Angular_Momentum(self)
   !!{
   Return the specific angular momentum of the standard dark core at a given radius.
   !!}
   implicit none
   class (nodeComponentDarkCoreStandard), intent(inout) :: self
   double precision                                     :: massDarkCore, velocityDarkCore
   
   !Get the velocity and the dunamical mass of the dark core.
   velocityDarkCore = self%velocity()
   massDarkCore     = self%massStellar()

   if (massDarkCore > 0.0d0) then
       Node_Component_Dark_Core_Standard_Angular_Momentum = massDarkCore * velocityDarkCore
   else
       Node_Component_Dark_Core_Standard_Angular_Momentum = 0.0d0
   end if
   return
end function Node_Component_Dark_Core_Standard_Angular_Momentum

function Node_Component_Dark_Core_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
  !!{
  Return the mass distribution for the standard nuclear star cluster component.
  !!}
  use :: Mass_Distributions                    , only : massDistributionClass    , massDistributionSphericalScaler, massDistributionComposite   , massDistributionList   , &
        &                                               massDistributionSpherical, massDistributionMatches_
  use :: Error                                 , only : Error_Report
  use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_ , kinematicDistribution_
  use :: Galactic_Structure_Options            , only : componentTypeDarkCore    , massTypeStellar                , massTypeGaseous             , enumerationWeightByType, &
       &                                                weightByMass             , weightByLuminosity             , enumerationComponentTypeType, enumerationMassTypeType
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_
  class           (nodeComponentDarkCoreStandard  ), intent(inout)           :: self
  type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
  type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
  type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
  integer                                          , intent(in   ), optional :: weightIndex
  type            (massDistributionSphericalScaler), pointer                 :: massDistributionStellar   
  type            (stellarLuminosities            ), save                    :: luminosities
  !$omp threadprivate(luminosities)
  double precision                                                           :: massStellar               , radiusScale
  logical                                                                    :: includeStars
  !![
  <optionalArgument name="weightBy" defaultsTo="weightByMass" />
  !!]
  
  ! Determine which components of the spheroid to include.
  includeStars=massDistributionMatches_(componentTypeDarkCore,massTypeStellar,componentType,massType) .and. (weightBy_ == weightByMass .or. weightBy_ == weightByLuminosity)  
  ! Get properties of the mass distribution and ensure they are physical.
  if      (weightBy_ == weightByMass      ) then
     massStellar       =  max (0.0d0,self        %massStellar        (           ))
  else if (weightBy_ == weightByLuminosity) then
     call Error_Report('Dark core mass distribution can not be weighted by Luminosity'//{introspection:location})
  else
     massDistribution_ => null()
     return
  end if
  ! Determine which components to build.
  radiusScale=self%radius()
  if (radiusScale <= 0.0d0 .or. .not.includeStars) then
     ! Disk has non-positive size, or no components matched. Return a null distribution.
     massDistribution_ => null()
  else
     ! Build the individual distributions.
     massDistributionStellar => null()
     if (includeStars) then
        allocate(massDistributionStellar)
        select type (massDistributionStellar_)
        class is (massDistributionSpherical)
           !![
           <referenceConstruct object="massDistributionStellar" constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_)"/>
           !!]
        end select
        call massDistributionStellar%setKinematicsDistribution(kinematicDistribution_)
     end if
     ! Combine the distributions as necessary.
     if (includeStars) then
        ! Return just the stellar component.
        massDistribution_ => massDistributionStellar
     end if
  end if
  return
end function Node_Component_Dark_Core_Standard_Mass_Distribution

double precision function Node_Component_Dark_Core_Standard_Mass_Baryonic(self) result(massBaryonic)
  !!{
  Return the baryonic mass for the standard dark core component.
  !!}
  implicit none
  class(nodeComponentDarkCoreStandard), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar())
  return
end function Node_Component_Dark_Core_Standard_Mass_Baryonic

double precision function Node_Component_Dark_Core_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{
  Return the half-mass radius of the standard nuclear star cluster.
  !!}
  use :: Error                                 , only : Error_Report
  use :: Mass_Distributions                    , only : massDistributionSpherical
  use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentDarkCoreStandard), intent(inout) :: self

 select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=0.0d0
     call Error_Report('Dark Core mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Dark_Core_Standard_Half_Mass_Radius
