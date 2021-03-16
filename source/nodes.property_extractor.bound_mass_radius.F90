!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a property extractor class that extracts the radius enclosing the current bound mass.

  !# <nodePropertyExtractor name="nodePropertyExtractorRadiusBoundMass">
  !#  <description>
  !#   A property extractor class that extracts the radius enclosing the current bound mass.
  !#  </description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusBoundMass
     !% A property extractor class that extracts the radius enclosing the current bound mass.
     private
   contains
     procedure :: extract     => radiusBoundMassExtract
     procedure :: name        => radiusBoundMassName
     procedure :: description => radiusBoundMassDescription
     procedure :: unitsInSI   => radiusBoundMassUnitsInSI
     procedure :: type        => radiusBoundMassType
  end type nodePropertyExtractorRadiusBoundMass

  interface nodePropertyExtractorRadiusBoundMass
     !% Constructors for the ``radiusBoundMass'' output analysis class.
     module procedure radiusBoundMassConstructorParameters
  end interface nodePropertyExtractorRadiusBoundMass

contains

  function radiusBoundMassConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily radiusBoundMass} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorRadiusBoundMass)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters
    
    self=nodePropertyExtractorRadiusBoundMass()
    return
  end function radiusBoundMassConstructorParameters

  double precision function radiusBoundMassExtract(self,node,instance)
    !% Implement a bound mass radius property extractor.
    use :: Galacticus_Nodes                  , only : nodeComponentSatellite
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Radius_Enclosing_Mass
    implicit none
    class(nodePropertyExtractorRadiusBoundMass), intent(inout)           :: self
    type (treeNode                            ), intent(inout), target   :: node
    type (multiCounter                        ), intent(inout), optional :: instance
    class(nodeComponentSatellite              )               , pointer  :: satellite
    !$GLC attributes unused :: instance

    satellite              => node%satellite()
    radiusBoundMassExtract =  Galactic_Structure_Radius_Enclosing_Mass(node,mass=satellite%boundMass())
    return
  end function radiusBoundMassExtract

  function radiusBoundMassName(self)
    !% Return the name of the bound mass radius property.
    implicit none
    type (varying_string                      )                :: radiusBoundMassName
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassName=var_str('satelliteRadiusBoundMass')
    return
  end function radiusBoundMassName

  function radiusBoundMassDescription(self)
    !% Return a description of the bound mass radius property.
    implicit none
    type (varying_string                      )                :: radiusBoundMassDescription
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassDescription=var_str('Radius enclosing the bound mass of the halo [Mpc].')
    return
  end function radiusBoundMassDescription

  double precision function radiusBoundMassUnitsInSI(self)
    !% Return the units of the bound mass radius property in the SI system.
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassUnitsInSI=megaParsec
    return
  end function radiusBoundMassUnitsInSI

  integer function radiusBoundMassType(self)
    !% Return the type of the bound mass radius property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassType=outputAnalysisPropertyTypeLinear
    return
  end function radiusBoundMassType
