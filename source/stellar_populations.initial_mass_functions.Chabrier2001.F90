!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a stellar initial mass function class based on \cite{chabrier_galactic_2001}.
  
  !# <initialMassFunction name="initialMassFunctionChabrier2001">
  !#  <description>A stellar initial mass function class based on \cite{chabrier_galactic_2001}.</description>
  !# </initialMassFunction>
  type, extends(initialMassFunctionClass) :: initialMassFunctionChabrier2001
     !% A stellar initial mass function class based on \cite{chabrier_galactic_2001}.
     private
     double precision :: massLower               , massTransition        , &
          &              massUpper               , exponent              , &
          &              massCharacteristic      , sigma                 , &
          &              normalizationExponential, normalizationLogNormal       
   contains
     procedure :: massMinimum => chabrier2001MassMinimum
     procedure :: massMaximum => chabrier2001MassMaximum
     procedure :: phi         => chabrier2001Phi
     procedure :: tabulate    => chabrier2001Tabulate
     procedure :: label       => chabrier2001Label
  end type initialMassFunctionChabrier2001

  interface initialMassFunctionChabrier2001
     !% Constructors for the {\normalfont \ttfamily chabrier2001} initial mass function class.
     module procedure chabrier2001ConstructorParameters
     module procedure chabrier2001ConstructorInternal
  end interface initialMassFunctionChabrier2001

contains

  function chabrier2001ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily chabrier2001} initial mass function class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type            (initialMassFunctionChabrier2001)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: massLower         , massTransition, &
          &                                                             massUpper         , exponent      , &
          &                                                             massCharacteristic, sigma
    
    !# <inputParameter>
    !#   <name>massUpper</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>125.0d0</defaultValue>
    !#   <description>The upper mass limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massLower</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The lower mass limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massTransition</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The transition limit for the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigma</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.69d0</defaultValue>
    !#   <description>The width of the lognormal part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponent</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-2.3d0</defaultValue>
    !#   <description>The exponent of the power law part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massCharacteristic</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.08d0</defaultValue>
    !#   <description>Characteristic mass of the lognormal part of the \cite{chabrier_galactic_2001} \gls{imf}.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=initialMassFunctionChabrier2001(massLower,massTransition,massUpper,exponent,massCharacteristic,sigma)
    !# <inputParametersValidate source="parameters"/>
    return
  end function chabrier2001ConstructorParameters

  function chabrier2001ConstructorInternal(massLower,massTransition,massUpper,exponent,massCharacteristic,sigma) result(self)
    !% Internal constructor for the {\normalfont \ttfamily chabrier2001} initial mass function.
    use Error_Functions
    use Numerical_Constants_Math
    implicit none
    type            (initialMassFunctionChabrier2001)                :: self
    double precision                                 , intent(in   ) :: massLower         , massTransition, &
         &                                                              massUpper         , exponent      , &
         &                                                              massCharacteristic, sigma
    double precision                                                 :: normalization
    !# <constructorAssign variables="massLower,massTransition,massUpper,exponent,massCharacteristic,sigma"/>
    
    self%normalizationLogNormal  =+sqrt(Pi/2.0d0)                                     &
         &                        *self%sigma                                         &
         &                        *self%massCharacteristic                            &
         &                        *log(10.0d0)                                        &
         &                        *exp(                                               &
         &                             +0.5d0                                         &
         &                             *self%sigma**2                                 &
         &                             *log(10.0d0)     **2                           &
         &                            )                                               &
         &                        *(                                                  &
         &                          -Error_Function(                                  &
         &                                          +self%sigma                       &
         &                                          *log (10.0d0)                     &
         &                                          /sqrt( 2.0d0)                     &
         &                                          -log10(                           &
         &                                                 +self%massTransition       &
         &                                                 /self%massCharacteristic   &
         &                                                )                           &
         &                                          /sqrt( 2.0d0)                     &
         &                                          /self%sigma                       &
         &                                         )                                  &
         &                          +Error_Function(                                  &
         &                                          +self%sigma                       &
         &                                          *log (10.0d0)                     &
         &                                          /sqrt( 2.0d0)                     &
         &                                          -log10(                           &
         &                                                 +self%massLower            &
         &                                                 /self%massCharacteristic   &
         &                                                )                           &
         &                                          /sqrt( 2.0d0)                     &
         &                                          /self%sigma                       &
         &                                         )                                  &
         &                         )
    self%normalizationExponential=+exp(                                               &
         &                              -0.50d0                                       &
         &                              *log10(                                       &
         &                                     +self%massTransition                   &
         &                                     /self%massCharacteristic               &
         &                                    )                        ** 2               &
         &                              /self%sigma                    ** 2               &
         &                             )                                                  &
         &                         *self%massTransition                                   &
         &                         /                                     (2.0d0+exponent) &
         &                         *(                                                     &
         &                           +(                                                   &
         &                             +self%massUpper                                    &
         &                             /self%massTransition                               &
         &                            )                                **(2.0d0+exponent) &
         &                           -1.0d0                                               &
         &                          )
    normalization                =+self%normalizationLogNormal   &
         &                        +self%normalizationExponential
    self%normalizationLogNormal  =1.0d0/normalization
    self%normalizationExponential=+exp(                                    &
         &                              -0.50d0                            &
         &                              *log10(                            &
         &                                     +self%massTransition        &
         &                                     /self%massCharacteristic    &
         &                                    )                        **2 &
         &                              /self%sigma                    **2 &
         &                             )                                   &
         &                         /self%massTransition                    &
         &                         /normalization
    return
  end function chabrier2001ConstructorInternal

  double precision function chabrier2001MassMinimum(self)
    !% Return the minimum mass of stars in the \cite{chabrier_galactic_2001} \gls{imf}.
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self

    chabrier2001MassMinimum=self%massLower
    return
  end function chabrier2001MassMinimum

  double precision function chabrier2001MassMaximum(self)
    !% Return the maximum mass of stars in the \cite{chabrier_galactic_2001} \gls{imf}.
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self

    chabrier2001MassMaximum=self%massUpper
    return
  end function chabrier2001MassMaximum
  
  double precision function chabrier2001Phi(self,massInitial)
    !% Evaluate the \cite{chabrier_galactic_2001} stellar initial mass function.
    implicit none
    class           (initialMassFunctionChabrier2001), intent(inout) :: self
    double precision                                 , intent(in   ) :: massInitial

    if      (                                    &
         &    massInitial >= self%massLower      &
         &   .and.                               &
         &    massInitial <  self%massTransition &
         &  ) then
       chabrier2001Phi=+self%normalizationLogNormal           &
            &          *exp(                                  &
            &               -0.5d0                            &
            &               *(                                &
            &                 +log10(                         &
            &                        +     massInitial        &
            &                        /self%massCharacteristic &
            &                       )                         &
            &                 /self%sigma                     &
            &                )**2                             &
            &              )                                  &
            &          /massInitial
    else if (                                    &
         &    massInitial >= self%massTransition &
         &   .and.                               &
         &    massInitial <  self%massUpper      &
         &  ) then
       chabrier2001Phi=+self%normalizationExponential &
            &          *massInitial**self%exponent
    else
       chabrier2001Phi=0.0d0
    end if
    return
  end function chabrier2001Phi
  
  subroutine chabrier2001Tabulate(self,imfTable)
    !% Construct and return a tabulation of the \cite{chabrier_galactic_2001} \gls{imf}.
    implicit none
    class  (initialMassFunctionChabrier2001)             , intent(inout) :: self
    class  (table1D                        ), allocatable, intent(inout) :: imfTable
    integer                                 , parameter                  :: countTable=100
    integer                                                              :: i

    allocate(table1DLogarithmicLinear :: imfTable)
    select type (imfTable)
    type is (table1DLogarithmicLinear)
       call imfTable%create(                 &
            &               self%massLower , &
            &               self%massUpper , &
            &                    countTable  &
            &              )
       do i=1,countTable
          call imfTable%populate(self%phi(imfTable%x(i)),i)
       end do
    end select
    return
  end subroutine chabrier2001Tabulate

  function chabrier2001Label(self)
    !% Return a label for this \gls{imf}.
    implicit none
    class(initialMassFunctionChabrier2001), intent(inout) :: self
    type (varying_string                    )             :: chabrier2001Label
    !GCC$ attributes unused :: self

    chabrier2001Label="Chabrier2001"
    return
  end function chabrier2001Label