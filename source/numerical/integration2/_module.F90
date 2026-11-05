!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements a variety of numerical integrators.
!!}

module Numerical_Integration2
  !!{RST
  Implements a variety of numerical integrators.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  private
  public :: integrand1D                                     , integrandVectorized1D                          , &
       &    integratorCompositeTrapezoidal1D                , integratorVectorizedCompositeTrapezoidal1D     , &
       &    integratorAdaptiveCompositeTrapezoidal1D        , integratorCompositeGaussKronrod1D              , &
       &    integratorVectorizedCompositeGaussKronrod1D     , integrator2                                    , &
       &    integrator1D                                    , integratorVectorized1D                         , &
       &    integratorMulti1D                               , integratorMultiVectorized1D                    , &
       &    integratorMultiVectorizedCompositeGaussKronrod1D, integrandMulti1D                               , &
       &    integratorMulti                                 , integratorMultiVectorizedCompositeTrapezoidal1D

  ! Interval types.
  type :: interval
     private
     integer                             :: depth
     double precision                    :: a              , b       , &
          &                                 fa             , fb      , &
          &                                 error          , integral
     type            (interval), pointer :: next  => null()
  end type interval

  type :: intervalMulti
     private
     integer                                                    :: depth
     double precision                                           ::  a             ,  b
     double precision               , dimension(:), allocatable :: fa             , fb       , &
          &                                                        error          , integral
     type            (intervalMulti), pointer                   :: next  => null()
  end type intervalMulti

  type :: intervalMultiList
     type(intervalMulti), pointer :: interval_ => null()
  end type intervalMultiList

  ! Generic integrator.
  type :: integrator2
     !!{RST
     Generic numerical integrator class.
     !!}
     double precision :: toleranceAbsolute, toleranceRelative
   contains
     !![
     <methods docformat="rst">
       <method description="Set tolerances to use in this integrator." method="toleranceSet" />
     </methods>
     !!]
     procedure :: toleranceSet => toleranceSetGeneric
  end type integrator2

  ! Generic one-dimensional integrator.
  type, abstract, extends(integrator2) :: integrator1D
     !!{RST
     Generic one-dimensional numerical integrator class.
     !!}
     private
     procedure       (integrand1D), pointer, nopass :: integrand
   contains
     !![
     <methods docformat="rst">
       <method description="Set the integrand function to be integrated." method="integrandSet" />
       <method description="Evaluate the integral." method="evaluate" />
     </methods>
     !!]
     procedure                       :: integrandSet => integrandSet1D
     procedure(evaluate1D), deferred :: evaluate
  end type integrator1D
  abstract interface
     double precision function evaluate1D(self,a,b)
       import integrator1D
       class           (integrator1D), intent(inout) :: self
       double precision              , intent(in   ) :: a,b
     end function evaluate1D
  end interface
  abstract interface
     double precision function integrand1D(x)
       double precision, intent(in   ) :: x
     end function integrand1D
  end interface

  ! Composite trapezoidal 1D integrator.
  type, extends(integrator1D) :: integratorCompositeTrapezoidal1D
     !!{RST
     One-dimensional numerical integrator class using a composite trapezoidal rule.
     !!}
     private
     integer :: iterationsMaximum
   contains
     procedure :: evaluate   => compositeTrapezoidalEvaluate1D
  end type integratorCompositeTrapezoidal1D

  interface integratorCompositeTrapezoidal1D
     !!{RST
     Constructors for the :galacticus-class:`integratorCompositeTrapezoidal1D` numerical integrator class.
     !!}
     module procedure compositeTrapezoidalConstructor1D
  end interface integratorCompositeTrapezoidal1D

  ! Adaptive-composite trapezoidal 1D integrator.
  type, extends(integratorCompositeTrapezoidal1D) :: integratorAdaptiveCompositeTrapezoidal1D
     !!{RST
     One-dimensional numerical integrator class using an adaptive composite trapezoidal rule.
     !!}
     private
   contains
     procedure :: evaluate => adaptiveCompositeTrapezoidalEvaluate1D
  end type integratorAdaptiveCompositeTrapezoidal1D

  interface integratorAdaptiveCompositeTrapezoidal1D
     !!{RST
     Constructors for the :galacticus-class:`integratorAdaptiveCompositeTrapezoidal1D` numerical integrator class.
     !!}
     module procedure adaptiveCompositeTrapezoidalConstructor1D
  end interface integratorAdaptiveCompositeTrapezoidal1D

  ! Composite Gauss-Kronrod 1D integrator.
  type, extends(integrator1D) :: integratorCompositeGaussKronrod1D
     !!{RST
     One-dimensional numerical integrator class using a composite Gauss-Kronrod rule.
     !!}
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: xKronrod         , wGauss, wKronrod
   contains
     !![
     <methods docformat="rst">
       <method description="Evaluate the integral over an interval and also return the error on the integral." method="evaluateInterval" />
     </methods>
     !!]
     procedure :: evaluate         => compositeGaussKronrod1DEvaluate
     procedure :: evaluateInterval => compositeGaussKronrod1DEvaluateInterval
  end type integratorCompositeGaussKronrod1D

  interface integratorCompositeGaussKronrod1D
     !!{RST
     Constructors for the :galacticus-class:`integratorCompositeGaussKronrod1D` numerical integrator class.
     !!}
     module procedure compositeGaussKronrod1DConstructor
  end interface integratorCompositeGaussKronrod1D

  ! Generic one-dimensional vectorized integrator.
  type, abstract, extends(integrator2) :: integratorVectorized1D
     !!{RST
     Generic one-dimensional vectorized numerical integrator class.
     !!}
     private
     procedure       (integrandVectorized1D), pointer, nopass :: integrand
   contains
     !![
     <methods docformat="rst">
       <method description="Set the integrand function to be integrated." method="integrandSet" />
     </methods>
     !!]
     procedure                                 :: integrandSet => integrandVectorizedSet1D
     procedure(evaluateVectorized1D), deferred :: evaluate
  end type integratorVectorized1D
  abstract interface
     double precision function evaluateVectorized1D(self,a,b)
       import integratorVectorized1D
       class           (integratorVectorized1D), intent(inout) :: self
       double precision                        , intent(in   ) :: a,b
     end function evaluateVectorized1D
  end interface
  abstract interface
     function integrandVectorized1D(x)
       double precision, intent(in   ), dimension(     : ) :: x
       double precision               , dimension(size(x)) :: integrandVectorized1D
     end function integrandVectorized1D
  end interface

  ! Vectorized composite trapezoidal 1D integrator.
  type, extends(integratorVectorized1D) :: integratorVectorizedCompositeTrapezoidal1D
     !!{RST
     One-dimensional numerical integrator class using a vectorized composite trapezoidal rule.
     !!}
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: d
   contains
     !![
     <methods docformat="rst">
       <method description="Evaluate the integral." method="evaluate" />
     </methods>
     !!]
     procedure :: evaluate   => vectorizedCompositeTrapezoidalEvaluate1D
  end type integratorVectorizedCompositeTrapezoidal1D

  interface integratorVectorizedCompositeTrapezoidal1D
     !!{RST
     Constructors for the :galacticus-class:`integratorVectorizedCompositeTrapezoidal1D` numerical integrator class.
     !!}
     module procedure vectorizedCompositeTrapezoidalConstructor1D
  end interface integratorVectorizedCompositeTrapezoidal1D

  ! Vectorized composite Gauss-Kronrod 1D integrator.
  type, extends(integratorVectorized1D) :: integratorVectorizedCompositeGaussKronrod1D
     !!{RST
     One-dimensional numerical integrator class using a vectorized composite Gauss-Kronrod rule.
     !!}
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: xKronrod         , wGauss, wKronrod
   contains
     !![
     <methods docformat="rst">
       <method description="Evaluate the integral." method="evaluate" />
       <method description="Evaluate the integral over an interval and also return the error on the integral." method="evaluateInterval" />
     </methods>
     !!]
     procedure :: evaluate         => vectorizedCompositeGaussKronrod1DEvaluate
     procedure :: evaluateInterval => vectorizedCompositeGaussKronrod1DEvaluateInterval
  end type integratorVectorizedCompositeGaussKronrod1D

  interface integratorVectorizedCompositeGaussKronrod1D
     !!{RST
     Constructors for the :galacticus-class:`integratorVectorizedCompositeGaussKronrod1D` numerical integrator class.
     !!}
     module procedure vectorizedCompositeGaussKronrod1DConstructor
  end interface integratorVectorizedCompositeGaussKronrod1D

  ! Generic multi-integrand integrator.
  type :: integratorMulti
     !!{RST
     Generic numerical integrator class.
     !!}
     double precision, allocatable, dimension(:) :: toleranceAbsolute, toleranceRelative
   contains
     !![
     <methods docformat="rst">
       <method description="Set tolerances to use in this integrator." method="tolerancesSet" />
     </methods>
     !!]
     procedure :: tolerancesSet => tolerancesSetGeneric
  end type integratorMulti

  ! Generic one-dimensional multi-integrand integrator.
  type, abstract, extends(integratorMulti) :: integratorMulti1D
     !!{RST
     Generic one-dimensional multi-integrand numerical integrator class.
     !!}
     private
     integer                                      :: integrandCount
     procedure(integrandMulti1D), pointer, nopass :: integrand
   contains
     !![
     <methods docformat="rst">
       <method description="Set the integrand function to be integrated." method="integrandSet" />
       <method description="Evaluate the integral." method="evaluate" />
     </methods>
     !!]
     procedure                            :: integrandSet => integrandMulti1DSet
     procedure(evaluateMulti1D), deferred :: evaluate
  end type integratorMulti1D
  abstract interface
     function evaluateMulti1D(self,a,b,status) result(integral)
       import integratorMulti1D
       class           (integratorMulti1D), intent(inout)                            :: self
       double precision                   , intent(in   )                            :: a         , b
       integer                            , intent(  out)                 , optional :: status
       double precision                   , dimension(self%integrandCount)           :: integral
     end function evaluateMulti1D
  end interface
  abstract interface
     function integrandMulti1D(n,x,e)
       integer         , intent(in   ) :: n
       double precision, intent(in   ) :: x
       logical         , intent(inout) :: e
       double precision, dimension(n)  :: integrandMulti1D
     end function integrandMulti1D
  end interface

  ! Generic one-dimensional multi-integrand vectorized integrator.
  type, abstract, extends(integratorMulti) :: integratorMultiVectorized1D
     !!{RST
     Generic one-dimensional multi-integrand vectorized numerical integrator class.
     !!}
     private
     integer                                                :: integrandCount
     procedure(integrandMultiVectorized1D), pointer, nopass :: integrand
   contains
     !![
     <methods docformat="rst">
       <method description="Set the integrand function to be integrated." method="integrandSet" />
     </methods>
     !!]
     procedure                                      :: integrandSet => integrandMultiVectorizedSet1D
     procedure(evaluateMultiVectorized1D), deferred :: evaluate
  end type integratorMultiVectorized1D
  abstract interface
     function evaluateMultiVectorized1D(self,a,b,status) result(integral)
       import integratorMultiVectorized1D
       class           (integratorMultiVectorized1D), intent(inout)                            :: self
       double precision                             , intent(in   )                            :: a       , b
       integer                                      , intent(  out)                 , optional :: status
       double precision                             , dimension(self%integrandCount)           :: integral
     end function evaluateMultiVectorized1D
  end interface
  abstract interface
     subroutine integrandMultiVectorized1D(n,x,e,integrand)
       integer         , intent(in   )                       :: n
       double precision, intent(in   ), dimension(       : ) :: x
       logical         , intent(inout), dimension(       : ) :: e
       double precision, intent(  out), dimension(n,size(x)) :: integrand
     end subroutine integrandMultiVectorized1D
  end interface

  ! Vectorized composite multi-integrand Gauss-Kronrod 1D integrator.
  type, extends(integratorMultiVectorized1D) :: integratorMultiVectorizedCompositeGaussKronrod1D
     !!{RST
     One-dimensional multi-integrand numerical integrator class using a vectorized composite Gauss-Kronrod rule.
     !!}
     private
     integer                                     :: intervalsMaximum
     double precision, allocatable, dimension(:) :: xKronrod         , wGauss, wKronrod
   contains
     !![
     <methods docformat="rst">
       <method description="Evaluate the integrals." method="evaluate" />
       <method description="Evaluate the integrals over an interval and also return errors on the integrals." method="evaluateInterval" />
     </methods>
     !!]
     procedure :: evaluate         => multiVectorizedCompositeGaussKronrod1DEvaluate
     procedure :: evaluateInterval => multiVectorizedCompositeGaussKronrod1DEvaluateInterval
  end type integratorMultiVectorizedCompositeGaussKronrod1D

  interface integratorMultiVectorizedCompositeGaussKronrod1D
     !!{RST
     Constructors for the :galacticus-class:`integratorMultiVectorizedCompositeGaussKronrod1D` numerical integrator class.
     !!}
     module procedure multiVectorizedCompositeGaussKronrod1DConstructor
  end interface integratorMultiVectorizedCompositeGaussKronrod1D

  ! Vectorized composite multi-integrand trapezoidal 1D integrator.
  type, extends(integratorMultiVectorized1D) :: integratorMultiVectorizedCompositeTrapezoidal1D
     !!{RST
     One-dimensional multi-integrand numerical integrator class using a vectorized composite trapezoidal rule.
     !!}
     private
     integer :: intervalsMaximum
   contains
     !![
     <methods docformat="rst">
       <method description="Evaluate the integrals." method="evaluate" />
     </methods>
     !!]
     procedure :: evaluate   => multiVectorizedCompositeTrapezoidal1DEvaluate
  end type integratorMultiVectorizedCompositeTrapezoidal1D

  interface integratorMultiVectorizedCompositeTrapezoidal1D
     !!{RST
     Constructors for the :galacticus-class:`integratorMultiVectorizedCompositeTrapezoidal1D` numerical integrator class.
     !!}
     module procedure multiVectorizedCompositeTrapezoidal1DConstructor
  end interface integratorMultiVectorizedCompositeTrapezoidal1D

  ! Implementations are in the submodules in this directory.
  interface
     module subroutine intervalInsert(head, newFirst, newLast)
       type(interval), intent(inout), pointer :: head
       type(interval), intent(in   ), pointer :: newFirst, newLast
     end subroutine intervalInsert

     module subroutine intervalListDestroy(head)
       type(interval), intent(inout), pointer :: head
     end subroutine intervalListDestroy

     module subroutine intervalMultiSearch(current, previous, converged, errorsMaximum)
       type            (intervalMulti), intent(inout), pointer :: current, previous
       logical                        , intent(in   ), dimension(:) :: converged
       double precision               , intent(in   ), dimension(:) :: errorsMaximum
     end subroutine intervalMultiSearch

     module subroutine intervalMultiLink(head, current, previous, newFirst, newLast, converged, errorsMaximum)
       type            (intervalMulti), intent(inout), pointer :: head
       type            (intervalMulti), intent(in   ), pointer :: current, previous, newFirst, newLast
       logical                        , intent(in   ), dimension(:) :: converged
       double precision               , intent(in   ), dimension(:) :: errorsMaximum
     end subroutine intervalMultiLink

     module function intervalMultiNew(a, b, integrandCount) result(interval_)
       type            (intervalMulti), pointer :: interval_
       double precision               , intent(in   ) :: a, b
       integer                        , intent(in   ) :: integrandCount
     end function intervalMultiNew

     module subroutine intervalMultiWorkspaceEnsure(list, listValue, listRank, intervalCount)
       type            (intervalMultiList), intent(inout), dimension(:), allocatable :: list
       double precision                   , intent(inout), dimension(:), allocatable :: listValue
       integer         (c_size_t         ), intent(inout), dimension(:), allocatable :: listRank
       integer         (c_size_t         ), intent(in   ) :: intervalCount
     end subroutine intervalMultiWorkspaceEnsure

     module subroutine intervalMultiResort(head, intervalCount, converged, errorScale, list, listValue, listRank)
       type            (intervalMulti    ), intent(inout), pointer :: head
       integer         (c_size_t         ), intent(in   ) :: intervalCount
       logical                            , intent(in   ), dimension(:) :: converged
       double precision                   , intent(in   ), dimension(:) :: errorScale
       type            (intervalMultiList), intent(inout), dimension(:), allocatable :: list
       double precision                   , intent(inout), dimension(:), allocatable :: listValue
       integer         (c_size_t         ), intent(inout), dimension(:), allocatable :: listRank
     end subroutine intervalMultiResort

     module subroutine intervalMultiListDestroy(head)
       type(intervalMulti), intent(inout), pointer :: head
     end subroutine intervalMultiListDestroy

     module subroutine intervalMultiFailureReport(a, b, current, head, integral, error, converged, statusPresent)
       double precision                , intent(in   ) :: a, b
       type            (intervalMulti ), intent(in   ), pointer :: current, head
       double precision                , intent(in   ), dimension(:) :: integral, error
       logical                         , intent(in   ), dimension(:) :: converged
       logical                         , intent(in   ) :: statusPresent
     end subroutine intervalMultiFailureReport

     module subroutine gaussKronrodRuleSet(order, xKronrod, wGauss, wKronrod)
       integer                                     , intent(in   ) :: order
       double precision, allocatable, dimension(:) , intent(inout) :: xKronrod, wGauss, wKronrod
     end subroutine gaussKronrodRuleSet

     elemental module function gaussKronrodError(error, integralAbsolute, integralAsc)
       double precision :: gaussKronrodError
       double precision, intent(in   ) :: error, integralAbsolute, integralAsc
     end function gaussKronrodError

     module function gaussKronrodEvaluate(self, a, b, iterationsMaximum) result(integral)
       double precision :: integral
       class           (integrator2                      ), intent(inout) :: self
       integer                                            , intent(in   ) :: iterationsMaximum
       double precision                                   , intent(in   ) :: a, b
     end function gaussKronrodEvaluate

     module subroutine gaussKronrodEvaluateInterval(self, a, b, integralKronrod, error)
       class           (integrator2), intent(inout) :: self
       double precision             , intent(in   ) :: a, b
       double precision             , intent(  out) :: integralKronrod, error
     end subroutine gaussKronrodEvaluateInterval

     module function compositeGaussKronrod1DConstructor(iterationsMaximum, order) result(self)
       type  (integratorCompositeGaussKronrod1D) :: self
       integer                                   , intent(in   ) :: iterationsMaximum, order
     end function compositeGaussKronrod1DConstructor

     module function compositeGaussKronrod1DEvaluate(self, a, b)
       double precision :: compositeGaussKronrod1DEvaluate
       class           (integratorCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                   , intent(in   ) :: a, b
     end function compositeGaussKronrod1DEvaluate

     module subroutine compositeGaussKronrod1DEvaluateInterval(self, a, b, integralKronrod, error)
       class           (integratorCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                   , intent(in   ) :: a, b
       double precision                                   , intent(  out) :: integralKronrod, error
     end subroutine compositeGaussKronrod1DEvaluateInterval

     module function vectorizedCompositeGaussKronrod1DConstructor(iterationsMaximum, order) result(self)
       type  (integratorVectorizedCompositeGaussKronrod1D) :: self
       integer                                             , intent(in   ) :: iterationsMaximum, order
     end function vectorizedCompositeGaussKronrod1DConstructor

     module function vectorizedCompositeGaussKronrod1DEvaluate(self, a, b)
       double precision :: vectorizedCompositeGaussKronrod1DEvaluate
       class           (integratorVectorizedCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                             , intent(in   ) :: a, b
     end function vectorizedCompositeGaussKronrod1DEvaluate

     module subroutine vectorizedCompositeGaussKronrod1DEvaluateInterval(self, a, b, integralKronrod, error)
       class           (integratorVectorizedCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                             , intent(in   ) :: a, b
       double precision                                             , intent(  out) :: integralKronrod, error
     end subroutine vectorizedCompositeGaussKronrod1DEvaluateInterval

     module function multiVectorizedCompositeGaussKronrod1DConstructor(intervalsMaximum, order) result(self)
       type  (integratorMultiVectorizedCompositeGaussKronrod1D) :: self
       integer                                                  , intent(in   ) :: intervalsMaximum, order
     end function multiVectorizedCompositeGaussKronrod1DConstructor

     module function multiVectorizedCompositeGaussKronrod1DEvaluate(self, a, b, status) result(integral)
       class           (integratorMultiVectorizedCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                                  , intent(in   ) :: a, b
       integer                                                           , intent(  out)                 , optional :: status
       double precision                                                  , dimension(self%integrandCount) :: integral
     end function multiVectorizedCompositeGaussKronrod1DEvaluate

     module subroutine multiVectorizedCompositeGaussKronrod1DEvaluateInterval(self, a, b, integralKronrod, error, mustEvaluate)
       class           (integratorMultiVectorizedCompositeGaussKronrod1D), intent(inout) :: self
       double precision                                                  , intent(in   ) :: a, b
       double precision                                                  , intent(  out), dimension(self%integrandCount                        ) :: integralKronrod, error
       logical                                                           , intent(inout), dimension(self%integrandCount                        ) :: mustEvaluate
     end subroutine multiVectorizedCompositeGaussKronrod1DEvaluateInterval

     module function compositeTrapezoidalConstructor1D(iterationsMaximum) result(self)
       type  (integratorCompositeTrapezoidal1D) :: self
       integer                                  , intent(in   ) :: iterationsMaximum
     end function compositeTrapezoidalConstructor1D

     module function compositeTrapezoidalEvaluate1D(self, a, b)
       double precision :: compositeTrapezoidalEvaluate1D
       class           (integratorCompositeTrapezoidal1D), intent(inout) :: self
       double precision                                  , intent(in   ) :: a, b
     end function compositeTrapezoidalEvaluate1D

     module function adaptiveCompositeTrapezoidalConstructor1D(iterationsMaximum) result(self)
       type   (integratorAdaptiveCompositeTrapezoidal1D) :: self
       integer                                          , intent(in   ) :: iterationsMaximum
     end function adaptiveCompositeTrapezoidalConstructor1D

     module function adaptiveCompositeTrapezoidalEvaluate1D(self, a, b)
       double precision :: adaptiveCompositeTrapezoidalEvaluate1D
       class           (integratorAdaptiveCompositeTrapezoidal1D), intent(inout) :: self
       double precision                                          , intent(in   ) :: a, b
     end function adaptiveCompositeTrapezoidalEvaluate1D

     module function vectorizedCompositeTrapezoidalConstructor1D(iterationsMaximum) result(self)
       type           (integratorVectorizedCompositeTrapezoidal1D) :: self
       integer                                                     , intent(in   ) :: iterationsMaximum
     end function vectorizedCompositeTrapezoidalConstructor1D

     module function vectorizedCompositeTrapezoidalEvaluate1D(self, a, b)
       double precision :: vectorizedCompositeTrapezoidalEvaluate1D
       class           (integratorVectorizedCompositeTrapezoidal1D), intent(inout) :: self
       double precision                                            , intent(in   ) :: a, b
     end function vectorizedCompositeTrapezoidalEvaluate1D

     module function multiVectorizedCompositeTrapezoidal1DConstructor(intervalsMaximum) result(self)
       type  (integratorMultiVectorizedCompositeTrapezoidal1D) :: self
       integer                                                 , intent(in   ) :: intervalsMaximum
     end function multiVectorizedCompositeTrapezoidal1DConstructor

     module function multiVectorizedCompositeTrapezoidal1DEvaluate(self, a, b, status) result(integral)
       class           (integratorMultiVectorizedCompositeTrapezoidal1D), intent(inout) :: self
       double precision                                                 , intent(in   ) :: a, b
       integer                                                          , intent(  out)                   , optional :: status
       double precision                                                 , dimension(self%integrandCount  ) :: integral
     end function multiVectorizedCompositeTrapezoidal1DEvaluate
  end interface

contains

  subroutine iterationsMaximumValidate(iterationsMaximum)
    !!{RST
    Validate the maximum number of iterations given to an integrator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer, intent(in   ) :: iterationsMaximum

    if (iterationsMaximum < 3)                                      &
         & call Error_Report(                                       &
         &                   'at least 3 iterations are required'// &
         &                    {introspection:location}              &
         &                  )
    return
  end subroutine iterationsMaximumValidate

  subroutine toleranceSetGeneric(self,toleranceAbsolute,toleranceRelative)
    !!{RST
    Initialize the tolerances for numerical integrators.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (integrator2), intent(inout)           :: self
    double precision             , intent(in   ), optional :: toleranceAbsolute,toleranceRelative

    if (.not.(present(toleranceAbsolute).or.present(toleranceRelative)))                                 &
         &  call Error_Report(                                                                &
         &                               'at least one of absolute and relative tolerance must be set'// &
         &                               {introspection:location}                                        &
         &                              )
    self%toleranceAbsolute=0.0d0
    self%toleranceRelative=0.0d0
    if (present(toleranceAbsolute)) then
       if (toleranceAbsolute < 0.0d0) call Error_Report('absolute tolerance must be non-negative'//{introspection:location})
       self%toleranceAbsolute=toleranceAbsolute
    end if
    if (present(toleranceRelative)) then
       if (toleranceRelative < 0.0d0) call Error_Report('relative tolerance must be non-negative'//{introspection:location})
       self%toleranceRelative=toleranceRelative
    end if
    return
  end subroutine toleranceSetGeneric

  subroutine integrandSet1D(self,integrand)
    !!{RST
    Initialize the integrand for one-dimensional numerical integrators.
    !!}
    implicit none
    class    (integrator1D), intent(inout) :: self
    procedure(integrand1D )                :: integrand

    self%integrand => integrand
    return
  end subroutine integrandSet1D

  subroutine integrandVectorizedSet1D(self,integrand)
    !!{RST
    Initialize the integrand for one-dimensional numerical integrators.
    !!}
    implicit none
    class    (integratorVectorized1D), intent(inout) :: self
    procedure(integrandVectorized1D )                :: integrand

    self%integrand => integrand
    return
  end subroutine integrandVectorizedSet1D

  subroutine tolerancesSetGeneric(self,toleranceAbsolute,toleranceRelative)
    !!{RST
    Initialize the tolerances for multi-integrand numerical integrators.
    !!}
    use :: Error , only : Error_Report
    implicit none
    class           (integratorMulti), intent(inout)                         :: self
    double precision                 , intent(in   ), dimension(:), optional :: toleranceAbsolute, toleranceRelative
    integer                                                                  :: toleranceCount

    if (.not.(present(toleranceAbsolute).or.present(toleranceRelative)))                                 &
         &  call Error_Report(                                                                &
         &                               'at least one of absolute and relative tolerance must be set'// &
         &                               {introspection:location}                                        &
         &                              )
    toleranceCount=0
    if (present(toleranceAbsolute)) then
       toleranceCount=size(toleranceAbsolute)
       if (present(toleranceRelative).and.size(toleranceRelative) /= size(toleranceAbsolute))                   &
         &  call Error_Report(                                                                       &
         &                               'absolute and relative tolerances must have same number of elements'// &
         &                               {introspection:location}                                               &
         &                              )
    else
       toleranceCount=size(toleranceRelative)
    end if
    if (allocated(self%toleranceAbsolute)) deallocate(self%toleranceAbsolute)
    if (allocated(self%toleranceRelative)) deallocate(self%toleranceRelative)
    allocate(self%toleranceAbsolute(toleranceCount))
    allocate(self%toleranceRelative(toleranceCount))
    self%toleranceAbsolute=0.0d0
    self%toleranceRelative=0.0d0
    if (present(toleranceAbsolute)) then
       if (any(toleranceAbsolute < 0.0d0)) call Error_Report('absolute tolerances must be non-negative'//{introspection:location})
       self%toleranceAbsolute=toleranceAbsolute
    end if
    if (present(toleranceRelative)) then
       if (any(toleranceRelative < 0.0d0)) call Error_Report('relative tolerances must be non-negative'//{introspection:location})
       self%toleranceRelative=toleranceRelative
    end if
    return
  end subroutine tolerancesSetGeneric

  subroutine integrandMulti1DSet(self,integrandCount,integrand)
    !!{RST
    Initialize the integrand for one-dimensional numerical integrators.
    !!}
    implicit none
    class    (integratorMulti1D), intent(inout) :: self
    integer                     , intent(in   ) :: integrandCount
    procedure(integrandMulti1D )                :: integrand

    self%integrandCount =  integrandCount
    self%integrand      => integrand
    return
  end subroutine integrandMulti1DSet

  subroutine integrandMultiVectorizedSet1D(self,integrandCount,integrand)
    !!{RST
    Initialize the integrand for one-dimensional numerical integrators.
    !!}
    implicit none
    class    (integratorMultiVectorized1D), intent(inout) :: self
    integer                               , intent(in   ) :: integrandCount
    procedure(integrandMultiVectorized1D )                :: integrand

    self%integrandCount =  integrandCount
    self%integrand      => integrand
    return
  end subroutine integrandMultiVectorizedSet1D

end module Numerical_Integration2
