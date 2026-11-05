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
Contains a submodule which provides composite Gauss-Kronrod numerical integrators.
!!}

submodule (Numerical_Integration2) Numerical_Integration2_Gauss_Kronrod
  !!{RST
  Provides composite Gauss-Kronrod numerical integrators.
  !!}
  implicit none

contains

  module procedure gaussKronrodRuleSet
    !!{RST
    Set the evaluation points and weights for a Gauss-Kronrod rule of the given ``order``. Evaluation points and weights are taken from those used in the :term:`GSL`.
    !!}
    use :: Error, only : Error_Report
    implicit none

    select case (order)
    case (15) ! 15-point Kronrod rule.
       allocate(xKronrod(8))
       allocate(wKronrod(8))
       allocate(wGauss  (4))
       xKronrod =[                                                  0.991455371120812639206854697526329d0,            0.949107912342758524526189684047851d0,            0.864864423359769072789712788640926d0,            0.741531185599394439863864773280788d0,            0.586087235467691130294144838258730d0,            0.405845151377397166906606412076961d0,            0.207784955007898467600689403773245d0,            0.000000000000000000000000000000000d0           ]
       wGauss  =[                                                   0.129484966168869693270611432679082d0,            0.279705391489276667901467771423780d0,            0.381830050505118944950369775488975d0,            0.417959183673469387755102040816327d0           ]
       wKronrod=[                                                   0.022935322010529224963732008058970d0,            0.063092092629978553290700663189204d0,            0.104790010322250183839876322541518d0,            0.140653259715525918745189590510238d0,            0.169004726639267902826583426598550d0,            0.190350578064785409913256402421014d0,            0.204432940075298892414161999234649d0,            0.209482141084727828012999174891714d0           ]
    case (61) ! 61-point Kronrod rule.
       allocate(xKronrod(31))
       allocate(wKronrod(31))
       allocate(wGauss  (15))
       xKronrod =[                                                  0.999484410050490637571325895705811d0,            0.996893484074649540271630050918695d0,            0.991630996870404594858628366109486d0,            0.983668123279747209970032581605663d0,            0.973116322501126268374693868423707d0,            0.960021864968307512216871025581798d0,            0.944374444748559979415831324037439d0,            0.926200047429274325879324277080474d0,            0.905573307699907798546522558925958d0,            0.882560535792052681543116462530226d0,            0.857205233546061098958658510658944d0,            0.829565762382768397442898119732502d0,            0.799727835821839083013668942322683d0,            0.767777432104826194917977340974503d0,            0.733790062453226804726171131369528d0,            0.697850494793315796932292388026640d0,            0.660061064126626961370053668149271d0,            0.620526182989242861140477556431189d0,            0.579345235826361691756024932172540d0,            0.536624148142019899264169793311073d0,            0.492480467861778574993693061207709d0,            0.447033769538089176780609900322854d0,            0.400401254830394392535476211542661d0,            0.352704725530878113471037207089374d0,            0.304073202273625077372677107199257d0,            0.254636926167889846439805129817805d0,            0.204525116682309891438957671002025d0,            0.153869913608583546963794672743256d0,            0.102806937966737030147096751318001d0,            0.051471842555317695833025213166723d0,            0.000000000000000000000000000000000d0            ]
       wGauss   =[                                                  0.007968192496166605615465883474674d0,            0.018466468311090959142302131912047d0,            0.028784707883323369349719179611292d0,            0.038799192569627049596801936446348d0,            0.048402672830594052902938140422808d0,            0.057493156217619066481721689402056d0,            0.065974229882180495128128515115962d0,            0.073755974737705206268243850022191d0,            0.080755895229420215354694938460530d0,            0.086899787201082979802387530715126d0,            0.092122522237786128717632707087619d0,            0.096368737174644259639468626351810d0,            0.099593420586795267062780282103569d0,            0.101762389748405504596428952168554d0,            0.102852652893558840341285636705415d0            ]
       wKronrod =[                                                  0.001389013698677007624551591226760d0,            0.003890461127099884051267201844516d0,            0.006630703915931292173319826369750d0,            0.009273279659517763428441146892024d0,            0.011823015253496341742232898853251d0,            0.014369729507045804812451432443580d0,            0.016920889189053272627572289420322d0,            0.019414141193942381173408951050128d0,            0.021828035821609192297167485738339d0,            0.024191162078080601365686370725232d0,            0.026509954882333101610601709335075d0,            0.028754048765041292843978785354334d0,            0.030907257562387762472884252943092d0,            0.032981447057483726031814191016854d0,            0.034979338028060024137499670731468d0,            0.036882364651821229223911065617136d0,            0.038678945624727592950348651532281d0,            0.040374538951535959111995279752468d0,            0.041969810215164246147147541285970d0,            0.043452539701356069316831728117073d0,            0.044814800133162663192355551616723d0,            0.046059238271006988116271735559374d0,            0.047185546569299153945261478181099d0,            0.048185861757087129140779492298305d0,            0.049055434555029778887528165367238d0,            0.049795683427074206357811569379942d0,            0.050405921402782346840893085653585d0,            0.050881795898749606492297473049805d0,            0.051221547849258772170656282604944d0,            0.051426128537459025933862879215781d0,            0.051494729429451567558340433647099d0            ]
    case default
       call Error_Report('unknown order'//{introspection:location})
    end select
    return
  end procedure gaussKronrodRuleSet

  module procedure gaussKronrodError
    !!{RST
    Compute the error estimate on an integral over an interval from the raw difference between the Kronrod and Gauss rules. This is based on the :term:`GSL` implementation. The logic is modified so that no power/sqrt is performed until after the evaluation of the conditional (no need to compute an expensive function if we don't really need to). Additionally, the expression for the error in the case where the conditional evaluates true is modified so that it involves a sqrt() instead of a ()**1.5 as this is better optimized.
    !!}
    implicit none
    double precision, parameter :: errorScaleFactor   =200.0d0
    double precision, parameter :: errorScaleFactorPow=errorScaleFactor**1.5d0

    gaussKronrodError=error
    if (integralAsc /= 0.0d0 .and. error /= 0.0d0) then
       if (errorScaleFactor*error < integralAsc) then
          gaussKronrodError=errorScaleFactorPow*error*sqrt(error/integralAsc)
       else
          gaussKronrodError=integralAsc
       end if
    end if
    if (integralAbsolute > tiny(1.0d0)/(50.0d0*epsilon(1.0d0))) gaussKronrodError=max(gaussKronrodError,50.0d0*epsilon(1.0d0)*integralAbsolute)
    return
  end procedure gaussKronrodError

  ! This illustrates how we might vectorize this. On each loop, you want to split the stack in two such that the sum of the errors below the split is less than the required error on the total integral. This means that we know that we need to bisect *every* interval above the split since in the best-case scenario the errors after bisecting these go to zero and we have then just sneaked in under the required error budget. Then, build an array of all midpoints in the intervals to be bisected - do a vector eval of the function at those points - then create the new intervals and populate them with the required function values. Could then look at optimizing the creation of the new intervals:
  !==> suppose we're making N new intervals
  !==> allocate a pointer array of intervals of size N instead of doing these one at a time
  !==> efficient means to sort them? use quickSort for example instead of insertion?

  module procedure gaussKronrodEvaluate
    !!{RST
    Evaluate a one-dimension integral by adaptive bisection using a numerical composite Gauss-Kronrod rule. The bisection logic is common to the scalar and vectorized integrators, which differ only in how they evaluate the rule over an interval.
    !!}
    use :: Error, only : Error_Report
    implicit none
    logical :: converged
    double precision :: error, midpoint
    type            (interval                         ), pointer :: head, newInterval1, newInterval2, current

    ! Create our first estimate of the integral, using a single interval.
    allocate(head)
    head%a        =  a
    head%b        =  b
    head%depth    =  1
    head%next     => null()
    call gaussKronrodEvaluateInterval(self,head%a,head%b,head%integral,head%error)
    ! Initialize current integral and error estimates, and flag that we're not converged.
    integral=head%integral
    error                          =head%error
    converged                      =                                                    abs(error) <= self%toleranceAbsolute                                       .or.                                                                         abs(error) <= self%toleranceRelative*abs(integral)
    ! Iterate until convergence is reached.
    do while (.not.converged)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       if (current%depth > iterationsMaximum)                     call Error_Report(                                                    'maximum iterations exceeded'//                      {introspection:location}                          )
       midpoint =  (current%b+current%a)/2.0d0
       ! Create two new subintervals.
       allocate(newInterval1)
       allocate(newInterval2)
       newInterval1%next  => null()
       newInterval2%next  => null()
       newInterval1%depth =  current%depth+1
       newInterval2%depth =  current%depth+1
       newInterval1%a    =  current%a
       newInterval1%b    =  midpoint
       newInterval2%a    =  midpoint
       newInterval2%b    =  current%b
       ! Compute the integral and error estimate in each new subinterval.
       call gaussKronrodEvaluateInterval(self,newInterval1%a,newInterval1%b,newInterval1%integral,newInterval1%error)
       call gaussKronrodEvaluateInterval(self,newInterval2%a,newInterval2%b,newInterval2%integral,newInterval2%error)
       ! Update the total integral and error estimate for the new subintervals.
       integral           =integral   -current     %integral             +newInterval1%integral             +newInterval2%integral
       error                                     =error                             -current     %error                +newInterval1%error                +newInterval2%error
       ! Test for convergence.
       converged=                                                                          abs(error) <= self%toleranceAbsolute                                       .or.                                                                         abs(error) <= self%toleranceRelative*abs(integral)
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       call intervalInsert(head,newInterval1,newInterval1)
       call intervalInsert(head,newInterval2,newInterval2)
    end do
    ! Destroy the stack.
    call intervalListDestroy(head)
    return
  end procedure gaussKronrodEvaluate

  module procedure gaussKronrodEvaluateInterval
    !!{RST
    Evaluate a Gauss-Kronrod rule over an interval using the method of the given integrator.
    !!}
    use :: Error, only : Error_Report
    implicit none

    select type (self)
    class is (integratorCompositeGaussKronrod1D          )
       call self%evaluateInterval(a,b,integralKronrod,error)
    class is (integratorVectorizedCompositeGaussKronrod1D)
       call self%evaluateInterval(a,b,integralKronrod,error)
    class default
       call Error_Report('unknown integrator class'//{introspection:location})
    end select
    return
  end procedure gaussKronrodEvaluateInterval

  module procedure compositeGaussKronrod1DConstructor
    !!{RST
    Constructor for a one-dimensional, composite Gauss-Kronrod numerical integrator. Evaluation points and weights are taken from those used in the :term:`GSL`.
    !!}
    implicit none

    call iterationsMaximumValidate(iterationsMaximum)
    self%iterationsMaximum=iterationsMaximum
    call gaussKronrodRuleSet(order,self%xKronrod,self%wGauss,self%wKronrod)
    return
  end procedure compositeGaussKronrod1DConstructor

  module procedure compositeGaussKronrod1DEvaluate
    !!{RST
    Evaluate a one-dimension integral using a numerical composite Gauss-Kronrod rule.
    !!}
    implicit none

    compositeGaussKronrod1DEvaluate=gaussKronrodEvaluate(self,a,b,self%iterationsMaximum)
    return
  end procedure compositeGaussKronrod1DEvaluate

  module procedure compositeGaussKronrod1DEvaluateInterval
    !!{RST
    Evaluate the integral over an interval using the Gauss-Kronrod method (also estimates the error). Specific implementation is based on that in the :term:`GSL`.
    !!}
    implicit none
    double precision                                   , dimension(size(self%xKronrod)-1) :: fValue1, fValue2
    double precision :: integralGauss, fCenter, xCenter, halfLength, fSum, x, halfLengthAbsolute, integralAbsolute, mean, integralAsc
    integer :: pointCountKronrod, pointCountGauss, jtw, jtwm1, j

    pointCountKronrod =size(self%wKronrod)
    pointCountGauss   =size(self%wGauss  )
    xCenter           =0.5d0*(b+a)
    halfLength        =0.5d0*(b-a)
    halfLengthAbsolute=abs(halfLength)
    fCenter           =self%integrand(xCenter)
    integralKronrod   =fCenter*self%wKronrod(pointCountKronrod)
    integralGauss     =0.0d0
    integralAbsolute  =abs(integralKronrod)
    if (mod(pointCountKronrod,2) == 0) integralGauss=fCenter*self%wGauss(pointCountKronrod/2)
    do j=0,(pointCountKronrod-1)/2-1
       jtw             =j*2+1
       x               =halfLength*self%xKronrod(jtw+1)
       fValue1(jtw+1)  =self%integrand(xCenter-x)
       fValue2(jtw+1)  =self%integrand(xCenter+x)
       fSum            =fValue1(jtw+1)+fValue2(jtw+1)
       integralGauss   =integralGauss   +self%wGauss  (j  +1)*fSum
       integralKronrod =integralKronrod +self%wKronrod(jtw+1)*fSum
       integralAbsolute=integralAbsolute+self%wKronrod(jtw+1)*(abs(fValue1(jtw+1))+abs(fValue2(jtw+1)))
    end do
    do j=0,pointCountKronrod/2-1
       jtwm1           =j*2
       x               =halfLength*self%xKronrod(jtwm1+1)
       fValue1(jtwm1+1)=self%integrand(xCenter-x)
       fValue2(jtwm1+1)=self%integrand(xCenter+x)
       integralKronrod =integralKronrod +self%wKronrod(jtwm1+1)*(    fValue1(jtwm1+1) +    fValue2(jtwm1+1))
       integralAbsolute=integralAbsolute+self%wKronrod(jtwm1+1)*(abs(fValue1(jtwm1+1))+abs(fValue2(jtwm1+1)))
    end do
    mean        =integralKronrod*0.5d0
    integralAsc=self%wKronrod(pointCountKronrod)*abs(fCenter-mean)
    do j=0,pointCountKronrod-2
       integralAsc=integralAsc+self%wKronrod(j+1)*(abs(fValue1(j+1)-mean)+abs(fValue2(j+1)-mean))
    end do
    error           =abs((integralKronrod-integralGauss)*halfLength)
    integralKronrod =integralKronrod *halfLength
    integralAbsolute=integralAbsolute*halfLengthAbsolute
    integralAsc     =integralAsc     *halfLengthAbsolute
    error           =gaussKronrodError(error,integralAbsolute,integralAsc)
    return
  end procedure compositeGaussKronrod1DEvaluateInterval

  module procedure vectorizedCompositeGaussKronrod1DConstructor
    !!{RST
    Constructor for a one-dimensional, vectorized composite Gauss-Kronrod numerical integrator. Evaluation points and weights are taken from those used in the :term:`GSL`.
    !!}
    implicit none

    call iterationsMaximumValidate(iterationsMaximum)
    self%iterationsMaximum=iterationsMaximum
    call gaussKronrodRuleSet(order,self%xKronrod,self%wGauss,self%wKronrod)
    return
  end procedure vectorizedCompositeGaussKronrod1DConstructor

  module procedure vectorizedCompositeGaussKronrod1DEvaluate
    !!{RST
    Evaluate a one-dimension integral using a numerical composite Gauss-Kronrod rule.
    !!}
    implicit none

    vectorizedCompositeGaussKronrod1DEvaluate=gaussKronrodEvaluate(self,a,b,self%iterationsMaximum)
    return
  end procedure vectorizedCompositeGaussKronrod1DEvaluate

  module procedure vectorizedCompositeGaussKronrod1DEvaluateInterval
    !!{RST
    Evaluate the integral over an interval using the Gauss-Kronrod method (also estimates the error). Specific implementation is based on that in the :term:`GSL`.
    !!}
    implicit none
    double precision                                             , pointer      , dimension(:                      ) :: fValue1, fValue2
    double precision                                             , target       , dimension(2*size(self%xKronrod)-1) :: fUnion
    double precision                                                            , dimension(2*size(self%xKronrod)-1) :: xUnion
    double precision :: integralGauss, xCenter, halfLength, halfLengthAbsolute, integralAbsolute, mean, integralAsc
    integer :: pointCountKronrod, pointCountGauss

    ! Establish point counts and interval.
    pointCountKronrod =size(self%wKronrod)
    pointCountGauss   =size(self%wGauss  )
    if (mod(pointCountKronrod,2) == 0) pointCountGauss=pointCountGauss-1
    xCenter           =0.5d0*(b+a)
    halfLength        =0.5d0*(b-a)
    halfLengthAbsolute=abs(halfLength)
    ! Evaluate function at all required points.
    xUnion(1                                            )=  xCenter
    xUnion(2                    :2+  pointCountKronrod-2)=  xCenter-self%xKronrod(1:pointCountKronrod-1)*halfLength
    xUnion(2+pointCountKronrod-1:2+2*pointCountKronrod-3)=  xCenter+self%xKronrod(1:pointCountKronrod-1)*halfLength
    fUnion                                               =  self%integrand(xUnion)
    fValue1                                              => fUnion(2                    :2+  pointCountKronrod-2)
    fValue2                                              => fUnion(2+pointCountKronrod-1:2+2*pointCountKronrod-3)
    ! Evaluate integrals.
    integralGauss=sum(                                                 +self%wGauss(1:  pointCountGauss  )              *(                                                 +fValue1  (2:2*pointCountGauss:2)                +fValue2  (2:2*pointCountGauss:2)               )                                              )
    if (mod(pointCountKronrod,2) == 0) integralGauss=integralGauss+fUnion(1)*self%wGauss(pointCountKronrod/2)
    integralKronrod=+    self%wKronrod  (  pointCountKronrod  )                  *       fUnion (1                    )             +sum(                                                        +self%wKronrod (1:pointCountKronrod-1)                  *(                                                        +     fValue1(1:pointCountKronrod-1)                    +     fValue2(1:pointCountKronrod-1)                   )                                                     )
    integralAbsolute=+self%wKronrod     (  pointCountKronrod  )                   *  abs(fUnion (1                    ))             +sum(                                                        +self%wKronrod(1:pointCountKronrod-1)                   *(                                                        +abs(fValue1(1:pointCountKronrod-1))                    +abs(fValue2(1:pointCountKronrod-1))                   )                                                     )
    mean       =integralKronrod*0.5d0
    integralAsc=+sum(                                                  +self%wKronrod(1:pointCountKronrod-1)             *(                                                  +abs(fValue1   -mean)                             +abs(fValue2   -mean)                            )                                               )                                             +    self%wKronrod(  pointCountKronrod  )         *       abs(fUnion (1)-mean)
    ! Evaluate error.
    error           =abs((integralKronrod-integralGauss)*halfLength)
    integralKronrod =integralKronrod *halfLength
    integralAbsolute=integralAbsolute*halfLengthAbsolute
    integralAsc     =integralAsc     *halfLengthAbsolute
    error           =gaussKronrodError(error,integralAbsolute,integralAsc)
    return
  end procedure vectorizedCompositeGaussKronrod1DEvaluateInterval

  module procedure multiVectorizedCompositeGaussKronrod1DConstructor
    !!{RST
    Constructor for a one-dimensional, multi-integrand, vectorized composite Gauss-Kronrod numerical integrator. Evaluation points and weights are taken from those used in the :term:`GSL`.
    !!}
    use :: Error , only : Error_Report
    implicit none

    if (intervalsMaximum < 1)                                                 call Error_Report(                                                                   'at least 1 interval is required'//                                 {introspection:location}                                         )
    self%intervalsMaximum=intervalsMaximum
    call gaussKronrodRuleSet(order,self%xKronrod,self%wGauss,self%wKronrod)
    return
  end procedure multiVectorizedCompositeGaussKronrod1DConstructor

  module procedure multiVectorizedCompositeGaussKronrod1DEvaluate
    !!{RST
    Evaluate a one-dimension integral using a numerical composite Gauss-Kronrod rule.
    !!}
    use            :: Error, only : Error_Report, errorStatusFail, errorStatusSuccess
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    double precision                                                  , dimension(self%integrandCount) :: error, errorScale, errorsMaximum
    logical                                                           , dimension(self%integrandCount) :: converged, mustEvaluate, convergedPrevious
    double precision :: midpoint, errorMaximum
    type            (intervalMultiList                               ), dimension(:)                 , allocatable :: list
    double precision                                                  , dimension(:)                 , allocatable :: listValue
    integer         (c_size_t                                        ), dimension(:)                 , allocatable :: listRank
    type            (intervalMulti                                   ), pointer :: head, newInterval1, newInterval2, current, previous, newInterval
    integer         (c_size_t) :: iInterval, intervalCount

    ! If the interval has zero size, return a zero result.
    if (a == b) then
       integral=0.0d0
       if (present(status)) status=errorStatusSuccess
       return
    end if
    ! Create our first estimate of the integral, using a single interval.
    intervalCount=1_c_size_t
    head => intervalMultiNew(a,b,self%integrandCount)
    mustEvaluate =  .true.
    call self%evaluateInterval(head%a,head%b,head%integral,head%error,mustEvaluate)
    ! Initialize current integral and error estimates, and flag that we're not converged.
    integral =head%integral
    error    =head%error
    converged=                                                   abs(error) <= self%toleranceAbsolute                .or.                                                  abs(error) <= self%toleranceRelative*abs(integral)
    ! Iterate until convergence is reached.
    do while (.not.all(converged) .and. intervalCount < self%intervalsMaximum)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       intervalCount=intervalCount-1_c_size_t
       midpoint =  (current%b+current%a)/2.0d0
       ! Create two new subintervals.
       newInterval1 => intervalMultiNew(current%a,midpoint ,self%integrandCount)
       newInterval2 => intervalMultiNew(midpoint ,current%b,self%integrandCount)
       ! Check for loss of precision in interval extent.
       if     (                         midpoint==current%a    .or.                     midpoint==current%b   ) call Error_Report("loss of precision in integration interval"//{introspection:location})
       ! Compute the integral and error estimate in each new subinterval.
       mustEvaluate=.not.converged
       call self%evaluateInterval(newInterval1%a,newInterval1%b,newInterval1%integral,newInterval1%error,mustEvaluate)
       call self%evaluateInterval(newInterval2%a,newInterval2%b,newInterval2%integral,newInterval2%error,mustEvaluate)
       ! Update the total integral and error estimate for the new subintervals, for integrands which were evaluated. For other
       ! integrands, the new interval integrals and errors are set to half of that of the parent interval.
       convergedPrevious=converged
       where (mustEvaluate)
          integral                        =integral                -current     %integral   +newInterval1%integral   +newInterval2%integral
          error                           =error                   -current     %error      +newInterval1%error      +newInterval2%error
          ! Test for convergence.
          converged=                                              error <= self%toleranceAbsolute                .or.                                             error <= self%toleranceRelative*abs(integral)
       elsewhere
          newInterval1%integral=current%integral/2.0d0
          newInterval2%integral=current%integral/2.0d0
          newInterval1%error   =current%error   /2.0d0
          newInterval2%error   =current%error   /2.0d0
       end where
       if (.not.all(converged) .and. any(converged .neqv. convergedPrevious).and.intervalCount > 1_c_size_t) then
          ! One or more integrals have converged. We must resort the linked list of intervals to ensure that it is in order of
          ! descending error measure.
          errorScale=max(                                                self%toleranceAbsolute              ,           self%toleranceRelative*abs(integral)           )
          call intervalMultiResort(head,intervalCount,converged,errorScale,list,listValue,listRank)
       end if
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       do iInterval=1,2
          if (iInterval == 1) then
             newInterval => newInterval1
          else
             newInterval => newInterval2
          end if
          ! Check if the stack is empty.
          if (associated(head)) then
             ! Stack is not empty. Perform an insertion sort to insert our new subinterval into the sorted stack. The sorting is
             ! such that the stack is ordered with the interval for which the scaled error (i.e. error divided by the larger of
             ! the absolute and relative tolerances) of non-converged integrals is descending.
             errorScale    =   max(                                                   self%toleranceAbsolute              ,                self%toleranceRelative*abs(integral)                )
             errorMaximum  =   max(0.0d0,maxval(newInterval%error/errorScale,mask=.not.converged))
             errorsMaximum =  +errorScale               *errorMaximum
             current       =>  head
             previous      =>  head
             call intervalMultiSearch(current,previous,converged,errorsMaximum)
             call intervalMultiLink  (head,current,previous,newInterval,newInterval,converged,errorsMaximum)
          else
             ! Stack is empty - simply point the head to our new subinterval.
             head => newInterval
          end if
       end do
       intervalCount=intervalCount+2_c_size_t
    end do
    ! Destroy the stack.
    call intervalMultiListDestroy(head)
    ! Report error if number of intervals was exceeded.
    if (present(status)) status=errorStatusSuccess
    if (intervalCount >= self%intervalsMaximum .and. .not.all(converged)) then
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('maximum number of intervals exceeded'//{introspection:location})
       end if
    end if
    return
  end procedure multiVectorizedCompositeGaussKronrod1DEvaluate

  module procedure multiVectorizedCompositeGaussKronrod1DEvaluateInterval
    !!{RST
    Evaluate the integral over an interval using the Gauss-Kronrod method (also estimates the error). Specific implementation is based on that in the :term:`GSL`.
    !!}
    implicit none
    double precision                                                                 , dimension(self%integrandCount                        ) :: integralGauss, integralAbsolute, mean, integralAsc
    double precision                                                  , pointer      , dimension(:                  ,:                      ) :: fValue1, fValue2
    double precision                                                  , target       , dimension(self%integrandCount,2*size(self%xKronrod)-1) :: fUnion
    double precision                                                                 , dimension(                    2*size(self%xKronrod)-1) :: xUnion
    double precision :: halfLength, halfLengthAbsolute, xCenter
    integer :: pointCountKronrod, pointCountGauss, i

    ! Establish point counts and interval.
    pointCountKronrod =size(self%wKronrod)
    pointCountGauss   =size(self%wGauss  )
    if (mod(pointCountKronrod,2) == 0) pointCountGauss=pointCountGauss-1
    xCenter           =0.5d0*(b+a)
    halfLength        =0.5d0*(b-a)
    halfLengthAbsolute=abs(halfLength)
    ! Evaluate function at all required points.
    xUnion(1                                            )=  xCenter
    xUnion(2                    :2+  pointCountKronrod-2)=  xCenter-self%xKronrod(1:pointCountKronrod-1)*halfLength
    xUnion(2+pointCountKronrod-1:2+2*pointCountKronrod-3)=  xCenter+self%xKronrod(1:pointCountKronrod-1)*halfLength
    call self%integrand(self%integrandCount,xUnion,mustEvaluate,fUnion)
    fValue1                                              => fUnion(:,2                    :2+  pointCountKronrod-2)
    fValue2                                              => fUnion(:,2+pointCountKronrod-1:2+2*pointCountKronrod-3)
    ! Evaluate integrals for which the integrand was computed.
    do i=1,self%integrandCount
       if (.not.mustEvaluate(i)) cycle
       integralGauss(i)=sum(                                                      +self%wGauss(  1:  pointCountGauss  )                 *(                                                      +fValue1  (i,2:2*pointCountGauss:2)                   +fValue2  (i,2:2*pointCountGauss:2)                  )                                                   )
       integralKronrod(i)=+    self%wKronrod  (    pointCountKronrod  )                *            fUnion (i,1                    )                +sum(                                                             +self%wKronrod (  1:pointCountKronrod-1)                     *(                                                             +     fvalue1(i,1:pointCountKronrod-1)                       +     fvalue2(i,1:pointCountKronrod-1)                      )                                                          )
       integralAbsolute(i)=+self%wKronrod     (    pointCountKronrod  )                 *       abs(fUnion (i,1                    ))                +sum(                                                             +self%wKronrod(  1:pointCountKronrod-1)                      *(                                                             +abs(fvalue1(i,1:pointCountKronrod-1))                       +abs(fvalue2(i,1:pointCountKronrod-1))                      )                                                          )
       mean       (i)=+0.5d0                                                                                                     *integralKronrod(i)
       integralAsc(i)=+sum(self%wKronrod(1:pointCountKronrod-1)*(abs(fValue1(i,:)-mean(i))+abs(fValue2(i,:)-mean(i))))           +    self%wKronrod(  pointCountKronrod  )* abs(fUnion (i,1)-mean(i))
       if (mod(pointCountKronrod,2) == 0) integralGauss(i)=integralGauss(i)+fUnion(i,1)*self%wGauss(pointCountKronrod/2)
       ! Evaluate error.
       error           (i)=abs((integralKronrod(i)-integralGauss(i))*halfLength)
       integralKronrod (i)=integralKronrod (i)*halfLength
       integralAbsolute(i)=integralAbsolute(i)*halfLengthAbsolute
       integralAsc     (i)=integralAsc     (i)*halfLengthAbsolute
       error           (i)=gaussKronrodError(error(i),integralAbsolute(i),integralAsc(i))
    end do
    return
  end procedure multiVectorizedCompositeGaussKronrod1DEvaluateInterval

end submodule Numerical_Integration2_Gauss_Kronrod
