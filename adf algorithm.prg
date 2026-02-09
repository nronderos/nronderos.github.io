'Nicolas Ronderos Pulido. November de 2014.
'--------------------------Inputs-------------------------------
%series="y" 'Series name
scalar significance=0.01 'Significance levels for the tests 0.1,0.05,0.025 or 0.01
'------------------------------------------------------------------
'Step 1
matrix results
freeze(uroot_ct){%series}.uroot(trend,info=aic,save=results)
if uroot_ct(7,5)<=significance then 'If the null is rejected (there is not unit root)
	@uiprompt("The series does not contain a unit root")
      delete results significance uroot_ct
	stop
endif

'Step 2

!lags=results(2,1)
for !i=1 to !lags
	%lag=@str(-!i)
	%lags=%lags+"@d("+%series+"("+%lag+"))"
next
equation restricted.ls @d({%series}) c {%lags}
!ssr_r=@ssr
equation unrestricted.ls @d({%series}) {%series}(-1) {%lags} c @trend 
freeze(unrestricted_tc) unrestricted
!ssr_ur=@ssr
scalar phi_3=((!ssr_r-!ssr_ur)/2)/(!ssr_ur/(results(1,1)-!lags-3))
scalar obs=@obs({%series})
vector critical_values
call cv(obs,significance,critical_values)
if phi_3>critical_values(1,3) then ' If the null is rejected (the trend and gamma are non-zero)
	equation t_test1.ls @d({%series}) c @trend {%lags}
	freeze(t_test_1) t_test1
	if t_test_1(10,5)<=significance then 'If the null is rejected (the trend is significative)
		if unrestricted_tc(9,5)<=significance then ' If the null is rejected (gamma is non-zero)
			@uiprompt("The series is trend-stationary")
			stop
		else
			@uiprompt("The series contain a unit root and a quadratic trend") '
                 delete critical_values obs phi_3 restricted results significance t_test1 t_test_1 unrestricted unrestricted_tc uroot_ct
			stop
		endif
	endif
endif

'Step 3

freeze(uroot_c){%series}.uroot(const,info=aic,save=results)
if uroot_c(7,5)<=significance then 'If the null is rejected (there is not unit root)
	@uiprompt("The series does not contain a unit root")
	stop
endif

if !lags=0 then
	genr first=@d({%series})
	!ssr_r=@sumsq(first)
else
	equation restricted.ls @d({%series}) {%lags}
	!ssr_r=@ssr
endif
equation unrestricted.ls @d({%series}) c {%series}(-1) {%lags}
freeze(unrestricted_c) unrestricted
!ssr_ur=@ssr
scalar phi_1=((!ssr_r-!ssr_ur)/2)/(!ssr_ur/(results(1,1)-!lags-2))
if phi_1>critical_values(1,1) then 'If the null is rejected (the constant and gamma are non-zero)
	equation t_test2.ls @d({%series}) c {%lags}
	freeze(t_test_2) t_test2
	if t_test_2(9,5)<=significance then 'If the null is rejected (the constant is significative)
		if unrestricted_c(10,5)<=significance then 'If the null is rejected (gamma is non-zero)
			@uiprompt("The series is stationary around a nonzero mean")
			stop
		else
			@uiprompt("The series contain a unit root and a drift") '
                delete critical_values obs phi_1 phi_3 restricted significance t_test2 t_test_2 unrestricted unrestricted_c unrestricted_tc uroot_c uroot_ct results
			stop
		endif
	endif
endif

'Step 4

freeze(uroot){%series}.uroot(none,info=aic,save=results)
if uroot(7,5)<=significance then  'If the null is rejected (there is not unit root)
	@uiprompt("The series does not contain a unit root")
	stop
else
	@uiprompt("The series contain a unit root") '
     'delete critical_values obs first phi_1 phi_3 restricted results significance unrestricted unrestricted_c unrestricted_tc uroot uroot_c uroot_ct
	stop
endif


'Critical values subroutine
subroutine local cv(scalar obs,scalar significance,vector critical_values)
	genr n_=NA
	n_.fill 25,50,100,250,500,1000
	stom(n_,n)
	vector(6) nearby=@abs(n-obs)
	
	for !i=1 to 6
		if nearby(!i,1)=@min(nearby) then
			!position=!i
		endif
	next
	
	genr criticalphi1=NA
	genr criticalphi2=NA
	genr criticalphi3=NA
	
	if significance=0.1 then
		criticalphi1.fill 4.12,3.94,3.86,3.81,3.79,3.78
		criticalphi2.fill 4.67,4.31,4.16,4.07,4.05,4.03
		criticalphi3.fill 5.91,5.61,5.47,5.39,5.36,5.34
		stom(criticalphi1,critical_phi1)
		stom(criticalphi2,critical_phi2)
		stom(criticalphi3,critical_phi3)
	endif
	
	if significance=0.05 then
		criticalphi1.fill 5.18,4.86,4.71,4.63,4.61,4.59
		criticalphi2.fill 5.68,5.13,4.88,4.75,4.71,4.68
		criticalphi3.fill 7.24,6.73,6.49,6.34,6.30,6.25
		stom(criticalphi1,critical_phi1)
		stom(criticalphi2,critical_phi2)
		stom(criticalphi3,critical_phi3)
	endif
	
	if significance=0.025 then
		criticalphi1.fill 6.30,5.80,5.57,5.45,5.41,5.38
		criticalphi2.fill 6.75,5.94,5.59,5.40,5.35,5.31
		criticalphi3.fill 8.65,7.81,7.44,7.25,7.20,7.16
		stom(criticalphi1,critical_phi1)
		stom(criticalphi2,critical_phi2)
		stom(criticalphi3,critical_phi3)
	endif
	
	if significance=0.01 then
		criticalphi1.fill 7.88,7.06,6.70,6.52,6.47,6.43
		criticalphi2.fill 8.21,7.02,6.50,6.22,6.15,6.09
		criticalphi3.fill 10.61,9.31,8.73,8.43,8.34,8.27
		stom(criticalphi1,critical_phi1)
		stom(criticalphi2,critical_phi2)
		stom(criticalphi3,critical_phi3)
	endif
	
	vector(1,3) critical_values
	for !i=1 to 3
		critical_values(1,!i)=critical_phi!i(!position)
	next
	
endsub
