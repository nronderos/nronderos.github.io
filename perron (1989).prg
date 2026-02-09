'Nicolas Ronderos Pulido
'-------------------------------------------------------Inputs---------------------------------------------------
!quiebre=300 'Break date
%series="ser02" 'Series object
!lags=20 'Maximim lag
%info_criteria="aic" 'schwarz, aic or hq      
!caso=2 'case 1,2 or 3
!significancia=0.1 '0.01,0.025,0.05 y 0.1
'------------------------------------------------------------------------------------------------------------------
genr d_t=0
genr du=0
genr d_tb=0
for !t=1 to @obsrange
	if !caso=1 or !caso=3 then
		if !t=!quiebre+1 then
			d_t(!t)=1 'Cambio en la tendencia en solo un momento del tiempo - D(TB)t
			d_tb(!t)=1
		endif
	endif
	if !t>!quiebre then
		if !caso=2 then
			d_t(!t)=!t-!quiebre 'Cambio en la pendiente - DT*t
			else if !caso=3 then
			d_t(!t)=!t ' Cambio en la pendiente y en el nivel - DTt
			endif
		endif
		du(!t)=1 ' DUt
	endif
next

for !i=1 to !lags
	%lag=@str(-!i)
	%lags=%lags+"@d("+%series+"("+%lag+"))"
	if !caso=3 then
		equation caso_!caso.ls {%series} c du @trend d_t {%series}(-1) d_tb {%lags}
	else
		equation caso_!caso.ls {%series} c du @trend d_t {%series}(-1) {%lags}
	endif
	vector(!lags) information_criteria(!i)=@{%info_criteria}
next
for !i=1 to !lags
	%lag_=@str(-!i)
	%lags_=%lags_+"@d("+%series+"("+%lag_+"))"
	if information_criteria(!i)=@min(information_criteria) then
		if !caso=3 then
			equation caso_!caso.ls {%series} c du @trend d_t {%series}(-1) d_tb {%lags_}
		else
			equation caso_!caso.ls {%series} c du @trend d_t {%series}(-1) {%lags_}
		endif
		exitloop
	endif
next

scalar t=(c(5)-1)/@stderrs(5)
scalar lamda=!quiebre/@obsrange
scalar caso=!caso
scalar significance=@floor(100*!significancia)
scalar cv
call cv(lamda,caso,significance,cv)

if !caso=1 then
	if cv<t then
		@uiprompt("Unit root with level change")
	else
		@uiprompt("No unit root")
	endif
endif

if !caso=2 then
	if cv<t then
		@uiprompt("Unit root with drift change")
	else
		@uiprompt("No unit root")
	endif
endif

if !caso=3 then
	if cv<t then
		@uiprompt("Unit root with level and drift change")
	else
		@uiprompt("No unit root")
	endif
endif

subroutine local cv(scalar lamda,scalar caso,scalar significance,scalar cv)
	genr lamdas_=NA
	lamdas_.fill 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
	stom(lamdas_,lamdas)
	vector(9) nearby=@abs(lamdas-lamda)
	
	for !i=1 to 9
		if nearby(!i,1)=@min(nearby) then
			!position=!i
		endif
	next
	
	genr critical1=NA
	genr critical2=NA
	genr critical5=NA
	genr critical10=NA
	
	if caso=1 then
		critical1.fill -4.30,-4.39,-4.39,-4.34,-4.32,-4.45,-4.42,-4.33,-4.27
		critical2.fill -3.93,-4.08,-4.03,-4.01,-4.01,-4.09,-4.07,-3.99,-3.97
		critical5.fill -3.68,-3.77,-3.76,-3.72,-3.76,-3.76,-3.80,-3.75,-3.69
		critical10.fill -3.40,-3.47,-3.46,-3.44,-3.46,-3.47,-3.51,-3.46,-3.38
		stom(critical1,critical_1)
		stom(critical2,critical_2)
		stom(critical5,critical_5)
		stom(critical10,critical_10)
	endif
	
	if caso=2 then
		critical1.fill -4.27,-4.41,-4.51,-4.55,-4.56,-4.57,-4.51,-4.38,-4.26
		critical2.fill -3.94,-4.08,-4.17,-4.2,-4.26,-4.2,-4.13,-4.07,-3.96
		critical5.fill -3.65,-3.8,-3.87,-3.94,-3.96,-3.95,-3.85,-3.82,-3.68
		critical10.fill -3.36,-3.49,-3.58,-3.66,-3.68,-3.66,-3.57,-3.5,-3.35
		stom(critical1,critical_1)
		stom(critical2,critical_2)
		stom(critical5,critical_5)
		stom(critical10,critical_10)
	endif
	
	if caso=3 then
		critical1.fill -4.38,-4.65,-4.78,-4.81,-4.9,-4.88,-4.75,-4.7,-4.41
		critical2.fill -4.01,-4.32,-4.46,-4.48,-4.53,-4.49,-4.44,-4.31,-4.1
		critical5.fill -3.75,-3.99,-4.17,-4.22,-4.24,-4.24,-4.18,-4.04,-3.8
		critical10.fill -3.45,-3.66,-3.87,-3.95,-3.96,-3.95,-3.86,-3.69,-3.46
		stom(critical1,critical_1)
		stom(critical2,critical_2)
		stom(critical5,critical_5)
		stom(critical10,critical_10)
	endif
	
	!significance_=significance
	scalar cv=critical_!significance_(!position)
	
endsub
'------------------------------------------------------------------------------------------------------------------


