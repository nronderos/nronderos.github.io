wfcreate u 1 200
genr ser01=nrnd
genr serie=ser01-@mean(ser01)
stom(serie,serie1)
scalar varianza=@var(serie1)
scalar lags=@rows(serie1)-1

group rezagos serie(1)
for !i=2 to lags
group rezagos rezagos serie(!i)
next

stomna(rezagos,matriz)

for !i=1 to lags+1
    for !j=1 to lags
        if matriz(!i,!j)=NA then
            matriz(!i,!j)=0
        endif
     next
next

%tabname = "Autocorrelacion"
if @isobject(%tabname) then 
	%tabname = @getnextname(%tabname)
endif
table {%tabname}
{%tabname}(1,1)="Rezagos"
{%tabname}(1,2)="ACF"
%tabname01 = "Autocorrelacion_Estrella"
if @isobject(%tabname01) then 
	%tabname01 = @getnextname(%tabname01)
endif
table {%tabname01}
{%tabname01}(1,1)="Rezagos"
{%tabname01}(1,2)="ACF"
for !i=1 to lags
	'Gorro
     vector r_!i=(1/(lags+1)*@transpose(serie1)*@columnextract(matriz,!i))/varianza
     {%tabname}(!i+1,1)=!i
     matrix(lags,1) provisional01(!i,1)=r_!i(1,1)
     {%tabname}(!i+1,2)=provisional01(!i,1)
     delete r_!i
	'Estrella
     vector r_!i=(1/(lags+1-!i)*@transpose(serie1)*@columnextract(matriz,!i))/varianza
	{%tabname01}(!i+1,1)=!i
     matrix(lags,1) provisional02(!i,1)=r_!i(1,1)
     {%tabname01}(!i+1,2)=provisional02(!i,1)
     delete r_!i
next
matrix diferencia=provisional01-provisional02
pagestruct(end=199)
mtos(diferencia,diferencial)
show ser02.line
'delete matriz
'delete lags
'delete rezagos
'delete serie 
'delete varianza
'delete serie1 
'delete provisional


