genr serie_1=x-@mean(x)
stom(serie_1,serie1)
scalar varianza1=@var(serie1)
scalar lags=@rows(serie1)-1

genr serie_2=y-@mean(y)
stom(serie_2,serie2)
scalar varianza2=@var(serie2)
scalar lags=@rows(serie2)-1

group rezagos1 serie_1(1)
group rezagos2 serie_2(1)
for !i=2 to lags
    group rezagos1 rezagos1 serie_1(!i)
    group rezagos2 rezagos2 serie_2(!i)
next

stomna(rezagos1,matriz1)
stomna(rezagos2,matriz2)

for !i=1 to lags+1
    for !j=1 to lags
        if matriz1(!i,!j)=NA then
            matriz1(!i,!j)=0
        endif
        if matriz2(!i,!j)=NA then
            matriz2(!i,!j)=0
        endif
     next
next


%tabname = "Autocorrelacion"
if @isobject(%tabname) then 
	%tabname = @getnextname(%tabname)
endif
table {%tabname}
{%tabname}(1,1)="Rezagos"
{%tabname}(1,2)="Rezagada"
{%tabname}(1,3)="Adelantada"
{%tabname}(2,1)=0
matrix zero1=(1/(lags+1)*@transpose(serie1)*serie2)/(@sqrt(varianza1*varianza2))
matrix zero2=(1/(lags+1)*@transpose(serie2)*serie1)/(@sqrt(varianza1*varianza2))
{%tabname}(2,2)=zero1(1,1)
{%tabname}(2,3)=zero2(1,1)
for !i=1to lags
     vector r_1!i=(1/(lags+1)*@transpose(serie1)*@columnextract(matriz2,!i))/(@sqrt(varianza1*varianza2))
     vector r_2!i=(1/(lags+1)*@transpose(serie2)*@columnextract(matriz1,!i))/(@sqrt(varianza1*varianza2))
     {%tabname}(!i+2,1)=!i
     matrix(lags,2) provisional(!i,2)=r_1!i(1,1)
     provisional(!i,1)=r_2!i(1,1)
     {%tabname}(!i+2,2)=provisional(!i,1)
     {%tabname}(!i+2,3)=provisional(!i,2)
     delete r_1!i
     delete r_2!i
next

delete matriz1
delete matriz2
delete lags
delete rezagos1
delete rezagos2
delete serie_1
delete serie_2
delete varianza1
delete varianza2
delete serie1
delete serie2
delete provisional
delete zero1
delete zero2
