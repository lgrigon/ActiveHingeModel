# Run this to get the MFPT estimations from SR~exp(-t/MFPT).
# To group all the results run:
# ls | while read p; do cat $p; done

R = '2.0'
D = '200'
rho = '0.75'
ang = '1.047'
F = '0 15 30 45 60 75 90'

aux = 0
f(x) = a*exp(-x/b)
set fit maxiter 100.
set fit quiet
set fit logfile '/dev/null'
set fit errorvariables

do for [i in rho] {
do for [l in ang] {
do for [j in F] {
do for [k in D] {

	fname = sprintf("mfpt/surv_R%s-pf%s-f%s-ang%s-D%s.dat",R,i,j,l,k)

	a = 1.0
	stats fname u 3 nooutput
	nabs = STATS_mean
	stats fname u 4 nooutput
	tini = STATS_mean
	stats fname u 5 nooutput
	m2 = STATS_mean
	b = m2

	t21 = -1.0
	t22 = -1.0
	t23 = -1.0
	t23_err = 1.0
	if (nabs > 5) {

		#fit f(x) "<awk '($2>0.0) {print $1,$2}'  ".fname u 1:2 via a,b
		#t21 = b
		#fit f(x) "<awk '($2<1.0 && $2>0.0) {print $1,$2}'  ".fname u 1:2 via a,b
		#t22 = b #+ tini

		stats "<awk '($2>0.0) {print $1,$2}'  ".fname u 1 nooutput
		tmax = STATS_max

		#fit f(x) "<awk '($2<1.0 && $1>".sprintf("%.1f",tmax*0.2)." && $1<".sprintf("%.1f",0.5*tmax).") {print $1,$2}'  ".fname u 1:2 via a,b
		fit f(x) "<awk '($2<1.0 && $1>".sprintf("%.1f",tmax*0.2)." && $1<".sprintf("%.1f",0.5*tmax).") {print $1,$2,$6}'  ".fname u 1:2:3 yerror via a,b
		t23 = b
		t23e = b_err

		p "<awk '($2<1.0 && $1>".sprintf("%.1f",tmax*0.2)." && $1<".sprintf("%.1f",0.5*tmax).") {print $1,$2,$6}'  ".fname u 1:2:3 w errorbars notitle, f(x) w l title sprintf("%f",b)
	}

	aux1 = (aux<10? "00".aux : aux<100 ? "0".aux : aux)

	#set print "temp/".aux1.".dat"; print sprintf("%s %s %s %s %s %f %f %f %f %f %f %f %f",i,l,j,R,k,0.0,0.0,0.0,t21,t22,t23,0.0,m2)
	set print "temp/".aux1.".dat"; print sprintf("%s %s %s %s %s %f %f",i,l,j,R,k,t23,t23e)
	aux = aux + 1

}
}
}
}