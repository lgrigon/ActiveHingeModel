# Run this to get the MFPT estimations from SR~exp(-t/MFPT).
# To group all the results run:
# ls | while read p; do cat $p; done

rho = '0.50 0.55 0.60 0.65 0.70 0.75'
F = '0 30 60 90'
D = '200'
ang = '0.785 1.047'
R = '2.0'

aux = 0
f(x) = a*exp(-x/b)
set fit maxiter 100.
set fit quiet
set fit logfile '/dev/null'

do for [i in rho] {
do for [l in ang] {
do for [j in F] {
do for [k in D] {

	fname = sprintf("mfpt/surv_R%s-pf%s-f%s-ang%s-D%s.dat",R,i,j,l,k)

	a = 1.0
	stats fname u 5 nooutput
	norm = STATS_mean
	stats fname u 6 nooutput
	nabs = STATS_mean
	stats fname u 8 nooutput
	tini = STATS_mean
	stats fname u 9 nooutput
	m1 = STATS_mean
	b = m1

	t11 = -1.0
	t12 = -1.0
	if (nabs > 5) {
		fit f(x) "<awk '($2>0.0) {print $1,$2}'  ".fname u 1:2 via a,b
		t11 = b
		fit f(x) "<awk '($2<1.0 && $2>0.0) {print $1,$2}'  ".fname u 1:2 via a,b
		t12 = b + tini
	}

	a = 1.0
	stats fname u 6 nooutput
	norm = STATS_mean
	stats fname u 7 nooutput
	nabs = STATS_mean
	stats fname u 10 nooutput
	tini = STATS_mean
	stats fname u 11 nooutput
	m2 = STATS_mean
	b = m2

	t21 = -1.0
	t22 = -1.0
	if (nabs > 5) {
		fit f(x) "<awk '($3>0.0) {print $1,$3}'  ".fname u 1:2 via a,b
		t21 = b
		fit f(x) "<awk '($3<1.0 && $3>0.0) {print $1,$3}'  ".fname u 1:2 via a,b
		t22 = b + tini
	}

	aux1 = (aux<10? "00".aux : aux<100 ? "0".aux : aux)

	set print "temp/".aux1.".dat"; print sprintf("%s %s %s %s %s %f %f %f %f %f %f",i,l,j,R,k,t11,t12,t21,t22,m1,m2)
	aux = aux + 1

}
}
}
}