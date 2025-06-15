# Run this to get the T_N estimations from SR(N)~exp(-t/T_N).
# To group all the results run:
# ls | while read p; do cat $p; done

R = '2.0'
D = '200'
rho = '0.75'
ang = '1.047'
F = '0 15 30 45 60 75 90'

dt = 0.00005*1000

aux = 0
f(x) = a*exp(-x/b)
set fit maxiter 100.
set fit quiet
set fit logfile '/dev/null'
set fit errorvariables

fname = "mfpt/mfpt_terms_sr.dat"

do for [i in rho] {
do for [l in ang] {
do for [j in F] {

	#"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2) {print $7,$8}'  ".fname

	t0 = 0.0
	t1 = 0.0
	t2 = 0.0
	t3 = 0.0
	t4 = 0.0

	#N = 0
	#stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==0&&$8>0.0) {print $7}'  ".fname u 1 nooutput
	#tmax = STATS_max
	stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==0) {print $8}'  ".fname using (column(1) != 0.0 ? 1 : 0) name 'non_zero' nooutput
	tmax = dt*(non_zero_sum-1)
	if (tmax>15*dt) {
		stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==0&&$8>0.0) {print $7}'  ".fname using 1 name 'data' nooutput
		b = data_mean
		a = 1.0

		b = 5

		fit f(x) "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==0&&$8<1.0 && $7>".sprintf("%.1f",tmax*0.1)." && $7<".sprintf("%.1f",0.6*tmax).") {print $7,$8}'  ".fname u 1:2 via a,b
		t0 = b

		p [0:tmax]"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==0&&$8<1.1) {print $7,$8}'  ".fname u 1:2 notitle, f(x) w l title sprintf("0 %f",b)
	}

	#N = 1
	#stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==1&&$8>0.0) {print $7}'  ".fname u 1 nooutput
	#tmax = STATS_max
	stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==1) {print $8}'  ".fname using (column(1) != 0.0 ? 1 : 0) name 'non_zero' nooutput
	tmax = dt*(non_zero_sum-1)
	if (tmax>15*dt) {
		stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==1&&$8>0.0) {print $7}'  ".fname using 1 name 'data' nooutput
		b = data_mean
		a = 1.0

		fit f(x) "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==1&&$8<1.0 && $7>".sprintf("%.1f",tmax*0.1)." && $7<".sprintf("%.1f",0.6*tmax).") {print $7,$8}'  ".fname u 1:2 via a,b
		t1 = b

		p [0:tmax]"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==1&&$8<1.1) {print $7,$8}'  ".fname u 1:2 notitle, f(x) w l title sprintf("1 %f",b)
	}

	#N = 2
	#stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2&&$8>0.0) {print $7}'  ".fname u 1 nooutput
	#tmax = STATS_max
	stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2) {print $8}'  ".fname using (column(1) != 0.0 ? 1 : 0) name 'non_zero' nooutput
	tmax = dt*(non_zero_sum-1)
	if (tmax>15*dt) {
		stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2&&$8>0.0) {print $7}'  ".fname using 1 name 'data' nooutput
		b = data_mean
		a = 1.0

		fit f(x) "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2&&$8<1.0 && $7>".sprintf("%.1f",tmax*0.1)." && $7<".sprintf("%.1f",0.6*tmax).") {print $7,$8}'  ".fname u 1:2 via a,b
		t2 = b

		p [0:tmax]"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==2&&$8<1.1) {print $7,$8}'  ".fname u 1:2 notitle, f(x) w l title sprintf("2 %f",b)
	}

	#N = 3
	#stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3&&$8>0.0) {print $7}'  ".fname u 1 nooutput
	#tmax = STATS_max
	stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3) {print $8}'  ".fname using (column(1) != 0.0 ? 1 : 0) name 'non_zero' nooutput
	tmax = dt*(non_zero_sum-1)
	if (tmax>15*dt) {
		stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3&&$8>0.0) {print $7}'  ".fname using 1 name 'data' nooutput
		b = data_mean
		a = 1.0

		fit f(x) "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3&&$8<1.0 && $7>".sprintf("%.1f",tmax*0.1)." && $7<".sprintf("%.1f",0.6*tmax).") {print $7,$8}'  ".fname u 1:2 via a,b
		t3 = b

		p [0:tmax]"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3&&$8<1.1) {print $7,$8}'  ".fname u 1:2 notitle, f(x) w l title sprintf("3 %f",b)
	}

	#N = all
	#stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==3&&$8>0.0) {print $7}'  ".fname u 1 nooutput
	#tmax = STATS_max
	stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==4) {print $8}'  ".fname using (column(1) != 0.0 ? 1 : 0) name 'non_zero' nooutput
	tmax = dt*(non_zero_sum-1)
	if (tmax>15*dt) {
		stats "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==4&&$8>0.0) {print $7}'  ".fname using 1 name 'data' nooutput
		b = data_mean
		a = 1.0

		fit f(x) "<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==4&&$8<1.0 && $7>".sprintf("%.1f",tmax*0.1)." && $7<".sprintf("%.1f",0.6*tmax).") {print $7,$8}'  ".fname u 1:2 via a,b
		t4 = b

		p [0:tmax]"<awk '($1==".i."&&$2==".l."&&$3==".j."&&$6==4&&$8<1.1) {print $7,$8}'  ".fname u 1:2 notitle, f(x) w l title sprintf("3 %f",b)
	}

	aux1 = (aux<10? "00".aux : aux<100 ? "0".aux : aux)

	set print "temp/".aux1.".dat"; print sprintf("%s %s %s %s %s %f %f %f %f %f",i,l,j,R,D,t0,t1,t2,t3,t4)
	aux = aux + 1

}
}
}
