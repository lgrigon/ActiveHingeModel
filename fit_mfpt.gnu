# Run this to get the MFPT estimations from SR~exp(-t/MFPT).
# To group all the results run:
# ls | while read p; do cat $p; done

R = '2.0' #'1.4 1.6 1.8 1.9 2.0'
mub = '0.10 0.25 0.50 0.75 1.00 2.00 3.00 4.00 6.00 8.00 10.00'
mup = '0.05'
rho = '0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70'
ang = '1.047' #'0.524 0.628 0.733 0.785 0.838 0.890 0.942 0.995 1.047'
F = '160.0' #'0.0 5.0 10.0 15.0 20.0 25.0 50.0 75.0 100.0 140.0 160.0 170.0'
KBT = '10.0'

Smax_t0 = 0.995
Smax_fit = 0.90
Smin_fit = 0.30

Nw = 10          # number of additional windows (tune if you like)
min_pts = 10     # min points required for a fit in a window
max_shift = 0.5    # do not shift Tmin beyond 50% of the window

aux = 0
f(x) = a*exp(-x/b)
set fit maxiter 100.
set fit quiet
set fit logfile '/dev/null'
set fit errorvariables

file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0

do for [kb in KBT] {
do for [k in mub] {
do for [k2 in mup] {
do for [m in R] {

do for [l in ang] {
do for [j in F] {
do for [i in rho] {


	fname = sprintf("mfpt/surv_L%s-pf%s-F%s-ang%s-mub%s-mup%s-T%s.dat",m,i,j,l,k,k2,kb)




	if ( file_exists(fname) ) {

		t21 = -1.0
		t22 = -1.0
		t22e = 0.0
		t23 = -1.0
		t23e = 0.0
		t0f = 0.0
		t0f1 = 0.0

		if (1==0) {
		### ---- first rod ---- ###

		# check if enough points
		stats fname u 9 nooutput name "E"

		if (E_mean > 0) {

			# t0
			cmd_t0 = sprintf("< awk '($6<%.6f){print $1}' %s", Smax_t0, fname)
			stats cmd_t0 u 1 nooutput name "T0"
			t0 = T0_min

			# setting range
			stats fname u 6 nooutput name "S"
			Smin = (S_min < Smin_fit ? Smin_fit : S_min)
			Smax = (S_max > Smax_fit ? Smax_fit : S_max)

			cmd_s = sprintf("< awk '( $6>%.6f && $6<%.6f ){print $1,$6,$7}' %s", Smin, Smax, fname)
			stats cmd_s u 1 nooutput name "T"
			Tmax = T_max
			Tmin = T_min
			Tmin = (Smax - Smin > 0.1 ? Tmin : t0)

			# initial guess for MFPT
			#stats fname u 3 nooutput name "B"
			b = Tmax - Tmin #T_mean
			a = 1.0
			
			# --- central fit on [Tmin, Tmax] ---
			cmd = sprintf("< awk '( $1>%.6f && $1<%.6f ){print $1,$6,$7}' %s", Tmin, Tmax, fname)
			stats cmd u 3 nooutput name "test"
			if (exists("test_records") && test_records > min_pts) {

				# central fit: this is your original behavior
				a = 1.0
				b = Tmax - Tmin
				fit f(x) cmd u 1:2:3 yerror via a,b

				a0    = a
				b0    = b
				b0err = b_err

				# --- scan over slightly different fit windows to estimate systematic error on b ---

				array bvals[Nw+1]

				nvalid = 0

				# store central fit as first sample
				nvalid = nvalid + 1
				bvals[nvalid] = b0

				# vary Tmin inside [Tmin, Tmax], keeping Tmax fixed
				do for [iw=1:Nw] {
					Tmin_w = Tmin + max_shift*iw*(Tmax - Tmin)/(Nw+1.)

					cmdw = sprintf("< awk '( $1>%.6f && $1<%.6f ){print $1,$6,$7}' %s", Tmin_w, Tmax, fname)
					stats cmdw u 3 nooutput name "TESTW"

					if (exists("TESTW_records") && TESTW_records > min_pts) {
						# use central fit as initial guess
						a = a0
						b = b0
						fit f(x) cmdw u 1:2:3 yerror via a,b

						nvalid = nvalid + 1
						bvals[nvalid] = b
					}
				}

				# compute mean(b) and stddev(b) over all valid windows
				bmean = 0.
				do for [k=1:nvalid] {
					bmean = bmean + bvals[k]
				}
				bmean = bmean / nvalid

				if (nvalid > 1) {
					bvar = 0.
					do for [k=1:nvalid] {
						bvar = bvar + (bvals[k] - bmean)**2
					}
					bstd = sqrt( bvar / (nvalid - 1) )
				} else {
					# fall back to gnuplot's covariance error if only one window worked
					bstd = b0err
				}

				# use central fit to define t0f (same as before)
				t0f1 = b0*log(a0)

				# MFPT estimate and its error
				t22  = bmean
				t22e = sqrt( b0err**2 + bstd**2 )

				# optional plot: central window + mean fit (they'll be very close)
				#p cmd w errorbars notitle, \
				#  f(x) w l title sprintf("%s -> %.1f", fname, t23)

				#pause 3.0

			} else {
				t22  = -Tmin
				t22e = 0. 
			}
		}
		}

		### ---- second rod ---- ###

		# check if enough points
		stats fname u 5 nooutput name "E"

		if (E_mean > 0) {

			# t0
			cmd_t0 = sprintf("< awk '($2<%.6f){print $1}' %s", Smax_t0, fname)
			stats cmd_t0 u 1 nooutput name "T0"
			t0 = T0_min

			# setting range
			stats fname u 2 nooutput name "S"
			Smin = (S_min < Smin_fit ? Smin_fit : S_min)
			Smax = (S_max > Smax_fit ? Smax_fit : S_max)

			cmd_s = sprintf("< awk '( $2>%.6f && $2<%.6f ){print $1,$2,$3}' %s", Smin, Smax, fname)
			stats cmd_s u 1 nooutput name "T"
			Tmax = T_max
			Tmin = T_min
			Tmin = (Smax - Smin > 0.1 ? Tmin : t0)

			# initial guess for MFPT
			#stats fname u 3 nooutput name "B"
			b = Tmax - Tmin #T_mean
			a = 1.0
			
			# --- central fit on [Tmin, Tmax] ---
			cmd = sprintf("< awk '( $1>%.6f && $1<%.6f ){print $1,$2,$3}' %s", Tmin, Tmax, fname)
			stats cmd u 3 nooutput name "test"
			if (exists("test_records") && test_records > min_pts) {

				# central fit: this is your original behavior
				a = 1.0
				b = Tmax - Tmin
				fit f(x) cmd u 1:2:3 yerror via a,b

				a0    = a
				b0    = b
				b0err = b_err

				# --- scan over slightly different fit windows to estimate systematic error on b ---

				array bvals[Nw+1]

				nvalid = 0

				# store central fit as first sample
				nvalid = nvalid + 1
				bvals[nvalid] = b0

				# vary Tmin inside [Tmin, Tmax], keeping Tmax fixed
				do for [iw=1:Nw] {
					Tmin_w = Tmin + max_shift*iw*(Tmax - Tmin)/(Nw+1.)

					cmdw = sprintf("< awk '( $1>%.6f && $1<%.6f ){print $1,$2,$3}' %s", Tmin_w, Tmax, fname)
					stats cmdw u 3 nooutput name "TESTW"

					if (exists("TESTW_records") && TESTW_records > min_pts) {
						# use central fit as initial guess
						a = a0
						b = b0
						fit f(x) cmdw u 1:2:3 yerror via a,b

						nvalid = nvalid + 1
						bvals[nvalid] = b
					}
				}

				# compute mean(b) and stddev(b) over all valid windows
				bmean = 0.
				do for [k=1:nvalid] {
					bmean = bmean + bvals[k]
				}
				bmean = bmean / nvalid

				if (nvalid > 1) {
					bvar = 0.
					do for [k=1:nvalid] {
						bvar = bvar + (bvals[k] - bmean)**2
					}
					bstd = sqrt( bvar / (nvalid - 1) )
				} else {
					# fall back to gnuplot's covariance error if only one window worked
					bstd = b0err
				}

				# use central fit to define t0f (same as before)
				t0f = b0*log(a0)

				# MFPT estimate and its error
				t23  = bmean
				t23e = sqrt( b0err**2 + bstd**2 )

				# optional plot: central window + mean fit (they'll be very close)
				p cmd w errorbars notitle, \
				  f(x) w l title sprintf("%s -> %.1f", fname, t23)

				#pause 3.0

			} else {
				t23  = -Tmin
				t23e = 0. 
			}
		}

		aux1 = (aux<10? "00".aux : aux<100 ? "0".aux : aux)

		set print "temp/".aux1.".dat"; print sprintf("%s %s %s %s %s %s %s %f %f %f %f %f %f",i,l,j,m,k,k2,kb,t23,t0f,t23e,t22,t0f1,t22e)
		aux = aux + 1

	}

}
}
}
}
}
}
}