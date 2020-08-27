import mysubsprograms as my

evolve_profile = 'testing'
inlist_evolve = 'inlist_evolve_10.0_0.005_0.24_0.02_0.25_7.4_0.1'
irrad_mod = 'irrad_10.0_0.005_0.24_0.02_0.25_7.4.mod'
evolve_mod = 'testing'
n_frac = 0.1
a = 1.0
ms = 1.0
orb_sep = 0.1
ec = 1e9
column_depth = 25
flux_dayside = 20357683.70173884
formation_time = 6e6
teq = 547.35
BA = 0.2
escape_regime = 0
diff_sep = 1
homopause_temp = 10000

my.run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod,
			n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
			formation_time, teq, BA, escape_regime, diff_sep, homopause_temp)