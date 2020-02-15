import mysubsprograms as my

evolve_profile = 'profile_evolve10.0_0.01_0.24_0.02_0.1_8.0_0.1'
inlist_evolve = 'inlist_evolve_10.0_0.01_0.24_0.02_0.1_8.0_0.1'
irrad_mod = 'irrad_10.0_0.01_0.24_0.02_0.1_8.0.mod'
evolve_mod = 'Test_BAD'
n_frac = 0.1
a = 1.0
ms = 1.0
orb_sep = 0.1
ec = 1e9
column_depth = 50
flux_dayside = 1164284558.89515 / 1000
formation_time = 6e6
teq = 1500
BA = 0.2
escape_regime = 0
diff_sep = 1

my.run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod,
			n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
			formation_time, teq, BA, escape_regime, diff_sep)