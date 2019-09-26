C There are a total of 70 entries in the algebraic variable array.
C There are a total of 19 entries in each of the rate and state variable arrays.
C There are a total of 53 entries in the constant variable array.
C
C
C VOI is time in component environment (millisecond).
C STATES(1) is V in component membrane (millivolt).
C CONSTS(1) is R in component membrane (joule_per_mole_kelvin).
C CONSTS(2) is T in component membrane (kelvin).
C CONSTS(3) is F in component membrane (coulomb_per_millimole).
C CONSTS(4) is Cm in component membrane (microF).
C CONSTS(5) is V_c in component membrane (micrometre3).
C ALGBRC(48) is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
C ALGBRC(55) is i_to in component transient_outward_current (picoA_per_picoF).
C ALGBRC(49) is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
C ALGBRC(50) is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
C ALGBRC(53) is i_CaL in component L_type_Ca_current (picoA_per_picoF).
C ALGBRC(56) is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
C ALGBRC(51) is i_Na in component fast_sodium_current (picoA_per_picoF).
C ALGBRC(52) is i_b_Na in component sodium_background_current (picoA_per_picoF).
C ALGBRC(57) is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
C ALGBRC(54) is i_b_Ca in component calcium_background_current (picoA_per_picoF).
C ALGBRC(59) is i_p_K in component potassium_pump_current (picoA_per_picoF).
C ALGBRC(58) is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
C ALGBRC(13) is i_Stim in component membrane (picoA_per_picoF).
C CONSTS(6) is stim_start in component membrane (millisecond).
C CONSTS(7) is stim_period in component membrane (millisecond).
C CONSTS(8) is stim_duration in component membrane (millisecond).
C CONSTS(9) is stim_amplitude in component membrane (picoA_per_picoF).
C ALGBRC(26) is E_Na in component reversal_potentials (millivolt).
C ALGBRC(34) is E_K in component reversal_potentials (millivolt).
C ALGBRC(42) is E_Ks in component reversal_potentials (millivolt).
C ALGBRC(44) is E_Ca in component reversal_potentials (millivolt).
C CONSTS(10) is P_kna in component reversal_potentials (dimensionless).
C CONSTS(11) is K_o in component potassium_dynamics (millimolar).
C CONSTS(12) is Na_o in component sodium_dynamics (millimolar).
C STATES(2) is K_i in component potassium_dynamics (millimolar).
C STATES(3) is Na_i in component sodium_dynamics (millimolar).
C CONSTS(13) is Ca_o in component calcium_dynamics (millimolar).
C STATES(4) is Ca_i in component calcium_dynamics (millimolar).
C CONSTS(14) is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
C ALGBRC(47) is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
C ALGBRC(45) is alpha_K1 in component inward_rectifier_potassium_current (dimensionless).
C ALGBRC(46) is beta_K1 in component inward_rectifier_potassium_current (dimensionless).
C CONSTS(15) is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
C STATES(5) is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
C STATES(6) is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
C ALGBRC(1) is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
C ALGBRC(14) is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
C ALGBRC(27) is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
C ALGBRC(35) is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
C ALGBRC(2) is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
C ALGBRC(15) is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
C ALGBRC(28) is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
C ALGBRC(36) is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
C CONSTS(16) is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
C STATES(7) is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
C ALGBRC(3) is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
C ALGBRC(16) is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
C ALGBRC(29) is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
C ALGBRC(37) is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
C CONSTS(17) is g_Na in component fast_sodium_current (nanoS_per_picoF).
C STATES(8) is m in component fast_sodium_current_m_gate (dimensionless).
C STATES(9) is h in component fast_sodium_current_h_gate (dimensionless).
C STATES(10) is j in component fast_sodium_current_j_gate (dimensionless).
C ALGBRC(4) is m_inf in component fast_sodium_current_m_gate (dimensionless).
C ALGBRC(17) is alpha_m in component fast_sodium_current_m_gate (dimensionless).
C ALGBRC(30) is beta_m in component fast_sodium_current_m_gate (dimensionless).
C ALGBRC(38) is tau_m in component fast_sodium_current_m_gate (millisecond).
C ALGBRC(5) is h_inf in component fast_sodium_current_h_gate (dimensionless).
C ALGBRC(18) is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
C ALGBRC(31) is beta_h in component fast_sodium_current_h_gate (per_millisecond).
C ALGBRC(39) is tau_h in component fast_sodium_current_h_gate (millisecond).
C ALGBRC(6) is j_inf in component fast_sodium_current_j_gate (dimensionless).
C ALGBRC(19) is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
C ALGBRC(32) is beta_j in component fast_sodium_current_j_gate (per_millisecond).
C ALGBRC(40) is tau_j in component fast_sodium_current_j_gate (millisecond).
C CONSTS(18) is g_bna in component sodium_background_current (nanoS_per_picoF).
C CONSTS(19) is g_CaL in component L_type_Ca_current (litre_per_farad_second).
C STATES(11) is Ca_ss in component calcium_dynamics (millimolar).
C STATES(12) is d in component L_type_Ca_current_d_gate (dimensionless).
C STATES(13) is f in component L_type_Ca_current_f_gate (dimensionless).
C STATES(14) is f2 in component L_type_Ca_current_f2_gate (dimensionless).
C STATES(15) is fCass in component L_type_Ca_current_fCass_gate (dimensionless).
C ALGBRC(7) is d_inf in component L_type_Ca_current_d_gate (dimensionless).
C ALGBRC(20) is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
C ALGBRC(33) is beta_d in component L_type_Ca_current_d_gate (dimensionless).
C ALGBRC(41) is gamma_d in component L_type_Ca_current_d_gate (millisecond).
C ALGBRC(43) is tau_d in component L_type_Ca_current_d_gate (millisecond).
C ALGBRC(8) is f_inf in component L_type_Ca_current_f_gate (dimensionless).
C ALGBRC(21) is tau_f in component L_type_Ca_current_f_gate (millisecond).
C ALGBRC(9) is f2_inf in component L_type_Ca_current_f2_gate (dimensionless).
C ALGBRC(22) is tau_f2 in component L_type_Ca_current_f2_gate (millisecond).
C ALGBRC(10) is fCass_inf in component L_type_Ca_current_fCass_gate (dimensionless).
C ALGBRC(23) is tau_fCass in component L_type_Ca_current_fCass_gate (millisecond).
C CONSTS(20) is g_bca in component calcium_background_current (nanoS_per_picoF).
C CONSTS(21) is g_to in component transient_outward_current (nanoS_per_picoF).
C STATES(16) is s in component transient_outward_current_s_gate (dimensionless).
C STATES(17) is r in component transient_outward_current_r_gate (dimensionless).
C ALGBRC(11) is s_inf in component transient_outward_current_s_gate (dimensionless).
C ALGBRC(24) is tau_s in component transient_outward_current_s_gate (millisecond).
C ALGBRC(12) is r_inf in component transient_outward_current_r_gate (dimensionless).
C ALGBRC(25) is tau_r in component transient_outward_current_r_gate (millisecond).
C CONSTS(22) is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
C CONSTS(23) is K_mk in component sodium_potassium_pump_current (millimolar).
C CONSTS(24) is K_mNa in component sodium_potassium_pump_current (millimolar).
C CONSTS(25) is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
C CONSTS(26) is K_sat in component sodium_calcium_exchanger_current (dimensionless).
C CONSTS(27) is alpha in component sodium_calcium_exchanger_current (dimensionless).
C CONSTS(28) is gamma in component sodium_calcium_exchanger_current (dimensionless).
C CONSTS(29) is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
C CONSTS(30) is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
C CONSTS(31) is g_pCa in component calcium_pump_current (picoA_per_picoF).
C CONSTS(32) is K_pCa in component calcium_pump_current (millimolar).
C CONSTS(33) is g_pK in component potassium_pump_current (nanoS_per_picoF).
C STATES(18) is Ca_SR in component calcium_dynamics (millimolar).
C ALGBRC(68) is i_rel in component calcium_dynamics (millimolar_per_millisecond).
C ALGBRC(60) is i_up in component calcium_dynamics (millimolar_per_millisecond).
C ALGBRC(61) is i_leak in component calcium_dynamics (millimolar_per_millisecond).
C ALGBRC(62) is i_xfer in component calcium_dynamics (millimolar_per_millisecond).
C ALGBRC(67) is O in component calcium_dynamics (dimensionless).
C STATES(19) is R_prime in component calcium_dynamics (dimensionless).
C ALGBRC(65) is k1 in component calcium_dynamics (per_millimolar2_per_millisecond).
C ALGBRC(66) is k2 in component calcium_dynamics (per_millimolar_per_millisecond).
C CONSTS(34) is k1_prime in component calcium_dynamics (per_millimolar2_per_millisecond).
C CONSTS(35) is k2_prime in component calcium_dynamics (per_millimolar_per_millisecond).
C CONSTS(36) is k3 in component calcium_dynamics (per_millisecond).
C CONSTS(37) is k4 in component calcium_dynamics (per_millisecond).
C CONSTS(38) is EC in component calcium_dynamics (millimolar).
C CONSTS(39) is max_sr in component calcium_dynamics (dimensionless).
C CONSTS(40) is min_sr in component calcium_dynamics (dimensionless).
C ALGBRC(63) is kcasr in component calcium_dynamics (dimensionless).
C CONSTS(41) is V_rel in component calcium_dynamics (per_millisecond).
C CONSTS(42) is V_xfer in component calcium_dynamics (per_millisecond).
C CONSTS(43) is K_up in component calcium_dynamics (millimolar).
C CONSTS(44) is V_leak in component calcium_dynamics (per_millisecond).
C CONSTS(45) is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
C ALGBRC(64) is Ca_i_bufc in component calcium_dynamics (dimensionless).
C ALGBRC(69) is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
C ALGBRC(70) is Ca_ss_bufss in component calcium_dynamics (dimensionless).
C CONSTS(46) is Buf_c in component calcium_dynamics (millimolar).
C CONSTS(47) is K_buf_c in component calcium_dynamics (millimolar).
C CONSTS(48) is Buf_sr in component calcium_dynamics (millimolar).
C CONSTS(49) is K_buf_sr in component calcium_dynamics (millimolar).
C CONSTS(50) is Buf_ss in component calcium_dynamics (millimolar).
C CONSTS(51) is K_buf_ss in component calcium_dynamics (millimolar).
C CONSTS(52) is V_sr in component calcium_dynamics (micrometre3).
C RATES(1) is d/dt V in component membrane (millivolt).
C RATES(5) is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
C RATES(6) is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
C RATES(7) is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
C RATES(8) is d/dt m in component fast_sodium_current_m_gate (dimensionless).
C RATES(9) is d/dt h in component fast_sodium_current_h_gate (dimensionless).
C RATES(10) is d/dt j in component fast_sodium_current_j_gate (dimensionless).
C RATES(12) is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
C RATES(13) is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
C RATES(14) is d/dt f2 in component L_type_Ca_current_f2_gate (dimensionless).
C RATES(15) is d/dt fCass in component L_type_Ca_current_fCass_gate (dimensionless).
C RATES(16) is d/dt s in component transient_outward_current_s_gate (dimensionless).
C RATES(17) is d/dt r in component transient_outward_current_r_gate (dimensionless).
C RATES(19) is d/dt R_prime in component calcium_dynamics (dimensionless).
C RATES(4) is d/dt Ca_i in component calcium_dynamics (millimolar).
C RATES(18) is d/dt Ca_SR in component calcium_dynamics (millimolar).
C RATES(11) is d/dt Ca_ss in component calcium_dynamics (millimolar).
C RATES(3) is d/dt Na_i in component sodium_dynamics (millimolar).
C RATES(2) is d/dt K_i in component potassium_dynamics (millimolar).
C
      SUBROUTINE initConsts(CONSTS, STATES)
c the input and output of each subroutine is derived form the matlab code
c inputs: nothing
c outputs: CONSTS, STATES
c >>  I removed the dummy argument RATES   <<
      double precision CONSTS(*), STATES(*)
      STATES(1) = -85.23
      CONSTS(1) = 8314.472
      CONSTS(2) = 310
      CONSTS(3) = 96485.3415
      CONSTS(4) = 0.185
      CONSTS(5) = 0.016404
      CONSTS(6) = 10
      CONSTS(7) = 1000
      CONSTS(8) = 1
      CONSTS(9) = 52
      CONSTS(10) = 0.03
      CONSTS(11) = 5.4
      CONSTS(12) = 140
      STATES(2) = 136.89
      STATES(3) = 8.604
      CONSTS(13) = 2
      STATES(4) = 0.000126
      CONSTS(14) = 5.405
      CONSTS(15) = 0.153
      STATES(5) = 0.00621
      STATES(6) = 0.4712
      CONSTS(16) = 0.392
      STATES(7) = 0.0095
      CONSTS(17) = 14.838
      STATES(8) = 0.00172
      STATES(9) = 0.7444
      STATES(10) = 0.7045
      CONSTS(18) = 0.00029
      CONSTS(19) = 0.0000398
      STATES(11) = 0.00036
      STATES(12) = 3.373D-5
      STATES(13) = 0.7888
      STATES(14) = 0.9755
      STATES(15) = 0.9953
      CONSTS(20) = 0.000592
      CONSTS(21) = 0.294
      STATES(16) = 0.999998
      STATES(17) = 2.42D-8
      CONSTS(22) = 2.724
      CONSTS(23) = 1
      CONSTS(24) = 40
      CONSTS(25) = 1000
      CONSTS(26) = 0.1
      CONSTS(27) = 2.5
      CONSTS(28) = 0.35
      CONSTS(29) = 1.38
      CONSTS(30) = 87.5
      CONSTS(31) = 0.1238
      CONSTS(32) = 0.0005
      CONSTS(33) = 0.0146
      STATES(18) = 3.64
      STATES(19) = 0.9073
      CONSTS(34) = 0.15
      CONSTS(35) = 0.045
      CONSTS(36) = 0.06
      CONSTS(37) = 0.005
      CONSTS(38) = 1.5
      CONSTS(39) = 2.5
      CONSTS(40) = 1
      CONSTS(41) = 0.102
      CONSTS(42) = 0.0038
      CONSTS(43) = 0.00025
      CONSTS(44) = 0.00036
      CONSTS(45) = 0.006375
      CONSTS(46) = 0.2
      CONSTS(47) = 0.001
      CONSTS(48) = 10
      CONSTS(49) = 0.3
      CONSTS(50) = 0.4
      CONSTS(51) = 0.00025
      CONSTS(52) = 0.001094
      CONSTS(53) = 0.00005468
      RETURN
      END
      SUBROUTINE computeRates(VOI, CONSTS,  RATES, STATES, ALGBRC)
c inputs: VOI, STATES, CONSTS
c outputs: RATES, ALGBRC
      double precision VOI, CONSTS(*), RATES(*), STATES(*), ALGBRC(*)
      ALGBRC(8) = 1.00000/(1.00000+EXP((STATES(1)+20.0000)/7.00000))
      ALGBRC(21) = 1102.50*EXP(-(STATES(1)+27.0000) ** 2.00000/225.000)
     & +200.000/(1.00000+EXP((13.0000 - STATES(1))/10.0000))+180.000/
     & (1.00000+EXP((STATES(1)+30.0000)/10.0000))+20.0000
      RATES(13) = (ALGBRC(8) - STATES(13))/ALGBRC(21)
      ALGBRC(9) = 0.670000/(1.00000+EXP((STATES(1)+35.0000)/7.00000))
     & +0.330000
      ALGBRC(22) =  562.000*EXP(- (STATES(1)+27.0000) ** 2.00000/
     & 240.000)+31.0000/(1.00000+EXP((25.0000 - STATES(1))/
     & 10.0000))+80.0000/(1.00000+EXP((STATES(1)+30.0000)/10.0000))
      RATES(14) = (ALGBRC(9) - STATES(14))/ALGBRC(22)
      ALGBRC(10) = 0.600000/(1.00000+(STATES(11)/0.0500000) ** 2.00000)
     & + 0.400000
      ALGBRC(23) = 80.0000/(1.00000+(STATES(11)/0.0500000) ** 2.00000)
     & + 2.00000
      RATES(15) = (ALGBRC(10) - STATES(15))/ALGBRC(23)
      ALGBRC(11) = 1.00000/(1.00000+EXP((STATES(1)+20.0000)/5.00000))
      ALGBRC(24) = 85.0000*EXP(- (STATES(1)+45.0000) ** 2.00000/
     & 320.000)+5.00000/(1.00000+EXP((STATES(1) - 20.0000)/5.00000))
     & +3.00000
      RATES(16) = (ALGBRC(11) - STATES(16))/ALGBRC(24)
      ALGBRC(12)=1.00000/(1.00000+EXP((20.0000 - STATES(1))/6.00000))
      ALGBRC(25) = 9.50000*EXP(- (STATES(1)+40.0000) ** 2.00000/
     & 1800.00)+0.800000
      RATES(17) = (ALGBRC(12) - STATES(17))/ALGBRC(25)
      ALGBRC(1)=1.00000/(1.00000+EXP((-26.0000 - STATES(1))/7.00000))
      ALGBRC(14) = 450.000/(1.00000+EXP((- 45.0000 - STATES(1))
     & /10.0000))
      ALGBRC(27) = 6.00000/(1.00000+EXP((STATES(1)+30.0000)/11.5000))
      ALGBRC(35) =  1.00000*ALGBRC(14)*ALGBRC(27)
      RATES(5) = (ALGBRC(1) - STATES(5))/ALGBRC(35)
      ALGBRC(2) = 1.00000/(1.00000+EXP((STATES(1)+88.0000)/24.0000))
      ALGBRC(15) = 3.00000/(1.00000+EXP((- 60.0000 - STATES(1))/
     & 20.0000))
      ALGBRC(28)=1.12000/(1.00000+EXP((STATES(1) - 60.0000)/20.0000))
      ALGBRC(36) =  1.00000*ALGBRC(15)*ALGBRC(28)
      RATES(6) = (ALGBRC(2) - STATES(6))/ALGBRC(36)
      ALGBRC(3) = 1.00000/(1.00000+EXP((- 5.00000 - STATES(1))/
     & 14.0000))
      ALGBRC(16) = 1400.00/ (1.00000+EXP((5.00000 - STATES(1))/
     & 6.00000)) ** (1.0 / 2)
      ALGBRC(29)=1.00000/(1.00000+EXP((STATES(1) - 35.0000)/15.0000))
      ALGBRC(37) =  1.00000*ALGBRC(16)*ALGBRC(29)+80.0000
      RATES(7) = (ALGBRC(3) - STATES(7))/ALGBRC(37)
      ALGBRC(4) = 1.00000/(1.00000+EXP((-56.8600 - STATES(1))/
     & 9.03000)) ** 2.00000
      ALGBRC(17) = 1.00000/(1.00000+EXP((- 60.0000 - STATES(1))/
     & 5.00000))
      ALGBRC(30) = 0.100000/(1.00000+EXP((STATES(1)+35.0000)/
     & 5.00000)) + 0.100000/(1.00000+EXP((STATES(1) - 50.0000)/
     & 200.000))
      ALGBRC(38) =  1.00000*ALGBRC(17)*ALGBRC(30)
      RATES(8) = (ALGBRC(4) - STATES(8))/ALGBRC(38)
      ALGBRC(5) = 1.00000/(1.00000+EXP(STATES(1)+71.5500)/7.43000)
     & ** 2.00000
c      ALGBRC(18) = TERNRY(STATES(1).LT. -40.0000,  0.0570000*
c     & EXP(-(STATES(1)+80.0000)/6.80000), 0.00000D0)
      IF (STATES(1) .LT. -40.0000d0) THEN
        algbrc(18) = 0.0570000*EXP(-(STATES(1)+80.0000)/6.80000)
      ELSE
        algbrc(18) = 0.00000D0
      END IF
c      ALGBRC(31) = TERNRY(STATES(1).LT.- 40.0000, 2.70000*EXP(
c     & 0.0790000*STATES(1))+ 310000.*EXP(0.348500*STATES(1)),
c     & 0.770000/(0.130000*(1.00000+EXP((STATES(1)+10.6600)/ 
c     & (-11.1000)))))
      IF (STATES(1) .LT. -40.0000d0) THEN
        algbrc(31) = 2.70000*EXP(0.0790000*STATES(1))+ 310000.*
     & EXP(0.348500*STATES(1))
      ELSE
        algbrc(31) = 0.770000/(0.130000*(1.00000+EXP((STATES(1)+ 
     & 10.6600)/(-11.1000))))
      END IF
      ALGBRC(39) = 1.00000/(ALGBRC(18)+ALGBRC(31))
      RATES(9) = (ALGBRC(5) - STATES(9))/ALGBRC(39)
      ALGBRC(6) = 1.00000/(1.00000 +
     & EXP(STATES(1) + 71.5500)/7.43000) ** 2.00000
c      ALGBRC(19) = TERNRY(STATES(1).LT.-40.0000, (((-25428.0*EXP(
c     & 0.244400*STATES(1)) - 6.94800D-06*EXP(-0.0439100*STATES(1)))*
c     & (STATES(1)+37.7800))/1.00000)/(1.00000+EXP( 0.311000*
c     & (STATES(1)+79.2300))), 0.00000D0)
      IF (STATES(1) .LT. -40.0000d0) THEN
       algbrc(19)=(((-25428.0*EXP(0.244400*STATES(1)) - 6.94800D-06*
     & EXP(-0.0439100*STATES(1)))*(STATES(1)+37.7800))/1.00000)/
     & (1.00000+EXP( 0.311000*(STATES(1)+ 79.2300)))
      ELSE
        algbrc(19) =  0.00000D0
      END IF
c      ALGBRC(32) = TERNRY(STATES(1).LT.-40.0000, (0.0242400*EXP(
c     & -0.0105200*STATES(1)))/(1.00000+EXP(-0.137800*(STATES(1)+
c     & 40.1400))), (0.600000*EXP(0.0570000*STATES(1)))/(1.00000+
c     & EXP(-0.100000*(STATES(1) + 32.0000))))
      IF (STATES(1) .LT. -40.0000d0) THEN
       algbrc(32)=(0.0242400*EXP(-0.0105200*STATES(1)))/(1.00000+
     & EXP(-0.137800*(STATES(1)+40.1400)))
      ELSE
        algbrc(32) = (0.600000*EXP(0.0570000*STATES(1)))/(1.00000+
     & EXP(-0.100000*(STATES(1) + 32.0000)))
      END IF
      ALGBRC(40) = 1.00000/(ALGBRC(19)+ALGBRC(32))
      RATES(10) = (ALGBRC(6) - STATES(10))/ALGBRC(40)
      ALGBRC(7) = 1.00000/(1.00000+EXP((- 8.00000 - STATES(1))/
     & 7.50000))
      ALGBRC(20) = 1.40000/(1.00000+EXP((- 35.0000 - STATES(1))/
     & 13.0000)) + 0.250000
      ALGBRC(33) = 1.40000/(1.00000+EXP((STATES(1)+5.00000)/5.00000))
      ALGBRC(41) = 1.00000/(1.00000+EXP((50.0000 - STATES(1))/
     & 20.0000))
      ALGBRC(43) = 1.00000*ALGBRC(20)*ALGBRC(33)+ALGBRC(41)
      RATES(12) = (ALGBRC(7) - STATES(12))/ALGBRC(43)
      ALGBRC(56) =((((CONSTS(22)*CONSTS(11))/(CONSTS(11)+CONSTS(23)))
     & *STATES(3))/(STATES(3)+CONSTS(24)))/(1.00000+ 0.124500*EXP((
     & -0.100000*STATES(1)*CONSTS(3))/( CONSTS(1)*CONSTS(2)))+ 
     & 0.0353000*EXP((-STATES(1)*CONSTS(3))/(CONSTS(1)*CONSTS(2))))
      ALGBRC(26) =  ((CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(12)/STATES(3))
      ALGBRC(51) =  CONSTS(17)*STATES(8) ** 3.00000*STATES(9)*
     & STATES(10)*(STATES(1) - ALGBRC(26))
      ALGBRC(52) =  CONSTS(18)*(STATES(1) - ALGBRC(26))
      ALGBRC(57) = (CONSTS(25)*(EXP((CONSTS(28)*STATES(1)*CONSTS(3))
     & /( CONSTS(1)*CONSTS(2)))*STATES(3) ** 3.00000*CONSTS(13) -  
     & EXP(((CONSTS(28) - 1.00000)*STATES(1)*CONSTS(3))/(CONSTS(1)
     & *CONSTS(2)))*CONSTS(12) ** 3.00000*STATES(4)*CONSTS(27)))/((
     & CONSTS(30) ** 3.00000+CONSTS(12) ** 3.00000)*(CONSTS(29)+
     & CONSTS(13))*(1.00000+ CONSTS(26)*EXP(((CONSTS(28) - 1.00000)
     & *STATES(1)*CONSTS(3))/( CONSTS(1)*CONSTS(2)))))
      RATES(3) =  ((-1.00000*(ALGBRC(51)+ALGBRC(52)+ 3.00000*
     & ALGBRC(56)+ 3.00000*ALGBRC(57)))/( 1.00000*CONSTS(5)*
     & CONSTS(3)))*CONSTS(4)
      ALGBRC(34) =  (( CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(11)/STATES(2))
      ALGBRC(45) = 0.100000/(1.00000+EXP(0.0600000*((STATES(1)
     & - ALGBRC(34)) - 200.000)))
      ALGBRC(46) = (3.00000*EXP(0.000200000*((STATES(1) - ALGBRC(34))
     & + 100.000))+EXP(0.100000*((STATES(1) - ALGBRC(34)) - 
     & 10.0000)))/(1.00000+EXP( - 0.500000*(STATES(1) - ALGBRC(34))))
      ALGBRC(47) = ALGBRC(45)/(ALGBRC(45)+ALGBRC(46))
      ALGBRC(48) =  CONSTS(14)*ALGBRC(47)* (CONSTS(11)/5.40000) ** 
     & (1.0 / 2)*(STATES(1) - ALGBRC(34))
      ALGBRC(55) =  CONSTS(21)*STATES(17)*STATES(16)*
     & (STATES(1) - ALGBRC(34))
      ALGBRC(49) =  CONSTS(15)* (CONSTS(11)/5.40000) ** (1.0 / 2)*
     & STATES(5)*STATES(6)*(STATES(1) - ALGBRC(34))
      ALGBRC(42) =  (( CONSTS(1)*CONSTS(2))/CONSTS(3))*log((
     & CONSTS(11)+ CONSTS(10)*CONSTS(12))/(STATES(2)+ 
     & CONSTS(10)*STATES(3)))
      ALGBRC(50) =  CONSTS(16)*STATES(7) ** 2.00000*(STATES(1) - 
     & ALGBRC(42))
      ALGBRC(53) = (((CONSTS(19)*STATES(12)*STATES(13)*STATES(14)*
     & STATES(15)*4.00000*(STATES(1) - 15.0000)*(CONSTS(3) ** 2.00000))
     & /(CONSTS(1)*CONSTS(2)))*( 0.250000*STATES(11)*EXP((2.00000*(
     & STATES(1) - 15.0000)*CONSTS(3))/(CONSTS(1)*CONSTS(2))) - 
     & CONSTS(13)))/(EXP(( 2.00000*(STATES(1) - 15.0000)*CONSTS(3))
     & /(CONSTS(1)*CONSTS(2))) - 1.00000)
      ALGBRC(44) =  ((0.500000*CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(13)/STATES(4))
      ALGBRC(54) =  CONSTS(20)*(STATES(1) - ALGBRC(44))
      ALGBRC(59) = ( CONSTS(33)*(STATES(1) - ALGBRC(34)))/
     & (1.00000+EXP((25.0000 - STATES(1))/5.98000))
      ALGBRC(58) = ( CONSTS(31)*STATES(4))/(STATES(4)+CONSTS(32))
c      ALGBRC(13) = TERNRY(VOI - INT(VOI/CONSTS(7))*CONSTS(7).GE.
c     & CONSTS(6).AND.VOI - INT(VOI/CONSTS(7))*CONSTS(7).LE.CONSTS(6)
c     & +CONSTS(8), - CONSTS(9), 0.00000D0)
      IF (VOI - INT(VOI/CONSTS(7))*CONSTS(7).GE.CONSTS(6).AND.VOI - 
     & INT(VOI/CONSTS(7))*CONSTS(7).LE.CONSTS(6)+CONSTS(8)) THEN
       algbrc(13)= - CONSTS(9)
      ELSE
        algbrc(13) = 0.00000D0
      END IF
      RATES(1) = (-1.00000/1.00000)*(ALGBRC(48)+ALGBRC(55)+ALGBRC(49)
     & +ALGBRC(50)+ALGBRC(53)+ALGBRC(56)+ALGBRC(51)+ALGBRC(52)+
     & ALGBRC(57)+ALGBRC(54)+ALGBRC(59)+ALGBRC(58)+ALGBRC(13))
      RATES(2) = ((-1.00000*((ALGBRC(48)+ALGBRC(55)+ALGBRC(49)+
     & ALGBRC(50)+ALGBRC(59)+ALGBRC(13)) - 2.00000*ALGBRC(56)))/
     & (1.00000*CONSTS(5)*CONSTS(3)))*CONSTS(4)
      ALGBRC(60) = CONSTS(45)/(1.00000+CONSTS(43) ** 2.00000/
     & STATES(4) ** 2.00000)
      ALGBRC(61) =  CONSTS(44)*(STATES(18) - STATES(4))
      ALGBRC(62) =  CONSTS(42)*(STATES(11) - STATES(4))
      ALGBRC(64) = 1.00000/(1.00000+(CONSTS(46)*CONSTS(47))/(STATES(4)
     & +CONSTS(47)) ** 2.00000)
      RATES(4) = ALGBRC(64)*((((ALGBRC(61) - ALGBRC(60))*CONSTS(52))/
     & CONSTS(5)+ALGBRC(62)) - (1.00000*((ALGBRC(54)+ALGBRC(58)) - 
     & 2.00000*ALGBRC(57))*CONSTS(4))/(2.00000*1.00000*CONSTS(5)*
     & CONSTS(3)))
      ALGBRC(63) = CONSTS(39) - (CONSTS(39) - CONSTS(40))/
     & (1.00000+(CONSTS(38)/STATES(18)) ** 2.00000)
      ALGBRC(66) =  CONSTS(35)*ALGBRC(63)
      RATES(19) = - ALGBRC(66)*STATES(11)*STATES(19)+ CONSTS(37)*
     & (1.00000 - STATES(19))
      ALGBRC(65) = CONSTS(34)/ALGBRC(63)
      ALGBRC(67) = ( ALGBRC(65)*STATES(11) ** 2.00000*STATES(19))/
     & (CONSTS(36)+ ALGBRC(65)*STATES(11) ** 2.00000)
      ALGBRC(68) =  CONSTS(41)*ALGBRC(67)*(STATES(18) - STATES(11))
      ALGBRC(69) = 1.00000/(1.00000+( CONSTS(48)*CONSTS(49))/
     & (STATES(18)+CONSTS(49)) ** 2.00000)
      RATES(18) =  ALGBRC(69)*(ALGBRC(60) - (ALGBRC(68)+ALGBRC(61)))
      ALGBRC(70) = 1.00000/(1.00000+( CONSTS(50)*CONSTS(51))/
     & (STATES(11)+CONSTS(51)) ** 2.00000)
      RATES(11) =  ALGBRC(70)*(((-1.00000*ALGBRC(53)*CONSTS(4))/
     & (2.00000*1.00000*CONSTS(53)*CONSTS(3))+(ALGBRC(68)*CONSTS(52))
     & /CONSTS(53)) - (ALGBRC(62)*CONSTS(5))/CONSTS(53))
      RETURN
      END
      SUBROUTINE computeVariables(VOI, CONSTS, STATES, ALGBRC)
c it is the same as function compute Algebraic in matlab
c inputs: ALGBRC, CONSTS, STATES, VOI
c outputs: ALGBRC
c >>  I removed the dummy argument RATES   <<
      double precision VOI, CONSTS(*), STATES(*), ALGBRC(*)
      ALGBRC(8) = 1.00000/(1.00000+EXP((STATES(1)+20.0000)/7.00000))
      ALGBRC(21) =  1102.50*EXP(- (STATES(1)+27.0000) ** 2.00000/
     & 225.000)+200.000/(1.00000+EXP((13.0000 - STATES(1))/10.0000))
     & + 180.000/(1.00000+EXP((STATES(1)+30.0000)/10.0000))+ 20.0000
      ALGBRC(9) = 0.670000/(1.00000+EXP((STATES(1)+35.0000)/7.00000))
     & + 0.330000
      ALGBRC(22) = 562.000*EXP(-(STATES(1)+27.0000) ** 2.00000/240.000)
     & +31.0000/(1.00000+EXP((25.0000 - STATES(1))/10.0000))+80.0000/
     & (1.00000+EXP((STATES(1)+30.0000)/10.0000))
      ALGBRC(10) = 0.600000/(1.00000+(STATES(11)/0.0500000) ** 2.00000)
     & + 0.400000
      ALGBRC(23) = 80.0000/(1.00000+(STATES(11)/0.0500000) ** 2.00000)
     & + 2.00000
      ALGBRC(11) = 1.00000/(1.00000+EXP((STATES(1)+20.0000)/5.00000))
      ALGBRC(24) = 85.0000*EXP(- (STATES(1)+45.0000) ** 2.00000/
     & 320.000)+5.00000/(1.00000+EXP((STATES(1) - 20.0000)/5.00000))
     & +3.00000
      ALGBRC(12)=1.00000/(1.00000+EXP((20.0000 - STATES(1))/6.00000))
      ALGBRC(25) = 9.50000*EXP(- (STATES(1)+40.0000) ** 2.00000/
     & 1800.00)+0.800000
      ALGBRC(1)=1.00000/(1.00000+EXP((-26.0000 - STATES(1))/7.00000))
      ALGBRC(14) = 450.000/(1.00000+EXP((-45.0000 - STATES(1))/
     & 10.0000))
      ALGBRC(27) = 6.00000/(1.00000+EXP((STATES(1)+30.0000)/11.5000))
      ALGBRC(35) =  1.00000*ALGBRC(14)*ALGBRC(27)
      ALGBRC(2) = 1.00000/(1.00000+EXP((STATES(1)+88.0000)/24.0000))
      ALGBRC(15) = 3.00000/(1.00000+EXP((-60.0000 - STATES(1))/
     & 20.0000))
      ALGBRC(28) = 1.12000/(1.00000+EXP((STATES(1) - 60.0000)/
     & 20.0000))
      ALGBRC(36) =  1.00000*ALGBRC(15)*ALGBRC(28)
      ALGBRC(3)=1.00000/(1.00000+EXP((-5.00000 - STATES(1))/14.0000))
      ALGBRC(16) = 1400.00/ (1.00000+EXP((5.00000 - STATES(1))/
     & 6.00000)) ** (1.0 / 2)
      ALGBRC(29)=1.00000/(1.00000+EXP((STATES(1) - 35.0000)/15.0000))
      ALGBRC(37) = 1.00000*ALGBRC(16)*ALGBRC(29)+80.0000
      ALGBRC(4) = 1.00000/(1.00000+EXP((- 56.8600 - STATES(1))/
     & 9.03000)) ** 2.00000
      ALGBRC(17) = 1.00000/(1.00000+EXP((-60.0000 - STATES(1))
     & /5.00000))
      ALGBRC(30) = 0.100000/(1.00000+EXP((STATES(1)+35.0000)/5.00000)
     & )+0.100000/(1.00000+EXP((STATES(1) - 50.0000)/200.000))
      ALGBRC(38) =  1.00000*ALGBRC(17)*ALGBRC(30)
      ALGBRC(5) = 1.00000/(1.00000+EXP(STATES(1)+71.5500)/7.43000) 
     & ** 2.00000
c      ALGBRC(18) = TERNRY(STATES(1).LT.- 40.0000, 0.0570000*EXP(-(
c     & STATES(1)+80.0000)/6.80000), 0.00000D0)
      IF (STATES(1).LT.- 40.0000) THEN
       algbrc(18)= 0.0570000*EXP(-(STATES(1)+80.0000)/6.80000)
      ELSE
        algbrc(18) = 0.00000D0
      END IF
c      ALGBRC(31) = TERNRY(STATES(1).LT.-40.0000,2.70000*EXP(0.0790000
c     & *STATES(1))+ 310000.*EXP( 0.348500*STATES(1)), 0.770000/
c     & (0.130000*(1.00000+EXP((STATES(1)+10.6600)/(-11.1000)))))
      IF (STATES(1).LT.- 40.0000) THEN
       algbrc(31)=2.70000*EXP(0.0790000*STATES(1))+ 310000.*
     & EXP( 0.348500*STATES(1))
      ELSE
        algbrc(31)=0.770000/(0.130000*(1.00000+
     & EXP((STATES(1)+10.6600)/(-11.1000))))
      END IF
      ALGBRC(39) = 1.00000/(ALGBRC(18)+ALGBRC(31))
      ALGBRC(6) = 1.00000/(1.00000 +
     & EXP(STATES(1) + 71.5500)/7.43000) ** 2.00000
c      ALGBRC(19) = TERNRY(STATES(1).LT.-40.0000,(((-25428.0*EXP(
c     & 0.244400*STATES(1)) - 6.94800D-06*EXP(-0.0439100*STATES(1)))*
c     & (STATES(1)+37.7800))/1.00000)/(1.00000+EXP( 0.311000*(
c     & STATES(1)+79.2300))), 0.00000D0)
      IF (STATES(1).LT.- 40.0000) THEN
       algbrc(19)=(((-25428.0*EXP(
     & 0.244400*STATES(1)) - 6.94800D-06*EXP(-0.0439100*STATES(1)))*
     & (STATES(1)+37.7800))/1.00000)/(1.00000+EXP( 0.311000*(STATES(1)+ 
     & 79.2300)))
      ELSE
        algbrc(19)=0.00000D0
      END IF
c      ALGBRC(32) = TERNRY(STATES(1).LT.-40.0000,(0.0242400*EXP(
c     & -0.0105200*STATES(1)))/(1.00000+EXP(-0.137800*(STATES(1) + 
c     & 40.1400))),(0.600000*EXP(0.0570000*STATES(1)))/(1.00000+EXP
c     & (-0.100000*(STATES(1)+32.0000))))
      IF (STATES(1).LT.- 40.0000) THEN
       algbrc(32)=(0.0242400*EXP(-0.0105200*STATES(1)))/(1.00000+
     & EXP(-0.137800*(STATES(1) + 40.1400)))
      ELSE
        algbrc(32)=(0.600000*EXP(0.0570000*STATES(1)))/(1.00000+EXP
     & (-0.100000*(STATES(1) + 32.0000)))
      END IF
      ALGBRC(40) = 1.00000/(ALGBRC(19)+ALGBRC(32))
      ALGBRC(7) = 1.00000/(1.00000+EXP((- 8.00000 - STATES(1))/
     & 7.50000))
      ALGBRC(20) = 1.40000/(1.00000+EXP((-35.0000 - STATES(1))/
     & 13.0000)) + 0.250000
      ALGBRC(33) = 1.40000/(1.00000+EXP((STATES(1)+5.00000)/5.00000))
      ALGBRC(41)=1.00000/(1.00000+EXP((50.0000 - STATES(1))/20.0000))
      ALGBRC(43) = 1.00000*ALGBRC(20)*ALGBRC(33)+ALGBRC(41)
      ALGBRC(56) = (( (( CONSTS(22)*CONSTS(11))/(CONSTS(11)+CONSTS(
     & 23)))*STATES(3))/(STATES(3)+CONSTS(24)))/(1.00000+ 0.124500*
     & EXP((-0.100000*STATES(1)*CONSTS(3))/(CONSTS(1)*CONSTS(2)))+ 
     & 0.0353000*EXP((-STATES(1)*CONSTS(3))/(CONSTS(1)*CONSTS(2))))
      ALGBRC(26) =  (( CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(12)/STATES(3))
      ALGBRC(51) =  CONSTS(17)*STATES(8) ** 3.00000*STATES(9)*
     & STATES(10)*(STATES(1) - ALGBRC(26))
      ALGBRC(52) =  CONSTS(18)*(STATES(1) - ALGBRC(26))
      ALGBRC(57) = (CONSTS(25)*(EXP((CONSTS(28)*STATES(1)*CONSTS(3))/
     & (CONSTS(1)*CONSTS(2)))*STATES(3) ** 3.00000*CONSTS(13) - 
     & EXP(((CONSTS(28) - 1.00000)*STATES(1)*CONSTS(3))/(CONSTS(1)
     & *CONSTS(2)))*CONSTS(12) ** 3.00000*STATES(4)*CONSTS(27)))/
     & ((CONSTS(30) ** 3.00000+CONSTS(12) ** 3.00000)*(CONSTS(29)+
     & CONSTS(13))*(1.00000+ CONSTS(26)*EXP(((CONSTS(28) - 1.00000)*
     & STATES(1)*CONSTS(3))/( CONSTS(1)*CONSTS(2)))))
      ALGBRC(34) =  (( CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(11)/STATES(2))
      ALGBRC(45) = 0.100000/(1.00000+EXP(0.0600000*((STATES(1) -
     & ALGBRC(34)) - 200.000)))
      ALGBRC(46) = (3.00000*EXP(0.000200000*((STATES(1) - ALGBRC(34))
     & + 100.000))+EXP(0.100000*((STATES(1) - ALGBRC(34)) - 10.0000))
     & )/(1.00000+EXP( - 0.500000*(STATES(1) - ALGBRC(34))))
      ALGBRC(47) = ALGBRC(45)/(ALGBRC(45)+ALGBRC(46))
      ALGBRC(48) = CONSTS(14)*ALGBRC(47)* (CONSTS(11)/5.40000) ** 
     & (1.0 / 2)*(STATES(1) - ALGBRC(34))
      ALGBRC(55) =  CONSTS(21)*STATES(17)*STATES(16)*(STATES(1) - 
     & ALGBRC(34))
      ALGBRC(49) = CONSTS(15)* (CONSTS(11)/5.40000) ** (1.0 / 2)*
     & STATES(5)*STATES(6)*(STATES(1) - ALGBRC(34))
      ALGBRC(42) =  ((CONSTS(1)*CONSTS(2))/CONSTS(3))*log((CONSTS(11)
     & + CONSTS(10)*CONSTS(12))/(STATES(2)+ CONSTS(10)*STATES(3)))
      ALGBRC(50) =  CONSTS(16)*STATES(7) ** 2.00000*(STATES(1) - 
     & ALGBRC(42))
      ALGBRC(53) = ( (( CONSTS(19)*STATES(12)*STATES(13)*STATES(14)*
     & STATES(15)*4.00000*(STATES(1) - 15.0000)*(CONSTS(3) ** 2.00000))
     & /(CONSTS(1)*CONSTS(2)))*( 0.250000*STATES(11)*EXP((2.00000*
     & (STATES(1) - 15.0000)*CONSTS(3))/(CONSTS(1)*CONSTS(2))) - 
     & CONSTS(13)))/(EXP((2.00000*(STATES(1) - 15.0000)*CONSTS(3))/
     & (CONSTS(1)*CONSTS(2))) - 1.00000)
      ALGBRC(44) =  (( 0.500000*CONSTS(1)*CONSTS(2))/CONSTS(3))*
     & log(CONSTS(13)/STATES(4))
      ALGBRC(54) =  CONSTS(20)*(STATES(1) - ALGBRC(44))
      ALGBRC(59) = ( CONSTS(33)*(STATES(1) - ALGBRC(34)))/
     & (1.00000+EXP((25.0000 - STATES(1))/5.98000))
      ALGBRC(58) = ( CONSTS(31)*STATES(4))/(STATES(4)+CONSTS(32))
c      ALGBRC(13) = TERNRY(VOI -  INT(VOI/CONSTS(7))*CONSTS(7).GE.
c     & CONSTS(6).AND.VOI - INT(VOI/CONSTS(7))*CONSTS(7).LE.
c     & CONSTS(6)+CONSTS(8), - CONSTS(9), 0.00000D0)
      IF (VOI -  INT(VOI/CONSTS(7))*CONSTS(7).GE.CONSTS(6).AND.VOI - 
     & INT(VOI/CONSTS(7))*CONSTS(7).LE.CONSTS(6)+CONSTS(8)) THEN
       algbrc(13) = - CONSTS(9)
      ELSE
        algbrc(13) = 0.00000D0
      END IF
      ALGBRC(60) = CONSTS(45)/(1.00000+CONSTS(43) ** 2.00000/
     & STATES(4) ** 2.00000)
      ALGBRC(61) =  CONSTS(44)*(STATES(18) - STATES(4))
      ALGBRC(62) =  CONSTS(42)*(STATES(11) - STATES(4))
      ALGBRC(64) = 1.00000/(1.00000+(CONSTS(46)*CONSTS(47))/(STATES(4)
     & +CONSTS(47)) ** 2.00000)
      ALGBRC(63) = CONSTS(39) - (CONSTS(39) - CONSTS(40))/(1.00000+
     & (CONSTS(38)/STATES(18)) ** 2.00000)
      ALGBRC(66) =  CONSTS(35)*ALGBRC(63)
      ALGBRC(65) = CONSTS(34)/ALGBRC(63)
      ALGBRC(67) = ( ALGBRC(65)*STATES(11) ** 2.00000*STATES(19))/
     & (CONSTS(36)+ ALGBRC(65)*STATES(11) ** 2.00000)
      ALGBRC(68) = CONSTS(41)*ALGBRC(67)*(STATES(18) - STATES(11))
      ALGBRC(69) = 1.00000/(1.00000+( CONSTS(48)*CONSTS(49))/
     & (STATES(18)+CONSTS(49)) ** 2.00000)
      ALGBRC(70) = 1.00000/(1.00000+( CONSTS(50)*CONSTS(51))/
     & (STATES(11)+CONSTS(51)) ** 2.00000)
      RETURN
      END
c      real FUNCTION TERNRY(TEST, VALA, VALB)
c      LOGICAL TEST
c      double precision VALA, VALB
c      IF (TEST) THEN
c        TERNRY = VALA
c      ELSE
c        TERNRY = VALB
c      ENDIF
c      RETURN
c      END
