*------------------------------*
|Simuation section of the paper|
*------------------------------*;

/************************/
/*macros for simulations*/
/************************/;
%macro simulation(
  nsim    =1000,  /*number of simulations runs*/
  ncls    =48,    /*number of clusters*/
  nrep    =5,     /*number of periods*/
  nsub    =50,    /*number of subjects per period per cluster*/
  vc      =,      /*variance components, vs > 0*/
  b       =0.2,   /*slope of period effects*/
  theta   =-0.4,  /*overall treatment effect*/
  theta_j =,      /*treatment period interactions*/
  mu      =1,     /*intercept*/
  design  =,      /*USD design: pgd, dsd, dsdopt, usd, usdopt*/
  c0      =0,     /*integer, number of clusters allocated to control patterns*/
  c1      =0,     /*integer, number of clusters allocated to 2nd patterns*/
  c2      =0,     /*integer, number of clusters allocated to 3rd and inner patterns*/
  seed1   =19890619, 
  seed2   =19600623, 
  seed3   =19881110
);
  
/*  options mprint mlogic nocenter nodate nonumber;*/
  options nodate nonumber nonotes;

  proc iml;
    start exch(dim, vc, cv);
      return(j(dim,dim,cv) + diag(j(dim,1,vc) ) );
    finish;

    start ar1(dim, vc, rho);
      u = cuprod(j(1,dim-1,rho)); /* cumulative product */
      return( vc # toeplitz(1 || u) );
    finish;

    call randseed(&seed1);
    u = repeat(&mu, &nrep, 1);
    a = RandNormal(&ncls*&nsim, &mu, &vc);
    aij= colvec(repeat(colvec(a), 1, &nrep*&nsub));
    sim_idx = colvec(repeat(t(1:&nsim), 1, &ncls * &nrep * &nsub));
    cls_idx = colvec(repeat(t(1:&ncls), &nsim, &nrep * &nsub));
    rep_idx = colvec(repeat(t(1:&nrep), &nsim * &ncls, &nsub));
    sub_idx = repeat(t(1:&nsub), &nsim * &ncls * &nrep, 1);
    combn = sim_idx || cls_idx || rep_idx || sub_idx || aij;
    create individual from combn[c={"sim", "cluster", "period", "id", "intcpt"}];
    append from combn;
    close individual;
  quit;

  ods select none;

  data individual;
    set individual;
    %if &design = pgd %then %do;
      xijk = (cluster <= (0.5 * &ncls));
    %end;
    %if &design = usd %then %do;
      cls_per_switch = &ncls / (&nrep + 1);
      xijk = ((cluster / cls_per_switch) <= period);
    %end;
    %if &design = dsd %then %do;
      cls_per_switch = &ncls / 3;
      select;
        when ((cluster/cls_per_switch) <= 1) xijk = 1;
        when ((cluster/cls_per_switch) > 2) xijk = 0;
        when (1 < (cluster/cls_per_switch) <= 2) xijk = (period >=3);
        otherwise;
      end;
    %end;
    %if &design = dsdopt %then %do;
      select;
        when (cluster <= &c0) xijk = 0;
        when (&c0 < cluster <= &c0+&c1) xijk = (period >=3);
        when (cluster > &c0+&c1) xijk = 1;
        otherwise;
      end; 
    %end;

    %if &design = usdopt %then %do;
      select;
        when (cluster <= &c0) xijk = 0;
        when (&c0  < cluster <= &c0+&c1) xijk = (period >= 5);
        when (&c0+&c1 < cluster <= &c0+&c1+&c2) xijk = (period >= 4);
        when (&c0+&c1+&c2 < cluster <= &c0+&c1+2*&c2) xijk = (period >= 3);
        when (&c0+&c1+2*&c2 < cluster <= &c0+2*&c1+2*&c2) xijk = (period >= 2);
        when (cluster > &c0+2*&c1+2*&c2) xijk = 1;
        otherwise;
      end;
    %end;

    pijk = probnorm(intcpt + (&theta + &theta_j * (period - 3)) * xijk + &b * (period - 3));
    yijk = ranbin(&seed3, 1, pijk);
    trt = xijk;
  run;

  proc means data=individual sum n noprint;
    by sim cluster period;
    var yijk; 
    output out=aggregate sum=event n=total;
  run;

  data aggregate;
    set aggregate;
    drop _FREQ_ _TYPE_;
    %if &design = pgd %then %do;
      xij = (cluster <= (0.5 * &ncls));
    %end;
    %if &design = usd %then %do;
      cls_per_switch = &ncls / (&nrep + 1);
      xij = ((cluster / cls_per_switch) <= period);
    %end;
    %if &design = swd %then %do;
      cls_per_switch = &ncls / (&nrep - 1);
      xij = ((cluster / cls_per_switch) <= (period - 1));
    %end;
    %if &design = dsd %then %do;
      cls_per_switch = &ncls / 3;
      select;
        when ((cluster/cls_per_switch) <= 1) xij = 1;
        when ((cluster/cls_per_switch) > 2) xij = 0;
        when (1 < (cluster/cls_per_switch) <= 2) xij = (period >=3);
        otherwise;
      end;
    %end;
    %if &design = dsdopt %then %do;
      select;
        when (cluster <= &c0) xij = 0;
        when (&c0 < cluster <= &c0+&c1) xij = (period >= 3);
        when (cluster > &c0+&c1) xij = 1;
        otherwise;
      end; 
    %end;

    %if &design = usdopt %then %do;
      select;
        when (cluster <= &c0) xij = 0;
        when (&c0  < cluster <= &c0+&c1) xij = (period >= 5);
        when (&c0+&c1 < cluster <= &c0+&c1+&c2) xij = (period >= 4);
        when (&c0+&c1+&c2 < cluster <= &c0+&c1+2*&c2) xij = (period >= 3);
        when (&c0+&c1+2*&c2 < cluster <= &c0+2*&c1+2*&c2) xij = (period >= 2);
        when (cluster > &c0+2*&c1+2*&c2) xij = 1;
        otherwise;
      end;
    %end;

    trt = xij;
  run;

/*************************/;
/*Marginal model with GEE*/;
/*************************/;
  ods select none;
  ods output GEEEmpPEst = result_mm;
  ods output Type3 = power_mm;
  proc genmod data=aggregate;
    class cluster period(ref='3');
    model event/total = period trt period*trt/link = probit dist = binomial type3;
    repeated subject = cluster/type = exch corrw modelse;
    by sim;
  run;

  data result_mm;
    set result_mm;
    length varname $ 15;
    varname = cats(Parm, Level1);
    if stderr = . then delete;
  run;

/********************/;
/*GLMM Cluster level*/;
/********************/;
  ods select none;
  ods output ParameterEstimates = result_mix;
  ods output tests3 = power_mix;
  proc glimmix data = aggregate method=laplace;
    class cluster period(ref = '3');
    model event/total = period trt trt*period/link = probit dist = binomial solution cl chisq ddfm=none;
    random intercept/subject = cluster;
    by sim;
  run;

  data result_mix;
    set result_mix;
    length varname $ 15;
    varname = compress(cats(Effect, period), "_");
  run;

/***************************/;
/*Individual model with GEE*/;
/***************************/;
  ods select none;
  ods output GEEEmpPEst = result_mm2;
  ods output Type3 = power_mm2;
  proc genmod data = individual;
    class cluster period(ref = '3');
    model yijk(desc) = period trt period*trt/link = probit dist = bin type3;
    repeated subject = cluster/type = EXCH corrw modelse;
    by sim;
  run;

  data result_mm2;
    set result_mm2;
    length varname $ 15;
    varname = cats(Parm, Level1);
    if stderr = . then delete;
  run;

/***********************/;
/*GLMM Individual level*/;
/***********************/;
  ods select none;
  ods output ParameterEstimates = result_mix2;
  ods output tests3 = power_mix2;
  proc glimmix data = individual;
    class cluster period(ref = '3');
    model yijk(desc) = period trt trt*period/link = probit dist = binary solution cl chisq ddfm=none;
    random intercept/subject = cluster;
    by sim;
  run;

  data result_mix2;
    set result_mix2;
    length varname $ 15;
    varname = compress(cats(Effect, period), "_");
  run;

/*****************************/;
/*Aggregating all sim results*/;
/*****************************/;
  proc sql;
    create table result_agg as
    select COALESCE(result_mix.sim, result_mm.sim) as sim,
           COALESCE(result_mix.varname, result_mm.varname) as parms,
           (case when COALESCE(result_mix.varname, result_mm.varname) = "Intercept" then &mu
                 when COALESCE(result_mix.varname, result_mm.varname) = "period1" then &b * (1-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "period2" then &b * (2-3) 
                 when COALESCE(result_mix.varname, result_mm.varname) = "period3" then 0
                 when COALESCE(result_mix.varname, result_mm.varname) = "period4" then &b * (4-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "period5" then &b * (5-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt" then &theta
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt*period1" then &theta_j * (1-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt*period2" then &theta_j * (2-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt*period3" then 0
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt*period4" then &theta_j * (4-3)
                 when COALESCE(result_mix.varname, result_mm.varname) = "trt*period5" then &theta_j * (5-3) end) as true_mixed,
           (calculated true_mixed / sqrt(1 + &vc)) as true_mm,
           result_mix.estimate  as glmm_agg_est,
           result_mm.estimate   as gee_agg_est,
           result_mix.lower as glmm_agg_est_lower,
           result_mix.upper as glmm_agg_est_upper,
           result_mm.lowercl as gee_agg_est_lower,
           result_mm.uppercl as gee_agg_est_upper,
           (glmm_agg_est - calculated true_mixed) as glmm_agg_bias,
           (gee_agg_est - calculated true_mm) as gee_agg_bias,
           (calculated glmm_agg_bias)**2 as glmm_agg_mse,
           (calculated gee_agg_bias)**2 as gee_agg_mse
    from result_mix full join result_mm on (result_mix.sim = result_mm.sim and result_mix.varname = result_mm.varname);
  quit;


    proc sql;
    create table result_ind as
    select COALESCE(result_mix2.sim, result_mm2.sim) as sim,
           COALESCE(result_mix2.varname, result_mm2.varname) as parms,
           (case when COALESCE(result_mix2.varname, result_mm2.varname) = "Intercept" then &mu
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "period1" then &b * (1-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "period2" then &b * (2-3) 
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "period3" then 0
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "period4" then &b * (4-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "period5" then &b * (5-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt" then &theta
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt*period1" then &theta_j * (1-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt*period2" then &theta_j * (2-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt*period3" then 0
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt*period4" then &theta_j * (4-3)
                 when COALESCE(result_mix2.varname, result_mm2.varname) = "trt*period5" then &theta_j * (5-3) end) as true_mixed,
           (calculated true_mixed / sqrt(1 + &vc)) as true_mm,
           result_mix2.estimate   as glmm_ind_est,
           result_mm2.estimate  as gee_ind_est,
           result_mix2.lower as glmm_ind_est_lower,
           result_mix2.upper as glmm_ind_est_upper,
           result_mm2.lowercl as gee_ind_est_lower,
           result_mm2.uppercl as gee_ind_est_upper,
           (glmm_ind_est - calculated true_mixed) as glmm_ind_bias,
           (gee_ind_est - calculated true_mm) as gee_ind_bias,
           (calculated glmm_ind_bias)**2 as glmm_ind_mse,
           (calculated gee_ind_bias)**2 as gee_ind_mse
    from result_mix2 full join result_mm2 on (result_mix2.sim = result_mm2.sim and result_mix2.varname = result_mm2.varname);
  quit;

  data results;
    merge result_agg result_ind;
    by sim parms true_mixed true_mm;
  run;

  data results;
    set results;
    where parms ~= "Scale";
    glmm_agg_cp = (glmm_agg_est_lower <= true_mixed <= glmm_agg_est_upper);
    gee_agg_cp = (gee_agg_est_lower <= true_mm <= gee_agg_est_upper);
    if (glmm_agg_est_lower = . or glmm_agg_est_upper = .) then glmm_agg_cp = .;
    if (gee_agg_est_lower = . or gee_agg_est_upper = .) then gee_agg_cp = .;

    glmm_ind_cp = (glmm_ind_est_lower <= true_mixed <= glmm_ind_est_upper);
    gee_ind_cp = (gee_ind_est_lower <= true_mm <= gee_ind_est_upper);
    if (glmm_ind_est_lower = . or glmm_ind_est_upper = .) then glmm_ind_cp = .;
    if (gee_ind_est_lower = . or gee_ind_est_upper = .) then gee_ind_cp = .;
  run;

  proc sql;
    create table power_agg as  
    select COALESCE(power_mix.sim, power_mm.sim) as sim,
           COALESCE(power_mix.effect, power_mm.source) as parms, 
           power_mix.ProbChiSq as glmm_agg_chisq,
           power_mm.probchisq  as gee_agg_chisq
    from power_mix full join power_mm on (power_mix.sim = power_mm.sim and power_mix.effect = power_mm.source);
  quit;

  proc sql;
    create table power_ind as  
    select COALESCE(power_mix2.sim, power_mm2.sim) as sim,
           COALESCE(power_mix2.effect, power_mm2.source) as parms,
           power_mix2.ProbChiSq as glmm_ind_chisq,
           power_mm2.probchisq  as gee_ind_chisq
    from power_mix2 full join power_mm2 on (power_mix2.sim = power_mm2.sim and power_mix2.effect =power_mm2.source);
  quit;

  data powers;
    merge power_agg power_ind;
    by sim parms;
  run;

  data powers;
    set powers;
    glmm_agg_pow = (glmm_agg_chisq <= 0.05);
    gee_agg_pow = (gee_agg_chisq <= 0.05);
    glmm_ind_pow = (glmm_ind_chisq <= 0.05);
    gee_ind_pow = (gee_ind_chisq <= 0.05);
    if glmm_agg_chisq = . then glmm_agg_pow =.;
    if gee_agg_chisq = . then gee_agg_pow = .;
    if glmm_ind_chisq = . then glmm_ind_pow =.;
    if gee_ind_chisq = . then gee_ind_pow = .;
  run;

  proc means data = results n mean;
    var glmm_agg_bias gee_agg_bias glmm_ind_bias gee_ind_bias glmm_agg_cp gee_agg_cp glmm_ind_cp gee_ind_cp;
    class parms;
    output out=res_est;
  run;

  data res_est1;
    set res_est;
    where _TYPE_=1 and _STAT_="MEAN";
    drop _TYPE_ _FREQ_;
  run;
  
  data res_est2;
    set res_est;
    where _TYPE_=1 and _STAT_="N";
    drop _TYPE_ _FREQ_;
  run;

  proc transpose data=res_est1 out=long_est1(rename=(COL1=bias)) name=group;
    by parms;
    var glmm_agg_bias -- gee_ind_bias;
  run;
  
  data long_est1;
    set long_est1;
    method=catx("_", scan(group,1, "_"), scan(group,2,"_"));
    drop group;
  run;

  proc sort data=long_est1; by method parms; run;

  proc transpose data=res_est1 out=long_est2(rename=(COL1=cp)) name=group;
    by parms;
    var glmm_agg_cp -- gee_ind_cp;
  run;

  data long_est2;
    set long_est2;
    method=catx("_", scan(group,1, "_"), scan(group,2,"_"));
    drop group;
  run;

  proc sort data=long_est2; by method parms; run;
  
  proc transpose data=res_est2 out=long_est3(rename=(COL1=N)) name=group;
    by parms;
    var glmm_agg_bias -- gee_ind_bias;
  run;

  data long_est3;
    set long_est3;
    method=catx("_", scan(group,1, "_"), scan(group,2,"_"));
    drop group;
  run;

  proc sort data=long_est3; by method parms; run;

  proc means data = results n std;
    var glmm_agg_est glmm_ind_est gee_agg_est gee_ind_est;
    class parms;
    output out=res_se;
  run;

  data res_se;
    set res_se;
    where _TYPE_=1 and _STAT_="STD";
    drop _TYPE_ _FREQ_;
  run;

  proc transpose data=res_se out=long_se(rename=(COL1=SE)) name=group;
    by parms;
    var glmm_agg_est -- gee_ind_est;
  run;

  data long_se;
    set long_se;
    method=catx("_", scan(group,1, "_"), scan(group,2,"_"));
    drop group;
  run;
  
  proc sort data=long_se; by method parms; run;

  data final_est;
    merge long_est1 long_est2 long_est3 long_se;
    by method parms;
  run;
  
  data final_est;
    set final_est;
    vc = &vc;
    theta_j=&theta_j;
    design = "&design";
    c0=&c0; 
    c1=&c1; 
    c2=&c2;
  run;
  
  data final_est;
    retain design c0 c1 c2 vc theta_j method parms N bias SE cp;
    set final_est;
  run;

  proc export data=final_est
    outfile="results_est.xlsx"
    dbms=xlsx
    replace;
    sheet="&design;&vc;&c0;&c1;&c2;&theta_j";
  run;

  proc means data = powers n mean;
    var glmm_agg_pow -- gee_ind_pow;
    class parms;
    output out=res_pwr;
  run;
  
  data res_pwr;
    set res_pwr;
    where _TYPE_=1 and _STAT_="MEAN";
    drop _TYPE_ _FREQ_;
  run;

  proc transpose data=res_pwr out=long_pwr(rename=(COL1=power)) name=group;
    by parms;
    var glmm_agg_pow -- gee_ind_pow;
  run;

  data long_pwr;
    set long_pwr;
    method=catx("_", scan(group,1, "_"), scan(group,2,"_"));
    drop group;
  run;
  
  proc sort data=long_pwr; by method parms; run;

  data final_pwr;
    set long_pwr;
    vc = &vc;
    theta_j=&theta_j;
    design = "&design";
    c0=&c0; 
    c1=&c1; 
    c2=&c2;
  run;
  
  data final_pwr;
    retain design c0 c1 c2 vc theta_j method parms power;
    set final_pwr;
  run;
  
  proc export data=final_pwr
    outfile="results_pwr.xlsx"
    dbms=xlsx
    replace;
    sheet="&design;&vc;&c0;&c1;&c2;&theta_j";
  run;
quit;
%mend simulation;

/*************/;
/*Simulation  /;
/*************/;

/*number of cluster=48*/
%simulation(vc=0.2, theta_j=0,     design=pgd);
%simulation(vc=0.2, theta_j=-0.1,  design=pgd);
%simulation(vc=0.5, theta_j=0,     design=pgd);
%simulation(vc=0.5, theta_j=-0.1,  design=pgd);
%simulation(vc=0.8, theta_j=0,     design=pgd);
%simulation(vc=0.8, theta_j=-0.1,  design=pgd);

%simulation(vc=0.2, theta_j=0,     design=dsd);
%simulation(vc=0.2, theta_j=-0.1,  design=dsd);
%simulation(vc=0.5, theta_j=0,     design=dsd);
%simulation(vc=0.5, theta_j=-0.1,  design=dsd);
%simulation(vc=0.8, theta_j=0,     design=dsd);
%simulation(vc=0.8, theta_j=-0.1,  design=dsd);

%simulation(vc=0.2, theta_j=0,     design=usd);
%simulation(vc=0.2, theta_j=-0.1,  design=usd);
%simulation(vc=0.5, theta_j=0,     design=usd);
%simulation(vc=0.5, theta_j=-0.1,  design=usd);
%simulation(vc=0.8, theta_j=0,     design=usd);
%simulation(vc=0.8, theta_j=-0.1,  design=usd);

%simulation(vc=0.2, theta_j=0,     design=dsdopt, c0=14, c1=24);
%simulation(vc=0.2, theta_j=-0.1,  design=dsdopt, c0=14, c1=24);
%simulation(vc=0.5, theta_j=0,     design=dsdopt, c0=14, c1=24);
%simulation(vc=0.5, theta_j=-0.1,  design=dsdopt, c0=14, c1=24);
%simulation(vc=0.8, theta_j=0,     design=dsdopt, c0=14, c1=24);
%simulation(vc=0.8, theta_j=-0.1,  design=dsdopt, c0=14, c1=24);

%simulation(vc=0.2, theta_j=0,     design=usdopt, c0=6, c1=9, c2=9);
%simulation(vc=0.2, theta_j=-0.1,  design=usdopt, c0=6, c1=9, c2=9);
%simulation(vc=0.5, theta_j=0,     design=usdopt, c0=6, c1=9, c2=9);
%simulation(vc=0.5, theta_j=-0.1,  design=usdopt, c0=6, c1=9, c2=9);
%simulation(vc=0.8, theta_j=0,     design=usdopt, c0=6, c1=9, c2=9);
%simulation(vc=0.8, theta_j=-0.1,  design=usdopt, c0=6, c1=9, c2=9);

%simulation(vc=0.2, theta_j=0,     design=usdopt, c0=5, c1=10, c2=9);
%simulation(vc=0.2, theta_j=-0.1,  design=usdopt, c0=5, c1=10, c2=9);
%simulation(vc=0.5, theta_j=0,     design=usdopt, c0=5, c1=10, c2=9);
%simulation(vc=0.5, theta_j=-0.1,  design=usdopt, c0=5, c1=10, c2=9);
%simulation(vc=0.8, theta_j=0,     design=usdopt, c0=5, c1=10, c2=9);
%simulation(vc=0.8, theta_j=-0.1,  design=usdopt, c0=5, c1=10, c2=9);
