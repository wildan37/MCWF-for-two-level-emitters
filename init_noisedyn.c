#include "header_noise.h"

void init_params(void){
   
    d = 1;
    atomiclevel = 2;
    Nstates = (int)pow(atomiclevel, Nspins);
    Nbasis = 2*Nstates;
    Lx = (int)sqrt(Nspins);
    Ly = (int)sqrt(Nspins);
//    Nreal = 1e0;
    
    R3 = V0*V0*V0;
    R6 = V0*V0*V0*V0*V0*V0;
    delta = 2.20;
    ///Time propagation
    tmax = 2e1;
    dt = 1e-3;        //second
    notimesteps = (int)(tmax/dt);
    t = 0.0;
    Niter = 1e3;
//    if (p0==0)Gamma=0;
//    else Gamma = sqrt((1-2*p0)/p0);

    ///Interaction range
    epsilon = 1.0e0;
    alpha = 6;
    cut_off = 1.0/pow ( epsilon, 1.0/(alpha - 1.0));
    cut_off = cut_off >= Nspins ? Nspins-1 : cut_off;
    printf("%.2lf\n", cos(0.5*M_PI));
//    printf("#%u\t%lf\t%lf\t%u\n", realisation, Gamma, p0, realisation);
    
    //store vector
    excited_num = malloc (Nstates*sizeof(unsigned));
    k = malloc (Nstates*sizeof(unsigned*));
    re = malloc (Nstates*sizeof(double));
    im = malloc (Nstates*sizeof(double));
    localnoise = malloc(Nspins *sizeof(double*));
    
//    basis = malloc(Nstates*sizeof(double));
    y = malloc (Nbasis * sizeof(double));
    xold = malloc (Nbasis * sizeof(double));
//    collect = malloc(Ncol * sizeof(double));
    
    for (state = 0; state < Nstates; state+=1 ){
        k[state] = malloc (Nspins*sizeof(unsigned));
        
        if (k[state]==NULL){
            
            fprintf(stderr, "out of memory kstate\n");
            exit(-5);
        }
    }
    
//    sprintf (filename, "%datoms_V0_%.2lf_xi%.2lf_gbar%.2lf_delta%.2lf_noise%d_re%d.dat", Nspins, V0, xi, frictionbar, delta, deph, realisation);

//    sprintf (filename, "%datoms_V0%.2lf_Om%.2lf_po%.2lf_fric%.2lf_width%.2lf_mu%.2lf_dim%d_noise%d_int%d_re%d.dat", Nspins, V0, Omega, p0, friction, specwidth, mu, d, deph, typeint, realisation);

    sprintf (filename, "%datoms_V0%.2lf_Om%.2lf_gdec%.2lf_fric%.2lf_width%.2lf_mu%.2lf_dim%d_noise%d_int%d_re%d.dat", Nspins, V0, Omega, Gamma, friction, specwidth, mu, d, deph, typeint, realisation);
    //sprintf (filename, "wienernoise.dat");
//    sprintf(filenamestate, "%datoms_V0_%.2lf_xi%.2lf_gbar%.2lf_delta%.2lf_CBLR_re%d_checkstate.dat", Nspins, V0, xi, frictionbar, delta, realisation);

//    sprintf(filenamestate, "%datoms_V0_%.2lf_delta%.2lf_xi%.2lf_OBLR_re%d_checkstate.dat", Nspins, V0, delta, xi, realisation);

//    sprintf (filestat, "counting%d_V%.2lf_xi%.2lf_gbar%.2lf_delta%.2lf_OBLRglo_re%d.dat", Nspins, V0, xi, frictionbar, delta, realisation);
//
    if ((outfile = fopen(filename, "w"))==NULL){
        fprintf (stderr, "error writing file %s", filename);
        exit(-7);
    }
    
//    if ((outfile1 = fopen(filestat, "w"))==NULL){
//        fprintf (stderr, "error writing file %s", filestat);
//        exit(-10);
//    }
    
}

void free_alloc(void){
    
    free(excited_num);
    free(re);
    free(im);
    for (state = 0; state < Nstates; state+=1 ){
        free(k[state]);
    }
    free(k);
    
//    for (col = 0; col < Ncol; col+=1){
//        collect[col] = malloc(Nrow*sizeof(double));
//    }
//    free(collect);

    free(localnoise);
    fclose(outfile);
//    fclose(outfile1);
//    free(basis);
    free(y);
    free(xold);
//    printf("\n#%u\t%lf\t%u\t%u\n", realisation, delta, Nstates, realisation);
}

void init_ground_noise(void){
    
    for (base =0; base < Nbasis; base+=1){
        y[base] = 0;
    }
    
    y[0] =1;
}

void init_inverted_noise(void){
    
    for (base = 0; base < Nbasis; base+=1){
        y[base] = 0;
    }
    y[Nbasis-2] = 1;
}

void init_configuration(void){
    
    ///initializing global noise
    
    gettimeofday(&tv, NULL);
    nowtime = tv.tv_sec;
    init_genrand(tv.tv_usec);
    initnoise = 0.0;
    if (deph==1)noise  = initnoise;
    
    ///initializing local noise & basis
    for (state=0; state < Nstates; state+=1){
        
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            k[state][spin_i] = 0.0;
        }
        excited_num[state] = 0;
        
        re[state] = 0.0;
        im[state] = 0.0;
    }
    
    if (deph==0){
        for (spin_i = 0; spin_i < Nspins; spin_i+=1)localnoise [spin_i] = initnoise;
    }
    
    re[0] = 1.0;
    dW = sqrt (dt) * rand_gauss();
    
    init_ground_noise();
//    init_inverted_noise();
}

double max (double data[], long unsigned Ndata){
    
    unsigned i;
    double max_value;
    
    max_value = data[0];
    
    for (i = 1; i < Ndata; i+=1){
        
        if (data[i] > max_value) {
            max_value = data[i];
        }
    }
    
    return max_value;
}

