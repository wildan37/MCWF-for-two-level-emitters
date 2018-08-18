#include "header_noise.h"

double rand_gauss (void){

   static double V1, V2, S;
   static int phase = 0;
   double X;

   if(phase == 0) {
        do {
            double U1 = genrand_real2();
            double U2 = genrand_real2();
                                                     
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
                                                                                                    
        X = V1 * sqrt(-2 * log(S) / S);
   } else
   
   X = V2 * sqrt(-2 * log(S) / S);
   phase = 1 - phase;
   return X;     
}

void counting_excited_state(void){
    
    for (state = 0; state < Nstates; state+=1){
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            if ((1<<spin_i)&state)excited_num[state]+=1;
        }
    }

}
//sigmax basis
void flip_state (void){
    
    for (state = 0; state < Nstates; state+=1){
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            k[state][spin_i] = state^(1<<spin_i);
            
        }
    }
    
}

void steady_noise(void){
    
    for (spin_i = 0; spin_i < Nspins; spin_i+=1){
        localnoise[spin_i] = rand_gauss()*mu;
        printf ("%lf ", localnoise[spin_i]);
    }
    
}

void local_noise(void){
    ///igor parameter
    /*
    for (spin_i = 0; spin_i < Nspins; spin_i+=1){
        dW = sqrt (dt) * rand_gauss();
        localnoise [spin_i] -= dt * localnoise[spin_i] * frictionbar;
        localnoise [spin_i] += 2* frictionbar * dW;
    }
    */
    
    ///general parameter
    for (spin_i = 0; spin_i < Nspins; spin_i+=1){
        dW = sqrt (dt) * rand_gauss();
        localnoise [spin_i] -= dt * localnoise[spin_i] * friction;
        localnoise [spin_i] += friction * sqrt(specwidth) * dW;
    }
    
}

void globalnoise(void){
    
    dW = sqrt (dt) * rand_gauss();
    ///igor parameter
//    noise -= dt * frictionbar * noise;
//    noise += 2 * frictionbar * dW;
    
    ///general parameter
    noise -= dt * friction * noise;
    noise += friction * sqrt(specwidth) * dW;
}

void eulermeruyama_global(int nruns){
    int step;
    double dW_sum;
    double dt_small;
    
    dt_small = dt/nruns;

    dW_sum = 0.0;
    for (step = 0; step < nruns; step+=1){
        dW_sum += rand_gauss();
    }
    
    noise -= dt * frictionbar * noise;
    noise += 2 * frictionbar * dW_sum*sqrt(dt_small);
}

void eulermeruyama_local(int nruns){
    
    int step;
    double dW_sum;
    double dt_small;
    
    dt_small = dt/nruns;
    
    dW_sum = 0.0;
    for (step = 0; step < nruns; step+=1){
        dW_sum += rand_gauss();
    }

    dW_sum*=dW_sum*sqrt(dt_small);
    
    for (spin_i = 0; spin_i < Nspins; spin_i+=1){
        localnoise [spin_i] -= dt * localnoise[spin_i] * frictionbar;
        localnoise [spin_i] += 2* frictionbar * dW_sum*sqrt(dt_small);
    }
}

double globalnoise_rk (double x, double gammabar){
    
    double k1, k2;
    
    dW = sqrt(dt)*rand_gauss();
    k1 = -gammabar*dt*x + 2*gammabar * dW;
    k2 = -gammabar*dt*(x+k1) + 2*gammabar * dW;
    x += 0.5*(k1+k2);
    
    return x;
}

void spectrum_symmetric (int Nspins){
    
//    char filenamespec[100];
//    FILE *outspec;
    double totalevec_i;
    double rsquare;
    double rad;
    double tingali;
    int i, ix, iy, iz;
    int j, jx, jy, jz;
    int nlargest;
    size_t *p;
    p = malloc(Nspins*sizeof(int));
    nlargest = 2;
    /*
    sprintf (filenamespec, "spec%datoms_V0%.2lf_Om%.2lf_gdec%.2lf_fric%.2lf_width%.2lf_delta%.2lf_noise%d_int%d_re%d.dat", Nspins, V0, Omega, Gamma, friction, specwidth, delta, deph, typeint, realisation);
    
    if ((outspec=fopen(filenamespec, "w"))==NULL){
        fprintf(stderr, "error writing file %s", filenamespec);
        exit(-1);
    }
    */
    gsl_matrix *Hsym;
    gsl_matrix *evec;
    gsl_vector *eval;
    gsl_vector *sumevec;
    double *xi;
    
    Hsym = gsl_matrix_alloc(Nspins, Nspins);
    evec = gsl_matrix_alloc(Nspins, Nspins);
    eval  = gsl_vector_alloc(Nspins);
    sumevec = gsl_vector_alloc(Nspins);
    xi = malloc((Nspins)*sizeof(double));
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(Nspins);
    gsl_matrix_set_zero(Hsym);

    /*
     for (spin_i = 0; spin_i < Nspins; spin_i+=1){
     for (spin_j = 0; spin_j < Nspins; spin_j+=1){
     yeig[spin_i*Nspins + spin_j] = 0.0;
     yeigim[spin_i*Nspins + spin_j] = 0.0;
     }
     }
     */
    
    for (i = 0; i < Nspins; i+=1){
        tingali = 0.0;
        rad = 0.0;
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            //interactions
            
            if (d==1){
                for (spin_j = 0; spin_j < Nspins; spin_j+=1) {
                    if (spin_i!=spin_j){
                        rad = abs((int)spin_i - (int)spin_j);
                        if ((i!=spin_i&&i==spin_j)){
                            gsl_matrix_set (Hsym, spin_i, spin_j, (V0/2.0)/(rad*rad*rad));
                        }
                    }
                }
            }
            
            if (d==2){
                iy = spin_i/Lx;
                ix = spin_i - Lx*iy;
                
                for (j = 0; j < Nspins; j+=1){
                    if (spin_i!=j){
                        jy = j/Lx;
                        jx = j - Lx*jy;
                        
                        spin_j = (jx + Lx*jy);
                        rsquare = (jx-ix)*(jx-ix) + (jy-iy)*(jy-iy);
                        rad = sqrt(rsquare);
//                        rad = 1;
//                        printf("(%d,%d)=%.2lf\t",jx,jy, rad);
//                        printf("(%d)=%.2lf\t",spin_j, rad);
                        
                        //forster interactions
                        if ((i!=spin_i&&i==spin_j)){
                            gsl_matrix_set (Hsym, spin_i, spin_j, (V0/2.0)/(rad*rad*rad));
                        }
                    }
                }
                
                //                printf("\n");
            }
            
            if (d==3) {
                iz = spin_i/(Ly*Ly);
                iy = spin_i/Lx - iz*Ly;
                ix = spin_i - iy*Lx - iz*(Lz*Lz);
                for (j = 0; j < Nspins; j+=1) {
                    if (spin_i!=j){
                        jz = j/(Ly*Ly);
                        jy = j/Lx - jz*Ly;
                        jx = j - jy*Lx - jz*(Lz*Lz);
                        
                        spin_j = (jx + Lx*jy + jz*Ly*Ly);
                        rsquare = (jx-ix)*(jx-ix) + (jy-iy)*(jy-iy) + (jz-iz)*(jz-iz);
                        rad = sqrt(rsquare);
//                        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf\n",
//                               spin_i, spin_j, ix, iy, iz, jx, jy, jz, rad);
                        
                        ///forster interactions
                        
                        if (i!=spin_i&&i==spin_j){
                            gsl_matrix_set(Hsym, spin_i, spin_j, (V0/2)/(rad*rad*rad));
                        }
                        
                    }
                }
            }
            
        }
    }

    /*
     for (i = 0; i < (Nspins-rubah); i+=1){
        for (j = 0; j < (Nspins-rubah); j+=1){
     //     if (i==j)gsl_matrix_set(Hsym, i, j, delta);
            printf("%.2lf\t", gsl_matrix_get(Hsym, i, j));
        }
        printf("\n");
     }
     */
    
    gsl_eigen_symmv (Hsym, eval, evec, w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_symmv_free(w);
    
    for (i = 0; i < Nspins; i+=1){
        //        printf("%.6lf\t", gsl_vector_get(eval, i));
        totalevec_i = 0.0;
        for (j = 0; j < Nspins; j+=1){
            totalevec_i += gsl_matrix_get(evec, j, i);
            /*
            totalevec_i += GSL_REAL(gsl_matrix_complex_get(evec, i, j));
            yeig[i*(Nspins-rubah) + j] = gsl_matrix_get(evec, j, i);
            yeig[i*(Nspins-rubah) + j] = GSL_REAL(gsl_matrix_complex_get(evec, j, i));
            yeigim[i*(Nspins-rubah) + j] = GSL_IMAG(gsl_matrix_complex_get(evec, j, i));
            printf("%.6lf\t", gsl_matrix_get(evec, i, j));
                        printf("%.6lf\t", gsl_matrix_get(Hsym, i, j));
            printf("%.6e\t", yeig[i*(Nspins-rubah)+j]);
            */
        }
        //        printf("\n");
        gsl_vector_set(sumevec, i, totalevec_i);
        totalevec_i = gsl_vector_get(sumevec, i)*gsl_vector_get(sumevec, i);
        totalevec_i /= Nspins;
        xi[i] = totalevec_i;
        //        printf("%.6lf\t%.6lf\n", gsl_vector_get(eval, i), totalevec_i);
//        fprintf(outspec, "%.6lf\t%.6lf\n", gsl_vector_get(eval, i), totalevec_i);
    }
    
    gsl_sort_largest_index(p, nlargest, xi, 1, Nspins);
    
    /*
     for (i = 0; i < nlargest;i+=1){
     printf("%d\t%.3lf\t%.3lf\n", p[i], gsl_vector_get(eval, p[i]), xi[p[i]]);
     }
     */
    
    Delta_E = gsl_vector_get(eval, p[0]) - gsl_vector_get(eval, p[1]);
    delta = - gsl_vector_get(eval, p[0]);
    printf("%.6e\t%.6e\t", tingali, delta);
//    printf("%.6e\t%.6e\n",Delta_E, delta);

    gsl_vector_free(sumevec);
    gsl_matrix_free(Hsym);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    free(xi);
    free(p);
    
//    fclose(outspec);
    
}

int
func_wave(double t, const double y[], double f[],
          void *params){
    double sumreal;
    double sumimag;
    int jx, ix;
    int jy, iy;
    int jz, iz;
    int dx, dy;
    unsigned j;
    double rsquare;
    int sigma_up;
    int sigma_down;
    double rad;
    
    struct ODE_params *p = (struct ODE_params *)params;
    double Delta0 = p->Delta0;
    
    
    for (state = 0; state < Nstates; state+=1){
        sumimag = 0.0;
        sumreal = 0.0;
        
        ///rabi frequency
        
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            sumimag+= Omega/2 * y[2*(k[state][spin_i]) + 1];
            sumreal-= Omega/2 * y[2*(k[state][spin_i])];
        }
        
        ///laser detuning
        
        sumimag +=  delta * y[2*state + 1] * excited_num[state];
        sumreal -=  delta * y[2*state] * excited_num[state];
        
        //global noise
        if (deph==1){
            sumimag -= Delta0 * y[2*state + 1] * excited_num[state];
            sumreal += Delta0 * y[2*state] * excited_num[state];
        }
        
        if (deph==0){
            //local noise
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if ((1<<spin_i)&state){
                    sumimag -=  localnoise[spin_i] * y[2*state + 1];
                    sumreal +=  localnoise[spin_i] * y[2*state];
                }
            }
        }
        
        ///spontaneous decay
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            if ((1<<spin_i)&state){
                sumimag -= Gamma/2 * y[2*state];
                sumreal -= Gamma/2 * y[2*state + 1];
            }
        }
    
        ///interactions
        if (V0 > 0){
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if (d==1){
                    /*
                    for (dx = -(int)cut_off; dx <= (int)cut_off; dx+=1){
                         if (dx!=0){
                             jx = spin_i + dx;
                             if (jx < 0)jx+=Nspins;
                             if (jx >= Nspins)jx-=Nspins;
                                 spin_j = jx;
                                 rsquare = dx*dx;
                                 rad = sqrt(rsquare);
                         
                             if (typeint ==0){
                                 ///vdW interaction
                                 if (rad < 1.1){
                                     // printf("(%d,%d)=%.2lf\t",jx,jy, rad);
                                     if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                         sumimag += (V0/2.0) * y[2 * state +1]
                                         /(rsquare*rsquare*rsquare);
                                         sumreal -= (V0/2.0)  *y[2 * state]
                                         /(rsquare*rsquare*rsquare);
                                     }
                                 }
                             }
                             if (typeint ==1){
                                 ///dipole dipole interaction
                                 if (rad < 1.1){
                                     if ((1<<spin_j)&state){
                                         sigma_down = (1<<spin_j)^state;
                                         if (!((1<<spin_i)&sigma_down)){
                                             sigma_up = (1<<spin_i)^sigma_down;
                                             
                                             sumimag += (V0/2.0) * y[2 * sigma_up+1]
                                             /(rad*rad*rad);
                                             sumreal -= (V0/2.0) * y[2 * sigma_up]
                                             /(rad*rad*rad);
                                         }
                                     }
                                 }
                             }
                         }
                     }
                    */
                    for (spin_j = 0; spin_j < Nspins; spin_j+=1){
                        if (spin_i!=spin_j){
                            rad = abs((int)(spin_i - spin_j));
                            rsquare = rad*rad;
//                            printf("%d\n", rad1);
                            if (typeint == 0){
                                ///vdW interaction
                                if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                    sumimag += (V0/2.0) * y[2 * state +1]
                                    /(rsquare*rsquare*rsquare);
                                    sumreal -= (V0/2.0)  *y[2 * state]
                                    /(rsquare*rsquare*rsquare);
                                }
                            }
                            else{
                                ///dipole dipole interaction
                                if ((1<<spin_j)&state){
                                    sigma_down = (1<<spin_j)^state;
                                    if (!((1<<spin_i)&sigma_down)){
                                        sigma_up = (1<<spin_i)^sigma_down;
                                        
                                        sumimag += (V0/2.0) * y[2 * sigma_up+1]
                                        /(rad*rad*rad);
                                        sumreal -= (V0/2.0) * y[2 * sigma_up]
                                        /(rad*rad*rad);
                                    }
                                }
                            }
                        }
                    }
                    
                }
                if (d==2){
                    
                    iy = spin_i/Lx;
                    ix = spin_i - Lx*iy;
//                    printf("\n(%d,%d),\t",ix,iy);
                    
                    for (dy = -(int)cut_off; dy <= (int)cut_off; dy+=1){
                        for (dx = -(int)cut_off; dx <= (int)cut_off; dx+=1){
                            if (dx!=0||dy!=0){
                                jx = ix + dx;
                                jy = iy + dy;
                                
                                if (jx<0)jx+=Lx;
                                if (jx>=Lx)jx-=Lx;
                                if (jy<0)jy+=Ly;
                                if (jy>=Ly)jy-=Ly;
                                
                                spin_j = (jx + Lx*jy);
                                rsquare = (jx-ix)*(jx-ix) + (jy-iy)*(jy-iy);
                                rad = sqrt(rsquare);
//                                rad = 1.0;
//                                printf("(%d,%d)=%.2lf\t",jx,jy, rad);
                            
                            if (typeint ==0){
                                ///vdW interaction
                                if (rad < 1.1){
//                                    printf("(%d,%d)=%.2lf\t",jx,jy, rad);
                                    if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                        sumimag += (V0/2.0) * y[2 * state +1]
                                        /(rsquare*rsquare*rsquare);
                                        sumreal -= (V0/2.0)  *y[2 * state]
                                        /(rsquare*rsquare*rsquare);
                                    }
                                }
                            }
                            if (typeint ==1){
                                ///dipole dipole interaction
                                if (rad < 1.1){
                                    if ((1<<spin_j)&state){
                                        sigma_down = (1<<spin_j)^state;
                                        if (!((1<<spin_i)&sigma_down)){
                                            sigma_up = (1<<spin_i)^sigma_down;
                                            
                                            sumimag += (V0/2.0) * y[2 * sigma_up+1]
                                            /(rad*rad*rad);
                                            sumreal -= (V0/2.0) * y[2 * sigma_up]
                                            /(rad*rad*rad);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                     
                    /*
                    for (j = 0; j < Nspins; j+=1){
                        if (spin_i!=j){
                            jy = j/Lx;
                            jx = j - Lx*jy;
                            
                            spin_j = (jx + Lx*jy);
                            rsquare = (jx-ix)*(jx-ix) + (jy-iy)*(jy-iy);
                            rad = sqrt(rsquare);
//                            rad = 1.0;
//                            printf("(%d,%d)=%.2lf\t",jx,jy, rad);
                            
                            if (typeint ==0){
                                ///vdW interaction
                                if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                    sumimag += (V0/2.0) * y[2 * state +1]
                                    /(rsquare*rsquare*rsquare);
                                    sumreal -= (V0/2.0)  *y[2 * state]
                                    /(rsquare*rsquare*rsquare);
                                }
                            }
                            if (typeint ==1){
                                ///dipole dipole interaction
                                if ((1<<spin_j)&state){
                                    sigma_down = (1<<spin_j)^state;
                                    if (!((1<<spin_i)&sigma_down)){
                                        sigma_up = (1<<spin_i)^sigma_down;
                                        
                                        sumimag += (V0/2.0) * y[2 * sigma_up+1]
                                        /(rad*rad*rad);
                                        sumreal -= (V0/2.0) * y[2 * sigma_up]
                                        /(rad*rad*rad);
                                    }
                                }
                            }
                        }
                    }
                    */
                    
//                    printf("\n");
                }
                if (d==3){
                    iz = spin_i /(Ly*Ly);
                    iy = spin_i /Lx - iz*Ly;
                    ix = spin_i  - iy*Lx - iz*(Lz*Lz);
                    for (j = 0; j < Nspins; j+=1) {
                        if (spin_i!=j){
                            jz = j/(Ly*Ly);
                            jy = j/Lx - jz*Ly;
                            jx = j - jy*Lx - jz*(Lz*Lz);
                            
                            spin_j = (jx + Lx*jy + jz*Ly*Ly);
                            rsquare = (jx-ix)*(jx-ix) + (jy-iy)*(jy-iy) + (jz-iz)*(jz-iz);
                            rad = sqrt(rsquare);
                            
                            if (typeint ==0){
                                ///vdW interaction
                                if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                    sumimag += (V0/2.0) * y[2 * state +1]
                                    /(rsquare*rsquare*rsquare);
                                    sumreal -= (V0/2.0)  *y[2 * state]
                                    /(rsquare*rsquare*rsquare);
                                }
                            }
                            if (typeint ==1){
                                ///dipole dipole interaction
                                if ((1<<spin_j)&state){
                                    sigma_down = (1<<spin_j)^state;
                                    if (!((1<<spin_i)&sigma_down)){
                                        sigma_up = (1<<spin_i)^sigma_down;
                                        
                                        sumimag += (V0/2.0) * y[2 * sigma_up+1]
                                        /(rad*rad*rad);
                                        sumreal -= (V0/2.0) * y[2 * sigma_up]
                                        /(rad*rad*rad);
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
        }
        
        f[2* state] = sumimag;
        f[2* state + 1] = sumreal;
    }
    
    return GSL_SUCCESS;

}

int
funcnoise (double t, const double y[], double f[],
          void *params){
    //    int col;
    double sumreal;
    double sumimag;
//    int jx, dx;
//    double rsquare;
//    int sigma_up;
//    int sigma_down;
//    int a,b;
    double rad;
    
    struct ODE_params *p = (struct ODE_params *)params;
    double Delta0 = p->Delta0;
    
    
    for (state = 0; state < Nstates; state+=1){
        sumimag = 0.0;
        sumreal = 0.0;
        
        ///rabi frequency
        
        for (spin_i = 0; spin_i < Nspins; spin_i+=1){
            sumimag+= (xi/8) * y[2*(k[state][spin_i]) + 1];
            sumreal-= (xi/8) * y[2*(k[state][spin_i])];
        }

        ///laser detuning
        
        sumimag +=  delta * xi*xi/8  * y[2*state + 1] * excited_num[state];
        sumreal -=  delta * xi*xi/8 * y[2*state] * excited_num[state];

        //global noise
        if (deph==1){
            sumimag -= Delta0 * xi/4 * y[2*state + 1] * excited_num[state];
            sumreal += Delta0 * xi/4 * y[2*state] * excited_num[state];
        }
        
        if (deph==0){
            //local noise
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if ((1<<spin_i)&state){
                    sumimag -=  localnoise[spin_i]  * xi/4 * y[2*state + 1];
                    sumreal +=  localnoise[spin_i] * xi/4 * y[2*state];
                }
            }
        }
        
        ///interactions
        if (V0 > 0){
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                /*
                for (dx = -(int)cut_off; dx <= (int)cut_off; dx+=1){
                    if (dx!=0){
                        jx = spin_i + dx;
                        if (jx>=0&&jx<Nspins){
                            
                            if (jx < 0)jx+=Nspins;
                            if (jx >= Nspins)jx-=Nspins;
                            
                            spin_j = jx;
                            rsquare = dx*dx;
                            rad = sqrt(rsquare);
                            
                            ///vdW interaction
                            if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                                
                                sumimag += xi * xi/8 * R6 * y[2 * state +1]
                                /(rsquare*rsquare*rsquare);
                                
                                sumreal -= xi * xi/8 * R6  *y[2 * state]
                                /(rsquare*rsquare*rsquare);
                            }
                            
                            
                            if ((1<<spin_j)&state){
                                sigma_down = (1<<spin_j)^state;
                                if (!((1<<spin_i)&sigma_down)){
                                    sigma_up = (1<<spin_i)^sigma_down;
                                   
                                    sumimag += xi * xi/8 * R3 * y[2 * sigma_up+1]
                                    /(rad*rad*rad);
                                    sumreal -= xi * xi/8 * R3 * y[2 * sigma_up]
                                    /(rad*rad*rad);
                                    
                                }
                            }
                            
                        }
                    }
                }
                 */
                
                for (spin_j = spin_i+1; spin_j < Nspins; spin_j+=1){
                    
                    rad = abs(spin_i - spin_j);
                    
                    if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                        
                        sumimag += xi*xi/4 * R6 * y[2 * state +1]
                        /(rad*rad*rad*rad*rad*rad);
                        
                        sumreal -= xi*xi/4 * R6  *y[2 * state]
                        /(rad*rad*rad*rad*rad*rad);
                    }
                }
                
            }
        }
        
        f[2* state] = sumimag;
        f[2* state + 1] = sumreal;
    }
    
    return GSL_SUCCESS;
    
}

void trajectory_noise (void){
    
    double norm;
    double expectation;
    double expectation2;
    double dptot;
    unsigned int jdown;
    int w;
    double ground;
//    double *autocor;
    dp = malloc(Nspins*sizeof(double));
//    autocor = malloc((notimesteps+1)*sizeof(double));
    const gsl_odeiv2_step_type * T
    = gsl_odeiv2_step_rk8pd;
    
    //    write_data_traj();
    
    gsl_odeiv2_step * s
    = gsl_odeiv2_step_alloc (T, Nbasis);
    gsl_odeiv2_control * c
    = gsl_odeiv2_control_y_new (1e-12, 1e-6);
    gsl_odeiv2_evolve * e
    = gsl_odeiv2_evolve_alloc (Nbasis);

    w = 0;
    steady_noise();
    while (t < tmax){
        for (spin_i=0; spin_i < Nspins; spin_i+=1) {
            dp[spin_i] = 0.0;
        }
        for (state = 0; state <Nstates; state+=1) {
            xold[2*state] = y[2*state];
            xold[2*state+1] = y[2*state+1];
        }
        
        struct ODE_params rod_param={noise};//ran26
//        gsl_odeiv2_system sys = {funcnoise, 0, Nbasis, &rod_param};
        gsl_odeiv2_system sys = {func_wave, 0, Nbasis, &rod_param};
        
        int status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &t, dt, y);

        if (status != GSL_SUCCESS){
            break;
        }
//        autocor[w] = localnoise[0];
        dptot = 0.0;
        for (spin_i =0; spin_i < Nspins; spin_i+=1) {
            for (state = 0; state <Nstates; state+=1) {
                if ((1<<spin_i)&state){
                    dp[spin_i]+=dt*Gamma*(xold[2*state]*xold[2*state] + xold[2*state + 1]*xold[2*state + 1]);
                }
            }
            dptot+=dp[spin_i];
        }
        w+=1;
        
        expectation = 0.0;
        expectation2 = 0.0;
        norm = 0.0;
        ground = 0.0;
        ground+=(y[2*0]*y[2*0] + y[2*0 + 1]*y[2*0 + 1]);
        
        for (state = 0; state < Nstates; state+=1){

            norm += (y[2*state]*y[2*state] + y[2*state + 1]*y[2*state + 1]);
            expectation += (y[2*state]*y[2*state] + y[2*state + 1]*y[2*state + 1])
            *excited_num[state];
            expectation2 += (y[2*state]*y[2*state] + y[2*state + 1]*y[2*state + 1])
            *excited_num[state]*excited_num[state];            
        }
//        printf("%.3lf\t%.6e\t%.6e\n", t, norm,dptot);


        if (dptot < genrand_real2()){
            for (state = 0; state < Nstates; state+=1){
                y[2*state]/=sqrt(norm);
                y[2*state+1]/=sqrt(norm);
            }
        }
        else{
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if (dp[spin_i]/dptot > genrand_real2()){
                    for (state=0; state < Nstates; state+=1){
                        if ((1<<spin_i)&state){
                            jdown = (1<<spin_i)^state;
                            y[2*jdown] = sqrt(Gamma)*y[2*state]/sqrt(dp[spin_i]/dt);
                            y[2*jdown+1] = sqrt(Gamma)*y[2*state+1]/sqrt(dp[spin_i]/dt);
                            y[2*state] = 0;
                            y[2*state+1] = 0;
                        }
                    }
                }
            }
        }
        
        
        if (t < 1.0){
//            dt = 1e-3;
            if ((w%5)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else if (t < 10.0){
//            dt = 1e-3;
            if ((w%10)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else if (t < 100.0){
//            dt = 1e-3;
            if ((w%20)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else if (t < 1000.0){
//            dt = 1.0e-3;
//            dt = 7.5e-4;
            if ((w%50)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else if (t < 10000.0){
//            dt = 4.0e-4;
            if ((w%25000)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else if (t < 100000){
//            dt = 3e-4;
            if ((w%50000)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
                //                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        else{
//            dt = 2.5e-4;
            if ((w%100000)==0){
                fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
//                printf("%.3lf\t%.6e\t%.6e\t%.6e\t%.6e\n", t, norm, ground, expectation, expectation2);
            }
        }
        
        if (deph==1)globalnoise();
        if (deph==0)local_noise();


    }
//    FILE *corelat;
//    char filename_cor[80];
//    sprintf(filename_cor, "noiseauto.txt");
//    corelat = fopen(filename_cor, "w");
//    double total;
//    for (i = 0; i < notimesteps; i+=1){
//        total = autocor[i]*autocor[i+1];
//        fprintf(corelat, "%d\t%.3lf\t%.3lf\n", i, autocor[i], total);
//    }

    free(dp);
//    fclose(corelat);
//    free(autocor);
    printf ("%d\n", realisation);
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    
}

void euller_igor (void){
    
    long int i;
    int n;
//    long int j;
    
//    int dx;
//    int jx;
//    double count_prob_state;
//   int sigma_up, sigma_down; 
    //double sum;
    //double Xtrue;
    //sum = 0.0;
    
    i = 0;
    t = 0.0;
    n = 16;
    while (t <= tmax){
        
        for ( state = 0; state < Nstates; state+=1){
        
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re[state]+= dt/8 * (xi/4) * im[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im[state] -= dt/4 * (xi/4) * re[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re[state] += dt/8 * (xi/4) * im[k[state][spin_i]];
            }
			
            re [state] += dt/4.0 * delta * xi*xi/8 * im[state] * excited_num [state];
            im [state] -= dt/2.0 * delta * xi*xi/8 * re[state] * excited_num [state];
            re [state] += dt/4.0 * delta * xi*xi/8 * im[state] * excited_num [state];
            
            if (deph==1){
                re [state] -= dt/4.0 * noise * xi/4 * im[state] * excited_num [state];
                im [state] += dt/2.0 * noise * xi/4 * re[state] * excited_num [state];
                re [state] -= dt/4.0 * noise * xi/4 * im[state] * excited_num [state];
            }
            if (deph == 0){
                for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if (state & (1<<spin_i)){
                    re [state] -= dt/4.0 * localnoise[spin_i] * xi/4 * im[state];
                    im [state] += dt/2.0 * localnoise[spin_i] * xi/4 * re[state];
                    re [state] -= dt/4.0 * localnoise[spin_i] * xi/4 * im[state];
                }
            }
            }
            
            if (V0 > 0.0){
                for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                    /*
                    for (dx = -(int)cut_off; dx <= (int)cut_off; dx+=1){
                        if (dx!=0){
                            jx = ilattice + dx;
                            if (jx>=0&&jx<Nspins){
                            
                                if (jx < 0)jx+=Nspins;
                                if (jx >= Nspins)jx-=Nspins;
                                
                                jlattice = jx;
                                rsquare = dx*dx;
                                rad = sqrt(rsquare);   
    //                            printf ("%.2lf\t%d\t%d\n", rad, ilattice, jlattice); 
                                //if (rsquare < 1.1){
    //                            if ( (state&(1<<ilattice))&&(state&(1<<jlattice))){
                                    im [state] -= dt/4.0 * re[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
                                    re [state] += dt/2.0 * im[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
                                    im [state] -= dt/4.0 * re[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
                                    
                                    re [state] += dt/4.0 * im[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
                                    im [state] -= dt/2.0 * re[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
                                    re [state] += dt/4.0 * im[state]*xi*xi/8 * R6/(rsquare*rsquare*rsquare);
    //                            }

                                if ((1<<jlattice)&state){
                                    sigma_down = (1<<jlattice)^state;
                                    if (!((1<<ilattice)&sigma_down)){
                                        sigma_up = (1<<ilattice)^sigma_down;
    //                                    printf ("%d\t%d\n", state, sigma_up);                                   
                                        
                                        re[state] += dt/4.0 * xi * xi/8 * R3 * im[sigma_up]
                                        /(rad*rad*rad);
                                        im[state] -= dt/2.0 *xi * xi/8 * R3 * re[sigma_up]
                                        /(rad*rad*rad);
                                        re[state] += dt/4.0 * xi * xi/8 * R3 * im[sigma_up]
                                        /(rad*rad*rad);

                                        im[state] -= dt/4.0 *xi * xi/8 * R3 * re[sigma_up]
                                        /(rad*rad*rad);
                                        re[state] += dt/2.0 * xi * xi/8 * R3 * im[sigma_up]
                                        /(rad*rad*rad);
                                        im[state] -= dt/4.0 *xi * xi/8 * R3 * re[sigma_up]
                                        /(rad*rad*rad);

                                    }
                                }
                            }
                        }
                    }
                     */
                    for (spin_j = spin_i+1; spin_j < Nspins; spin_j+=1){
                        
                        rad = abs(spin_i - spin_j);
                        
                        if ( (state&(1<<spin_i))&&(state&(1<<spin_j))){
                            
                            im [state] -= dt/4.0 * re[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);
                            re [state] += dt/2.0 * im[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);
                            im [state] -= dt/4.0 * re[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);
                            
                            re [state] += dt/4.0 * im[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);
                            im [state] -= dt/2.0 * re[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);
                            re [state] += dt/4.0 * im[state]*xi*xi/4 * R6/(rad*rad*rad*rad*rad*rad);

                        }
                    }
                }
            }

            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im [state] -=dt/8 * (xi/4) * re[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re [state] += dt/4 * (xi/4) * im[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im [state] -= dt/8 * (xi/4) * re[k[state][spin_i]];
            }
            
            im [state] -= dt/4.0 * delta * xi*xi/8 * re[state] * excited_num [state];
            re [state] += dt/2.0 * delta * xi*xi/8 * im[state] * excited_num [state];
            im [state] -= dt/4.0 * delta * xi*xi/8 * re[state] * excited_num [state];
            
            if (deph==1){
                im [state] += dt/4.0 * noise * xi/4 * re[state] * excited_num [state];
                re [state] -= dt/2.0 * noise * xi/4 * im[state] * excited_num [state];
                im [state] += dt/4.0 * noise * xi/4 * re[state] * excited_num [state];
            }
            
            if (deph == 0){
                for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                if (state & (1<<spin_i)){
                    im [state] += dt/4.0 * localnoise[spin_i] * xi/4 * re[state];
                    re [state] -= dt/2.0 * localnoise[spin_i] * xi/4 * im[state];
                    im [state] += dt/4.0 * localnoise[spin_i] * xi/4 * re[state];
                }
            }
            }
        }
        
        prob = 0.0;
        Nryd = 0.0;
        Nryd2 = 0.0;
        
		numerator = 0.0;
        denom1 = 0.0;
        denom2 = 0.0;
        /*
        for (state = 0; state < Nstates; state+=1){
            basis[state] = 0.0;
        }
        */
//        if((i%100)==0)fprintf (outfile, "%.3lf\t", t);
        
        for (state = 0; state < Nstates; state+=1){
            prob+=(re[state]*re[state] + im[state]*im[state]);
//            basis[state]+=(re[state]*re[state] + im[state]*im[state]);
        }
        
        for (state = 0; state < Nstates; state+=1){
            Nryd +=(re[state]*re[state] + im[state]*im[state])*excited_num[state];
			Nryd2 += (re[state]*re[state] + im[state]*im[state])*excited_num[state]* excited_num[state];
//            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
//                numerator+=(re[k[state][spin_i]]*re[k[state][spin_i]]+im[k[state][spin_i]]*im[k[state][spin_i]]);
//            }

			//if ((i%100)==0)printf ("%.12e\t", re[state]*re[state] + im[state]*im[state]);
//		    if((i%100)==0)fprintf (outfile, "%.6e\t", re[state]*re[state] + im[state]*im[state]);
	    }
        
//		if((i%100)==0)fprintf (outfile, "\n");
        //DeltaNryd = sqrt ( Nryd2 - (Nryd * Nryd) );
		//Q = ((Nryd2 - (Nryd * Nryd))/Nryd)-1;
        //g2 = 2*numerator/sqrt(denom1*denom2) ;
        
        if (t < 1.0){
            if ((i%10)==0){
                
    //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
    //                    i*dt, prob, Nryd, Nryd2, numerator);
                
                fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", i*dt, prob, Nryd, Nryd2);
            
            }
        }
        else if (t < 10.0){
            
            if ((i%75)==0){
                
    //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
    //                   i*dt, prob, Nryd, Nryd2, numerator);
                 fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            }
        }
        else if (t < 100.0){
            
            if ((i%250)==0){
                
                //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
                //                   i*dt, prob, Nryd, Nryd2, numerator);
                fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            }
        }
        else if (t < 1000.0){
            dt = 5e-4;
            if ((i%1000)==0){
                //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
                //                   i*dt, prob, Nryd, Nryd2, numerator);
                fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            }
        }
        else if (t < 10000.0){
            dt = 5e-4;
            if ((i%5000)==0){
                
                //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
                //                   i*dt, prob, Nryd, Nryd2, numerator);
                fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            }
        }
        else if (t < 100000.0){
            
            if ((i%10000)==0){
                
                //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
                //                   i*dt, prob, Nryd, Nryd2, numerator);
                fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            }
        }
        else{
            
            if ((i%50000)==0){
                
    //            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\n",
    //                   i*dt, prob, Nryd, Nryd2, numerator);
                
                  fprintf(outfile, "%.6lf\t%.6e\t%.6e\t%.6e\n", t, prob, Nryd, Nryd2);
            
            }
     
        }
        
        i+=1;
        t+=dt;
    
        ///global noise
        if (deph==1){
            globalnoise();
//            eulermeruyama_global(n);
        }
        
        ///local noise
        if (deph==0){
            local_noise();
//            eulermeruyama_local(n);
        }
    }
    
    /*
    for (spin_i = 0; spin_i < Nspins+1; spin_i+=1){
        count_prob_state = 0.0;
        for (state = 0; state < Nstates; state+=1){
            if (excited_num[state]==spin_i){
                count_prob_state += (re[state]*re[state]+im[state]*im[state]);
            }
        }
        printf("%d\t%.6e\n", spin_i, count_prob_state);
        fprintf(outfile1, "%d\t%.6e\n", spin_i, count_prob_state);
    }
     */
    
    printf ("%d\n", realisation);        
}

void euller (void){
    
    int dx;
    int jx;
    int i;
    //double sum;
    //double Xtrue;
    //sum = 0.0;

    i = 0;
    while (i < notimesteps){
        
        t = i*dt;
        
        for ( state = 0; state < Nstates; state+=1){
            //printf("%.12e\t%.12e\n", re[state], im[state]);
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re[state]+= dt/8 *  im[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im[state] -= dt/4 * re[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re[state] += dt/8 * im[k[state][spin_i]];
            }
			
            re [state] += dt * Delta * im[state] * excited_num [state];
            //re [state] -= dt * noise * im[state] * excited_num [state];
            im [state] -= dt * Delta * re[state] * excited_num [state];
            //im [state] += dt * noise * re[state] * excited_num [state];
            
//            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
//                if (state & (1<<spin_i)){
//                    re [state] -= dt * localnoise[spin_i]  * im[state];
//                    im [state] += dt * localnoise[spin_i]  * re[state];
//                }
//            }
            
   
 	    //printf ("%lf\n", noise);
	    if (V0 > 0.0){ 
            
            for (ilattice = 0; ilattice < Nspins; ilattice+=1){
                
//				for (jlattice = ilattice+1; jlattice < Nspins; jlattice+=1){
//					rsquare = (ilattice - jlattice)*(ilattice - jlattice)*rad*rad;
//					rsquare*=(ilattice - jlattice)*(ilattice - jlattice)*rad*rad;
//					rsquare*=(ilattice - jlattice)*(ilattice - jlattice)*rad*rad;
//					if ( (state&(1<<ilattice))&&(state&(1<<jlattice))){
//						im [state] = im [state] - dt * re[state]*V0/rsquare;
//						re [state] = re [state] + dt * im[state]*V0/rsquare;
//					}
//				}
                

    //                CB
                  
                for (dx = -(int)cut_off; dx <= (int)cut_off; dx+=1){
                   
                    if (dx!=0){
                        jx = ilattice + dx;
                        if (jx>=0&&jx<Nspins){
                        if (jx < 0)jx+=Nspins;
                        if (jx >= Nspins)jx-=Nspins;
                        
                        jlattice = jx;
                        rsquare = dx*dx;
                           
                        //printf ("%d\t%d\t%d\n", dx, ilattice, jlattice); 
                        //if (rsquare < 1.1){
                            if ( (state&(1<<ilattice))&&(state&(1<<jlattice))){
                            im [state] = im [state] - dt * re[state]*V0/(rsquare*rsquare*rsquare*2);
                            re [state] = re [state] + dt * im[state]*V0/(rsquare*rsquare*rsquare*2);
                            }
                        }
                    }
                }
				
            }
        }

            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im [state] -=dt/8 * re[k[state][spin_i]];
            }
            
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                re [state] += dt/4 * im[k[state][spin_i]];
            }
            
            for (spin_i = 0; spin_i < Nspins; spin_i+=1){
                im [state] -= dt/8 * re[k[state][spin_i]];
            }
        }
        
        prob = 0.0;
        Nryd = 0.0;
        Nryd2 = 0.0;
		numerator = 0.0;
		denom1 = 0.0;
		denom2 = 0.0;
//        ree = 0.0;
        for (state = 0; state < Nstates; state+=1){
            prob+=(re[state]*re[state] + im[state]*im[state]);
        }
//        ree += (re[Nstates-1]*re[Nstates-1] + im[Nstates-1]*im[Nstates-1]);
        
        norm = sqrt (prob); 
        for (state = 0; state < Nstates; state+=1){
            re[state]/=norm;
            im[state]/=norm;
            Nryd +=(re[state]*re[state] + im[state]*im[state])*excited_num[state];
			Nryd2 += (re[state]*re[state] + im[state]*im[state])*excited_num[state]* excited_num[state];
			
			
			for (spin_i = 0; spin_i < Nspins; spin_i+=1){
				
				for (spin_j = spin_i + 1; spin_j < Nspins; spin_j+=1){

					//if ((i%100)==0)printf ("%d\t%d\t%d\n",state, spin_i, spin_j );
						
                    //if (spin_i!=spin_j){
							
						if (state & (1<<spin_i)){

							denom1+= re[state] * re[state];
							denom1+= im[state] * im[state];

                            	//printf ("%d\n", state);					
						}

						if (state & (1<<spin_j)){

                            denom2+= re[state] * re[state];
                            denom2+= im[state] * im[state];

                            	//printf ("%d\n", state);					
						}
							
                        if ( (state & (1<<spin_i))&&(state & (1<<spin_j))){

                            numerator += re[state] * re[state];
                            numerator += im[state] * im[state];

                            //printf ("%d\n", state);					
                        }	
				    //}
					
				}
			
			}

	    }
		
        if ((i%50)==0){
            
//            printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
//                    i*dt, prob, Nryd, Nryd2, numerator, denom1, denom2);
            //printf("%.2e\n", Nryd);
            //printf("%.2e\n", i*dt);
            fprintf(outfile, "%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                    i*dt, prob, Nryd, Nryd2, numerator, denom1, denom2);
            //fprintf(outfile, "%.12e\n", Nryd);
            //fprintf(outfile, "%.12e\n", i*dt);
            
        }


	i+=1;

	//sum += dW;
	//dW = sqrt (dt) * rand_gauss();
    //noise -= dt * noise * friction;
    //noise += friction * sqrt(Gamma) * dW;
    
	for (spin_i = 0; spin_i < Nspins; spin_i+=1){
        
	    dW = sqrt (dt) * rand_gauss();
        localnoise [spin_i] -= dt * localnoise[spin_i] * friction;
        localnoise [spin_i] += friction * sqrt(Gamma) * dW;
        //localnoise [state] -= 0.5 * (dW*dW - dt);

    }

	//Xtrue = initnoise * exp ( (-friction - 0.5 * D * D)*(i*dt) + D * sum );
	//Xtrue = D * sum;

	//printf("%.2lf\t%.6e\n", i*dt, noise[i]);
        //fprintf(outfile, "%.6lf\t%.12e\n", i*dt, noise[i]);

    }
    
//    printf("%u\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
//           realisation, prob, Nryd, Nryd2, numerator, denom1, denom2);
//    
}




