//
//  main.c
//  noisedynamics.c
//
//  Created by Wildan Abdussalam on 4/24/14.
//  Copyright (c) 2014 Wildan Abdussalam. All rights reserved.
//

#include "header_noise.h"

void steady_state_run (void){
    

    rad = 1.0;
    while (rad > 0.99){
        
        init_params();
        ///init sum
        sumNryd = 0.0;
        sumprob = 0.0;
        sumNryd2 = 0.0;
        sumnumerator = 0.0;
        sumdenom1 = 0.0;
        sumdenom2 = 0.0;
        sumree = 0.0;
        for (iter = 0; iter < Nreal; iter+=1){
    
            init_configuration();
            //printf ("%lf\t%lf\n", rand_gauss(), genrand_real2());
            counting_excited_state();
            flip_state();
            trajectory_noise();
            sumNryd+=Nryd;
            sumNryd2+=Nryd2;
            sumnumerator+=numerator;
            sumprob+=prob;
            sumdenom1+=denom1;
            sumdenom2+=denom2;
            sumree += ree;
            printf("#%u\n", iter);
        }
    
        printf("%d\n", realisation);
        aveprob = sumprob/Nreal;
        aveNryd = sumNryd/Nreal;
        aveNryd2 = sumNryd2/Nreal;
        avenumerator = sumnumerator/Nreal;
        avedenom1 = sumdenom1/Nreal;
        avedenom2 = sumdenom2/Nreal;
        averee =  sumree/Nreal;

        printf("%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
               rad, averee, aveNryd, aveNryd2, avenumerator, avedenom1, avedenom2);
        //printf("%.2e\n", Nryd);
        //printf("%.2e\n", i*dt);
//              fprintf(outfile, "%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
//                    friction, averee, aveNryd, aveNryd2, avenumerator, avedenom1, avedenom2);
        fprintf(outfile, "%.6lf\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
                        frictionbar, prob, Nryd, Nryd2, numerator, denom1, denom2);

//        fprintf(outfile1, "%.6lf\t", delta);

//           for (state =0; state < Nstates; state+=1){
//               fprintf(outfile1, "%.12e\t", basis[state]);
//           }

//           fprintf(outfile1, "\n");
        
        free_alloc();
        rad-=0.01;
    }
}

void run_normal (void){

    init_params();
//    if (V0 > 0)spectrum_symmetric(Nspins);
    init_configuration();
    counting_excited_state();
    flip_state();
    trajectory_noise();
    free_alloc();
}

int main(int argc, const char * argv[])
{
   /*
   if(argc != 8){
      printf("Incorrect number of args\n");
      exit(-1);
   }
   if (sscanf(argv[1],"%d", &Nspins) != 1){
      printf("failed to insert Nspins.\n");
      exit(-2);
   }
   if (sscanf(argv[2],"%lf", &frictionbar) != 1){
      printf("failed to insert time correlation.\n");
      exit(-3);
   }
    if (sscanf(argv[3],"%lf", &xi) != 1){
        printf("failed to insert dephasing rate.\n");
        exit(-4);
    }
    if (sscanf(argv[4],"%lf", &delta) != 1){
        printf("failed to insert laser detuning.\n");
        exit(-5);
    }
    if (sscanf(argv[5],"%lf", &V0) != 1){
        printf("failed to insert interactions.\n");
        exit(-6);
    }
    if (sscanf(argv[6],"%deph", &deph) != 1){
      printf("failed to insert noise type 0 for local 1 for global.\n");
      exit(-7);
    }
    if (sscanf(argv[7],"%d", &realisation) != 1){
        printf("failed to insert realisation.\n");
        exit(-8);
    }
    */
   
    if(argc != 11){
        printf("Incorrect number of args\n");
        exit(-1);
    }
    if (sscanf(argv[1],"%d", &Nspins) != 1){
        printf("failed to insert Nspins.\n");
        exit(-2);
    }
    if (sscanf(argv[2],"%lf", &friction) != 1){
        printf("failed to insert time correlation.\n");
        exit(-3);
    }
    if (sscanf(argv[3],"%lf", &Omega) != 1){
        printf("failed to insert dephasing rate.\n");
        exit(-4);
    }
    if (sscanf(argv[4],"%lf", &mu) != 1){
        printf("failed to insert laser detuning.\n");
        exit(-5);
    }
    if (sscanf(argv[5],"%lf", &V0) != 1){
        printf("failed to insert interactions.\n");
        exit(-6);
    }
    if (sscanf(argv[6],"%d", &deph) != 1){
        printf("failed to insert noise type 0 for local 1 for global.\n");
        exit(-7);
    }
    if (sscanf(argv[7],"%lf", &Gamma) != 1){
        printf("failed to insert Gamma.\n");
        exit(-8);
    }
    if (sscanf(argv[8],"%lf", &specwidth) != 1){
        printf("failed to insert spectral width.\n");
        exit(-9);
    }
    if (sscanf(argv[9],"%d", &typeint) != 1){
        printf("failed to insert type of interaction.\n");
        exit(-10);
    }
    if (sscanf(argv[10],"%d", &realisation) != 1){
        printf("failed to insert realisation.\n");
        exit(-11);
    }
    
//    steady_state_run();
    run_normal();
    
    return 0;
}

