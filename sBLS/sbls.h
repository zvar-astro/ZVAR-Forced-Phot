/* C implementation of the sparse BLS algorithm
 * 
 * P. Mroz @ Caltech, Mar 18, 2021
 */

#ifndef SBLS_H_
#define SBLS_H_

#define MIN_POINTS_IN_TRANSIT 3
#define MIN_NUMBER_OF_POINTS 50
#define MAX_TRANSIT_DURATION 0.2

int prepare_data (double*, double*, double*, double*, int);
    
void heap_sort_pm (double*, int*, int);

void insertion_sort_pm (double*, int*, int);

void fold_and_sort_phases (double*, int*, double*, double, int);
    
void sparse_bls (double*, double*, int, double*, double*, double*, int);
    
#endif
