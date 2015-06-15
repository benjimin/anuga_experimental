// Python - C extension module for shallow_water.py
// for the case of subgrid topgraphy
//
// Gareth Davies, GA 2014 - 2015
//
// (Based on earlier shallow water domain ext code with authors including)
// Ole Nielsen, GA 2004
// Stephen Roberts, ANU 2009
// Gareth Davies, GA 2011 - 2015
//



#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
//#include "numpy_shim.h"

// Shared code snippets
#include "util_ext.h"
#include "sw_SG_domain.h"


const double pi = 3.14159265358979;

// Trick to compute n modulo d when d is a power of 2
inline unsigned int Mod_of_power_2(unsigned int n, unsigned int d)
{
  return ( n & (d-1) );
} 


// Computational function for rotation
inline int _rotate(double *q, double n1, double n2) {
  /*Rotate the last  2 coordinates of q (q[1], q[2])
    from x,y coordinates to coordinates based on normal vector (n1, n2).

    Result is returned in array 2x1 r
    To rotate in opposite direction, call rotate with (q, n1, -n2)

    Contents of q are changed by this function */


  double q1, q2;

  // Shorthands
  q1 = q[1];  // x coordinate
  q2 = q[2];  // y coordinate

  // Rotate
  q[1] =  n1*q1 + n2*q2;
  q[2] = -n2*q1 + n1*q2;

  return 0;
}


inline double eval_subgrid_function(double value, int k, int edge,
                             int var, int inverse, struct domain *D,
                             int use_last_lookup_info,
                             int use_last_cube_root 
                             //double *lookup_table
                             ){

    // Function to evaluate the subgrid functions.
    //
    // Mostly we are given a 'value' corresponding to a stage (reference stage)
    // and need to interpolate a particular subgrid function at this stage.
    // The function is stored in a lookup table along with many others, and is
    // coded by 'var', the triangle index k, and the 'edge' integer. 
    //
    // If inverse==1, then we are given a value of the subgrid function,
    // and need to return the depth.
    //
    /////////////////////////////////////////////////////////////////////////
    // INPUT DATA
    //
    // value = value used in interpolation
    //
    // k = triangle centroid index for interpolation function
    //
    // edge = triangle edge index 0,1, or 2 (or -1 to lookup cell-integrated
    //        quantities) for interpolation function
    //
    // var = code for the interpolation function to look up 
    //       0 - wetted area or length, over volume or edge respectively
    //       1 - integral (h), over volume or edge
    //       2 - integral 1/n h^(5/3), over volume or edge
    //       3 - integral 1/n^2 h^(7/3), over edge [not supported for volume]
    //       4 - integral (h^2), over edge [not supported for volume] 
    //
    // inverse = 0 or 1
    //       If 0, then 'value' is a stage value, and we lookup var 
    //       If 1, then 'value' is a 'var' value, and we lookup stage
    //
    // D = The domain object. We need the subgrid look tables and related info
    //
    // use_last_lookup_info = 0 or 1
    //       If the last call to eval_subgrid_function had the same 'value','k','edge','inverse=0','D'
    //       then set this to 1 to avoid recalculating various values when interpolating a different 'var'
    //       Otherwise set it to 0 
    //
    // use_last_cube_root = 0 or 1
    //       If the last call to eval_subgrid_function had the same 'value','k','edge','inverse=0','D'
    //       AND it involved computation of cube-root weights for the
    //       interpolation, then set this to 1 to avoid recalculating the cube
    //       roots.
    //
    // FIXME: Are the above optimizations worthwhile? Misuse will introduce errors
    //       Otherwise set it to zero
    //

    // To interpolate we use power-law transformations for increased accuracy,
    // see comments in the code
    int value_ind ;
    int var_row, found_bound, interp_ind, lb, ub;
    double w0, w1, return_value, pow_trans, base_value;
    int counter;
   
    // Variables defining where we look-up in the lookup table
    // Can sometimes be reused, so make them static
    static int ki ;
    static double *lookup_table ; // Pointer to appropriate lookup table 
    static int starting_index, subtable_ncol, subtable_lastcol;

    // These variables hold cube roots of various numbers from table lookup
    // Make them static to avoid recalculating in some cases
    static double ct1, ct2, ct3, t1, t2, t3;

    //
    // To make the interpolation more accurate add a power-law transformation
    // to the x-data first, and then linearly interpolate
    //
    // The interpolation will be exact for a power law of the given power, and
    // good so long as the interpolation curve grows approximately like this
    //
    // Define: h = (reference_stage-min_possible_reference_stage)
    static double var_power_trans[6] = { 1.0, // stage curve grows as h^1
                                         1.0, // area curve grows as ~ h^1
                                         1.0, // volume curve grows as ~ h^1
                                         5./3., // int 1/n h^(5/3) curve grows as ~ h^(5/3)
                                         7./3., // int 1/n^2 h^(7/3) curve grows as ~ h^(7/3)
                                         2.   };// int h^2 curve grows as ~ h^(2)

    // Store 1/var_power_trans, to avoid repeated division
    static double var_power_trans_inv[6] = { 1.0, // stage curve grows as h^1
                                             1.0, // area curve grows as ~ h^1
                                             1.0, // volume curve grows as ~ h^1
                                             3./5., // int 1/n h^(5/3) curve grows as ~ h^(5/3)
                                             3./7., // int 1/n^2 h^(7/3) curve grows as ~ h^(7/3)
                                             1./2.};// int h^2 curve grows as ~ h^(2)


    //////////////////////////////////////////////////////////////////////////////
    //
    // CHECK INPUT ARGS (consider removing)
    //
    /////////////////////////////////////////////////////////////////////////////
    
    // // Check that var has a valid value 
    // if(var < 0 | var >= 5){
    //     report_python_error(AT, "Invalid lookup var code: Must be in 0,1,2,3,4");
    //     return -1.0;
    // }

    // // Check that if we are interpolating from a volume curve, then var <= 2
    // // (Since we only store area/volume/integral(1/n h^(5/3)))
    // if ( edge == -1){
    //     if ( var > 2 ){
    //         report_python_error(AT, "Invalid lookup var code for centroid: Must be 0 or 1 or 2");
    //         return -1.0;
    //     }
    // }else{
    //     if ( edge > 2 | edge < -1){
    //         report_python_error(AT, "Invalid edge code: Must be -1 (for the volume) or 0 or 1 or 2 (for corresponding edges)");
    //         return -1.0;
    //     }
    // }

    // // Check that value is not NAN
    // if(value != value){
    //     printf("var : %d, inverse: %d \n", var, inverse);
    //     report_python_error(AT, "NAN lookup value");
    //     return -1.0;
    // }

    //////////////////////////////////////////////////////////////////////////////////
    //
    // EXTRACT REQUIRED DATA
    //
    //////////////////////////////////////////////////////////////////////////////////

    if(use_last_lookup_info == 0){
        // Make a pointer to the correct subgrid_table + other important constants
        // These variables are static, so we can skip this if use_last_lookup_info==1

        if (edge == -1){
            //lookup_table = D->subgrid_centroid_table ;
            //starting_index = D->subgrid_centroid_starting_index[k] ;

            // Pointer to start of the part of the array we can use
            lookup_table = &(D->subgrid_centroid_table[ D->subgrid_centroid_starting_index[k] ]);
            starting_index = 0;
            subtable_ncol = D->subgrid_centroid_i_ncol[k] ;
            subtable_lastcol = D->subgrid_centroid_last_lookup_col[k] ;
        }else{
            ki = 3*k+edge; // Linear index for edge 
            //lookup_table = D->subgrid_edge_table ;
            //starting_index = D->subgrid_edge_starting_index[ki] ;

            // Pointer to start of the part of the array we can use
            lookup_table = &(D->subgrid_edge_table[ D->subgrid_edge_starting_index[ki] ]) ;
            starting_index = 0;
            subtable_ncol = D->subgrid_edge_i_ncol[ki] ;
            subtable_lastcol = D->subgrid_edge_last_lookup_col[ki] ;
        }
    }
   

    //
    // Figure out which row to look up 
    //
    // area curve in row 1, 
    // integral(h) curve in row 2, 
    // integral(1/n h^(5/3)) curve in row 3
    // integral(1/n^2 h^(7/3)) curve in row 4
    // integral( h^(2)) curve in row 5
    var_row = var + 1;


    //////////////////////////////////////////////////////////////////////////
    // 
    // Find lower/upper bounds on indices containing the variables we want to
    // use to look-up, and find column index just below 'value' ( named 'value_ind'
    // below)
    //////////////////////////////////////////////////////////////////////////

    // Get lower/ upper possible bounds on the index we are searching for 
    if(inverse==0){
        // We are in the top row
        lb = starting_index;
    }else{
        // We are in another row
        lb = starting_index + var_row*subtable_ncol;
    }
    ub = lb + (subtable_ncol - 1);

    //// Debugging statement    
    // printf(" k: %d, edge: %d , inverse: %d, var: %d, si: %d, snc: %d, slc: %d, v_lb:, %e, v_ub: %e, value: %e \n", 
    //        k, edge, inverse, var, starting_index, subtable_ncol, subtable_lastcol, lookup_table[lb], lookup_table[ub], value);

    if(use_last_lookup_info == 0){ 
        // Get the 'value_ind' = lower index defining the 2 table values
        // used for interpolation

        // Check that value >= minimum allowed value 
        if(value <=  lookup_table[lb]){
            value_ind = lb;
        }else{
            // First guess of the value_ind to start search
            // This should give the lookup table column just below value
            value_ind = lb + subtable_lastcol;

            // Loop to find the correct value_ind
            found_bound = 0;
            while(found_bound == 0){

                if(lookup_table[value_ind] <= value && 
                   lookup_table[value_ind+1] >= value){
                    // The current value_ind can be used to interpolate. It is
                    // quick if we get here soon
                    found_bound = 1;
                }else{

                    if(lookup_table[value_ind] > value){
                        value_ind -= 1;
                    }else if(lookup_table[value_ind] < value){
                        value_ind += 1;
                        // If we get below here, then
                        // lookup_table[value_ind]==value
                        // So we should have satisfied the previous lookup
                    }else{
                        // Something bad has happened....
                        if(value == value){
                            printf("no NaN error\n");
                        }else{
                            printf("NaN error\n");
                        }
                        report_python_error(AT, "Interpolation search error");
                        return -1.0;
                    }

                    // Ensure we have not exceeded the upper bound of our table
                    if(value_ind >= ub){
                         printf("%d, %d, %d, %d \n", value_ind, lb, ub, var);
                         printf("... %e, %e, %e \n", value, 
                                lookup_table[lb],lookup_table[ub]);
                         report_python_error(AT, 
                             "Interpolation searched outside range of table values");
                         return -1.0;
                    }
                }
            }
        }

        // Update the last lower ind
        subtable_lastcol = value_ind - lb;
        if (edge == -1){
            D->subgrid_centroid_last_lookup_col[k] = subtable_lastcol;
        }else{
            D->subgrid_edge_last_lookup_col[ki] = subtable_lastcol;
        }

    }else{
        // In this case we just need to update value_ind
        value_ind = subtable_lastcol + lb;
    }    

    //////////////////////////////////////////////////////////////////////
    //
    // Compute interpolation
    //
    //////////////////////////////////////////////////////////////////////

    // 
    // Suppose [x0,y0], [x1,y1] are the tabulated values
    // and we interpolate at x [x0<=x<=x1] 
    // 
    // Then we interpolate as:
    // y = (y1*w0 + y0*w1)/(w0+w1)
    // where:
    // w0 = ( t1 - t2)
    // w1 = ( t3 - t1)
    // and 
    // t1 = (x-base_value)^pow_trans
    // t2 = (x0-base_value)^pow_trans
    // t3 = (x1 - base_value)^pow_trans
    // base_value = smallest allowed value of 'x' [bottom of the lookup table]
    //
    // This would be exact for a power law of the form
    //    y = (x-base_value)^(pow_trans)
    // If the latter relation is approximately true then it will be ok
    //
    // Note if inverse==1, we have to invert the pow_trans

    if(inverse == 0){
        // h_ref is given, need to get another variable
        interp_ind= starting_index + subtable_lastcol + var_row*subtable_ncol;
        pow_trans = var_power_trans[var_row];
    }else{
        // Need to back-calculate h_ref
        interp_ind = starting_index + subtable_lastcol;
        // Invert the power transformation
        pow_trans = var_power_trans_inv[var_row]; 
    }

    // base_value = smallest allowed value of h_ref [or the other variable]
    base_value = lookup_table[ value_ind - subtable_lastcol ];

    // Compute weights for interpolation, making quick exits if possible
    if(use_last_lookup_info == 0){
        // x - base_value -- power transforms are performed later
        t1 = value - base_value;
        t2 = lookup_table[value_ind] - base_value;
        t3 = lookup_table[value_ind+1] - base_value;
    }

    // Try various quick-exit options
    if(t1 <= 0.){
        // Quick exit
        return_value = lookup_table[interp_ind];
        return return_value;

    }else if(t1 == t2){
        // Quick exit -- value direct from table
        return_value = lookup_table[interp_ind];
        return return_value;

    }else if(t1 == t3){
        // Quick exit -- value direct table
        return_value = lookup_table[interp_ind+1];
        return return_value;

    }

    // If we got to here, 
    // value is between 2 table entries

    // Try to avoid using 'pow' if possible, it really slows things down
    if(pow_trans != 1.0){
        if(pow_trans != 2.0){
            // Efficient computation of powers 5/3 and 7/3

            // Compute cube roots efficiently 
            // (use old ones if user says)
            if(use_last_cube_root == 0){
                ct1 = cbrt(t1);
                ct2 = cbrt(t2);
                ct3 = cbrt(t3);
            }

            
            if( fabs(pow_trans*3.-5.) < 1.0e-06 ){
                // pow_trans = 5./3.
                // Compute cube roots efficiently 
                // (use old ones if user says)
                if(use_last_cube_root == 0){
                    ct1 = cbrt(t1);
                    ct2 = cbrt(t2);
                    ct3 = cbrt(t3);
                }
                t1 *= ct1*ct1; 
                t2 *= ct2*ct2;
                t3 *= ct3*ct3;

            }else if( fabs(pow_trans*3.-7.) < 1.0e-06 ){

                // pow_trans = 7./3.
                // Compute cube roots efficiently 
                // (use old ones if user says)
                if(use_last_cube_root == 0){
                    ct1 = cbrt(t1);
                    ct2 = cbrt(t2);
                    ct3 = cbrt(t3);
                }
                t1 *= t1*ct1;
                t2 *= t2*ct2;
                t3 *= t3*ct3;
            }else if( fabs(pow_trans*2.0 - 1.0) < 1.0e-06 ){
                // pow_trans = 0.5
                // Currently not used but support in-case of change
                t1 = sqrt(t1);
                t2 = sqrt(t2);
                t3 = sqrt(t3);

            }else{
                // Could support 3/5, 3/7 -- but not used at the moment
                report_python_error(AT,"power transformation value not supported");

            }

        }else{
            //pow_trans == 2.0
            t1 *= t1;
            t2 *= t2;
            t3 *= t3;

        }
    }

    // Weights from transformed 't' values
    w0 = (t1 - t2);
    w1 = (t3 - t1);

    // // Only positive weights are reasonable
    if(w0<0.0 | w1 < 0.0){
        printf("%e, %e, %e, %e,%e,%e \n", t1, t2, t3, w0, w1, base_value);
        report_python_error(AT,"Interpolation Error, negative weights");
        return -1.0;

    }

    // if(w0>0.0| w1>0.0){
    //     // Compute interpolated value here
    return_value = (w1*lookup_table[interp_ind] + 
                    w0*lookup_table[interp_ind+1])/(w0+w1);
    // }else{
    //     printf("%d, %e, %e, %e, %e,%e,%e,%e \n", inverse, t1, t2, t3, w0, w1, base_value, value);
    //     report_python_error(AT,"Interpolation Error, all zero weights");
    //     return -1.0;
    // }
   

    // // Check for NAN return value
    if(return_value != return_value){
        printf("var: %d, inverse: %d, value: %e, lb: %d, ub: %d, value_ind: %d, lower_val: %e, upper_val: %e, w0: %e, w1:%e, bv: %e\n",
                var, inverse, value, lb, ub, value_ind, 
                lookup_table[value_ind], lookup_table[value_ind+1], 
                w0, w1, base_value);
        report_python_error(AT, "Interpolated value is NAN");
        return -1.0;
    }

    return return_value;

}


//// Function to obtain speed from momentum and depth.
//// This is used by flux functions
//// Input parameters uh and h may be modified by this function.
//// Tried to inline, but no speedup was achieved 27th May 2009 (Ole)
////static inline double _compute_speed(double *uh, 
//double _compute_speed(double *uh, 
//		      double *h, 
//		      double epsilon, 
//		      double h0,
//		      double limiting_threshold) {
//  
//  double u;
//
//  if (*h < limiting_threshold) {   
//    // Apply limiting of speeds according to the ANUGA manual
//    if (*h < epsilon) {
//      //*h = max(0.0,*h);  // Could have been negative
//      u = 0.0;
//    } else {
//      u = *uh/(*h + h0/ *h);    
//    }
//  
//
//    // Adjust momentum to be consistent with speed
//    *uh = u * *h;
//  } else {
//    // We are in deep water - no need for limiting
//    u = *uh/ *h;
//  }
//  
//  return u;
//}
//
//// minmod limiter
//int _minmod(double a, double b){
//    // Compute minmod
//
//    if(sign(a)!=sign(b)){
//        return 0.0;
//    }else{
//        return min(fabs(a), fabs(b))*sign(a);
//    }
//
//
//}

// Innermost flux function (using stage w=z+h)
inline int _flux_function_central(double* ql, double* qr,
                           double* edgeflux, 
                           double n1, double n2,
                           double epsilon, 
                           double g,
                           double* max_speed, 
                           double* pressure_flux, 
                           double length, 
                           double u_left, double u_right, 
                           double hl, double hr) 
{

  /*Compute fluxes between volumes for the shallow water wave equation
    cast in terms of the 'stage', w = h+z using
    the 'central scheme' as described in

    Kurganov, Noelle, Petrova. 'Semidiscrete Central-Upwind Schemes For
    Hyperbolic Conservation Laws and Hamilton-Jacobi Equations'.
    Siam J. Sci. Comput. Vol. 23, No. 3, pp. 707-740.

    The implemented formula is given in equation (3.15) on page 714

    FIXME: Several variables in this interface are no longer used, clean up
  */

  int i, method;

  double w_left,  uh_left, vh_left; //, u_left;
  double w_right, uh_right, vh_right;//, u_right;
  double s_min, s_max, soundspeed_left, soundspeed_right;
  double denom, inverse_denominator;
  double uint, t1, t2, t3, min_speed, tmp;
  // Workspace (allocate once, use many)
  static double q_left_rotated[3], q_right_rotated[3], flux_right[3], flux_left[3];

  //method=1;
  //if(method==1){
  //    // Compute left / right soundspeeds using
  //    // appropriate 'mean-values' of integrated quantities
  //    if(fabs(ql[0]/length)>1.0e-12){
  //        // Velocity ~ uh_integral / h_integral
  //        u_left = ql[1]/(ql[0]);
  //    }else{
  //        u_left = 0.;
  //    }

  //    if(fabs(qr[0]/length)>1.0e-12){
  //        // Velocity ~ uh_integral / h_integral
  //        u_right = qr[1]/(qr[0]);
  //    }else{
  //        u_right = 0;
  //    }
  //}else if(method==2){

  //    if(fabs(ql[1])>1.0e-12){
  //       // u^2h integral / uh integral
  //       u_left = ql[3]/ql[1];
  //    }else{
  //       u_left = 0.;
  //    }
  //    
  //    if(fabs(qr[1])>1.0e-12){
  //       // u^2h integral / uh integral
  //       u_right = qr[3]/qr[1];
  //    }else{
  //       u_right = 0.;
  //    }
  //}else{
  //   report_python_error(AT, "UNKNOWN METHOD");
  //   return -1.0;
  //}

  // Quick Exit
  if((hl == 0.) && (hr == 0.)){
    // Both wave speeds are very small
    memset(edgeflux, 0, 3*sizeof(double));

    *max_speed = 0.0;
    *pressure_flux = 0.0;

    return 0;

  }

  // Maximal and minimal gravity wave speeds
  // Use edge-mean-depth = h-integral/edgelength
  //soundspeed_left  = sqrt(g*ql[0]/length);
  //soundspeed_right = sqrt(g*qr[0]/length);  
  soundspeed_left  = sqrt(g*hl);
  soundspeed_right = sqrt(g*hr);  

  //if(ql[0]>1.0e-12){
  //    soundspeed_left  = sqrt(g*ql[6]/ql[0]);
  //}else{
  //    soundspeed_left = 0.;
  //}

  //soundspeed_left = max(soundspeed_left, sqrt(g*ql[0]/length));

  //if(qr[0]>1.0e-12){
  //    soundspeed_right = sqrt(g*qr[6]/qr[0]);  
  //}else{
  //    soundspeed_right = 0.;  
  //}
  ////  
  //soundspeed_right = max(soundspeed_right, sqrt(g*qr[0]/length));
  


  // WAVESPEEDS

  s_max = max(u_left + soundspeed_left, u_right + soundspeed_right);
  if (s_max < 0.0) 
  {
    s_max = 0.0;
  }
  

  //if( hc < 1.0e-03){
  //  s_max = 0.0;
  //}


  s_min = min(u_left - soundspeed_left, u_right - soundspeed_right);
  if (s_min > 0.0)
  {
    s_min = 0.0;
  }
 
  // CHECK FOR NAN
  if( (s_min != s_min) || (s_max != s_max) ){
    printf("%e, %e, %e, %e", u_left, u_right, soundspeed_left, soundspeed_right);
    report_python_error(AT, "NAN WAVESPEEDS");
    return -1;
  }
 
  //if( hc_n < 1.0e-03){
  //  s_min = 0.0;
  //}

  method = 2; 
  if(method==1){ 
      // Flux formulas without gravity terms
      // integral [uh] *n1 + integral [vh] *n2
      // u_x * h
      flux_left[0] = ql[1];
      //flux_left[0] = u_left*ql[0]; //ql[1];
      // u_x^2 * h
      flux_left[1] = u_left*ql[1]; //ql[3]; 
      // u_x* v_x * h
      flux_left[2] = u_left*ql[2]; //ql[4];
      
      // u_x * h
      flux_right[0] = qr[1];
      //flux_right[0] = u_right*qr[0]; //qr[1];
      // u_x^2 * h
      flux_right[1] = u_right*qr[1]; //qr[3]; 
      // u_x* v_x * h
      flux_right[2] = u_right*qr[2]; //qr[4];
  }else if(method==2){
      // Flux formulas without gravity terms
      // integral [uh] *n1 + integral [vh] *n2
      // u_x * h
      flux_left[0] = ql[1];
      // u_x^2 * h
      flux_left[1] = ql[3]; 
      // u_x* v_x * h
      flux_left[2] = ql[4];
      
      // u_x * h
      flux_right[0] = qr[1];
      // u_x^2 * h
      flux_right[1] = qr[3]; 
      // u_x* v_x * h
      flux_right[2] = qr[4];

  }else{
     report_python_error(AT, "UNKNOWN METHOD");
     return -1.0;

  }
  

  // Flux computation
  denom = s_max - s_min;
  if (denom < epsilon) 
  { 
    // Both wave speeds are very small
    memset(edgeflux, 0, 3*sizeof(double));

    *max_speed = 0.0;
    *pressure_flux = 0.0;
    //*pressure_flux = 0.5*g*0.5*(h_left*h_left+h_right*h_right);
  } 
  else 
  {
    // Maximal wavespeed
    *max_speed = max(s_max, -s_min);
  
    inverse_denominator = 1.0/max(denom,1.0e-100);
    for (i = 0; i < 3; i++) 
    {
      edgeflux[i] = s_max*flux_left[i] - s_min*flux_right[i];

      // Standard smoothing term
      // edgeflux[i] += 1.0*(s_max*s_min)*(q_right_rotated[i] - q_left_rotated[i]);
      // Smoothing by stage alone can cause high velocities / slow draining for nearly dry cells
      //if(i==0) edgeflux[0] += (s_max*s_min)*(qr[0] - ql[0]);
      //if(i==1) edgeflux[1] += (s_max*s_min)*(qr[1] - ql[1]);
      //if(i==2) edgeflux[2] += (s_max*s_min)*(qr[2] - ql[2]);
      edgeflux[i] += (s_max*s_min)*(qr[i] - ql[i]);
      edgeflux[i] *= inverse_denominator;

    }

   // Rotate edgeflux back to x,y coordinate system
   _rotate(edgeflux,n1, -n2);

    // Separate pressure flux, so we can apply different wet-dry hacks to it
    *pressure_flux = 0.5*g*(s_max*ql[6] -s_min*qr[6])*inverse_denominator;
    

  }
  
  return 0;
}

////////////////////////////////////////////////////////////////

inline int _compute_flux_update_frequency(struct domain *D, double timestep){
    //
    // Compute the 'flux_update_frequency' for each edge.
    //
    // This determines how regularly we need
    // to update the flux computation (not every timestep)
    //
    // Allowed values are 1,2,4,8,... max_flux_update_frequency.
    //
    // For example, an edge with flux_update_frequency = 4 would
    // only have the flux updated every 4 timesteps
    //
    //
    // Local variables
    int k, i, k3, ki, m, n, nm, ii, j, ii2;
    long fuf;
    double notSoFast=1.0;
    static int cyclic_number_of_steps=-1;

    // QUICK EXIT
    if(D->max_flux_update_frequency==1){
        return 0;
    }
    
    // Count the steps
    cyclic_number_of_steps++;
    if(cyclic_number_of_steps==D->max_flux_update_frequency){
        // The flux was just updated in every cell
        cyclic_number_of_steps=0;
    }


    // PART 1: ONLY OCCURS FOLLOWING FLUX UPDATE

    for ( k = 0; k < D->number_of_elements; k++){
        for ( i = 0; i < 3; i++){
            ki = k*3 + i;
            if((Mod_of_power_2(cyclic_number_of_steps, D->flux_update_frequency[ki])==0)){
                // The flux was just updated, along with the edge_timestep
                // So we better recompute the flux_update_frequency
                n=D->neighbours[ki];
                if(n>=0){
                    m = D->neighbour_edges[ki];
                    nm = n * 3 + m; // Linear index (triangle n, edge m)
                }

                // Check if we have already done this edge
                // (Multiply already_computed_flux by -1 on the first update,
                // and again on the 2nd)
                if(D->already_computed_flux[ki] > 0 ){
                    // We have not fixed this flux value yet
                    D->already_computed_flux[ki] *=-1;
                    if(n>=0){
                        D->already_computed_flux[nm] *=-1;
                    }
                }else{
                    // We have fixed this flux value already 
                    D->already_computed_flux[ki] *=-1;
                    if(n>=0){
                        D->already_computed_flux[nm] *=-1;
                    }
                    continue;
                }

                // Basically int( edge_ki_timestep/timestep ) with upper limit + tweaks
                // notSoFast is ideally = 1.0, but in practice values < 1.0 can enhance stability
                // NOTE: edge_timestep[ki]/timestep can be very large [so int overflows].
                //       Do not pull the (int) inside the min term
                fuf = (int)min((D->edge_timestep[ki]/timestep)*notSoFast,D->max_flux_update_frequency*1.);
                // Account for neighbour
                if(n>=0){
                    fuf = min( (int)min(D->edge_timestep[nm]/timestep*notSoFast, D->max_flux_update_frequency*1.), fuf);
                }

                // Deal with notSoFast<1.0
                if(fuf<1){
                    fuf=1;
                }
                // Deal with large fuf 
                if(fuf> D->max_flux_update_frequency){
                    fuf = D->max_flux_update_frequency;
                }
                //// Deal with intermediate cases
                ii=2;
                while(ii< D->max_flux_update_frequency){
                    // Set it to 1,2,4, 8, ...
                    ii2=2*ii;
                    if((fuf>ii) && (fuf<ii2)){
                        fuf = ii;
                        continue;
                    }
                    ii=ii2;
                }

                // Set the values
                D->flux_update_frequency[ki]=fuf;
                if(n>=0){
                    D->flux_update_frequency[nm]= fuf;
                }
                
            }
        }
    } 

    //// PART 2 -- occcurs every timestep

    // At this point, both edges have the same flux_update_frequency.
    // Next, ensure that flux_update_frequency varies within a constant over each triangle
    // Experiences suggests this is numerically important
    // (But, it can result in the same edge having different flux_update_freq)
    for( k=0; k< D->number_of_elements; k++){
        k3=3*k;
        ii = 1*min(D->flux_update_frequency[k3],
                 min(D->flux_update_frequency[k3+1],
                     D->flux_update_frequency[k3+2]));
        
        D->flux_update_frequency[k3]=min(ii, D->flux_update_frequency[k3]);
        D->flux_update_frequency[k3+1]=min(ii, D->flux_update_frequency[k3+1]);
        D->flux_update_frequency[k3+2]=min(ii,D->flux_update_frequency[k3+2]);
            
    }
            
    // Now enforce the same flux_update_frequency on each edge
    // (Could have been broken above when we limited the variation on each triangle)
    // This seems to have nice behaviour. Notice how an edge
    // with a large flux_update_frequency, near an edge with a small flux_update_frequency,
    // will have its flux_update_frequency updated after a few timesteps (i.e. before max_flux_update_frequency timesteps)
    // OTOH, could this cause oscillations in flux_update_frequency?
    for( k = 0; k < D->number_of_elements; k++){
        D->update_extrapolation[k]=0;
        for( i = 0; i< 3; i++){
            ki=3*k+i;
            // Account for neighbour
            n=D->neighbours[ki];
            if(n>=0){
                m = D->neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
                D->flux_update_frequency[ki]=min(D->flux_update_frequency[ki], D->flux_update_frequency[nm]);
            }
            // Do we need to update the extrapolation?
            // (We do if the next flux computation will actually compute a flux!)
            if(Mod_of_power_2((cyclic_number_of_steps+1),D->flux_update_frequency[ki])==0){
                D->update_next_flux[ki]=1;
                D->update_extrapolation[k]=1;
            }else{
                D->update_next_flux[ki]=0;
            }
        }
    }

    // Check whether the timestep can be increased in the next compute_fluxes call
    if(cyclic_number_of_steps+1==D->max_flux_update_frequency){
        // All fluxes will be updated on the next timestep
        // We also allow the timestep to increase then
        D->allow_timestep_increase[0]=1;
    }else{
        D->allow_timestep_increase[0]=0;
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////////

double adjust_edgeflux_with_weir(double *edgeflux,
                                 double h_left, double h_right, 
                                 double g, double weir_height,
                                 double Qfactor, 
                                 double s1, double s2, 
                                 double h1, double h2,
                                 double *max_speed_local
                                ){
    // Adjust the edgeflux to agree with a weir relation [including
    // subergence], but smoothly vary to shallow water solution when
    // the flow over the weir is much deeper than the weir, or the
    // upstream/downstream water elevations are too similar
    double rw, rw2; // 'Raw' weir fluxes
    double rwRat, hdRat,hdWrRat, scaleFlux, scaleFluxS, minhd, maxhd;
    double w1, w2; // Weights for averaging
    double newFlux;
    double twothirds = (2.0/3.0);
    //
    // Following constants control the 'blending' with the shallow water solution
    // They are now user-defined
    //
    //double s1=0.9; // At this submergence ratio, begin blending with shallow water solution
    //double s2=0.95; // At this submergence ratio, completely use shallow water solution
    //double h1=1.0; // At this (tailwater height above weir) / (weir height) ratio, begin blending with shallow water solution
    //double h2=1.5; // At this (tailwater height above weir) / (weir height) ratio, completely use the shallow water solution

    minhd = min(h_left, h_right);
    maxhd = max(h_left, h_right);
    // 'Raw' weir discharge = Qfactor*2/3*H*(2/3*g*H)**0.5
    rw = Qfactor * twothirds * maxhd * sqrt(twothirds * g * maxhd);
    // Factor for villemonte correction
    rw2 = Qfactor * twothirds * minhd * sqrt(twothirds * g * minhd);
    // Useful ratios
    rwRat = rw2 / max(rw, 1.0e-100);
    hdRat = minhd / max(maxhd, 1.0e-100);

    // (tailwater height above weir)/weir_height ratio
    hdWrRat = minhd / max(weir_height, 1.0e-100);
        
    // Villemonte (1947) corrected weir flow with submergence
    // Q = Q1*(1-Q2/Q1)**0.385
    rw = rw*pow(1.0 - rwRat, 0.385);

    if(h_right > h_left){
        rw *= -1.0;
    }

    if( (hdRat<s2) & (hdWrRat< h2) ){
        // Rescale the edge fluxes so that the mass flux = desired flux
        // Linearly shift to shallow water solution between hdRat = s1 and s2  
        // and between hdWrRat = h1 and h2

        //
        // WEIGHT WITH RAW SHALLOW WATER FLUX BELOW
        // This ensures that as the weir gets very submerged, the 
        // standard shallow water equations smoothly take over
        //

        // Weighted average constants to transition to shallow water eqn flow
        w1 = min( max(hdRat-s1, 0.) / (s2-s1), 1.0);
        
        // Adjust again when the head is too deep relative to the weir height
        w2 = min( max(hdWrRat-h1,0.) / (h2-h1), 1.0);

        newFlux = (rw*(1.0-w1)+w1*edgeflux[0])*(1.0-w2) + w2*edgeflux[0];

        if(fabs(edgeflux[0]) > 1.0e-100){
            scaleFlux = newFlux/edgeflux[0];
        }else{
            scaleFlux = 0.;
        }

        scaleFlux = max(scaleFlux, 0.);

        edgeflux[0] = newFlux;

        // FIXME: Do this in a cleaner way
        // IDEA: Compute momentum flux implied by weir relations, and use
        //       those in a weighted average (rather than the rescaling trick here)
        // If we allow the scaling to momentum to be unbounded,
        // velocity spikes can arise for very-shallow-flooded walls
        edgeflux[1] *= min(scaleFlux, 10.);
        edgeflux[2] *= min(scaleFlux, 10.);
    }

    // Adjust the max speed
    if (fabs(edgeflux[0]) > 0.){
        *max_speed_local = sqrt(g*(maxhd+weir_height)) + abs(edgeflux[0]/(maxhd + 1.0e-12));
    }
    //*max_speed_local += abs(edgeflux[0])/(maxhd+1.0e-100);
    //*max_speed_local *= max(scaleFlux, 1.0);

    return 0;
}


//////////////////////////////////////////////////////////////////////


 
inline int _get_rotated_subgrid_edge_quantities(
        double * ql, int ck, int ei, int eki, struct domain *D, 
        double stage_left, double stage_le, double n1, double n2
        ){
    //
    // Function to get important subgrid edge quantities
    // in the rotated coordinates system defined by the edge normal
    //
    // Approach:
    // Fills out a vector ql of length 8 with the quantities we need.
    //
    // First Step:
    //
    // Rotate (alphax, alphay) to (u_scale, v_scale) where u_scale is normal to the
    // edge, and v_scale is parallel to it [as this coordinate system is used by
    // the flux function]
    //
    // u_scale = n1 * alphax_edge_values + n2 * alphay_edge_values
    // v_scale = - n2 * alphax_edge_values + n1 * alphay_edge_values
    //
    // Second Step: Get the edge variables we need, already rotated
    // 
    // ql[0] = integral h
    // ql[1] = integral(1/n h^(5/3)) * u_scale         = integrated version of uh
    // ql[2] = integral(1/n h^(5/3)) * v_scale         = integrated version of vh
    // ql[3] = integral(1/n h^(7/3)) * u_scale^2       = integrated version of u^2h
    // ql[4] = integral(1/n h^(7/3)) * u_scale*v_scale = integrated version of uvh
    // ql[5] = integral(1/n h^(7/3)) * v_scale^2       = integrated version of v^2h
    // ql[6] = integral(h^2)                           = integrated version of h^2
    // ql[7] = 0             (was useful in another context)
    // ql[8] = edge_length

    // Local vars
    int inverse, var, use_last_lookup_info, use_last_cube_root, ki, i;
    double h_5on3n_integral, h_7on3n2_integral,u_scale, v_scale;
    //double* lookup_table;

    // Check that interpolation x values are not NAN
    if( (stage_left != stage_left) || (stage_le != stage_le) ){
        printf("%e, %e \n", stage_left, stage_le);
        report_python_error(AT, "NAN stage_left or stage_le");
        return -1;
    }

    // Rotate the alphax, alphay terms. This allows us to easily compute the
    // rotated quantities
    u_scale =  n1*D->alphax_edge_values[eki] + n2*D->alphay_edge_values[eki];
    v_scale = -n2*D->alphax_edge_values[eki] + n1*D->alphay_edge_values[eki];
 
    inverse = 0;

    ki = 3*ck + ei ;

    // h integral
    var = 1;
    use_last_lookup_info = 0;
    use_last_cube_root = 0;
    //ql[0] = eval_subgrid_function(stage_left, ck, ei, var, inverse, D);
    ql[0] = eval_subgrid_function(stage_le, ck, ei, var, inverse, D, 
                                  use_last_lookup_info, use_last_cube_root) ;//,
                                  //lookup_table);

    // If h integral is zero, then everything is zero
    if(ql[0]==0.){
        // Quick exit
        ql[1] = 0.;
        ql[2] = 0.;
        ql[3] = 0.;
        ql[4] = 0.;
        ql[5] = 0.;
        ql[6] = 0.;
        ql[7] = 0.;
    }else{

        // 1/n h^(5/3) integral
        var = 2;
        //h_5on3n_integral = eval_subgrid_function(stage_left, ck, ei, var, inverse, D);
        use_last_lookup_info = 1; // Can reuse the lookup indices
        h_5on3n_integral = eval_subgrid_function(stage_le, ck, ei, var, inverse, 
                                                D, use_last_lookup_info, use_last_cube_root); //, lookup_table);
        // uh integral
        ql[1]=h_5on3n_integral*u_scale;
        // vh integral
        ql[2]=h_5on3n_integral*v_scale;

        // 1/n^2 h^(7/3) integral
        // Terms u2h, v2h, uvh
        var = 3;
        //h_7on3n2_integral = eval_subgrid_function(stage_left, ck, ei, var, inverse, D);
        use_last_cube_root = 1; // reuse cube-roots in interpolation weights
        h_7on3n2_integral = eval_subgrid_function(stage_le, ck, ei, var, inverse, D, 
                                                use_last_lookup_info, use_last_cube_root) ; //, lookup_table);
        use_last_cube_root = 0; // set back for safety

        // u^2h integral
        ql[3]=h_7on3n2_integral*u_scale*u_scale;
        // uvh integral
        ql[4]=h_7on3n2_integral*u_scale*v_scale;
        // v^2h integral
        ql[5]=h_7on3n2_integral*v_scale*v_scale;

        // h^2 integral
        var = 4;
        //ql[6] = eval_subgrid_function(stage_left, ck, ei, var, inverse, D);
        ql[6] = eval_subgrid_function(stage_le, ck, ei, var, inverse, D, 
                                      use_last_lookup_info, use_last_cube_root) ; //, lookup_table);
        // h^2 integral value
        //ql[7] = eval_subgrid_function(stage_le, ck, ei, var, inverse, D);
        //ql[7] = eval_subgrid_function(D->stage_centroid_values[ck], 
        //                             ck, ei, var, inverse, D, use_last_lookup_info);
        ql[7] = 0.;
    }

    // Wet length
    var = 0;
    ql[8] = eval_subgrid_function(stage_le, ck, ei, var, inverse, D, 
                                use_last_lookup_info, use_last_cube_root) ; //, lookup_table);

    return 0;
}

/////////////////////////////////////////////////////////////////////
//
// Computational function for flux computation
inline double _compute_fluxes_central(struct domain *D, double timestep){

    // Local variables
    double max_speed_local, length, inv_area, zl, zr;
    double h_left, h_right, z_half ;  // For andusse scheme
    // FIXME: limiting_threshold is not used for DE
    double limiting_threshold = 10*D->H0;
    //
    int k, i, m, n,j, ii;
    int ki,k3, nm = 0, ki2,ki3, nm3, ck, ei,eki; // Index shorthands
    // Workspace (making them static actually made function slightly slower (Ole))
    double ql[9], qr[9], edgeflux[3]; // Work array for summing up fluxes
    int inverse, var, useint;
    double bedslope_work;
    static double local_timestep;
    long RiverWall_count, substep_count;
    double hle, hre, sle, sre, zc, zc_n, Qfactor, s1, s2, h1, h2; 
    double stage_edge_lim, outgoing_mass_edges, pressure_flux, hc, hc_n, tmp, tmp2;
    double h_left_tmp, h_right_tmp,sfx,sfy, max_stage, cross_section_area_l, cross_section_area_r;
    //double* lookup_table;
    static long call = 1; // Static local variable flagging already computed flux
    double speed_max_last, vol, weir_height, hl, hr, u_left, u_right, alpha_norm;

    call++; // Flag 'id' of flux calculation for this timestep

    // Set explicit_update to zero for all conserved_quantities.
    // This assumes compute_fluxes called before forcing terms
    memset((char*) D->vol_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->u_vol_explicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->v_vol_explicit_update, 0, D->number_of_elements * sizeof (double));

    memset((char*) D->u_vol_semi_implicit_update, 0, D->number_of_elements * sizeof (double));
    memset((char*) D->v_vol_semi_implicit_update, 0, D->number_of_elements * sizeof (double));


    // Counter for riverwall edges
    RiverWall_count = 0;
    // Which substep of the timestepping method are we on?
    substep_count = (call-2)%D->timestep_fluxcalls;

    // Fluxes are not updated every timestep,
    // but all fluxes ARE updated when the following condition holds
    if(D->allow_timestep_increase[0] == 1){
        // We can only increase the timestep if all fluxes are allowed to be updated
        // If this is not done the timestep can't increase (since local_timestep is static)
        local_timestep = 1.0e+100;
    }

    // For all triangles
    for (k = 0; k < D->number_of_elements; k++) {
        speed_max_last = 0.0;

        // Loop through neighbours and compute edge flux for each
        for (i = 0; i < 3; i++) {
            ki = k * 3 + i; // Linear index to edge i of triangle k
            ki2 = 2 * ki; //k*6 + i*2
            ki3 = 3*ki; 

            if ((D->already_computed_flux[ki] == call) || (D->update_next_flux[ki]!=1)) {
                // We've already computed the flux across this edge
                // Check if it is a riverwall
                if(D->edge_flux_type[ki] == 1){
                    // Update counter of riverwall edges == index of
                    // riverwall_elevation + riverwall_rowIndex
                    RiverWall_count += 1;
                }
                continue;
            }

            // Get left hand side values from triangle k, edge i
            zl = D->bed_edge_values[ki];
            zc = D->bed_centroid_values[k];

            // stage on left edge used for sub-grid lookup table computation
            sle = D->stage_edge_values[ki];

            // Fully-submerged edge length
            length = D->edgelengths[ki];

            // Get right hand side values either from neighbouring triangle
            // or from boundary array (Quantities at neighbour on nearest face).
            n = D->neighbours[ki];
            if (n < 0) {
                // Neighbour is a boundary condition
                m = -n - 1; // Convert negative flag to boundary index
                // Force zero flux boundaries for now
                zr = zl; //1.0e+06; // Extend bed elevation to boundary

                sre = D->stage_edge_values[ki];
                //sre = D->stage_boundary_values[m];
                // FIXME: Permit momentum type boundary conditions
                zc_n = zc;
            } else {
                // Neighbour is a real triangle
                zc_n = D->bed_centroid_values[n];
                m = D->neighbour_edges[ki];
                nm = n * 3 + m; // Linear index (triangle n, edge m)
                nm3 = nm*3;

                zr = D->bed_edge_values[nm];

                // stage on right edge used for sub-grid lookup table computation
                sre = D->stage_edge_values[nm];
            }

            // Audusse magic
            //z_half=max(zl,zr);
            // 
            //// Account for riverwalls
            //if(D->edge_flux_type[ki]==1){
            //    // Update counter of riverwall edges == index of
            //    // riverwall_elevation + riverwall_rowIndex
            //    RiverWall_count+=1;
            //    
            //    // Set central bed to riverwall elevation
            //    z_half=max(D->riverwall_elevation[RiverWall_count-1], z_half) ;

            //}

            // Define h left/right for Audusse flux method
            //h_left= max(hle+zl-z_half,0.);
            //h_right= max(hre+zr-z_half,0.);
            //h_left= hle+min(zl-z_half, 0.);
            //h_right= hre+min(zr-z_half,0.);

            // Quick exit
            if( n >= 0 ){
                if(D->vol_centroid_values[k] == 0.0 && D->vol_centroid_values[n] == 0.0){
                    // Neither cell will have any flux
                    D->already_computed_flux[ki] = call ;
                    D->pressuregrad_work[ki] = 0.0;
                    D->edge_flux_work[ki3] = 0.0;
                    D->edge_flux_work[ki3+1] = 0.0;
                    D->edge_flux_work[ki3+2] = 0.0;
                    D->edge_timestep[ki] = local_timestep; // Arbitrary value but will not slow down code

                    D->already_computed_flux[nm] = call ;
                    D->pressuregrad_work[nm] = 0.0;
                    D->edge_flux_work[nm3] = 0.0;
                    D->edge_flux_work[nm3+1] = 0.0;
                    D->edge_flux_work[nm3+2] = 0.0;
                    D->edge_timestep[nm] = local_timestep;// Arbitrary value but will not slow down code

                    continue;
                }
            }

            // SUBGRID VERSION
            // The 'interpolation stage' for the edge is evaluated in a funky
            // way consistent with the Audusse approach.
            //
            // The most naive approach (first order) would be:
            // interpolation_stage_for_edge = stage_c = h_c + z_c
            // This could be passed to either the left or right interpolation
            // functions.
            //
            // In general we want higher order accuracy, which comes
            // from extrapolating stage or depth to the edge. Stage extrapolation
            // is a bad idea because of its troublesome wet-dry behaviour. Depth
            // extrapolation (modified based on the edge bed elevations) is
            // like what we would use for the Audusse method.
            //
            // Consider the value of stage used at the left edge for subgrid
            // lookup table interpolation.
            // If the left lookup table is used, we have:
            // interpolation_stage_for_edge = h_e + z_c 
            //                              = (stage_e - z_l) + z_c
            //
            // If the right lookup table is used, we have:
            // interpolation_stage_for_edge = h_e_above_zr + zc_n
            //                              = h_e + min(z_l - z_r, 0.) + zc_n [ This is applying the Audusse approach ]
            //                              = stage_e   - z_r + zc_n [as long as z_l > z_r]
            //
            // The choice of whether to use the left or right interpolation
            // functions, depending on which produces the minimum wet
            // cross-sectional area at the edge
            //
            // Enquire which is smaller here -- the one with the smallest area will always
            // be the one to use
            //
            max_stage = max(sle, sre);
            //lookup_table = &(D->subgrid_edge_table[D->subgrid_edge_starting_index[ki]]) ;
            cross_section_area_l = eval_subgrid_function( max_stage - zl + zc, k, i, 1, 0, D, 0, 0);
            if(n >= 0){
                //lookup_table = &(D->subgrid_edge_table[D->subgrid_edge_starting_index[nm]]) ;
                cross_section_area_r = eval_subgrid_function( max_stage - zr + zc_n, n, m, 1, 0, D, 0, 0) ;
            }else{
                // Ensure left side will be used (cross_section_area_l < cross_section_area_r)
                cross_section_area_r = cross_section_area_l + 99999.0;
            }

            //// Another chance for a quick exit
            // FIXME: Can do this but need to coordination for bedslope term
            //        since the grad(w) part may still matter
            //if( (cross_section_area_l == 0.) && (cross_section_area_r == 0.)){
            //    // Neither cell will have any flux
            //    D->already_computed_flux[ki] = call ;
            //    D->pressuregrad_work[ki] = 0.;
            //    D->edge_flux_work[ki3] = 0.;
            //    D->edge_flux_work[ki3+1] = 0.;
            //    D->edge_flux_work[ki3+2] = 0.;
            //    D->edge_timestep[ki] = local_timestep; // Arbitrary value but will not slow down code

            //    if(n >= 0){
            //        D->already_computed_flux[nm] = call ;
            //        D->pressuregrad_work[nm] = 0.;
            //        D->edge_flux_work[nm3] = 0.;
            //        D->edge_flux_work[nm3+1] = 0.;
            //        D->edge_flux_work[nm3+2] = 0.;
            //        D->edge_timestep[nm] = local_timestep;// Arbitrary value but will not slow down code
            //    }

            //    continue;
            //}


            if( (cross_section_area_l < cross_section_area_r) || (n < 0)){
                // Use left hand side elevation / lookup function

                // Get 'ql'
                _get_rotated_subgrid_edge_quantities(ql, k, i, ki, D,
                    sle - zl + zc, sle - zl + zc, D->normals[ki2], D->normals[ki2+1]);

                // Get 'qr'
                if( n >= 0){
                    // Lookup from triangle k -- stage only from right side.
                    _get_rotated_subgrid_edge_quantities(qr, k, i, nm, D,
                        sre - zl + zc, sre - zl + zc, D->normals[ki2], D->normals[ki2+1]);
                }else{
                    // Use 'left' values [qr = ql]
                    _get_rotated_subgrid_edge_quantities(qr, k, i, ki, D,
                        sre - zl + zc, sre - zl + zc, D->normals[ki2], D->normals[ki2+1]);
                }
            }else{
                // Use right hand side elevation / lookup function

                // Get 'ql'
                // Lookup table from triangle n -- stage only from left side.
                _get_rotated_subgrid_edge_quantities(ql, n, m, ki, D, 
                    sle - zr + zc_n, sle - zr + zc_n, D->normals[ki2], D->normals[ki2+1]);

                // Get 'qr'
                _get_rotated_subgrid_edge_quantities(qr, n, m, nm, D, 
                    sre - zr + zc_n, sre - zr + zc_n, D->normals[ki2], D->normals[ki2+1]);
            } 
            //////////////////////////////////////////////////////////
           
            // CHECK for NAN
            //for(ii = 0; ii< 9; ii++){
            //    if(ql[ii]!=ql[ii] | qr[ii]!=qr[ii]){
            //        printf("Var %d is NAN, %e, %e, %e, %e\n", ii, ql[ii], qr[ii], h_left, h_right);
            //        report_python_error(AT, "qvar is Nan");
            //        return -1.0;
            //    }
            //}

            // Variables for wave-speed computation
            u_left  = ql[1]/(ql[0]+1.0e-10); // [Integral uh] / [Integral h]
            u_right = qr[1]/(qr[0]+1.0e-10);// 
            hl = ql[0]/ql[8]; // Integral [h] / edge_length 
            hr = qr[0]/qr[8];

            // Edge flux computation (triangle k, edge i)
            _flux_function_central(ql, qr, edgeflux, 
                    D->normals[ki2],D->normals[ki2 + 1],
                    D->epsilon, D->g,
                    &max_speed_local, &pressure_flux, 
                    length, 
                    u_left, u_right, 
                    hl, hr);

            if(n < 0){
                // ENFORCE REFLECTIVE BOUNDARY FOR NOW
                edgeflux[0] = 0.;
                edgeflux[1] = 0.;
                edgeflux[2] = 0.;
            }

            // Check for NAN
            if( (edgeflux[0] != edgeflux[0]) || 
                (edgeflux[1] != edgeflux[1]) || 
                (edgeflux[2] != edgeflux[2]) ){
                    printf("------\n");
                    for(ii=0; ii<9; ii++){
                        printf("%e, %e\n", ql[ii], qr[ii]);
                    }
                    printf("%e, %e, %e, %e, %e,%e\n", 
                            edgeflux[0], edgeflux[1],edgeflux[2], 
                            max_speed_local, 
                            D->alphax_centroid_values[k], D->alphay_centroid_values[k]);
                    report_python_error(AT, "NAN edgeflux");
                    return -1.0;
            }

            //// Force weir discharge to match weir theory
            //if(D->edge_flux_type[ki]==1){
            //    report_python_error(AT,"Weir not supported");
            //    return -1.0;
            //    // Reference weir height  
            //    weir_height=max(D->riverwall_elevation[RiverWall_count-1]-min(zl,zr), 0.);
            //    // If the weir height is zero, avoid the weir computation
            //    // entirely
            //    if(weir_height>0.){
            //        ///////////////////////////////////////////////////////////
            //        // Use first-order h's for weir -- as the
            //        // 'upstream/downstream' heads are measured away from the
            //        // weir itself
            //        h_left_tmp= max(D->stage_centroid_values[k]-z_half,0.);
            //        if(n>=0){
            //            h_right_tmp= max(D->stage_centroid_values[n]-z_half,0.);
            //        }else{
            //            h_right_tmp= max(hc_n+zr-z_half,0.);
            //        }

            //        ///////////////////////////////////////////////////////////
            //        // Get Qfactor index - multiply the idealised weir
            //        // discharge by this constant factor
            //        ii = D->riverwall_rowIndex[RiverWall_count-1] * 
            //             D->ncol_riverwall_hydraulic_properties;
            //        Qfactor = D->riverwall_hydraulic_properties[ii];
            //        // Get s1, submergence ratio at which we start blending
            //        // with the shallow water solution 
            //        ii+=1;
            //        s1= D->riverwall_hydraulic_properties[ii];
            //        // Get s2, submergence ratio at which we entirely use the
            //        // shallow water solution 
            //        ii+=1;
            //        s2= D->riverwall_hydraulic_properties[ii];
            //        // Get h1, tailwater head / weir height at which we start
            //        // blending with the shallow water solution 
            //        ii+=1;
            //        h1= D->riverwall_hydraulic_properties[ii];
            //        // Get h2, tailwater head / weir height at which we
            //        // entirely use the shallow water solution 
            //        ii+=1;
            //        h2= D->riverwall_hydraulic_properties[ii];
            //        
            //        // Weir flux adjustment 
            //        adjust_edgeflux_with_weir(
            //            edgeflux, h_left_tmp, h_right_tmp, D->g, 
            //            weir_height, Qfactor, 
            //            s1, s2, h1, h2);
            //        // NOTE: Should perhaps also adjust the wave-speed here??
            //        // Small chance of instability??
            //    }
            //}
            

            //// Don't allow an outward advective flux if the cell centroid
            ////   stage is < the edge value. Is this important (??). Seems not
            ////   to be with DE algorithms
            //if((hc<H0) && edgeflux[0] > 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed_local=0.;
            //    //pressure_flux=0.;
            //}
            ////
            //if((hc_n<H0) && edgeflux[0] < 0.){
            //    edgeflux[0] = 0.;
            //    edgeflux[1] = 0.;
            //    edgeflux[2] = 0.;
            //    //max_speed_local=0.;
            //    //pressure_flux=0.;
            //}

            D->edge_flux_work[ki3 + 0 ] = -edgeflux[0];
            D->edge_flux_work[ki3 + 1 ] = -edgeflux[1];
            D->edge_flux_work[ki3 + 2 ] = -edgeflux[2];

            // bedslope_work contains all gravity related terms
            //bedslope_work=length*(- D->g *0.5*(h_left*h_left - hle*hle -(hle+hc)*(zl-zc))+pressure_flux);

            // Approach with subgrid topography
            if( D->vol_centroid_values[k] > 0.0 ){
                //
                // g d grad(z) = g d grad(w) - 0.5g grad(h^2)
                //
                // Notice this formula us true for ANY grad(w) so long as z, w,
                // and h are consistent-- it does not have to be the physical
                // water surface gradient. However, we compute it in such a way
                // that for smooth flows it IS the physical water surface
                // gradient, whereas for rapidly varying flows the limiters
                // kick in and it will be closer to 0.
                //
                // Integrate this to get the terms we need
                // Evaluate grad(w) using (limited) cell edge values (assuming
                // a constant within-cell grad(w)), grad(h^2) using the cell
                // audusse 'pressure gradient integrated over edges'.
                ql[7] = (D->vol_centroid_values[k] / D->areas[k]) * sle * length  ;
                bedslope_work = (D->g * (-0.5 * ql[6] + ql[7]) + pressure_flux);

            }else{
                bedslope_work = 0.;

            }

            //// Valiani + Begnudelli approach
            //if(D->vol_centroid_values[k]>0.){
            //    bedslope_work=(- D->g*0.5*(ql[7])+pressure_flux);
            //}else{
            //    bedslope_work = 0.;
            //}

            D->pressuregrad_work[ki]=bedslope_work;
            
            D->already_computed_flux[ki] = call; // #k Done

            // Update neighbour n with same flux but reversed sign
            if (n >= 0) {

                D->edge_flux_work[nm3 + 0 ] = edgeflux[0];
                D->edge_flux_work[nm3 + 1 ] = edgeflux[1];
                D->edge_flux_work[nm3 + 2 ] = edgeflux[2];
                //bedslope_work=length*(- D->g*0.5*(h_right*h_right-hre*hre-(hre+hc_n)*(zr-zc_n))+pressure_flux);

                // Approach with subgrid topography
                if(D->vol_centroid_values[n] > 0.){
                    qr[7] = (D->vol_centroid_values[n] / D->areas[n]) * sre * length ;
                    bedslope_work = ( D->g * ( -0.5 * qr[6] + qr[7]) + pressure_flux);
                }else{
                    bedslope_work = 0.;
                }

                //
                //// Valiani + Begnudelli approach
                //if(D->vol_centroid_values[n]>0.){
                //    bedslope_work=(- D->g*0.5*(qr[7])+pressure_flux);
                //}else{
                //    bedslope_work=0.;
                //}
                D->pressuregrad_work[nm]=bedslope_work;

                D->already_computed_flux[nm] = call; // #n Done
            }

            // Update timestep based on edge i and possibly neighbour n
            // NOTE: We should only change the timestep on the 'first substep'
            //  of the timestepping method [substep_count==0]
            if(substep_count==0){

                // Compute the 'edge-timesteps' (useful for setting flux_update_frequency)
                //D->edge_timestep[ki] = D->radii[k] / max(max_speed_local, D->epsilon);
                // CFL condition based on Sanders et al (2008)
                D->edge_timestep[ki] = D->subgrid_wet_area[k] /(ql[8] * max(max_speed_local, D->epsilon));

                if (n >= 0) {
                    //D->edge_timestep[nm] = D->radii[n]/ max(max_speed_local, D->epsilon);
                    // CFL condition based on Sanders et al (2008)
                    D->edge_timestep[nm] = D->subgrid_wet_area[n] /(qr[8] * max(max_speed_local, D->epsilon));
                }

                // Update the timestep
                if ((D->tri_full_flag[k] == 1)) {

                    speed_max_last=max(speed_max_last, max_speed_local);

                    if (max_speed_local > D->epsilon) {
                        // Apply CFL condition for triangles joining this edge (triangle k and triangle n)

                        // CFL for triangle k
                        local_timestep = min(local_timestep, D->edge_timestep[ki]);

                        if (n >= 0) {
                            // Apply CFL condition for neigbour n (which is on the ith edge of triangle k)
                            local_timestep = min(local_timestep, D->edge_timestep[nm]);
                        }
                    }
                }
            }

        } // End edge i (and neighbour n)
        // Keep track of maximal speeds
        if(substep_count==0) D->max_speed[k] = speed_max_last; //max_speed;


    } // End triangle k

    ////// Limit edgefluxes, for mass conservation near wet/dry cells
    //for(k=0; k< D->number_of_elements; k++){
    //    // Loop over every edge
    //    for(i = 0; i<3; i++){
    //        if(i==0){
    //            // Add up the outgoing flux through the cell -- only do this once (i==0)
    //            outgoing_mass_edges=0.0;
    //            for(useint=0; useint<3; useint++){
    //                if(D->edge_flux_work[3*(3*k+useint)]< 0.){
    //                    //outgoing_mass_edges+=1.0;
    //                    outgoing_mass_edges+=(D->edge_flux_work[3*(3*k+useint)]);
    //                }
    //            }
    //            outgoing_mass_edges*=local_timestep;
    //        }

    //        ki=3*k+i;   
    //        ki2=ki*2;
    //        ki3 = ki*3;
    //        
    //        // Prevent outflow from 'seriously' dry cells
    //        // Idea: The cell will not go dry if:
    //        // total_outgoing_flux <= cell volume = Area_triangle*hc
    //        vol=D->vol_centroid_values[k];
    //        if((D->edge_flux_work[ki3]< 0.0) && (-outgoing_mass_edges> vol)){
    //            
    //            // This bound could be improved (e.g. we could actually sum the
    //            // + and - fluxes and check if they are too large).  However,
    //            // the advantage of this method is that we don't have to worry
    //            // about subsequent changes to the + edgeflux caused by
    //            // constraints associated with neighbouring triangles.
    //            tmp = vol/(-(outgoing_mass_edges)) ;
    //            if(tmp< 1.0){
    //                D->edge_flux_work[ki3+0]*=tmp;
    //                D->edge_flux_work[ki3+1]*=tmp;
    //                D->edge_flux_work[ki3+2]*=tmp;

    //                // Compute neighbour edge index
    //                n = D->neighbours[ki];
    //                if(n>=0){
    //                    nm = 3*n + D->neighbour_edges[ki];
    //                    nm3 = nm*3;
    //                    D->edge_flux_work[nm3+0]*=tmp;
    //                    D->edge_flux_work[nm3+1]*=tmp;
    //                    D->edge_flux_work[nm3+2]*=tmp;
    //                }
    //            }
    //        }
    //    }
    // }

    // Now add up stage, xmom, ymom explicit updates
    for(k=0; k < D->number_of_elements; k++){
        hc = max(D->stage_centroid_values[k] - D->bed_centroid_values[k],0.);

        for(i=0;i<3;i++){
            // FIXME: Make use of neighbours to efficiently set things
            ki=3*k+i;
            ki2=ki*2;
            ki3 = ki*3;
            n=D->neighbours[ki];

            // GD HACK
            // Option to limit advective fluxes
            //if(hc > H0){
                D->vol_explicit_update[k] += D->edge_flux_work[ki3+0];
                D->u_vol_explicit_update[k] += D->edge_flux_work[ki3+1];
                D->v_vol_explicit_update[k] += D->edge_flux_work[ki3+2];
            //}else{
            //    stage_explicit_update[k] += edge_flux_work[ki3+0];
            //}


            // If this cell is not a ghost, and the neighbour is a boundary
            // condition OR a ghost cell, then add the flux to the
            // boundary_flux_integral
            if( (n<0 & D->tri_full_flag[k]==1) | ( n>=0 && (D->tri_full_flag[k]==1 & D->tri_full_flag[n]==0)) ){
                // boundary_flux_sum is an array with length = timestep_fluxcalls
                // For each sub-step, we put the boundary flux sum in.
                D->boundary_flux_sum[substep_count] += D->edge_flux_work[ki3];
            }

            // GD HACK
            // Compute bed slope term
            //if(hc > H0){
            D->u_vol_explicit_update[k] -= D->normals[ki2]*D->pressuregrad_work[ki];
            D->v_vol_explicit_update[k] -= D->normals[ki2+1]*D->pressuregrad_work[ki];
            //}else{
            //    xmom_explicit_update[k] *= 0.;
            //    ymom_explicit_update[k] *= 0.;
            //}

        } // end edge i

        // FRICTION -- only once per cell
        alpha_norm = sqrt(pow(D->alphax_centroid_values[k],2)+pow(D->alphay_centroid_values[k],2));
        sfx = D->alphax_centroid_values[k]*alpha_norm ;
        sfy = D->alphay_centroid_values[k]*alpha_norm ; 
        //D->u_vol_explicit_update[k] -= D->g * D->vol_centroid_values[k]*sfx ;
        //D->v_vol_explicit_update[k] -= D->g * D->vol_centroid_values[k]*sfy ;
        D->u_vol_semi_implicit_update[k] -= D->g * D->vol_centroid_values[k]*sfx; 
        D->v_vol_semi_implicit_update[k] -= D->g * D->vol_centroid_values[k]*sfy; 

    }  // end cell k

    // Ensure we only update the timestep on the first call within each rk2/rk3 step
    if(substep_count==0){
        timestep=local_timestep; 
        //D->timestep = timestep;
    } 
    return timestep;
}


// Protect against the water elevation falling below the triangle bed
inline double  _protect(int N,
         double minimum_allowed_height,
         double maximum_allowed_speed,
         double epsilon,
         double* wc,
         double* wv,
         double* zc,
         double* zv,
         double* xmomc,
         double* ymomc,
         double* areas,
         double* xc,
         double* yc) {

  int k;
  double hc, bmin, bmax;
  double u, v, reduced_speed;
  double mass_error = 0.;
  // This acts like minimum_allowed height, but scales with the vertical
  // distance between the bed_centroid_value and the max bed_edge_value of
  // every triangle.
  //double minimum_relative_height=0.05; 
  int mass_added = 0;

  // Protect against inifintesimal and negative heights  
  //if (maximum_allowed_speed < epsilon) {
    for (k=0; k<N; k++) {
      hc = wc[k] - zc[k];
      if (hc < minimum_allowed_height*1.0 ){
            // Set momentum to zero and ensure h is non negative
            xmomc[k] = 0.;
            ymomc[k] = 0.;
        if (hc <= 0.0){
             bmin = zc[k];
             // Minimum allowed stage = bmin

             // WARNING: ADDING MASS if wc[k]<bmin
             if(wc[k] < bmin){
                 mass_error += (bmin-wc[k])*areas[k];
                 mass_added = 1; //Flag to warn of added mass                

                 wc[k] = bmin; 

                 // FIXME: Set vertex values as well. Seems that this shouldn't be
                 // needed. However, from memory this is important at the first
                 // time step, for 'dry' areas where the designated stage is
                 // less than the bed centroid value
                 wv[3*k] = min(bmin, wc[k]); //zv[3*k]-minimum_allowed_height);
                 wv[3*k+1] = min(bmin, wc[k]); //zv[3*k+1]-minimum_allowed_height);
                 wv[3*k+2] = min(bmin, wc[k]); //zv[3*k+2]-minimum_allowed_height);
            }
        }
      }
    }

  //if(mass_added==1){
  //  printf("Cumulative mass protection: %f m^3 \n", mass_error);
  //}
  return mass_error;
}

inline int find_qmin_and_qmax(double dq0, double dq1, double dq2, 
               double *qmin, double *qmax){
  // Considering the centroid of an FV triangle and the vertices of its 
  // auxiliary triangle, find 
  // qmin=min(q)-qc and qmax=max(q)-qc, 
  // where min(q) and max(q) are respectively min and max over the
  // four values (at the centroid of the FV triangle and the auxiliary 
  // triangle vertices),
  // and qc is the centroid
  // dq0=q(vertex0)-q(centroid of FV triangle)
  // dq1=q(vertex1)-q(vertex0)
  // dq2=q(vertex2)-q(vertex0)

  // This is a simple implementation 
  *qmax = max(max(dq0, max(dq0+dq1, dq0+dq2)), 0.0) ;
  *qmin = min(min(dq0, min(dq0+dq1, dq0+dq2)), 0.0) ;
 
  return 0;
}

inline int limit_gradient(double *dqv, double qmin, double qmax, double beta_w){
  // Given provisional jumps dqv from the FV triangle centroid to its
  // vertices/edges, and jumps qmin (qmax) between the centroid of the FV
  // triangle and the minimum (maximum) of the values at the auxiliary triangle
  // vertices (which are centroids of neighbour mesh triangles), calculate a
  // multiplicative factor phi by which the provisional vertex jumps are to be
  // limited
  
  int i;
  double r=1000.0, r0=1.0, phi=1.0;
  static double TINY = 1.0e-100; // to avoid machine accuracy problems.
  // FIXME: Perhaps use the epsilon used elsewhere.
  
  // Any provisional jump with magnitude < TINY does not contribute to 
  // the limiting process.
  //return 0;
  
  for (i=0;i<3;i++){
    if (dqv[i]<-TINY)
      r0=qmin/dqv[i];
      
    if (dqv[i]>TINY)
      r0=qmax/dqv[i];
      
    r=min(r0,r);
  }
  
  phi=min(r*beta_w,1.0);
  //phi=1.;
  dqv[0]=dqv[0]*phi;
  dqv[1]=dqv[1]*phi;
  dqv[2]=dqv[2]*phi;

  return 0;
}

// Computational routine
inline int _extrapolate_second_order_edge_sw(struct domain *D){
                  
  // Local variables
  double a, b; // Gradient vector used to calculate edge values from centroids
  int k, k0, k1, k2, k3, k6, coord_index, i, ii, ktmp, k_wetdry;
  double x, y, x0, y0, x1, y1, x2, y2, xv0, yv0, xv1, yv1, xv2, yv2; // Vertices of the auxiliary triangle
  double dx1, dx2, dy1, dy2, dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dq0, dq1, dq2, area2, inv_area2, dpth,momnorm;
  double dqv[3], qmin, qmax, hmin, hmax, bedmax,bedmin, stagemin;
  double hc, h0, h1, h2, beta_tmp, hfactor, xtmp, ytmp, weight, tmp, wet_frac_min, dt;
  double dk, dk_inv,dv0, dv1, dv2, de[3], demin, dcmax, r0scale, vel_norm, l1, l2, a_tmp, b_tmp, c_tmp,d_tmp;
  

  memset((char*) D->x_centroid_work, 0, D->number_of_elements * sizeof (double));
  memset((char*) D->y_centroid_work, 0, D->number_of_elements * sizeof (double));
 
  // Some numbers we need 
  a_tmp=0.2; // Highest depth ratio with hfactor=1
  b_tmp=0.1; // Highest depth ratio with hfactor=0
  c_tmp=1.0/(a_tmp-b_tmp); 
  d_tmp= 1.0-(c_tmp*a_tmp);
 
  dt = D->timestep;
  //printf("dt: %e \n", dt);
  //if(D->extrapolate_velocity_second_order==1){

  //    // Replace momentum centroid with velocity centroid to allow velocity
  //    // extrapolation This will be changed back at the end of the routine
  //    for (k=0; k< D->number_of_elements; k++){
  //        
  //        D->height_centroid_values[k] = max(D->stage_centroid_values[k] - D->bed_centroid_values[k], 0.);

  //        dk = D->height_centroid_values[k]; 
  //        if(dk> D->minimum_allowed_height){
  //            dk_inv=1.0/dk;
  //            D->x_centroid_work[k] = D->alphax_centroid_values[k];
  //            D->alphax_centroid_values[k] = D->alphax_centroid_values[k]*dk_inv;

  //            D->y_centroid_work[k] = D->alphay_centroid_values[k];
  //            D->alphay_centroid_values[k] = D->alphay_centroid_values[k]*dk_inv;
  //        }else{
  //            D->x_centroid_work[k] = 0.;
  //            D->alphax_centroid_values[k] = 0.;
  //            D->y_centroid_work[k] = 0.;
  //            D->alphay_centroid_values[k] = 0.;

  //       }
  //    }
  //}

  //// If a triangle is surrounded by dry cells (or dry cells + boundary
  //// condition) set its momentum to zero too. This prevents 'pits' of
  //// of water being trapped and unable to lose momentum, which can occur in
  //// some situations
  //for (k=0; k< D->number_of_elements;k++){
  //    
  //    k3=k*3;
  //    k0 = D->surrogate_neighbours[k3];
  //    k1 = D->surrogate_neighbours[k3 + 1];
  //    k2 = D->surrogate_neighbours[k3 + 2];

  //    if((D->height_centroid_values[k0] < D->minimum_allowed_height | k0==k) &
  //       (D->height_centroid_values[k1] < D->minimum_allowed_height | k1==k) &
  //       (D->height_centroid_values[k2] < D->minimum_allowed_height | k2==k)){
  //  	  	  //printf("Surrounded by dry cells\n");
  //            D->x_centroid_work[k] = 0.;
  //            D->alphax_centroid_values[k] = 0.;
  //            D->y_centroid_work[k] = 0.;
  //            D->alphay_centroid_values[k] = 0.;

  //    }


  //}

  // Begin extrapolation routine
  for (k = 0; k < D->number_of_elements; k++) 
  {

    // Don't update the extrapolation if the flux will not be computed on the
    // next timestep
    if(D->update_extrapolation[k]==0){
       continue;
    }


    // Useful indices
    k3=k*3;
    k6=k*6;

    if (D->number_of_boundaries[k]==3)
    {
      // No neighbours, set gradient on the triangle to zero

      D->stage_edge_values[k3]   = D->stage_centroid_values[k];
      D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
      D->stage_edge_values[k3+2] = D->stage_centroid_values[k];

      D->alphax_edge_values[k3]    = D->alphax_centroid_values[k];
      D->alphax_edge_values[k3+1]  = D->alphax_centroid_values[k];
      D->alphax_edge_values[k3+2]  = D->alphax_centroid_values[k];
      D->alphay_edge_values[k3]    = D->alphay_centroid_values[k];
      D->alphay_edge_values[k3+1]  = D->alphay_centroid_values[k];
      D->alphay_edge_values[k3+2]  = D->alphay_centroid_values[k];

      //dk = D->height_centroid_values[k];
      D->height_edge_values[k3] = D->stage_edge_values[k3] - D->bed_edge_values[k3];
      D->height_edge_values[k3+1] = D->stage_edge_values[k3+1] - D->bed_edge_values[k3+1];
      D->height_edge_values[k3+2] = D->stage_edge_values[k3+2] - D->bed_edge_values[k3+2];
      //D->height_edge_values[k3+1] = dk;
      //D->height_edge_values[k3+2] = dk;
      
      continue;
    }
    else 
    {
      // Triangle k has one or more neighbours. 
      // Get centroid and edge coordinates of the triangle

      // Get the edge coordinates
      xv0 = D->edge_coordinates[k6];
      yv0 = D->edge_coordinates[k6+1];
      xv1 = D->edge_coordinates[k6+2];
      yv1 = D->edge_coordinates[k6+3];
      xv2 = D->edge_coordinates[k6+4];
      yv2 = D->edge_coordinates[k6+5];

      // Get the centroid coordinates
      coord_index = 2*k;
      x = D->centroid_coordinates[coord_index];
      y = D->centroid_coordinates[coord_index+1];

      // Store x- and y- differentials for the edges of 
      // triangle k relative to the centroid
      dxv0 = xv0 - x;
      dxv1 = xv1 - x;
      dxv2 = xv2 - x;
      dyv0 = yv0 - y;
      dyv1 = yv1 - y;
      dyv2 = yv2 - y;

    }


    if (D->number_of_boundaries[k]<=1)
    {
      //==============================================
      // Number of boundaries <= 1
      // 'Typical case'
      //==============================================    


      // If no boundaries, auxiliary triangle is formed 
      // from the centroids of the three neighbours
      // If one boundary, auxiliary triangle is formed 
      // from this centroid and its two neighbours

      k0 = D->surrogate_neighbours[k3];
      k1 = D->surrogate_neighbours[k3 + 1];
      k2 = D->surrogate_neighbours[k3 + 2];

      // Get the auxiliary triangle's vertex coordinates 
      // (really the centroids of neighbouring triangles)
      coord_index = 2*k0;
      x0 = D->centroid_coordinates[coord_index];
      y0 = D->centroid_coordinates[coord_index+1];

      coord_index = 2*k1;
      x1 = D->centroid_coordinates[coord_index];
      y1 = D->centroid_coordinates[coord_index+1];

      coord_index = 2*k2;
      x2 = D->centroid_coordinates[coord_index];
      y2 = D->centroid_coordinates[coord_index+1];

      // Store x- and y- differentials for the vertices
      // of the auxiliary triangle
      dx1 = x1 - x0;
      dx2 = x2 - x0;
      dy1 = y1 - y0;
      dy2 = y2 - y0;

      // Calculate 2*area of the auxiliary triangle
      // The triangle is guaranteed to be counter-clockwise
      area2 = dy2*dx1 - dy1*dx2;
      //if(k==54) printf("K=54\n");
       
      //// Treat triangles with no neighbours (area2 <=0.)
      if ((area2 <= 0.))
      {


          // Isolated wet cell -- constant stage/depth extrapolation
          D->stage_edge_values[k3]   = D->stage_centroid_values[k];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k];

          //dk= D->height_centroid_values[k]; //max(stage_centroid_values[k]-bed_centroid_values[k],0.);
          //D->height_edge_values[k3] = dk;
          //D->height_edge_values[k3+1] = dk;
          //D->height_edge_values[k3+2] = dk;
          D->height_edge_values[k3] = D->stage_edge_values[k3] - D->bed_edge_values[k3];
          D->height_edge_values[k3+1] = D->stage_edge_values[k3+1] - D->bed_edge_values[k3+1];
          D->height_edge_values[k3+2] = D->stage_edge_values[k3+2] - D->bed_edge_values[k3+2];

          D->alphax_edge_values[k3]    = D->alphax_centroid_values[k];
          D->alphax_edge_values[k3+1]  = D->alphax_centroid_values[k];
          D->alphax_edge_values[k3+2]  = D->alphax_centroid_values[k];

          D->alphay_edge_values[k3]    = D->alphay_centroid_values[k];
          D->alphay_edge_values[k3+1]  = D->alphay_centroid_values[k];
          D->alphay_edge_values[k3+2]  = D->alphay_centroid_values[k];

          continue;
      }

      // // Calculate heights of neighbouring cells
      hc = D->vol_centroid_values[k] / D->areas[k]; //D->height_centroid_values[k];
      h0 = D->vol_centroid_values[k0] / D->areas[k0]; ////D->height_centroid_values[k0];
      h1 = D->vol_centroid_values[k1] / D->areas[k1]; //D->height_centroid_values[k1];
      h2 = D->vol_centroid_values[k2] / D->areas[k2]; //D->height_centroid_values[k2];
      hmin = max(min(min(h0, min(h1, h2)), hc), 0.);
      hmax = max(max(max(h0, max(h1, h2)), hc), 0.);
      // Look for strong changes in cell depth as an indicator of near-wet-dry
      // Reduce hfactor linearly from 1-0 between depth ratio (hmin/hc) of [a_tmp , b_tmp]
      // NOTE: If we have a more 'second order' treatment in near dry areas (e.g. with b_tmp being negative), then
      //       the water tends to dry more rapidly (which is in agreement with analytical results),
      //       but is also more 'artefacty' in important cases (tendency for high velocities, etc).
      //       
      // So hfactor = depth_ratio*(c_tmp) + d_tmp, but is clipped between 0 and 1.
      hfactor= max(0.,
                   min(c_tmp*max(hmin,0.0)/max(hc,1.0e-06)+d_tmp,
                       min(c_tmp*max(hc,0.)/max(hmax,1.0e-06)+d_tmp, 1.0))
                  );

      //hfactor = 1.0;

      //// Try applying lower beta as wet fraction varies from 0.3 to 0.1
      //hc = D->subgrid_wet_area[k]/D->areas[k]; //D->height_centroid_values[k];
      //h0 = D->subgrid_wet_area[k0]/D->areas[k0]; ////D->height_centroid_values[k0];
      //h1 = D->subgrid_wet_area[k1]/D->areas[k1]; //D->height_centroid_values[k1];
      //h2 = D->subgrid_wet_area[k2]/D->areas[k2]; //D->height_centroid_values[k2];
      //wet_frac_min =  min(hc, min(h0, min(h1, h2)));
      //hfactor = max(0., min(c_tmp*wet_frac_min + d_tmp, 1.0));

      //// Try applying lower beta as wet fraction varies from 0.3 to 0.1
      //hc = sqrt( 
      //          (D->u_vol_centroid_values[k]*D->u_vol_centroid_values[k] +
      //           D->v_vol_centroid_values[k]*D->v_vol_centroid_values[k])*dt*dt /
      //          ((D->vol_centroid_values[k] + 1.0e-6)*(D->vol_centroid_values[k] + 1.0e-06))
      //         ); //D->height_centroid_values[k];
      //h0 = sqrt( 
      //          (D->u_vol_centroid_values[k0]*D->u_vol_centroid_values[k0] +
      //           D->v_vol_centroid_values[k0]*D->v_vol_centroid_values[k0])*dt*dt /
      //          ( (D->vol_centroid_values[k0] + 1.0e-06)*(D->vol_centroid_values[k0] + 1.0e-06))
      //         ); //D->height_centroid_values[k];
      //h1 = sqrt( 
      //          (D->u_vol_centroid_values[k1]*D->u_vol_centroid_values[k1] +
      //           D->v_vol_centroid_values[k1]*D->v_vol_centroid_values[k1])*dt*dt /
      //          ((D->vol_centroid_values[k1] + 1.0e-06)*(D->vol_centroid_values[k1] + 1.0e-06))
      //         ); //D->height_centroid_values[k];
      //h2 = sqrt( 
      //          (D->u_vol_centroid_values[k2]*D->u_vol_centroid_values[k2] +
      //           D->v_vol_centroid_values[k2]*D->v_vol_centroid_values[k2])*dt*dt /
      //          ((D->vol_centroid_values[k2] + 1.0e-6)*(D->vol_centroid_values[k2] + 1.0e-06))
      //         ); //D->height_centroid_values[k];

      //wet_frac_min =  max(hc, min(h0, min(h1, h2)));

      //if(wet_frac_min > 0.9){
      //    hfactor = 0.;
      //}
      //h0 = D->subgrid_wet_area[k0]/D->areas[k0]; ////D->height_centroid_values[k0];
      //h1 = D->subgrid_wet_area[k1]/D->areas[k1]; //D->height_centroid_values[k1];
      //h2 = D->subgrid_wet_area[k2]/D->areas[k2]; //D->height_centroid_values[k2];

      //wet_frac_min =  min(hc, min(h0, min(h1, h2)));
      //hfactor = max(0., min(c_tmp*wet_frac_min + d_tmp, 1.0));
      // Set hfactor to zero smothly as hmin--> minimum_allowed_height. This
      // // avoids some 'chatter' for very shallow flows 
      // //hfactor=min( 1.2*max(hmin- D->minimum_allowed_height,0.)/(max(hmin,0.)+1.* D->minimum_allowed_height), hfactor);
      //hfactor = 1.0;

      inv_area2 = 1.0/area2;
      //-----------------------------------
      // stage
      //-----------------------------------

      beta_tmp = D->beta_w_dry + (D->beta_w - D->beta_w_dry) * hfactor;

      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary 
          // triangle and the centroid of triangle k
          dq0 = D->stage_centroid_values[k0] - D->stage_centroid_values[k];
          
          // Calculate differentials between the vertices 
          // of the auxiliary triangle (centroids of neighbouring triangles)
          dq1 = D->stage_centroid_values[k1] - D->stage_centroid_values[k0];
          dq2 = D->stage_centroid_values[k2] - D->stage_centroid_values[k0];
         
          // Calculate the gradient of stage on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;
          // Calculate provisional jumps in stage from the centroid 
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0 + b*dyv0;
          dqv[1] = a*dxv1 + b*dyv1;
          dqv[2] = a*dxv2 + b*dyv2;
        
          // Now we want to find min and max of the centroid and the 
          // vertices of the auxiliary triangle and compute jumps 
          // from the centroid to the min and max
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
          
          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          D->stage_edge_values[k3+0] = D->stage_centroid_values[k] + dqv[0];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k] + dqv[1];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k] + dqv[2];
      }else{
          // Fast alternative when beta_tmp==0
          D->stage_edge_values[k3+0] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+1] = D->stage_centroid_values[k];
          D->stage_edge_values[k3+2] = D->stage_centroid_values[k];
      }


      //-----------------------------------
      // height
      //-----------------------------------

      //D->height_edge_values[k3] = D->stage_edge_values[k3] - D->bed_edge_values[k3];
      //D->height_edge_values[k3+1] = D->stage_edge_values[k3+1] - D->bed_edge_values[k3+1];
      //D->height_edge_values[k3+2] = D->stage_edge_values[k3+2] - D->bed_edge_values[k3+2];
      if(beta_tmp>0.){ 
      //if(0==1){
          // Calculate the difference between vertex 0 of the auxiliary 
          // triangle and the centroid of triangle k
          dq0 = D->height_centroid_values[k0] - D->height_centroid_values[k];
          
          // Calculate differentials between the vertices 
          // of the auxiliary triangle (centroids of neighbouring triangles)
          dq1 = D->height_centroid_values[k1] - D->height_centroid_values[k0];
          dq2 = D->height_centroid_values[k2] - D->height_centroid_values[k0];
         
          // Calculate the gradient of height on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;
          // Calculate provisional jumps in height from the centroid 
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0 + b*dyv0;
          dqv[1] = a*dxv1 + b*dyv1;
          dqv[2] = a*dxv2 + b*dyv2;
        
          // Now we want to find min and max of the centroid and the 
          // vertices of the auxiliary triangle and compute jumps 
          // from the centroid to the min and max
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);
       
          // Limit the gradient
          // Same beta_tmp as for stage
          //beta_tmp = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          //beta_tmp = 0. + (beta_w - 0.) * hfactor;

          D->height_edge_values[k3+0] = D->height_centroid_values[k] + dqv[0];
          D->height_edge_values[k3+1] = D->height_centroid_values[k] + dqv[1];
          D->height_edge_values[k3+2] = D->height_centroid_values[k] + dqv[2];
      }else{
          // Fast alternative when beta_tmp==0
          D->height_edge_values[k3+0] = D->height_centroid_values[k];
          D->height_edge_values[k3+1] = D->height_centroid_values[k];
          D->height_edge_values[k3+2] = D->height_centroid_values[k];
      }
      //-----------------------------------
      // alphaxentum
      //-----------------------------------

      beta_tmp = D->beta_uh_dry + (D->beta_uh - D->beta_uh_dry) * hfactor;
      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->alphax_centroid_values[k0] - D->alphax_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle
          dq1 = D->alphax_centroid_values[k1] - D->alphax_centroid_values[k0];
          dq2 = D->alphax_centroid_values[k2] - D->alphax_centroid_values[k0];

          // Calculate the gradient of alphax on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;

          // Calculate provisional jumps in stage from the centroid
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0+b*dyv0;
          dqv[1] = a*dxv1+b*dyv1;
          dqv[2] = a*dxv2+b*dyv2;

          // Now we want to find min and max of the centroid and the 
          // vertices of the auxiliary triangle and compute jumps 
          // from the centroid to the min and max
          //
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          for (i=0; i < 3; i++)
          {
            D->alphax_edge_values[k3+i] = D->alphax_centroid_values[k] + dqv[i];
          }
      }else{
          // Fast alternative when beta_tmp==0
          for (i=0; i < 3; i++)
          {
            D->alphax_edge_values[k3+i] = D->alphax_centroid_values[k];
          }
      }

      //-----------------------------------
      // alphayentum
      //-----------------------------------

      beta_tmp = D->beta_vh_dry + (D->beta_vh - D->beta_vh_dry) * hfactor;

      if(beta_tmp>0.){
          // Calculate the difference between vertex 0 of the auxiliary
          // triangle and the centroid of triangle k
          dq0 = D->alphay_centroid_values[k0] - D->alphay_centroid_values[k];

          // Calculate differentials between the vertices
          // of the auxiliary triangle
          dq1 = D->alphay_centroid_values[k1] - D->alphay_centroid_values[k0];
          dq2 = D->alphay_centroid_values[k2] - D->alphay_centroid_values[k0];

          // Calculate the gradient of alphax on the auxiliary triangle
          a = dy2*dq1 - dy1*dq2;
          a *= inv_area2;
          b = dx1*dq2 - dx2*dq1;
          b *= inv_area2;

          // Calculate provisional jumps in stage from the centroid
          // of triangle k to its vertices, to be limited
          dqv[0] = a*dxv0 + b*dyv0;
          dqv[1] = a*dxv1 + b*dyv1;
          dqv[2] = a*dxv2 + b*dyv2;

          // Now we want to find min and max of the centroid and the
          // vertices of the auxiliary triangle and compute jumps
          // from the centroid to the min and max
          //
          find_qmin_and_qmax(dq0, dq1, dq2, &qmin, &qmax);

          // Limit the gradient
          limit_gradient(dqv, qmin, qmax, beta_tmp);

          for (i=0;i<3;i++)
          {
            D->alphay_edge_values[k3 + i] = D->alphay_centroid_values[k] + dqv[i];
          }
      }else{
          // Fast alternative when beta_tmp==0
          for (i=0;i<3;i++)
          {
            D->alphay_edge_values[k3 + i] = D->alphay_centroid_values[k];
          }

      }

    } // End number_of_boundaries <=1 
    else
    {

      //==============================================
      // Number of boundaries == 2
      //==============================================

      // One internal neighbour and gradient is in direction of the neighbour's
      // centroid

      // Find the only internal neighbour (k1?)
      for (k2 = k3; k2 < k3 + 3; k2++)
      {
      // Find internal neighbour of triangle k
      // k2 indexes the edges of triangle k

          if (D->surrogate_neighbours[k2] != k)
          {
             break;
          }
      }
      
      if ((k2 == k3 + 3)) 
      {
        // If we didn't find an internal neighbour
        report_python_error(AT, "Internal neighbour not found");
        return -1;
      }
      
      k1 = D->surrogate_neighbours[k2];
      
      // The coordinates of the triangle are already (x,y). 
      // Get centroid of the neighbour (x1,y1)
      coord_index = 2*k1;
      x1 = D->centroid_coordinates[coord_index];
      y1 = D->centroid_coordinates[coord_index + 1];
      
      // Compute x- and y- distances between the centroid of 
      // triangle k and that of its neighbour
      dx1 = x1 - x; 
      dy1 = y1 - y;
      
      // Set area2 as the square of the distance
      area2 = dx1*dx1 + dy1*dy1;
      
      // Set dx2=(x1-x0)/((x1-x0)^2+(y1-y0)^2) 
      // and dy2=(y1-y0)/((x1-x0)^2+(y1-y0)^2) which
      // respectively correspond to the x- and y- gradients 
      // of the conserved quantities
      dx2 = 1.0/area2;
      dy2 = dx2*dy1;
      dx2 *= dx1;
      
      
      //-----------------------------------
      // stage
      //-----------------------------------            

      // Compute differentials
      dq1 = D->stage_centroid_values[k1] - D->stage_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now limit the jumps
      if (dq1>=0.0)
      {
        qmin=0.0;
        qmax=dq1;
      }
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, D->beta_w);
      
      D->stage_edge_values[k3] = D->stage_centroid_values[k] + dqv[0];
      D->stage_edge_values[k3 + 1] = D->stage_centroid_values[k] + dqv[1];
      D->stage_edge_values[k3 + 2] = D->stage_centroid_values[k] + dqv[2];

      ////-----------------------------------
      //// height
      ////-----------------------------------
      //
      //// Compute differentials
      dq1 = D->height_centroid_values[k1] - D->height_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now limit the jumps
      if (dq1>=0.0)
      {
        qmin=0.0;
        qmax=dq1;
      }
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, D->beta_w);
      
      D->height_edge_values[k3] = D->height_centroid_values[k] + dqv[0];
      D->height_edge_values[k3 + 1] = D->height_centroid_values[k] + dqv[1];
      D->height_edge_values[k3 + 2] = D->height_centroid_values[k] + dqv[2];
      //D->height_edge_values[k3] = D->stage_edge_values[k3] - D->bed_edge_values[k3];
      //D->height_edge_values[k3+1] = D->stage_edge_values[k3+1] - D->bed_edge_values[k3+1];
      //D->height_edge_values[k3+2] = D->stage_edge_values[k3+2] - D->bed_edge_values[k3+2];

      //-----------------------------------
      // alphaxentum
      //-----------------------------------                        
      
      // Compute differentials
      dq1 = D->alphax_centroid_values[k1] - D->alphax_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0+b*dyv0;
      dqv[1] = a*dxv1+b*dyv1;
      dqv[2] = a*dxv2+b*dyv2;
      
      // Now limit the jumps
      if (dq1 >= 0.0)
      {
        qmin = 0.0;
        qmax = dq1;
      }
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient      
      limit_gradient(dqv, qmin, qmax, D->beta_w);
      
      for (i = 0; i < 3;i++)
      {
          D->alphax_edge_values[k3 + i] = D->alphax_centroid_values[k] + dqv[i];
      }
      
      //-----------------------------------
      // alphayentum
      //-----------------------------------                        

      // Compute differentials
      dq1 = D->alphay_centroid_values[k1] - D->alphay_centroid_values[k];
      
      // Calculate the gradient between the centroid of triangle k 
      // and that of its neighbour
      a = dq1*dx2;
      b = dq1*dy2;
      
      // Calculate provisional edge jumps, to be limited
      dqv[0] = a*dxv0 + b*dyv0;
      dqv[1] = a*dxv1 + b*dyv1;
      dqv[2] = a*dxv2 + b*dyv2;
      
      // Now limit the jumps
      if (dq1>=0.0)
      {
        qmin = 0.0;
        qmax = dq1;
      } 
      else
      {
        qmin = dq1;
        qmax = 0.0;
      }
      
      // Limit the gradient
      limit_gradient(dqv, qmin, qmax, D->beta_w);
      
      for (i=0;i<3;i++)
              {
              D->alphay_edge_values[k3 + i] = D->alphay_centroid_values[k] + dqv[i];
              }
    } // else [number_of_boundaries==2]
  } // for k=0 to number_of_elements-1


  // Compute vertex values of quantities
  for (k=0; k< D->number_of_elements; k++){
      //if(D->extrapolate_velocity_second_order==1){
      //    //Convert velocity back to momenta at centroids
      //    D->alphax_centroid_values[k] = D->x_centroid_work[k];
      //    D->alphay_centroid_values[k] = D->y_centroid_work[k];
      //}
     
      // Don't proceed if we didn't update the edge/vertex values
      if(D->update_extrapolation[k]==0){
         continue;
      }

      k3=3*k;
      
      // Compute stage vertex values 
      D->stage_vertex_values[k3] = D->stage_edge_values[k3+1] + D->stage_edge_values[k3+2] - D->stage_edge_values[k3] ; 
      D->stage_vertex_values[k3+1] =  D->stage_edge_values[k3] + D->stage_edge_values[k3+2]- D->stage_edge_values[k3+1]; 
      D->stage_vertex_values[k3+2] =  D->stage_edge_values[k3] + D->stage_edge_values[k3+1]- D->stage_edge_values[k3+2]; 
      
      // Compute height vertex values 
      D->height_vertex_values[k3] = D->height_edge_values[k3+1] + D->height_edge_values[k3+2] - D->height_edge_values[k3] ; 
      D->height_vertex_values[k3+1] =  D->height_edge_values[k3] + D->height_edge_values[k3+2]- D->height_edge_values[k3+1]; 
      D->height_vertex_values[k3+2] =  D->height_edge_values[k3] + D->height_edge_values[k3+1]- D->height_edge_values[k3+2]; 

      //// If needed, convert from velocity to momenta
      //if(D->extrapolate_velocity_second_order==1){
      //    // Re-compute momenta at edges
      //    for (i=0; i<3; i++){
      //        dk= D->height_edge_values[k3+i]; 
      //        D->alphax_edge_values[k3+i] = D->alphax_edge_values[k3+i]*dk;
      //        D->alphay_edge_values[k3+i] = D->alphay_edge_values[k3+i]*dk;
      //    }
      //}
      // Compute momenta at vertices -- 
      //D->alphax_vertex_values[k3]   =  D->alphax_edge_values[k3+1] + D->alphax_edge_values[k3+2] - D->alphax_edge_values[k3] ; 
      //D->alphax_vertex_values[k3+1] =  D->alphax_edge_values[k3] + D->alphax_edge_values[k3+2]- D->alphax_edge_values[k3+1]; 
      //D->alphax_vertex_values[k3+2] =  D->alphax_edge_values[k3] + D->alphax_edge_values[k3+1]- D->alphax_edge_values[k3+2]; 
      //D->alphay_vertex_values[k3]   =  D->alphay_edge_values[k3+1] + D->alphay_edge_values[k3+2] - D->alphay_edge_values[k3] ; 
      //D->alphay_vertex_values[k3+1] =  D->alphay_edge_values[k3] + D->alphay_edge_values[k3+2]- D->alphay_edge_values[k3+1]; 
      //D->alphay_vertex_values[k3+2] =  D->alphay_edge_values[k3] + D->alphay_edge_values[k3+1]- D->alphay_edge_values[k3+2]; 

      //// Compute new bed elevation
      D->bed_edge_values[k3]= D->stage_edge_values[k3]- D->height_edge_values[k3];
      D->bed_edge_values[k3+1]= D->stage_edge_values[k3+1]- D->height_edge_values[k3+1];
      D->bed_edge_values[k3+2]= D->stage_edge_values[k3+2]- D->height_edge_values[k3+2];
      D->bed_vertex_values[k3] = D->bed_edge_values[k3+1] + D->bed_edge_values[k3+2] - D->bed_edge_values[k3] ; 
      D->bed_vertex_values[k3+1] =  D->bed_edge_values[k3] + D->bed_edge_values[k3+2] - D->bed_edge_values[k3+1]; 
      D->bed_vertex_values[k3+2] =  D->bed_edge_values[k3] + D->bed_edge_values[k3+1] - D->bed_edge_values[k3+2]; 
  } 

  return 0;
}           

//=========================================================================
// Python Glue
//=========================================================================


//========================================================================
// Compute fluxes
//========================================================================

// Modified central flux function

PyObject *compute_fluxes_ext_central(PyObject *self, PyObject *args) {
  /*Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Compute total flux for each conserved quantity using "flux_function_central"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities
    stage, xmomentum and ymomentum

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.
  */
  struct domain D;
  PyObject *domain;

   
  double timestep;
  
  if (!PyArg_ParseTuple(args, "Od", &domain, &timestep)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }
    
  get_python_domain(&D,domain);

  timestep=_compute_fluxes_central(&D,timestep);

  // Return updated flux timestep
  return Py_BuildValue("d", timestep);
}


//PyObject *flux_function_central(PyObject *self, PyObject *args) {
//  //
//  // Gateway to innermost flux function.
//  // This is only used by the unit tests as the c implementation is
//  // normally called by compute_fluxes in this module.
//
//
//  PyArrayObject *normal, *ql, *qr,  *edgeflux;
//  double g, epsilon, max_speed, H0, zl, zr;
//  double h0, limiting_threshold, pressure_flux, smooth;
//  double hle, hre;
//
//  if (!PyArg_ParseTuple(args, "OOOddOddd",
//            &normal, &ql, &qr, &zl, &zr, &edgeflux,
//            &epsilon, &g, &H0)) {
//      report_python_error(AT, "could not parse input arguments");
//      return NULL;
//  }
//
//  printf("Warning: This interface is broken -- do not use \n");
//  
//  h0 = H0; 
//  hle=1.0; // Fake values to force this to compile
//  hre=1.0; // Fake values to force this to compile	      
//  limiting_threshold = 10*H0; // Avoid applying limiter below this
//                              // threshold for performance reasons.
//                              // See ANUGA manual under flux limiting  
// 
//  pressure_flux = 0.0; // Store the water pressure related component of the flux 
//  smooth = 1.0 ; // term to scale diffusion in riemann solver
//
//  _flux_function_central((double*) ql -> data, 
//			 (double*) qr -> data, 
//			 zl, 
//			 zr,
//             hle,
//             hre,                         
//			 ((double*) normal -> data)[0],
//			 ((double*) normal -> data)[1],          
//			 epsilon, h0, limiting_threshold,
//			 g,
//			 (double*) edgeflux -> data, 
//			 &max_speed,
//             &pressure_flux,
//			 ((double*) normal -> data)[0],
//			 ((double*) normal -> data)[1]
//             );
//  
//  return Py_BuildValue("d", max_speed);  
//}

//========================================================================
// Gravity
//========================================================================

PyObject *gravity(PyObject *self, PyObject *args) {
  //
  //  gravity(g, h, v, x, xmom, ymom)
  //
  
  
  PyArrayObject *h, *v, *x, *xmom, *ymom;
  int k, N, k3, k6;
  double g, avg_h, zx, zy;
  double x0, y0, x1, y1, x2, y2, z0, z1, z2;
  //double epsilon;
  
  if (!PyArg_ParseTuple(args, "dOOOOO",
			&g, &h, &v, &x,
			&xmom, &ymom)) {
    //&epsilon)) {
    PyErr_SetString(PyExc_RuntimeError, "shallow_water_ext.c: gravity could not parse input arguments");
    return NULL;
  }
  
  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(h);
  CHECK_C_CONTIG(v);
  CHECK_C_CONTIG(x);
  CHECK_C_CONTIG(xmom);
  CHECK_C_CONTIG(ymom);
  
  N = h -> dimensions[0];
  for (k=0; k<N; k++) {
    k3 = 3*k;  // base index
    
    // Get bathymetry
    z0 = ((double*) v -> data)[k3 + 0];
    z1 = ((double*) v -> data)[k3 + 1];
    z2 = ((double*) v -> data)[k3 + 2];
    
    // Optimise for flat bed
    // Note (Ole): This didn't produce measurable speed up.
    // Revisit later
    //if (fabs(z0-z1)<epsilon && fabs(z1-z2)<epsilon) {
    //  continue;
    //} 
    
    // Get average depth from centroid values
    avg_h = ((double *) h -> data)[k];
    
    // Compute bed slope
    k6 = 6*k;  // base index
    
    x0 = ((double*) x -> data)[k6 + 0];
    y0 = ((double*) x -> data)[k6 + 1];
    x1 = ((double*) x -> data)[k6 + 2];
    y1 = ((double*) x -> data)[k6 + 3];
    x2 = ((double*) x -> data)[k6 + 4];
    y2 = ((double*) x -> data)[k6 + 5];
    
    
    _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);
    
    // Update momentum
    ((double*) xmom -> data)[k] += -g*zx*avg_h;
    ((double*) ymom -> data)[k] += -g*zy*avg_h;
  }
  
  return Py_BuildValue("");
}

PyObject *compute_flux_update_frequency(PyObject *self, PyObject *args) {
  /*
    
    Compute how often we should update fluxes and extrapolations (perhaps not every timestep)

  */

  struct domain D;
  PyObject *domain;
  
    
  double timestep;
  int max_flux_update_frequency;
  
  if (!PyArg_ParseTuple(args, "Od", &domain, &timestep)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }
    
  get_python_domain(&D,domain);

  _compute_flux_update_frequency(&D, timestep);

  // Return 
  return Py_BuildValue("");
}


PyObject *extrapolate_second_order_edge_sw(PyObject *self, PyObject *args) {
  /*Compute the edge values based on a linear reconstruction 
    on each triangle
    
    Post conditions:
        The edges of each triangle have values from a 
        limited linear reconstruction
        based on centroid values

  */
 
  struct domain D; 
  PyObject *domain;

  int e;
  
  if (!PyArg_ParseTuple(args, "O", &domain)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }
  
  get_python_domain(&D, domain);

  // Call underlying flux computation routine and update
  // the explicit update arrays
  e = _extrapolate_second_order_edge_sw(&D);

  if (e == -1) {
    // Use error string set inside computational routine
    return NULL;
  }                   
  
  
  return Py_BuildValue("");
  
}// extrapolate_second-order_edge_sw

//========================================================================
// Protect -- to prevent the water level from falling below the minimum 
// bed_edge_value
//========================================================================

PyObject *protect(PyObject *self, PyObject *args) {
  //
  //    protect(minimum_allowed_height, maximum_allowed_speed, wc, zc, xmomc, ymomc)


  PyArrayObject
  *wc,            //Stage at centroids
  *wv,            //Stage at vertices
  *zc,            //Elevation at centroids
  *zv,            //Elevation at vertices
  *xmomc,         //Momentums at centroids
  *ymomc,
  *areas,         // Triangle areas
  *xc,
  *yc;

  int N;
  double mass_error;
  double minimum_allowed_height, maximum_allowed_speed, epsilon;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddOOOOOOOOO",
            &minimum_allowed_height,
            &maximum_allowed_speed,
            &epsilon,
            &wc, &wv, &zc, &zv, &xmomc, &ymomc, &areas, &xc, &yc)) {
    report_python_error(AT, "could not parse input arguments");
    return NULL;
  }  

  N = wc -> dimensions[0];

  mass_error = _protect(N,
       minimum_allowed_height,
       maximum_allowed_speed,
       epsilon,
       (double*) wc -> data,
       (double*) wv -> data,
       (double*) zc -> data,
       (double*) zv -> data,
       (double*) xmomc -> data,
       (double*) ymomc -> data,
       (double*) areas -> data,
       (double*) xc -> data,
       (double*) yc -> data );

  return Py_BuildValue("d", mass_error);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

PyObject *set_subgrid_volume_quantities_from_reference_quantities(PyObject *self, PyObject *args) {

  /*

    Update domain.quantities['vol'].centroid values
           domain.quantities['u_vol'].centroid values
           domain.quantities['v_vol'].centroid values

    from the reference quantities 'height', 'alphax', 'alphay'
  */

  struct domain D;
  PyObject *domain;

  int i,k,edge,var,inverse, number_of_elements;
  int use_last_lookup_info, use_last_cube_root;
  double h, uh;
  //double* lookup_table;

  if (!PyArg_ParseTuple(args, "O", &domain)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }
  
  get_python_domain(&D,domain);

  edge=-1; // Get centroid subgrid quantities
  inverse=0; // Go from reference_depth-->vol
   
  for(k=0; k< D.number_of_elements ; k++){

      // Reset the lookup: FIXME: Try to include this 'caching' inside eval_subgrid_function -- safer!
      use_last_lookup_info = 0;
      use_last_cube_root = 0;
    
      //lookup_table = &(D.subgrid_centroid_table[D.subgrid_centroid_starting_index[k]]);
      // Update height from stage/elevation, (just in case it is not yet updated)
      D.height_centroid_values[k] = D.stage_centroid_values[k] - D.bed_centroid_values[k];
      
      var = 0 ; // Get area 
      D.subgrid_wet_area[k] = eval_subgrid_function( D.stage_centroid_values[k],
                                                  k, edge, var, inverse, &D, 
                                                  use_last_lookup_info,
                                                  use_last_cube_root); //, lookup_table);
 
      // From here we can avoid recomputing lookup locations 
      use_last_lookup_info = 1;

      var = 1 ; // Get volume
      D.vol_centroid_values[k] = eval_subgrid_function( D.stage_centroid_values[k],
                                                        k, edge, var, inverse, &D,
                                                        use_last_lookup_info,
                                                        use_last_cube_root); //, lookup_table);

      var = 2 ; // Get 1/n h^(5/3) integral
      uh = eval_subgrid_function(D.stage_centroid_values[k],
                                 k, edge, var, inverse, &D,
                                 use_last_lookup_info,
                                 use_last_cube_root) ; //, lookup_table);

      // Multiply by 'alpha' values to get velocity volumes
      D.u_vol_centroid_values[k] = uh*(D.alphax_centroid_values[k]);
      D.v_vol_centroid_values[k] = uh*(D.alphay_centroid_values[k]);

  }

  return Py_BuildValue("");

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
PyObject *set_quantities_from_subgrid_volume_quantities(PyObject *self, PyObject *args) {

  /*

    Update centroids of 'stage', 'height', 'alphax', 'alphay' from the subgrid volume quantities

           domain.quantities['vol'].centroid values
           domain.quantities['u_vol'].centroid values
           domain.quantities['v_vol'].centroid values
  */

  struct domain D;
  PyObject *domain;

  int i,k,edge,var,inverse, number_of_elements;
  double h, uh, uh_inv, area_inv, zero_velocity_threshold;
  int use_last_lookup_info, use_last_cube_root;
  //double* lookup_table;
   
  if (!PyArg_ParseTuple(args, "O", &domain)) {
      report_python_error(AT, "could not parse input arguments");
      return NULL;
  }
  
  get_python_domain(&D,domain);

  edge=-1; // Get centroid subgrid quantities


  for(k=0; k< D.number_of_elements ; k++){

      // Reset the lookup
      use_last_lookup_info = 0;
      use_last_cube_root = 0;

      inverse=1; // Go from vol--> referencestage 
      var = 1 ; // Interpolate from volume to get stage

      D.stage_centroid_values[k] = eval_subgrid_function( D.vol_centroid_values[k],
                                                          k, edge, var, inverse, &D,
                                                          use_last_lookup_info,
                                                          use_last_cube_root);

      // Get height 
      D.height_centroid_values[k] = D.stage_centroid_values[k] - D.bed_centroid_values[k];

      // Need to update subgrid wet area -- this is not a 'conserved quantity'
      // so it still requires update
      var = 0; 
      inverse = 0; // This prevents us from reusing lookup info immediately
      D.subgrid_wet_area[k] = eval_subgrid_function( D.stage_centroid_values[k],
                                                  k, edge, var, inverse, &D, 
                                                 use_last_lookup_info, 
                                                 use_last_cube_root);

      if(D.vol_centroid_values[k] == 0.){
          // Quick exit without lookup
          uh = 0.;
          
      }else{

          // Get uh integral function to support alphax,alphay computation
          // u_vol = (1/n h^(5/3)_integral_function)*(alpha_x)
          var = 2; 
          inverse = 0; // This prevents us from reusing lookup info immediately
          use_last_lookup_info = 1; //From here we can avoid recomputing lookup index details
          uh = eval_subgrid_function(D.stage_centroid_values[k],
                                     k, edge, var, inverse, &D, 
                                     use_last_lookup_info, 
                                     use_last_cube_root);

      }

      // Back-calculate alphax, zero-ing if the 
      zero_velocity_threshold = 1.0e-10; //1.0e-03*D.areas[k];

      if(uh>(zero_velocity_threshold)){ 
          uh_inv=1./uh;
          //uh_inv=1./(uh+0.0e-06*D.areas[k]);
          D.alphax_centroid_values[k] = (D.u_vol_centroid_values[k])*uh_inv;
          D.alphay_centroid_values[k] = (D.v_vol_centroid_values[k])*uh_inv;

      }else{
          D.alphax_centroid_values[k] = 0.;
          D.alphay_centroid_values[k] = 0.;

      }

      // COSMETIC ONLY
      // Set xmom / ymom 
      if(D.subgrid_wet_area[k]>0.){
        area_inv = 1./D.subgrid_wet_area[k];
        D.xmom_centroid_values[k] = D.u_vol_centroid_values[k]*area_inv;
        D.ymom_centroid_values[k] = D.v_vol_centroid_values[k]*area_inv;

      }else{
        D.xmom_centroid_values[k] = 0.;
        D.ymom_centroid_values[k] = 0.;

      }
      for(i = 0; i < 3; i++){
         D.xmom_vertex_values[3*k+i] = D.xmom_centroid_values[k];
         D.ymom_vertex_values[3*k+i] = D.ymom_centroid_values[k];

      }

  }

  return Py_BuildValue("");

}


//========================================================================
// Method table for python module
//========================================================================

static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */
  //{"rotate", (PyCFunction)rotate, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {"compute_fluxes_ext_central", compute_fluxes_ext_central, METH_VARARGS, "Print out"},
  {"gravity_c",        gravity,            METH_VARARGS, "Print out"},
  //{"flux_function_central", flux_function_central, METH_VARARGS, "Print out"},  
  {"extrapolate_second_order_edge_sw", extrapolate_second_order_edge_sw, METH_VARARGS, "Print out"},
  {"compute_flux_update_frequency", compute_flux_update_frequency, METH_VARARGS, "Print out"},
  {"set_subgrid_volume_quantities_from_reference_quantities", set_subgrid_volume_quantities_from_reference_quantities, METH_VARARGS, "Print out"},
  {"set_quantities_from_subgrid_volume_quantities", set_quantities_from_subgrid_volume_quantities, METH_VARARGS, "Print out"},
  {"protect",          protect, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {NULL, NULL, 0, NULL}
};

// Module initialisation
void initswSG1_domain_ext(void){
  Py_InitModule("swSG1_domain_ext", MethodTable);

  import_array(); // Necessary for handling of NumPY structures
}
