/* This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license. */

// enable double precision (not enabled by default)

#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
    #error "IEEE-754 double precision not supported by OpenCL implementation."
#endif

kernel void MM(const global uint *A,const global uint *B,global bool *C)
{

  // Thread identifiers
    const uint globalRow = get_global_id(0); // Row ID of C (0..M)
    const uint globalCol = get_global_id(1); // Col ID of C (0..N)

	const uint num_rows=(uint)NR;
	const uint num_cols=(uint)NC;
	const uint num_i=(uint)NI;

    // Compute a single element (loop over K)
    int acc = 0;

uint n;

    for (uint k=0; k<num_i; k++) {
        n = A[globalRow*num_i + k]^B[globalCol*num_i + k];
    
        while(n) {
            if(n & 1) acc++;
            n >>= 1;
        }
       // acc += popcount(A[globalRow*num_i + k]^B[globalCol*num_i + k]);
        //if(A[globalRow*num_i + k]^B[globalCol*num_i + k])  acc++;
        

        if(acc>2) break;
    }

//        acc += A[k*num_rows + globalRow] * B[k*num_rows + globalCol];//B[globalCol*num_i + k];
// Store the result


C[globalCol*num_rows + globalRow] = false;
if(acc ==2) C[globalCol*num_rows + globalRow] = true;
//if (acc == suma-1) C[globalCol*num_rows + globalRow] = true;
//if (sumb == suma) C[globalCol*num_rows + globalRow] = true;

        



};
