#ifndef REFLECTIONSYMMETRY_H
#define REFLECTIONSYMMETRY_H

#include "cctk.h"
#include "cctk_Arguments.h"

void ReflectionSymmetry_Apply (CCTK_ARGUMENTS);

CCTK_INT
ReflectionSymmetry_Interpolate (CCTK_POINTER_TO_CONST restrict const cctkGH,
                                CCTK_INT const N_dims,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                CCTK_INT const coord_system_handle,
                                CCTK_INT const N_interp_points,
                                CCTK_INT const interp_coords_type,
                                CCTK_POINTER_TO_CONST restrict const interp_coords[],
                                CCTK_INT const N_input_arrays,
                                CCTK_INT const input_array_indices[],
                                CCTK_INT const N_output_arrays,
                                CCTK_INT const output_array_types[],
                                CCTK_POINTER restrict const output_arrays[],
                                CCTK_INT const faces);

void 
ReflectionSymmetry_GetManualParities(int table, int gi, int *parities);

#endif /* ! defined REFLECTIONSYMMETRY_H */
