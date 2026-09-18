/* =====================================================================
 * resasm_mat_abi.h  --  versioned C ABI for OTI-enabled replay materials
 * =====================================================================
 *
 * This header is the STABLE, compiler-independent boundary between a
 * collaborator's residual/sensitivity tool and a material-model developer's
 * (e.g. JHU's) closed-source constitutive binary.
 *
 * DESIGN RULES (do not violate; the whole privacy model depends on them):
 *   1. No OTI / dual / hypercomplex objects ever cross this boundary. The
 *      binary keeps OTI types INTERNAL and exposes only plain doubles: the
 *      real value AND, separately, the first-order derivative coefficient for
 *      each active seed direction.
 *   2. A `double*` NEVER "contains an OTI number". `stress[]` holds the real
 *      stress; `dstress_dseed[]` holds d(stress)/d(seeded parameter). They are
 *      distinct arrays with documented layout.
 *   3. All dimensions are passed explicitly in `resasm_mat_desc_t`; the callee
 *      must validate them and must not read past them.
 *   4. Model metadata (mat_describe_v1) is separate from evaluation
 *      (mat_eval_v1).
 *   5. Every array is caller-allocated. Ownership stays with the caller; the
 *      callee only writes into the provided `*_out` buffers.
 *   6. Failures are explicit: a non-zero return code means the outputs are
 *      NOT valid. The callee must not return partially-valid results.
 *   7. The exported symbol is versioned (mat_eval_v1). A future incompatible
 *      contract is mat_eval_v2, never a silent change to v1.
 *
 * ARRAY LAYOUT (row-major, 0-based flattening):
 *   stress[i]                 : i in [0,ntens)
 *   dstress_dseed[i*nseed + s]: d stress_i / d props[seed_indices[s]]
 *   ddsdde[i*ntens + j]       : d stress_i / d dstrain_j  (consistent tangent)
 *   state_out[k]              : k in [0,nstatev)
 *   dstate_dseed_out[k*nseed + s] : d state_k / d props[seed_indices[s]]
 *
 * KINEMATICS coding (desc.kinematics):
 *   RESASM_KIN_SMALL_STRAIN (0): kin holds the strain increment (or total
 *       small strain) as a length-`ntens` Voigt vector, order
 *       [11,22,33,12,13,23], engineering shear.
 *   RESASM_KIN_FINITE_STRAIN(1): kin holds [F0(9), F1(9)] row-major 3x3
 *       deformation gradients at increment start and end.
 *
 * SEEDING: seed_indices[s] is the 1-BASED index into props of the parameter
 *   that OTI direction s differentiates; nseed is how many are active this
 *   call. Multiple directions may be seeded simultaneously.
 *
 * PATH DEPENDENCE: for history-dependent models the caller replays increment
 *   by increment, passing the previous increment's state_out/dstate_dseed_out
 *   back in as state_in/dstate_dseed_in. Optional-input pointers that are not
 *   used (small-strain elasticity uses none) may be NULL; the callee MUST test
 *   for NULL before dereferencing.
 * ===================================================================== */
#ifndef RESASM_MAT_ABI_H
#define RESASM_MAT_ABI_H

#ifdef __cplusplus
extern "C" {
#endif

#define RESASM_MAT_ABI_VERSION 1

/* kinematics codes */
#define RESASM_KIN_SMALL_STRAIN 0
#define RESASM_KIN_FINITE_STRAIN 1

/* status / return codes (0 == success; outputs valid only on 0) */
#define RESASM_MAT_OK              0
#define RESASM_MAT_ERR_ABI_VERSION 1   /* desc.abi_version != RESASM_MAT_ABI_VERSION */
#define RESASM_MAT_ERR_DIMS        2   /* a dimension is non-positive / inconsistent */
#define RESASM_MAT_ERR_NULL        3   /* a required pointer was NULL */
#define RESASM_MAT_ERR_SEED        4   /* a seed index is out of [1,nprops] */
#define RESASM_MAT_ERR_KINEMATICS  5   /* desc.kinematics not supported by this binary */
#define RESASM_MAT_ERR_STATE       6   /* state layout / nstatev mismatch */
#define RESASM_MAT_ERR_CONVERGENCE 7   /* constitutive update failed to converge */
#define RESASM_MAT_ERR_INTERNAL    8   /* any other internal failure */

/* Model metadata / dimensions. Kept separate from evaluation data so a caller
 * can size its buffers and validate the twin BEFORE evaluating anything. */
typedef struct resasm_mat_desc {
    int abi_version;   /* must equal RESASM_MAT_ABI_VERSION */
    int ntens;         /* number of stress/strain components (6 for 3D) */
    int nprops;        /* number of material parameters in props[] */
    int nstatev;       /* number of state variables (0 for elasticity) */
    int kinematics;    /* one of RESASM_KIN_* */
    int order;         /* OTI order exposed (1 for first-order sensitivity) */
} resasm_mat_desc_t;

/* Fill desc_out with this binary's fixed dimensions/metadata, and (optionally)
 * copy up to model_id_cap-1 chars of the model id into model_id (NUL-terminated).
 * Returns RESASM_MAT_OK or an error code. Does NOT evaluate the model. */
int mat_describe_v1(resasm_mat_desc_t *desc_out,
                    char *model_id, int model_id_cap);

/* Evaluate the constitutive model at one material point for one increment,
 * returning the real response AND the first-order derivative coefficients for
 * every seeded parameter direction. Returns RESASM_MAT_OK on success; on any
 * non-zero return the *_out buffers are undefined. `status` receives the same
 * code (or a finer internal detail) and may be NULL. */
int mat_eval_v1(
    const resasm_mat_desc_t *desc,
    const double *props,             /* [nprops]                 (required) */
    const int    *seed_indices,      /* [nseed] 1-based          (required)  */
    int           nseed,
    const double *kin,               /* [ntens] or [18]          (required)  */
    const double *dkin_dseed,        /* [len(kin)*nseed] or NULL (optional)  */
    const double *state_in,          /* [nstatev] or NULL        (optional)  */
    const double *dstate_dseed_in,   /* [nstatev*nseed] or NULL  (optional)  */
    const double *time_data,         /* [4]={time,dtime,temp,dtemp} (required)*/
    double       *stress,            /* out [ntens]              (required)  */
    double       *dstress_dseed,     /* out [ntens*nseed]        (required)  */
    double       *state_out,         /* out [nstatev] or NULL    (optional)  */
    double       *dstate_dseed_out,  /* out [nstatev*nseed] or NULL (optional)*/
    double       *ddsdde,            /* out [ntens*ntens]        (required)  */
    int          *status);           /* out                      (optional)  */

#ifdef __cplusplus
}
#endif
#endif /* RESASM_MAT_ABI_H */
