#define TYPEDEPARGS 0, 1, 2, 7
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_maskedheapint
#define OCTFILEHELP "Computes masked heapint.\n\
Usage: c = comp_maskedheapint(s, itime, ifreq, mask, a, tol, do_timeinv, usephase);\n\
Yeah."


#include "ltfat_oct_template_helper.h"

static inline void
fwd_maskedheapint(const double *s, const double *tgrad, const double *fgrad,
                  const int* mask, const octave_idx_type a, const octave_idx_type M,
                  const octave_idx_type L, const octave_idx_type W,
                  const double tol, int phasetype, double *phase)
{
    if (phasetype == 2)
        maskedheapint_d( s, tgrad, fgrad, mask, a, M, L, W, tol, phase);
    else
        maskedheapint_relgrad_d( s, tgrad, fgrad, mask, a, M, L, W, tol,
                static_cast<dgt_phasetype>(phasetype), phase);
}

static inline void
fwd_maskedheapint(const float *s, const float *tgrad, const float *fgrad,
                  const int* mask, const octave_idx_type a, const octave_idx_type M,
                  const octave_idx_type L, const octave_idx_type W,
                  const float tol, int phasetype, float *phase)
{
    if (phasetype == 2)
        maskedheapint_s(s, tgrad, fgrad, mask, a, M, L, W, tol, phase);
    else
        maskedheapint_relgrad_s(s, tgrad, fgrad, mask, a, M, L, W, tol,
                static_cast<dgt_phasetype>(phasetype), phase);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    MArray<LTFAT_REAL> s = ltfatOctArray<LTFAT_REAL>(args(0));
    MArray<LTFAT_REAL> tgrad = ltfatOctArray<LTFAT_REAL>(args(1));
    MArray<LTFAT_REAL> fgrad = ltfatOctArray<LTFAT_REAL>(args(2));
    MArray<double> maskDouble  = ltfatOctArray<double>(args(3));
    const octave_idx_type a  = args(4).int_value();
    const double tol  = args(5).double_value();
    const int phasetype  = args(6).int_value();
    MArray<LTFAT_REAL> usephase = ltfatOctArray<LTFAT_REAL>(args(7));

    const octave_idx_type M = args(0).rows();
    const octave_idx_type N = args(0).columns();
    const octave_idx_type L = N * a;
    const octave_idx_type W = s.numel() / (M * N);

    MArray<LTFAT_REAL> phase(dim_vector(M, N, W));

    int* mask = new int[M * N * W];
    for (octave_idx_type w = 0; w < M * N * W;++w)
        mask[w] = (int) maskDouble.data()[w];

    memcpy(phase.fortran_vec(), usephase.data(), M * N * W * sizeof(LTFAT_REAL));

    fwd_maskedheapint(s.data(), tgrad.data(), fgrad.data(),
                      mask, a, M, L, W, static_cast<LTFAT_REAL>(tol),
                      phasetype, phase.fortran_vec());

    delete [] mask;
    return octave_value(phase);
}
