#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_isepdgt // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                    c=comp_idepdgt(coef,g,L,a,M,phasetype);\n\
                    Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer

static inline void
fwd_idgt_fb(const Complex *coef, const Complex *gf,
            const octave_idx_type L, const octave_idx_type gl,
            const octave_idx_type W, const octave_idx_type a,
            const octave_idx_type M, const octave_idx_type ptype,
            Complex *f)
{
    idgt_fb_d(reinterpret_cast<const fftw_complex*>(coef),
              reinterpret_cast<const fftw_complex*>(gf),
              L, gl, W, a, M,
              static_cast<dgt_phasetype>(ptype),
              reinterpret_cast<fftw_complex*>(f));
}

static inline void
fwd_idgt_fb(const FloatComplex *coef, const FloatComplex *gf,
            const octave_idx_type L, const octave_idx_type gl,
            const octave_idx_type W, const octave_idx_type a,
            const octave_idx_type M, const octave_idx_type ptype,
            FloatComplex *f)
{
    idgt_fb_s(reinterpret_cast<const fftwf_complex*>(coef),
              reinterpret_cast<const fftwf_complex*>(gf),
              L, gl, W, a, M,
              static_cast<dgt_phasetype>(ptype),
              reinterpret_cast<fftwf_complex*>(f));
}

static inline void
fwd_idgt_long(const Complex *coef, const Complex *gf,
              const octave_idx_type L, const octave_idx_type W,
              const octave_idx_type a, const octave_idx_type M,
              const octave_idx_type ptype, Complex *f)
{
    idgt_long_d(reinterpret_cast<const fftw_complex*>(coef),
                reinterpret_cast<const fftw_complex*>(gf),
                L, W, a, M, 
                static_cast<dgt_phasetype>(ptype),
                reinterpret_cast<fftw_complex*>(f));
}

static inline void
fwd_idgt_long(const FloatComplex *coef, const FloatComplex *gf,
              const octave_idx_type L, const octave_idx_type W,
              const octave_idx_type a, const octave_idx_type M,
              const octave_idx_type ptype, FloatComplex *f)
{
    idgt_long_s(reinterpret_cast<const fftwf_complex*>(coef),
                reinterpret_cast<const fftwf_complex*>(gf),
                L, W, a, M, 
                static_cast<dgt_phasetype>(ptype),
                reinterpret_cast<fftwf_complex*>(f));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    MArray<LTFAT_TYPE> coef = ltfatOctArray<LTFAT_TYPE>(args(0));
    MArray<LTFAT_TYPE> gf = ltfatOctArray<LTFAT_TYPE>(args(1));
    const octave_idx_type L = args(2).int_value();
    const octave_idx_type a = args(3).int_value();
    const octave_idx_type M = args(4).int_value();
    const octave_idx_type ptype = args(5).int_value();
    const octave_idx_type N = L / a;
    const octave_idx_type gl = gf.rows();

    const octave_idx_type W = coef.numel() / (M * N);

    MArray<LTFAT_COMPLEX> f(dim_vector(L, W));

    if (gl < L)
    {
        fwd_idgt_fb(coef.data(), gf.data(), L, gl, W, a, M, ptype,
                    f.fortran_vec());
    }
    else
    {
        fwd_idgt_long(coef.data(), gf.data(), L, W, a, M, ptype,
                      f.fortran_vec());
    }
    return octave_value(f);
}
