#define TYPEDEPARGS 0
#define MATCHEDARGS 1, 2
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_filterbankreassign // change to filename
#define OCTFILEHELP "Computes reassigned filterbank spectrogram.\n\
Usage: sr=comp_filterbankreassign(s,itime,ifreq,a,cfreq);\n\
Yeah."


#include "ltfat_oct_template_helper.h"
/*
   gabreassign forwarders
   */

    static inline void
fwd_filterbankreassign(const double *s[], const double *tgrad[], const double *fgrad[],
        const ltfatInt N[], const double a[], const double cfreq[],
        const ltfatInt M, double *sr[], fbreassOptOut* optout)
{
    filterbankreassign_d(s, tgrad, fgrad, N, a, cfreq, M, sr, REASS_DEFAULT, optout);
}

    static inline void
fwd_filterbankreassign(const float *s[], const float *tgrad[], const float *fgrad[],
        const ltfatInt N[], const double a[], const double cfreq[],
        const ltfatInt M, float *sr[], fbreassOptOut* optout)
{
    filterbankreassign_s(s, tgrad, fgrad, N, a, cfreq, M, sr, REASS_DEFAULT, optout);
}

    static inline void
fwd_filterbankreassign(const Complex *s[], const double *tgrad[], const double *fgrad[],
        const ltfatInt N[], const double a[], const double cfreq[],
        const ltfatInt M, Complex *sr[], fbreassOptOut* optout)
{
    filterbankreassign_cd(reinterpret_cast<const fftw_complex **>(s),
            tgrad, fgrad, N, a, cfreq, M,
            reinterpret_cast<fftw_complex **>(sr),
            REASS_DEFAULT, optout);
}

    static inline void
fwd_filterbankreassign(const FloatComplex *s[], const float *tgrad[], const float *fgrad[],
        const ltfatInt N[], const double a[], const double cfreq[],
        const ltfatInt M, FloatComplex *sr[], fbreassOptOut* optout)
{
    filterbankreassign_cs(reinterpret_cast<const fftwf_complex **>(s),
            tgrad, fgrad, N, a, cfreq, M,
            reinterpret_cast<fftwf_complex **>(sr),
            REASS_DEFAULT, optout);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
    octave_value_list
octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;
    Cell sCell      = args(0).cell_value();
    Cell tgradCell  = args(1).cell_value();
    Cell fgradCell  = args(2).cell_value();
    Matrix aMat     = args(3).matrix_value();
    Matrix cfreqMat = args(4).matrix_value();

    const ltfatInt M  = (ltfatInt) sCell.numel();


    OCTAVE_LOCAL_BUFFER (ltfatInt, NPtr, M);
    OCTAVE_LOCAL_BUFFER (double, aPtr, M);
    OCTAVE_LOCAL_BUFFER (const LTFAT_TYPE*, sPtr, M);
    OCTAVE_LOCAL_BUFFER (const LTFAT_REAL*, tgradPtr, M);
    OCTAVE_LOCAL_BUFFER (const LTFAT_REAL*, fgradPtr, M);
    OCTAVE_LOCAL_BUFFER (LTFAT_TYPE*, srPtr, M);

    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, s_elems, M);
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_REAL>, tgrad_elems, M);
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_REAL>, fgrad_elems, M);
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_TYPE>, sr_elems, M);

    for (octave_idx_type m = 0; m < M; ++m)
    {
        s_elems[m] = ltfatOctArray<LTFAT_TYPE>(sCell.elem(m));
        tgrad_elems[m] = ltfatOctArray<LTFAT_REAL>(tgradCell.elem(m));
        fgrad_elems[m] = ltfatOctArray<LTFAT_REAL>(fgradCell.elem(m));

        sPtr[m] = s_elems[m].data();
        tgradPtr[m] = tgrad_elems[m].data();
        fgradPtr[m] = fgrad_elems[m].data();
        NPtr[m] = s_elems[m].rows();

        sr_elems[m] = MArray<LTFAT_TYPE>(dim_vector(NPtr[m], 1));
        srPtr[m] = sr_elems[m].fortran_vec();
        aPtr[m] = aMat(m);
    }
    octave_value_list retlist;

    fbreassOptOut* optout = NULL;
    Cell repos;
    octave_idx_type sumN = 0;
    if (nargout > 1)
    {
        for (octave_idx_type ii = 0; ii < M; ii++) sumN += NPtr[ii];

        repos = Cell(dim_vector(sumN, 1));
        optout = fbreassOptOut_init(sumN, 16 );
    }

    // Adjust a
    if (aMat.columns() > 1)
    {
        for (octave_idx_type m = 0; m < M; ++m)
        {
            aPtr[m] /= aMat(m + M);
        }
    }

    fwd_filterbankreassign(sPtr, tgradPtr, fgradPtr, NPtr, aPtr, cfreqMat.data(), M,
            srPtr, optout);

    // Output cell
    Cell srCell(dim_vector(M, 1));
    for (octave_idx_type m = 0; m < M; ++m)
    {
        srCell.elem(m) = sr_elems[m];
    }
    retlist(0) = srCell;


    if (nargout > 1)
    {
        for (ltfatInt ii = 0; ii < sumN; ii++)
        {
            octave_idx_type l = optout->reposl[ii];
            MArray<double> cEl = MArray<double>(dim_vector(l, l > 0 ? 1 : 0));
            double* cElPtr = cEl.fortran_vec();
            for (octave_idx_type jj = 0; jj < optout->reposl[ii]; jj++)
            {
                cElPtr[jj] = static_cast<double>(optout->repos[ii][jj]) + 1.0;
            }
            repos.elem(ii) = cEl;
        }

        retlist(1) = repos;
        fbreassOptOut_destroy(optout);
    }
    return retlist;
}
