namespace nt.math.statistics
{
    using System;
    using System.Collections;
    using nt.math;

    public enum ARMA_FLAGS
    {
        X12A = 1 << 0,
        EXACT = 1 << 1,
        LS = 1 << 2,
        OLS = 1 << 3
    }

    public enum ARMA_PRIV_FLAGS
    {
        SEAS = 1 << 0,
        DSPEC = 1 << 1,
        XDIFF = 1 << 2,
        LBFGS = 1 << 3,
        VECH = 1 << 4,
        NAOK = 1 << 5,
        NAS = 1 << 6,
        LEV = 1 << 7,
        YDIFF = 1 << 8,
        AVGLL = 1 << 9
    }

    public struct BChecker
    {
        public int qmax;
        public double[] tmp;
        public double[] tmp2;
        public Complex[] roots;

        public BChecker(Arma p)
        {
            qmax = p.q + p.Q * p.pd;
            tmp = new double[qmax + 1];
            tmp2 = new double[qmax + 1];
            roots = new Complex[qmax];
        }
    }

    public class Arma
    {
        public ARMA_FLAGS flags;
        public ARMA_PRIV_FLAGS pflags;

        public char[] pmask;
        public char[] qmask;

        public int p;
        public int d;
        public int q;

        public int P;
        public int D;
        public int Q;

        public int np;
        public int nq;

        public int pd;

        public int ifc; // специальная константа

        public Arma()
        {
            flags = ARMA_FLAGS.X12A;
            pflags = ARMA_PRIV_FLAGS.SEAS;

            p = 0;
            d = 0;
            q = 0;
            P = 0;
            D = 0;
            Q = 0;

            np = 0;
            nq = 0;

            ifc = 0;
        }

        public unsafe int AddRoots(double[] coeff)
        {
            fixed (double* pcoeff = coeff)
            {
                double* phi = pcoeff + this.ifc;
                double* Phi = phi + this.np;
                double* theta = Phi + this.P;
                double* Theta = theta + this.nq;

                int nr = this.p + this.P + this.q + this.Q;

                int pmax = (this.p > this.P) ? this.p : this.P;
                int qmax = (this.q > this.Q) ? this.q : this.Q;
                int lmax = (pmax > qmax) ? pmax : qmax;

                if (pmax == 0 || qmax == 0) return 0;

                double* tmp = stackalloc double[lmax + 1];
                double* tmp2 = stackalloc double[lmax + 1];

                Complex* roots = stackalloc Complex[nr];
                Complex* rptr;

                int cerr = 0;

                tmp[0] = 1.0;
                rptr = roots;

                if (this.p > 0)
                {
                    /* A(L), non-seasonal */
                    int k = 0;

                    for (int i = 0; i < this.p; i++)
                    {
                        if (AR_Included(this, i))
                        {
                            tmp[i + 1] = -phi[k++];
                        }
                        else
                        {
                            tmp[i + 1] = 0.0;
                        }
                    }

                    cerr = Algebra.PolyRoots(tmp, tmp2, this.p, rptr);
                    rptr += this.p;
                }

                if (cerr != 0 && this.P > 0)
                {
                    /* B(L), seasonal */
                    for (int i = 0; i < this.P; i++)
                    {
                        tmp[i + 1] = -Phi[i];
                    }

                    cerr = Algebra.PolyRoots(tmp, tmp2, this.P, rptr);
                    rptr += this.P;
                }

                if (cerr != 0 && this.q > 0)
                {
                    /* C(L), non-seasonal */
                    int k = 0;

                    for (int i = 0; i < this.q; i++)
                    {
                        if (MA_Included(this, i))
                        {
                            tmp[i + 1] = theta[k++];
                        }
                        else
                        {
                            tmp[i + 1] = 0.0;
                        }
                    }

                    cerr = Algebra.PolyRoots(tmp, tmp2, this.q, rptr);
                    rptr += this.q;
                }

                if (cerr != 0 && this.Q > 0)
                {
                    /* D(L), seasonal */
                    for (int i = 0; i < this.Q; i++)
                    {
                        tmp[i + 1] = Theta[i];
                    }

                    cerr = Algebra.PolyRoots(tmp, tmp2, this.Q, rptr);
                }

                tmp = null;
                tmp2 = null;

                if (cerr == 0) 
                {
	                roots = null;
                } 
                else 
                {
	                //gretl_model_set_data(pmod, "roots", roots, GRETL_TYPE_CMPLX_ARRAY, nr * sizeof *roots);
                }

                return 0;
            }
        }

        #region Other function
        public int MAOutOfBounds(double[] theta, double[] Theta)
        {
            BChecker b = new BChecker(this);
            bool tzero = true;
            bool Tzero = true;

            int k = 0;

            for(int i = 0; i < this.q && tzero; i++)
            {
                if(MA_Included(this, i))
                {
                    if (theta[k++] != 0.0)
                        tzero = false;
                }
            }

            for (int i = 0; i < this.Q && Tzero; i++)
            {
                if (Theta[k++] != 0.0)
                    Tzero = false;
            }

            if (tzero && Tzero) return 0;

            b.tmp[0] = 1.0;

            // init non-seasonal MA or zero
            for(int i = 0; i < b.qmax; i++)
            {
                if(i < this.q && MA_Included(this, i))
                {
                    b.tmp[i + 1] = theta[k++];
                }
                else
                {
                    b.tmp[i + 1] = 0.0;
                }
            }

            int qtot = 0;

            if(Tzero)
            {
                qtot = this.q;
            }
            else
            {
                for(int i = 0; i < this.Q; i++)
                {
                    int si = (i + 1) * this.pd;

                    b.tmp[si] += Theta[i];
                    k = 0;
                    for (i = 0; i < this.q; i++)
                    {
                        if (MA_Included(this, i))
                        {
                            int m = si + (i + 1);

                            b.tmp[m] += Theta[i] * theta[k++];
                        }
                    }
                }
            }

            int cerr = Algebra.PolyRoots(b.tmp, b.tmp2, qtot, b.roots);
            int err = 0;

            for (int i = 0; i < qtot; i++)
            {
                double re = b.roots[i].Re;
                double im = b.roots[i].Im;
                double rt = re * re + im * im;

                if (rt > Constants.DoubleEpsilon && rt <= 1.0)
                {
                    //pprintf(ainfo->prn, _("MA root %d = %g\n"), i, rt);
                    err = 1;
                    break;
                }
            }

            return err;
        }
        #endregion

        #region Static functions
        public static ARMA_FLAGS Arma_By_X12A(Arma p)
        {
            return p.flags & ARMA_FLAGS.X12A;
        }

        public static ARMA_FLAGS Arma_Exact_ML(Arma p)
        {
            return p.flags & ARMA_FLAGS.EXACT;
        }

        public static ARMA_FLAGS Arma_Least_Squares(Arma p)
        {
            return p.flags & ARMA_FLAGS.LS;
        }

        public static void SetArma_Least_Squares(Arma p)
        {
            p.flags |= ARMA_FLAGS.LS;
        }

        public static ARMA_PRIV_FLAGS Arma_Has_Seasonal(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.SEAS;
        }

        public static ARMA_PRIV_FLAGS Arma_Is_Arima(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.DSPEC;
        }

        public static ARMA_PRIV_FLAGS Arma_Is_XDiff(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.XDIFF;
        }

        public static ARMA_PRIV_FLAGS Arma_Is_LBFGS(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.LBFGS;
        }

        public static ARMA_PRIV_FLAGS Arma_Using_Vech(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.VECH;
        }

        public static ARMA_PRIV_FLAGS Arma_Na_Ok(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.NAOK;
        }

        public static ARMA_PRIV_FLAGS Arma_Missvals(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.NAS;
        }

        public static ARMA_PRIV_FLAGS Arma_Levels(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.LEV;
        }

        public static ARMA_PRIV_FLAGS Arma_YDiff(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.YDIFF;
        }

        public static ARMA_PRIV_FLAGS Arma_AVGLL(Arma p)
        {
            return p.pflags & ARMA_PRIV_FLAGS.AVGLL;
        }

        public static bool AR_Included(Arma p, int i)
        {
            return (p.pmask == null || p.pmask[i] == '1');
        }

        public static bool MA_Included(Arma p, int i)
        {
            return (p.qmask == null || p.qmask[i] == '1');
        }
        #endregion
    }
}
