namespace nt.math
{
    using System;
    using nt.math;

    public static partial class Algebra
    {
        public static int PolyRoots(double[] xcof, double[] cof, int m, Complex[] root)
        {
            int n = m;
            bool final = false;

            if (n <= 0)
            {
                return 1;
            }

            if (n > 36)
            {
                return 2;
            }

            if (xcof[m] == 0.0)
            {
                return 4;
            }

            int n1 = n;
            int n2 = n;
            int nroot = 0;
            int nsav = n;

            for (int i = 0; i <= nsav; i++)
            {
                cof[n - i] = xcof[i];
            }

            Complex xsav = new Complex(0.0, 0.0);
            Complex x0 = new Complex();
            Complex x = new Complex();
            Complex u = new Complex();
            Complex ud = new Complex();
            Complex t = new Complex();
            Complex t1 = new Complex();
            Complex dx = new Complex();

            int retry, iter;
            double cofj, mag;

        nxtrut:
            x0.Re = 0.00500101;
            x0.Im = 0.01000101;

            retry = 0;

        tryagn:
            retry++;
            x.Re = x0.Re;

            x0.Re = -10.0 * x0.Im;
            x0.Im = -10.0 * x.Re;

            x.Re = x0.Re;
            x.Im = x0.Im;

        finitr:
            iter = 0;

            while (iter < 500)
            {
                u.Re = cof[n];

                if (u.Re == 0.0)
                {
                    x.Re = 0.0;
                    n1--;
                    n2--;
                    goto zerrut;
                }


                u.Im = 0.0;
                ud.Re = 0.0;
                ud.Im = 0.0;
                t.Re = 1.0;
                t.Im = 0.0;

                for (int i = 0; i < n; i++)
                {
                    t1.Re = x.Re * t.Re - x.Im * t.Im;
                    t1.Im = x.Re * t.Im + x.Im * t.Re;

                    cofj = cof[(n - 1) - i];

                    u.Re += cofj * t1.Re;
                    u.Im += cofj * t1.Im;

                    cofj = cofj * (i + 1); /* derivative */

                    ud.Re += cofj * t.Re;
                    ud.Im -= cofj * t.Im;

                    t.Re = t1.Re;
                    t.Im = t1.Im;
                }

                mag = ud.Re * ud.Re + ud.Im * ud.Im;

                if (mag == 0.0)
                {
                    if (!final)
                    {
                        if (retry < 50)
                        {

                            goto tryagn;
                        }
                        else
                        {
                            break;
                        }
                    }

                    x.Re = xsav.Re;
                    x.Im = xsav.Im;
                    goto findon;
                }

                dx.Re = (u.Im * ud.Im - u.Re * ud.Re) / mag;
                x.Re += dx.Re;

                dx.Im = -(u.Re * ud.Im + u.Im * ud.Re) / mag;
                x.Im += dx.Im;

                if ((Math.Abs(dx.Im) + Math.Abs(dx.Re)) < 1.0e-6)
                    goto lupdon;

                iter++;
            }

            if (final) goto lupdon;
            if (retry < 5) goto tryagn;

            return 3;

        lupdon:
            /* Swap original and reduced polynomials */

            for (int i = 0; i <= n2; i++)
            {
                cofj = xcof[nsav - i];
                xcof[(nsav - 1) - i] = cof[i];
                cof[i] = cofj;
            }

            int k = n;
            n = n1;
            n1 = k;

            if (!final)
            {
                final = true;

                if (Math.Abs(x.Im / x.Re) < 1.0e-4) x.Im = 0.0;

                xsav.Re = x.Re;
                xsav.Im = x.Im;
                goto finitr; /* do final iteration on original polynomial */
            }

        findon:
            final = false;

        zerrut:
            if (Math.Abs(x.Im / x.Re) >= 1.0e-5 && u.Re != 0.0)
            {
                cofj = x.Re + x.Re;
                mag = x.Re * x.Re + x.Im * x.Im;
                n -= 2;
            }
            else
            {
                x.Im = 0.0;
                cofj = x.Re;
                mag = 0.0;
                n--;
            }

            for (int i = 1; i < n; i++)
            {
                cof[i + 1] += cofj * cof[i] - mag * cof[i - 1];
            }

        setrut:
            root[nroot].Re = x.Re;
            root[nroot].Im = x.Im;
            nroot++;

            if (mag != 0.0)
            {
                x.Im = -x.Im;
                mag = 0;
                goto setrut;	/* fill in the complex conjugate root */
            }

            if (n > 0)
                goto nxtrut;

            return 0;
        }

        public static unsafe int PolyRoots(double* xcof, double* cof, int m, Complex* root)
        {
            int n = m;
            bool final = false;

            if (n <= 0)
            {
                return 1;
            }

            if (n > 36)
            {
                return 2;
            }

            if (xcof[m] == 0.0)
            {
                return 4;
            }

            int n1 = n;
            int n2 = n;
            int nroot = 0;
            int nsav = n;
            double* q = &xcof[0];
            double* p = &cof[n];

            for (int i = 0; i <= nsav; i++)
            {
                *p-- = *q++;
            }

            Complex xsav = new Complex(0.0, 0.0);
            Complex x0 = new Complex();
            Complex x = new Complex();
            Complex u = new Complex();
            Complex ud = new Complex();
            Complex t = new Complex();
            Complex t1 = new Complex();
            Complex dx = new Complex();

            int retry, iter;
            double cofj, mag;

        nxtrut:
            x0.Re = 0.00500101;
            x0.Im = 0.01000101;

            retry = 0;

        tryagn:
            retry += 1;
            x.Re = x0.Re;

            x0.Re = -10.0 * x0.Im;
            x0.Im = -10.0 * x.Re;

            x.Re = x0.Re;
            x.Im = x0.Im;

        finitr:
            iter = 0;

            while (iter < 500)
            {
                u.Re = cof[n];

                if (u.Re == 0.0)
                {
                    x.Re = 0.0;
                    n1 -= 1;
                    n2 -= 1;
                    goto zerrut;
                }


                u.Im = 0.0;
                ud.Re = 0.0;
                ud.Im = 0.0;
                t.Re = 1.0;
                t.Im = 0.0;

                p = &cof[n - 1];
 
                for (int i = 0; i < n; i++)
                {
                    t1.Re = x.Re * t.Re - x.Im * t.Im;
                    t1.Im = x.Re * t.Im + x.Im * t.Re;

                    cofj = *p--;

                    u.Re += cofj * t1.Re;
                    u.Im += cofj * t1.Im;

                    cofj = cofj * (i + 1); /* derivative */

                    ud.Re += cofj * t.Re;
                    ud.Im -= cofj * t.Im;

                    t.Re = t1.Re;
                    t.Im = t1.Im;
                }

                mag = ud.Re * ud.Re + ud.Im * ud.Im;

                if (mag == 0.0)
                {
                    if (!final)
                    {
                        if (retry < 50)
                        {

                            goto tryagn;
                        }
                        else
                        {
                            break;
                        }
                    }

                    x.Re = xsav.Re;
                    x.Im = xsav.Im;
                    goto findon;
                }

                dx.Re = (u.Im * ud.Im - u.Re * ud.Re) / mag;
                x.Re += dx.Re;

                dx.Im = -(u.Re * ud.Im + u.Im * ud.Re) / mag;
                x.Im += dx.Im;

                if ((Math.Abs(dx.Im) + Math.Abs(dx.Re)) < 1.0e-6)
                    goto lupdon;

                iter++;
            }

            if (final) goto lupdon;
            if (retry < 5) goto tryagn;

            return 3;

        lupdon:
            /* Swap original and reduced polynomials */
            q = &xcof[nsav];
            p = &cof[0];

            for (int i = 0; i <= n2; i++)
            {
                cofj = *q;
                *q-- = *p;
                *p++ = cofj;
            }

            int k = n;
            n = n1;
            n1 = k;

            if (!final)
            {
                final = true;

                if (Math.Abs(x.Im / x.Re) < 1.0e-4) x.Im = 0.0;

                xsav.Re = x.Re;
                xsav.Im = x.Im;
                goto finitr; /* do final iteration on original polynomial */
            }

        findon:
            final = false;

        zerrut:
            if (Math.Abs(x.Im / x.Re) >= 1.0e-5 && u.Re != 0.0)
            {
                cofj = x.Re + x.Re;
                mag = x.Re * x.Re + x.Im * x.Im;
                n -= 2;
            }
            else
            {
                x.Im = 0.0;
                cofj = x.Re;
                mag = 0.0;
                n -= 1;
            }

            p = &cof[1];
            *p += cofj * *(p - 1);

            for (int i = 1; i < n; i++)
            {
                *(p + 1) += cofj * *p - mag * *(p - 1);
                p++;
            }

        setrut:
            root[nroot].Re = x.Re;
            root[nroot].Im = x.Im;
            nroot += 1;

            if (mag != 0.0)
            {
                x.Im = -x.Im;
                mag = 0;
                goto setrut;	/* fill in the complex conjugate root */
            }

            if (n > 0)
                goto nxtrut;

            return 0;
        }
    }
}
