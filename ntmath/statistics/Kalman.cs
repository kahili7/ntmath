namespace nt.math.statistics
{
    using System;
    using nt.math;

    public enum KALMAN
    {
        ARMA_LL = 1 << 0,
        AVG_LL = 1 << 1,
        USER = 1 << 2,
        DIFFUSE = 1 << 3,
        FORWARD = 1 << 4,
        SMOOTH = 1 << 5,
        SIM = 1 << 6,
        CROSS = 1 << 7,
        ETT = 1 << 8,
        CHECK = 1 << 9
    }

    public enum KALMAN_CHECK
    {
        F = 0,
        A,
        H,
        Q,
        R,
        m,
        y,
        x,
        S,
        P,
        MMAX
    }

    public enum KALMAN_RESIZE
    {
        E,
        V,
        BIG_S,
        BIG_P,
        K,
        LL
    }

    public class Kalman
    {
        public int flags;
        public int fnlevel;

        public int r;
        public int n;
        public int k;
        public int p;
        public int T;
        public int okT;
        public int t;

        public int ifc;

        public double[] S0;
        public double[] S1;
        public double[,] P0;
        public double[,] P1;

        public double[] e;

        public static readonly double[,] F;
        public static readonly double[,] A;
        public static readonly double[,] H;
        public static readonly double[,] Q;
        public static readonly double[,] R;

        public static readonly double[] mu;
        public static readonly double[,] y;
        public static readonly double[,] x;

        public static readonly double[] Sini;
        public static readonly double[,] Pini;

        #region Static functions
        public static bool CheckMatrix(Kalman k, double[,] m, KALMAN_CHECK f)
        {
            int r = 0;
            int c = 0;
            bool symm = (f == KALMAN_CHECK.Q || f == KALMAN_CHECK.R);

            if (f == KALMAN_CHECK.F || f == KALMAN_CHECK.Q || f == KALMAN_CHECK.P) r = c = k.r;
            else if(f == KALMAN_CHECK.A)
            {
                r = k.k;
                c = k.n;
            }
            else if(f == KALMAN_CHECK.H)
            {
                r = k.r;
                c = k.n;
            }
            else if (f == KALMAN_CHECK.R) r = c = k.n;
            else if (f == KALMAN_CHECK.S || f == KALMAN_CHECK.m)
            {
                r = k.r;
                c = 1;
            }

            if (m.GetLength(0) != r || m.GetLength(1) != c) return false;
            else if (symm) return false;

            return true;
        }

        public static bool ResizeMatrix(Kalman k, double[,] m, KALMAN_RESIZE f)
        {
            int r = k.T;
            int c = 0;

            if (f == KALMAN_RESIZE.E) c = k.n;
            else if (f == KALMAN_RESIZE.V) c = (k.n * k.n + k.n) / 2;
            else if (f == KALMAN_RESIZE.BIG_S) c = k.r;
            else if (f == KALMAN_RESIZE.BIG_P) c = (k.r * k.r + k.r) / 2;
            else if (f == KALMAN_RESIZE.LL) c = 1;
            else if (f == KALMAN_RESIZE.K) c = k.r * k.n;
            else return false;

            m = null;

            if (m.GetLength(0) != r || m.GetLength(1) != c)
            {
                m = new double[r, c];
            }

            return m != null;
         }
        #endregion
    }
}
