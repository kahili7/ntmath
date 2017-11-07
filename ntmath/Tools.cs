namespace nt.math
{
    using System;

    public static class Tools
    {
        public static double Hypotenuse(double a, double b)
        {
            double r = 0.0;
            double absA = System.Math.Abs(a);
            double absB = System.Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * System.Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * System.Math.Sqrt(1 + r * r);
            }

            return r;
        }

        public static decimal Hypotenuse(decimal a, decimal b)
        {
            decimal r = 0;
            decimal absA = System.Math.Abs(a);
            decimal absB = System.Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * (decimal)System.Math.Sqrt((double)(1 + r * r));
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * (decimal)System.Math.Sqrt((double)(1 + r * r));
            }

            return r;
        }

        public static float Hypotenuse(float a, float b)
        {
            double r = 0;
            float absA = System.Math.Abs(a);
            float absB = System.Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * System.Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * System.Math.Sqrt(1 + r * r);
            }

            return (float)r;
        }

        public static double Interpolate1D(double value, double[] x, double[] y, double lower, double upper)
        {
            for (int i = 0; i < x.Length; i++)
            {
                if (value < x[i])
                {
                    if (i == 0)
                        return lower;

                    int start = i - 1;
                    int next = i;
                    double m = (value - x[start]) / (x[next] - x[start]);

                    return y[start] + (y[next] - y[start]) * m;
                }
            }

            return upper;
        }
    }
}
