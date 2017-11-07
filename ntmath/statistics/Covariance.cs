namespace nt.math.statistics
{
    using System;
    using System.Collections;
    using nt.math;

    public static class Covariance
    {
        public static unsafe double ComputeCovariance(double[] data1, double[] data2, double mean1, double mean2)
        {
            double res = 0.0;

            fixed(double *p1 = data1)
            fixed (double *p2 = data2)
            {
                for (int i = 0; i < data1.GetLength(0); i++)
                {
                    res += ((p1[i] - mean1) * (p2[i] - mean2) - res) / (i + 1);
                }
            }

            return res;
        }

        public static unsafe double CovarianceM(double[] data1, double[] data2, double mean1, double mean2)
        {
            int len = data1.GetLength(0);

            return ComputeCovariance(data1, data2, mean1, mean2) * ((double)len / (double)(len - 1));
        }

        public static unsafe double Mean(double[] data)
        {
            double res = 0.0;
            int len = data.GetLength(0);

            fixed (double* p = data)
            {
                for (int i = 0; i < len; i++)
                {
                    res += (p[i] - res) / (i + 1);
                }
            }

            return res;
        }

        public static double CovarianceFull(double[] data1, double[] data2)
        {
            return CovarianceM(data1, data2, Mean(data1), Mean(data2));
        }

        public static unsafe double Correlation(double[] data1, double[] data2)
        {
            double ratio;
            double delta_x, delta_y;
            double sum_xsq = 0.0;
            double sum_ysq = 0.0;
            double sum_cross = 0.0;
            int len = data1.GetLength(0);

            fixed (double* p1 = data1)
            fixed (double* p2 = data2)
            {
                double mean_x = p1[0];
                double mean_y = p2[0];

                for (int i = 1; i < len; i++)
                {
                    ratio = i / (i + 1.0);

                    delta_x = p1[i] - mean_x;
                    delta_y = p2[i] - mean_y;

                    sum_xsq += delta_x * delta_x * ratio;
                    sum_ysq += delta_y * delta_y * ratio;

                    sum_cross += delta_x * delta_y * ratio;
                    mean_x += delta_x / (i + 1);
                    mean_y += delta_y / (i + 1);
                }
            }

            return sum_cross / (Math.Sqrt(sum_xsq) * Math.Sqrt(sum_ysq));
        }

        public static unsafe double Correlation2(double[] data1, double[] data2)
        {
            int len = data1.GetLength(0);
            int n = 0;
            double r = 0.0;

            double mean_x = 0.0;
            double mean_y = 0.0;
            double var_x = 0.0;
            double var_y = 0.0;

            fixed (double* p1 = data1)
            fixed (double* p2 = data2)
            {
                for (int i = 0; i < len; i++)
                {
                    double current_x = p1[i];
                    double current_y = p2[i];

                    double delta_x = current_x - mean_x;
                    double scale_x = delta_x / ++n;

                    double delta_y = current_y - mean_y;
                    double scale_y = delta_y / n;

                    mean_x += scale_x;
                    mean_y += scale_y;

                    var_x += scale_x * delta_x * (n - 1);
                    var_y += scale_y * delta_y * (n - 1);

                    r += (delta_x * delta_y * (n - 1)) / n;
                }
            }

            return r / (Math.Sqrt(var_x * var_y));
        }

        public static unsafe void ComputeRank(this double[] data)
        {
            int n = data.GetLength(0);
            int i = 0;

            fixed (double* p = data)
            {
                while (i < n - 1)
                {
                    double v = p[i];

                    if (v == p[i + 1])
                    {
                        int j = i + 2;
                        int k = 0;
                        double rank = 0.0;

                        while (j < n && v == p[j]) ++j;

                        for (k = i; k < j; ++k) rank += k + 1.0;

                        rank /= (double)(j - i);

                        for (k = i; k < j; ++k) p[k] = rank;

                        i = j;
                    }
                    else
                    {
                        p[i] = i + 1.0;
                    }
                }

                if (i == n - 1) p[n - 1] = (double)n;
            }
        }

        public static unsafe double Stats_Spearman(double[] data1, double[] data2, int n, double[] work)
        {
            double[] rank1 = Matrix.Vector<double>(n, work[0]);
            double[] rank2 = Matrix.Vector<double>(n, work[n]);

            fixed (double* p1 = data1)
            fixed (double* p2 = data2)
            {
                for(int i = 0; i < n; i++)
                {
                    p1[i] = rank1[i];
                    p2[i] = rank2[i];
                }

                rank1.Sort2(rank2);
                data1.ComputeRank();

                rank2.Sort2(rank1);
                data2.ComputeRank();
            }

            return Correlation(data1, data2);
        }
    }
}
