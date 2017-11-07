namespace nt.math
{
    using System;
    using System.Collections.Generic;
    using nt.core;

    public static partial class Matrix
    {
        public static unsafe void DownHeap(this double[] data, int n, int k)
        {
            fixed(double* p = data)
            {
                double v = *(p + k);

                while(k <= n / 2)
                {
                    int j = 2 * k;

                    if (j < n && *(p + j) < *(p + (j + 1))) j++;
                    if (!(v < *(p + j))) break;

                    *(p + k) = *(p + j);
                    k = j;
                }

                *(p + k) = v;
            }
        }

        public static unsafe void DownHeap2(this double[] data1, double[] data2, int n, int k)
        {
            fixed (double* p1 = data1)
            fixed(double* p2 = data2)
            {
                double v1 = *(p1 + k);
                double v2 = *(p2 + k);

                while (k <= n / 2)
                {
                    int j = 2 * k;

                    if (j < n && *(p1 + j) < *(p1 + (j + 1))) j++;
                    if (!(v1 < *(p1 + j))) break;

                    *(p1 + k) = *(p1 + j);
                    *(p2 + k) = *(p2 + j);
                    k = j;
                }

                *(p1 + k) = v1;
                *(p2 + k) = v2;
            }
        }

        public static unsafe void Sort(this double[] data)
        {
            int cnt = data.GetLength(0);

            if(cnt == 0) throw new ArgumentException("Vector dimensions must match");

            int n = cnt - 1;
            int k = n / 2;

            k++;

            do{
                data.DownHeap(n, k);
            } while(k > 0);

            fixed(double* p = data)
            {
                while(n > 0)
                {
                    double tmp = *(p + 0);

                    *(p + 0) = *(p + n);
                    *(p + n) = tmp;
                    n--;

                    data.DownHeap(n, 0);
                }
            }
        }

        public static unsafe void Sort2(this double[] data1, double[] data2)
        {
            int cnt = data1.GetLength(0);

            if (cnt == 0) throw new ArgumentException("Vector dimensions must match");

            int n = cnt - 1;
            int k = n / 2;

            k++;

            do
            {
                data1.DownHeap2(data2, n, k);
            } while (k > 0);

            fixed (double* p1 = data1)
            fixed (double* p2 = data2)
            {
                while (n > 0)
                {
                    double tmp = *(p1 + 0);

                    *(p1 + 0) = *(p1 + n);
                    *(p1 + n) = tmp;

                    tmp = *(p2 + 0);
                    *(p2 + 0) = *(p2 + n);
                    *(p2 + n) = tmp;

                    n--;

                    data1.DownHeap2(data2, n, 0);
                }
            }
        }
    }
}
