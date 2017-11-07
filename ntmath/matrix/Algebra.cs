namespace nt.math
{
    using System;
    using nt.core;

    public static partial class Matrix
    {
        #region Matrix-Matrix Multiplication
        public static double[,] Multiply(this double[,] a, double[,] b)
        {
            double[,] r = new double[a.GetLength(0), b.GetLength(1)];

            Multiply(a, b, r);
            return r;
        }

        public static double[][] Multiply(this double[][] a, double[][] b)
        {
            int rows = a.Length;
            int cols = b[0].Length;
            double[][] r = new double[rows][];

            for (int i = 0; i < rows; i++)
                r[i] = new double[cols];

            Multiply(a, b, r);
            return r;
        }

        public static float[][] Multiply(this float[][] a, float[][] b)
        {
            int rows = a.Length;
            int cols = b[0].Length;
            float[][] r = new float[rows][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new float[cols];

            Multiply(a, b, r);
            return r;
        }

        public static double[][] Multiply(this float[][] a, double[][] b)
        {
            int rows = a.Length;
            int cols = b[0].Length;
            double[][] r = new double[rows][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new double[cols];

            Multiply(a, b, r);
            return r;
        }

        public static float[,] Multiply(this float[,] a, float[,] b)
        {
            float[,] r = new float[a.GetLength(0), b.GetLength(1)];

            Multiply(a, b, r);
            return r;
        }

        public static unsafe void Multiply(this double[,] a, double[,] b, double[,] result)
        {
            // TODO: enable argument checking
            if (a.GetLength(1) != b.GetLength(0))
                throw new ArgumentException("Matrix dimensions must match");

            int n = a.GetLength(1);
            int m = result.GetLength(0); //a.GetLength(0);
            int p = result.GetLength(1); //b.GetLength(1);

            fixed (double* ptrA = a)
            {
                double[] Bcolj = new double[n];

                for (int j = 0; j < p; j++)
                {
                    for (int k = 0; k < Bcolj.Length; k++)
                        Bcolj[k] = b[k, j];

                    double* Arowi = ptrA;

                    for (int i = 0; i < m; i++)
                    {
                        double s = 0;

                        for (int k = 0; k < Bcolj.Length; k++)
                            s += *(Arowi++) * Bcolj[k];

                        result[i, j] = s;
                    }
                }
            }
        }

        public static void Multiply(this double[][] a, double[][] b, double[][] result)
        {
            // TODO: enable argument checking
            if (a.GetLength(1) != b.GetLength(0))
                 throw new ArgumentException("Matrix dimensions must match");

            int n = a[0].Length;
            int m = a.Length;
            int p = b[0].Length;
            double[] Bcolj = new double[n];

            for (int j = 0; j < p; j++)
            {
                for (int k = 0; k < n; k++)
                    Bcolj[k] = b[k][j];

                for (int i = 0; i < m; i++)
                {
                    double[] Arowi = a[i];
                    double s = 0;

                    for (int k = 0; k < n; k++)
                        s += Arowi[k] * Bcolj[k];

                    result[i][j] = s;
                }
            }
        }

        public static unsafe void Multiply(this double[,] a, double[,] b, double[][] result)
        {
            // TODO: enable argument checking
            if (a.GetLength(1) != b.GetLength(0))
                throw new ArgumentException("Matrix dimensions must match");

            int n = a.GetLength(1);
            int m = result.GetLength(0); //a.GetLength(0);
            int p = result.GetLength(1); //b.GetLength(1);

            fixed (double* ptrA = a)
            {
                double[] Bcolj = new double[n];
                for (int j = 0; j < p; j++)
                {
                    for (int k = 0; k < Bcolj.Length; k++)
                        Bcolj[k] = b[k, j];

                    double* Arowi = ptrA;
                    for (int i = 0; i < m; i++)
                    {
                        double s = 0;
                        for (int k = 0; k < Bcolj.Length; k++)
                            s += *(Arowi++) * Bcolj[k];
                        result[i][j] = s;
                    }
                }
            }
        }

        public static void Multiply(this float[][] a, float[][] b, float[][] result)
        {
            // TODO: enable argument checking
            if (a.GetLength(1) != b.GetLength(0))
                throw new ArgumentException("Matrix dimensions must match");

            int n = a[0].Length;
            int m = a.Length;
            int p = b[0].Length;
            float[] Bcolj = new float[n];

            for (int j = 0; j < p; j++)
            {
                for (int k = 0; k < n; k++)
                    Bcolj[k] = b[k][j];

                for (int i = 0; i < m; i++)
                {
                    float[] Arowi = a[i];
                    float s = 0;

                    for (int k = 0; k < n; k++)
                        s += Arowi[k] * Bcolj[k];

                    result[i][j] = s;
                }
            }
        }

        public static void Multiply(this float[][] a, double[][] b, double[][] result)
        {
            // TODO: enable argument checking
            if (a.GetLength(1) != b.GetLength(0))
                throw new ArgumentException("Matrix dimensions must match");

            int n = a[0].Length;
            int m = a.Length;
            int p = b[0].Length;
            double[] Bcolj = new double[n];

            for (int j = 0; j < p; j++)
            {
                for (int k = 0; k < n; k++)
                    Bcolj[k] = b[k][j];

                for (int i = 0; i < m; i++)
                {
                    float[] Arowi = a[i];

                    double s = 0;
                    for (int k = 0; k < n; k++)
                        s += Arowi[k] * Bcolj[k];

                    result[i][j] = s;
                }
            }
        }

        public static unsafe void Multiply(this float[,] a, float[,] b, float[,] result)
        {
            int acols = a.GetLength(1);
            int arows = a.GetLength(0);
            int bcols = b.GetLength(1);

            fixed (float* ptrA = a)
            {
                float[] Bcolj = new float[acols];

                for (int j = 0; j < bcols; j++)
                {
                    for (int k = 0; k < acols; k++)
                        Bcolj[k] = b[k, j];

                    float* Arowi = ptrA;

                    for (int i = 0; i < arows; i++)
                    {
                        float s = 0;

                        for (int k = 0; k < acols; k++)
                            s += *(Arowi++) * Bcolj[k];

                        result[i, j] = s;
                    }
                }
            }
        }

        public static double[,] MultiplyByTranspose(this double[,] a, double[,] b)
        {
            double[,] r = new double[a.GetLength(0), b.GetLength(0)];

            MultiplyByTranspose(a, b, r);
            return r;
        }

        public static float[,] MultiplyByTranspose(this float[,] a, float[,] b)
        {
            float[,] r = new float[a.GetLength(0), b.GetLength(0)];
            MultiplyByTranspose(a, b, r);
            return r;
        }

        public static unsafe void MultiplyByTranspose(this double[,] a, double[,] b, double[,] result)
        {
            int n = a.GetLength(1);
            int m = a.GetLength(0);
            int p = b.GetLength(0);

            fixed (double* ptrA = a)
            fixed (double* ptrB = b)
            fixed (double* ptrR = result)
            {
                double* rc = ptrR;

                for (int i = 0; i < m; i++)
                {
                    double* bColj = ptrB;

                    for (int j = 0; j < p; j++)
                    {
                        double* aColi = ptrA + n * i;
                        double s = 0;

                        for (int k = 0; k < n; k++)
                            s += *(aColi++) * *(bColj++);

                        *(rc++) = s;
                    }
                }
            }
        }

        public static unsafe void MultiplyByTranspose(this float[,] a, float[,] b, float[,] result)
        {
            int n = a.GetLength(1);
            int m = a.GetLength(0);
            int p = b.GetLength(0);

            fixed (float* ptrA = a)
            fixed (float* ptrB = b)
            fixed (float* ptrR = result)
            {
                float* rc = ptrR;

                for (int i = 0; i < m; i++)
                {
                    float* bColj = ptrB;

                    for (int j = 0; j < p; j++)
                    {
                        float* aColi = ptrA + n * i;
                        float s = 0;

                        for (int k = 0; k < n; k++)
                            s += *(aColi++) * *(bColj++);

                        *(rc++) = s;
                    }
                }
            }
        }

        public static double[,] TransposeAndMultiply(this double[,] a, double[,] b)
        {
            double[,] r = new double[a.GetLength(1), b.GetLength(1)];

            TransposeAndMultiply(a, b, r);
            return r;
        }

        public static double[][] TransposeAndMultiply(this double[][] a, double[][] b)
        {
            int aCols = a[0].Length;
            int bCols = b[0].Length;
            double[][] r = new double[aCols][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new double[bCols];

            TransposeAndMultiply(a, b, r);
            return r;
        }

        public static unsafe void TransposeAndMultiply(this double[,] a, double[,] b, double[,] result)
        {
            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");
            if (result == null) throw new ArgumentNullException("result");

            // TODO: Check dimensions
            // TODO: Change result to be an "out" value

            int n = a.GetLength(0);
            int m = a.GetLength(1);
            int p = b.GetLength(1);
            double[] Bcolj = new double[n];

            for (int i = 0; i < p; i++)
            {
                for (int k = 0; k < n; k++)
                    Bcolj[k] = b[k, i];

                for (int j = 0; j < m; j++)
                {
                    double s = 0;

                    for (int k = 0; k < n; k++)
                        s += a[k, j] * Bcolj[k];

                    result[j, i] = s;
                }
            }
        }

        public static unsafe void TransposeAndMultiply(this double[][] a, double[][] b, double[][] result)
        {
            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");
            if (result == null) throw new ArgumentNullException("result");

            // TODO: Check dimensions
            // TODO: Change result to be an "out" value

            int n = a.Length;
            int m = a[0].Length;
            int p = b[0].Length;
            double[] Bcolj = new double[n];

            for (int i = 0; i < p; i++)
            {
                for (int k = 0; k < n; k++)
                    Bcolj[k] = b[k][i];

                for (int j = 0; j < m; j++)
                {
                    double s = 0;

                    for (int k = 0; k < n; k++)
                        s += a[k][j] * Bcolj[k];

                    result[j][i] = s;
                }
            }
        }

        public static double[] TransposeAndMultiply(this double[,] a, double[] b)
        {
            double[] r = new double[a.GetLength(1)];

            TransposeAndMultiply(a, b, r);
            return r;
        }

        public static unsafe void TransposeAndMultiply(this double[,] a, double[] b, double[] result)
        {
            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");
            if (result == null) throw new ArgumentNullException("result");

            for (int j = 0; j < result.Length; j++)
            {
                double s = 0;

                for (int k = 0; k < b.Length; k++)
                    s += a[k, j] * b[k];

                result[j] = s;
            }
        }

        public static double[,] TransposeAndMultiplyByDiagonal(this double[,] a, double[] b)
        {
            double[,] r = new double[a.GetLength(1), b.Length];

            TransposeAndMultiplyByDiagonal(a, b, r);
            return r;
        }

        public static unsafe void TransposeAndMultiplyByDiagonal(this double[,] a, double[] b, double[,] result)
        {
            if (a.GetLength(0) != b.Length)
                throw new ArgumentException("Matrix dimensions must match.");

            int m = a.GetLength(1);

            for (int i = 0; i < b.Length; i++)
                for (int j = 0; j < m; j++)
                    result[j, i] = a[i, j] * b[i];
        }

        public static double[][] MultiplyByDiagonal(this double[][] a, double[] b)
        {
            double[][] r = new double[a.Length][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new double[b.Length];

            MultiplyByDiagonal(a, b, r);
            return r;
        }

        public static float[][] MultiplyByDiagonal(this float[][] a, float[] b)
        {
            var r = new float[a.Length][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new float[b.Length];

            MultiplyByDiagonal(a, b, r);
            return r;
        }

        public static double[,] MultiplyByDiagonal(this double[,] a, double[] b)
        {
            double[,] r = new double[a.GetLength(0), b.Length];

            MultiplyByDiagonal(a, b, r);
            return r;
        }

        public static unsafe void MultiplyByDiagonal(this double[,] a, double[] b, double[,] result)
        {
            if (a.GetLength(1) != b.Length)
                throw new ArgumentException("Matrix dimensions must match.");


            int rows = a.GetLength(0);

            fixed (double* ptrA = a, ptrR = result)
            {
                double* A = ptrA;
                double* R = ptrR;

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < b.Length; j++)
                        *R++ = *A++ * b[j];
            }
        }

        public static unsafe void MultiplyByDiagonal(this double[][] a, double[] b, double[][] result)
        {
            int rows = a.Length;

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < b.Length; j++)
                    result[i][j] = a[i][j] * b[j];
        }

        public static unsafe void MultiplyByDiagonal(this float[][] a, float[] b, float[][] result)
        {
            int rows = a.Length;

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < b.Length; j++)
                    result[i][j] = a[i][j] * b[j];
        }

        public static float[,] MultiplyByDiagonal(this float[,] a, float[] b)
        {
            float[,] r = new float[a.GetLength(0), b.Length];

            MultiplyByDiagonal(a, b, r);
            return r;
        }

        public static unsafe void MultiplyByDiagonal(this float[,] a, float[] b, float[,] result)
        {
            if (a.GetLength(1) != b.Length)
                throw new ArgumentException("Matrix dimensions must match.");


            int rows = a.GetLength(0);

            fixed (float* ptrA = a, ptrR = result)
            {
                float* A = ptrA;
                float* R = ptrR;

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < b.Length; j++)
                        *R++ = *A++ * b[j];
            }
        }

        public static double[,] DivideByDiagonal(this double[,] a, double[] b)
        {
            double[,] r = new double[a.GetLength(0), b.Length];

            DivideByDiagonal(a, b, r);
            return r;
        }

        public static unsafe void DivideByDiagonal(this double[,] a, double[] b, double[,] result)
        {
            if (a.GetLength(1) != b.Length)
                throw new ArgumentException("Matrix dimensions must match.");

            int rows = a.GetLength(0);

            fixed (double* ptrA = a, ptrR = result)
            {
                double* A = ptrA;
                double* R = ptrR;

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < b.Length; j++)
                        *R++ = *A++ / b[j];
            }
        }
        #endregion

        #region Matrix-Vector Multiplication
        public static double[] Multiply(this double[] rowVector, double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (rows != rowVector.Length)
            {
                throw new DimensionMismatchException("matrix", "Matrix must have the same number of rows as the length of the vector.");
            }

            double[] r = new double[cols];

            for (int j = 0; j < cols; j++)
                for (int k = 0; k < rowVector.Length; k++)
                    r[j] += rowVector[k] * matrix[k, j];

            return r;
        }

        public static float[] Multiply(this float[] rowVector, float[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (rows != rowVector.Length)
                throw new DimensionMismatchException("matrix", "Matrix must have the same number of rows as the length of the vector.");

            float[] r = new float[cols];

            for (int j = 0; j < cols; j++)
                for (int k = 0; k < rowVector.Length; k++)
                    r[j] += rowVector[k] * matrix[k, j];

            return r;
        }

        public static double[] Multiply(this double[,] matrix, double[] columnVector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (cols != columnVector.Length)
                throw new DimensionMismatchException("columnVector", "Vector must have the same length as columns in the matrix.");

            double[] r = new double[rows];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columnVector.Length; j++)
                    r[i] += matrix[i, j] * columnVector[j];

            return r;
        }

        public static float[] Multiply(this float[][] matrix, float[] columnVector)
        {
            int rows = matrix.Length;
            float[] r = new float[rows];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columnVector.Length; j++)
                    r[i] += matrix[i][j] * columnVector[j];

            return r;
        }

        public static double[] Multiply(this double[][] matrix, double[] columnVector)
        {
            int rows = matrix.Length;
            int cols = matrix[0].Length;

            if (cols != columnVector.Length)
                throw new DimensionMismatchException("columnVector", "Vector must have the same length as columns in the matrix.");

            double[] r = new double[rows];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columnVector.Length; j++)
                    r[i] += matrix[i][j] * columnVector[j];

            return r;
        }

        public static float[] Multiply(this float[,] matrix, float[] columnVector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (cols != columnVector.Length)
                throw new DimensionMismatchException("columnVector", "Vector must have the same length as columns in the matrix.");

            var r = new float[rows];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columnVector.Length; j++)
                    r[i] += matrix[i, j] * columnVector[j];

            return r;
        }

        public static double[,] Multiply(this double[,] matrix, double x, bool inPlace = false)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[,] result = inPlace ? matrix : new double[rows, cols];

            Multiply(matrix, x, result);
            return result;
        }

        public static double[,] Multiply(this double[,] matrix, double x)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[,] result = new double[rows, cols];

            Multiply(matrix, x, result);
            return result;
        }

        public static float[,] Multiply(this float[,] matrix, float x)
        {
            float[,] result = new float[matrix.GetLength(0), matrix.GetLength(1)];

            Multiply(matrix, x, result);
            return result;
        }

        public unsafe static void Multiply(this double[,] matrix, double x, double[,] result)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int length = matrix.Length;

            fixed (double* ptrA = matrix, ptrR = result)
            {
                double* pa = ptrA, pr = ptrR;

                for (int i = 0; i < length; i++, pa++, pr++)
                    *pr = *pa * x;
            }
        }

        public unsafe static void Multiply(this float[,] matrix, float x, float[,] result)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int length = matrix.Length;

            fixed (float* ptrA = matrix, ptrR = result)
            {
                float* pa = ptrA, pr = ptrR;

                for (int i = 0; i < length; i++, pa++, pr++)
                    *pr = *pa * x;
            }
        }

        public static double[] Multiply(this double[] vector, double x)
        {
            double[] r = new double[vector.Length];

            for (int i = 0; i < vector.Length; i++)
                r[i] = vector[i] * x;

            return r;
        }

        public static double[] Multiply(this double[] vector, double x, bool inPlace)
        {
            double[] r = inPlace ? vector : new double[vector.Length];

            for (int i = 0; i < vector.Length; i++)
                r[i] = vector[i] * x;

            return r;
        }

        public static float[] Multiply(this float[] vector, float x)
        {
            float[] r = new float[vector.Length];

            for (int i = 0; i < vector.Length; i++)
                r[i] = vector[i] * x;

            return r;
        }

        public static double[,] Multiply(this double x, double[,] matrix)
        {
            return matrix.Multiply(x);
        }

        public static float[,] Multiply(this float x, float[,] matrix)
        {
            return matrix.Multiply(x);
        }

        public static double[] Multiply(this double x, double[] vector)
        {
            return vector.Multiply(x);
        }

        public static float[] Multiply(this float x, float[] vector)
        {
            return vector.Multiply(x);
        }

        public static double[] Multiply(this int x, double[] vector)
        {
            return vector.Multiply(x);
        }

        public static float[] Multiply(this int x, float[] vector)
        {
            return vector.Multiply(x);
        }
        #endregion

        public static double[,] Power(this double[,] matrix, int n)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            if (!matrix.IsSquare())
                throw new ArgumentException("Matrix must be square", "matrix");

            if (n == 0)
                return Matrix.Identity(matrix.GetLength(0));

            double[,] result = matrix;
            string bin = System.Convert.ToString(n, 2);

            for (int i = 1; i < bin.Length; i++)
            {
                result = Matrix.Multiply(result, result);

                if (bin[i] == '1')
                    result = Matrix.Multiply(result, matrix);
            }

            return result;
        }
    }
}
