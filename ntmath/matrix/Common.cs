namespace nt.math
{
    using System;
    using System.Collections;
    using nt.core;
    using nt.math.comparers;

    public static partial class Matrix
    {
        #region Transpose
        public static T[,] Transpose<T>(this T[,] matrix)
        {
            return Transpose(matrix, false);
        }

        public static T[][] Transpose<T>(this T[][] matrix)
        {
            return Transpose(matrix, false);
        }

        public static T[,] Transpose<T>(this T[,] matrix, bool inPlace)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (inPlace)
            {
                if (rows != cols)
                    throw new ArgumentException("Only square matrices can be transposed in place.", "matrix");

                for (int i = 0; i < rows; i++)
                {
                    for (int j = i; j < cols; j++)
                    {
                        T element = matrix[j, i];

                        matrix[j, i] = matrix[i, j];
                        matrix[i, j] = element;
                    }
                }

                return matrix;
            }
            else
            {
                T[,] result = new T[cols, rows];

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[j, i] = matrix[i, j];

                return result;
            }
        }

        public static T[][] Transpose<T>(this T[][] matrix, bool inPlace)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.Length;
            if (rows == 0) return new T[rows][];
            int cols = matrix[0].Length;

            if (inPlace)
            {
                if (rows != cols)
                    throw new ArgumentException("Only square matrices can be transposed in place.", "matrix");

                for (int i = 0; i < rows; i++)
                {
                    for (int j = i; j < cols; j++)
                    {
                        T element = matrix[j][i];

                        matrix[j][i] = matrix[i][j];
                        matrix[i][j] = element;
                    }
                }

                return matrix;
            }
            else
            {
                T[][] result = new T[cols][];

                for (int j = 0; j < cols; j++)
                {
                    result[j] = new T[rows];

                    for (int i = 0; i < rows; i++)
                        result[j][i] = matrix[i][j];
                }

                return result;
            }
        }

        public static T[,] Transpose<T>(this T[] rowVector)
        {
            if (rowVector == null) throw new ArgumentNullException("rowVector");

            T[,] trans = new T[rowVector.Length, 1];
            for (int i = 0; i < rowVector.Length; i++)
                trans[i, 0] = rowVector[i];

            return trans;
        }

        public static Array Transpose(this Array array, int[] order)
        {
            return transpose(array, order);
        }

        public static T Transpose<T>(this T array, int[] order) where T : class, ICloneable, IList
        {
            Array arr = array as Array;

            if (arr == null)
                throw new ArgumentException("The given object must inherit from System.Array.", "array");

            return transpose(arr, order) as T;
        }

        private static Array transpose(Array array, int[] order)
        {
            if (array.Length == 1 || array.Length == 0)
                return array;

            // Get the number of samples at each dimension
            int[] size = new int[array.Rank];

            for (int i = 0; i < size.Length; i++)
                size[i] = array.GetLength(i);

            Array r = Array.CreateInstance(array.GetType().GetElementType(), size.Submatrix(order));

            // Generate all indices for accessing the matrix 
            foreach (int[] pos in Combinatorics.Sequences(size, true))
            {
                int[] newPos = pos.Submatrix(order);
                object value = array.GetValue(pos);

                r.SetValue(value, newPos);
            }

            return r;
        }
        #endregion

        #region Matrix character
        public static bool IsSorted<T>(this T[] values, ComparerDirection direction) where T : IComparable<T>
        {
            if (direction == ComparerDirection.Descending)
            {
                for (int i = 1; i < values.Length; i++)
                    if (values[i - 1].CompareTo(values[i]) >= 0)
                        return false;
            }
            else
            {
                for (int i = 1; i < values.Length; i++)
                    if (values[i - 1].CompareTo(values[i]) <= 0)
                        return false;
            }

            return true;
        }

        public static bool IsSquare<T>(this T[,] matrix)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            return matrix.GetLength(0) == matrix.GetLength(1);
        }

        public static bool IsSymmetric<T>(this T[,] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            if (matrix.GetLength(0) == matrix.GetLength(1))
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        if (matrix[i, j].CompareTo(matrix[j, i]) != 0)
                            return false;
                    }
                }

                return true;
            }

            return false;
        }

        public static bool IsSymmetric<T>(this T[][] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            if (matrix.Length == matrix[0].Length)
            {
                for (int i = 0; i < matrix.Length; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        if (matrix[i][j].CompareTo(matrix[j][i]) != 0)
                            return false;
                    }
                }
                return true;
            }
            return false;
        }

        public static bool IsUpperTriangular<T>(this T[][] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.Length != matrix[0].Length)
                return false;

            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < i; j++)
                    if (matrix[i][j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static bool IsLowerTriangular<T>(this T[][] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.Length != matrix[0].Length)
                return false;

            for (int i = 0; i < matrix.Length; i++)
                for (int j = i + 1; j < matrix[i].Length; j++)
                    if (matrix[i][j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static bool IsUpperTriangular<T>(this T[,] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.GetLength(0) != matrix.GetLength(1))
                return false;

            for (int i = 0; i < matrix.GetLength(0); i++)
                for (int j = 0; j < i; j++)
                    if (matrix[i, j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static bool IsLowerTriangular<T>(this T[,] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.GetLength(0) != matrix.GetLength(1))
                return false;

            for (int i = 0; i < matrix.GetLength(0); i++)
                for (int j = i + 1; j < matrix.GetLength(1); j++)
                    if (matrix[i, j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static bool IsDiagonal<T>(this T[,] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.GetLength(0) != matrix.GetLength(1))
                return false;

            for (int i = 0; i < matrix.GetLength(0); i++)
                for (int j = 0; j < matrix.GetLength(1); j++)
                    if (i != j && matrix[i, j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static bool IsDiagonal<T>(this T[][] matrix) where T : IComparable
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            T zero = default(T);

            if (matrix.Length != matrix[0].Length)
                return false;

            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++)
                    if (i != j && matrix[i][j].CompareTo(zero) != 0)
                        return false;

            return true;
        }

        public static double Trace(this double[,] matrix)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            double trace = 0.0;

            for (int i = 0; i < rows; i++)
                trace += matrix[i, i];
            return trace;
        }

        public static unsafe double Trace(double[,] matrixA, double[,] matrixB)
        {
            if (matrixA.Length != matrixB.Length)
                throw new DimensionMismatchException("matrixB", "Matrices must have the same length.");

            int length = matrixA.Length;

            fixed (double* ptrA = matrixA)
            fixed (double* ptrB = matrixB)
            {
                double* a = ptrA;
                double* b = ptrB;
                double trace = 0.0;

                for (int i = 0; i < length; i++)
                    trace += (*a++) * (*b++);

                return trace;
            }
        }

        public static int Trace(this int[,] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int trace = 0;

            for (int i = 0; i < rows; i++)
                trace += matrix[i, i];

            return trace;
        }

        public static float Trace(this float[,] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            float trace = 0.0f;

            for (int i = 0; i < rows; i++)
                trace += matrix[i, i];

            return trace;
        }

        public static float Trace(this float[][] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            float trace = 0.0f;

            for (int i = 0; i < matrix.Length; i++)
                trace += matrix[i][i];

            return trace;
        }

        public static T[] Diagonal<T>(this T[][] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            T[] r = new T[matrix.Length];

            for (int i = 0; i < r.Length; i++)
                r[i] = matrix[i][i];

            return r;
        }

        public static T[] Diagonal<T>(this T[,] matrix)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            var r = new T[matrix.GetLength(0)];

            for (int i = 0; i < r.Length; i++)
                r[i] = matrix[i, i];

            return r;
        }

        /*public static double Determinant(this double[,] matrix)
        {
            return Determinant(matrix, false);
        }*/
        #endregion

        #region Summation
        public static float[] Sum(this float[,] matrix)
        {
            return Sum(matrix, 0);
        }

        public static double[] Sum(this double[,] matrix)
        {
            return Sum(matrix, 0);
        }

        public static double[] Sum(this double[,] matrix, int dimension)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[] sum;

            if (dimension == -1)
            {
                sum = new double[1];

                foreach (double a in matrix)
                    sum[0] += a;
            }
            else if (dimension == 0)
            {
                sum = new double[cols];

                for (int j = 0; j < cols; j++)
                {
                    double s = 0.0;

                    for (int i = 0; i < rows; i++)
                        s += matrix[i, j];

                    sum[j] = s;
                }
            }
            else if (dimension == 1)
            {
                sum = new double[rows];

                for (int j = 0; j < rows; j++)
                {
                    double s = 0.0;

                    for (int i = 0; i < cols; i++)
                        s += matrix[j, i];

                    sum[j] = s;
                }
            }
            else
            {
                throw new ArgumentException("Invalid dimension", "dimension");
            }

            return sum;
        }

        public static float[] Sum(this float[,] matrix, int dimension)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            float[] sum;

            if (dimension == -1)
            {
                sum = new float[1];

                foreach (float a in matrix)
                    sum[0] += a;
            }
            else if (dimension == 0)
            {
                sum = new float[cols];

                for (int j = 0; j < cols; j++)
                {
                    float s = 0.0f;

                    for (int i = 0; i < rows; i++)
                        s += matrix[i, j];

                    sum[j] = s;
                }
            }
            else if (dimension == 1)
            {
                sum = new float[rows];

                for (int j = 0; j < rows; j++)
                {
                    float s = 0.0f;

                    for (int i = 0; i < cols; i++)
                        s += matrix[j, i];

                    sum[j] = s;
                }
            }
            else
            {
                throw new ArgumentException("Invalid dimension", "dimension");
            }

            return sum;
        }

        public static double[] Sum(this double[][] matrix)
        {
            return Sum(matrix, 0);
        }

        public static double[] Sum(this double[][] matrix, int dimension)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.Length;
            int cols = matrix[0].Length;
            double[] sum;

            if (dimension == 0)
            {
                sum = new double[cols];

                for (int i = 0; i < matrix.Length; i++)
                {
                    double[] row = matrix[i];

                    for (int j = 0; j < row.Length; j++)
                        sum[j] += row[j];
                }
            }
            else if (dimension == 1)
            {
                sum = new double[rows];

                for (int j = 0; j < matrix.Length; j++)
                {
                    double[] row = matrix[j];
                    double s = 0.0;

                    for (int i = 0; i < row.Length; i++)
                        s += row[i];

                    sum[j] = s;
                }
            }
            else
            {
                throw new ArgumentException("Invalid dimension", "dimension");
            }

            return sum;
        }

        public static int[] Sum(this int[,] matrix)
        {
            return Sum(matrix, 0);
        }

        public static int[] Sum(this int[,] matrix, int dimension)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int[] sum;

            if (dimension == 0)
            {
                sum = new int[cols];

                for (int j = 0; j < cols; j++)
                {
                    int s = 0;

                    for (int i = 0; i < rows; i++)
                        s += matrix[i, j];

                    sum[j] = s;
                }
            }
            else if (dimension == 1)
            {
                sum = new int[rows];

                for (int j = 0; j < rows; j++)
                {
                    int s = 0;

                    for (int i = 0; i < cols; i++)
                        s += matrix[j, i];

                    sum[j] = s;
                }
            }
            else
            {
                throw new ArgumentException("Invalid dimension", "dimension");
            }

            return sum;
        }

        public static double Sum(this double[] vector)
        {
            if (vector == null) throw new ArgumentNullException("vector");

            double sum = 0.0;

            for (int i = 0; i < vector.Length; i++)
                sum += vector[i];

            return sum;
        }

        public static float Sum(this float[] vector)
        {
            if (vector == null) throw new ArgumentNullException("vector");

            float sum = 0.0f;

            for (int i = 0; i < vector.Length; i++)
                sum += vector[i];

            return sum;
        }

        public static int Sum(this int[] vector)
        {
            if (vector == null) throw new ArgumentNullException("vector");

            int sum = 0;

            for (int i = 0; i < vector.Length; i++)
                sum += vector[i];

            return sum;
        }

        public static double[] CumulativeSum(this double[] vector)
        {
            if (vector == null) throw new ArgumentNullException("vector");
            if (vector.Length == 0) return new double[0];

            double[] sum = new double[vector.Length];

            sum[0] = vector[0];

            for (int i = 1; i < vector.Length; i++)
                sum[i] += sum[i - 1] + vector[i];

            return sum;
        }

        public static double[][] CumulativeSum(this double[,] matrix, int dimension)
        {
            if (matrix == null) throw new ArgumentNullException("matrix");

            double[][] sum;

            if (dimension == 1)
            {
                sum = new double[matrix.GetLength(0)][];
                sum[0] = matrix.GetRow(0);

                // for each row
                for (int i = 1; i < matrix.GetLength(0); i++)
                {
                    sum[i] = new double[matrix.GetLength(1)];

                    // for each column
                    for (int j = 0; j < matrix.GetLength(1); j++)
                        sum[i][j] += sum[i - 1][j] + matrix[i, j];
                }
            }
            else if (dimension == 0)
            {
                sum = new double[matrix.GetLength(1)][];
                sum[0] = matrix.GetColumn(0);

                // for each column
                for (int i = 1; i < matrix.GetLength(1); i++)
                {
                    sum[i] = new double[matrix.GetLength(0)];

                    // for each row
                    for (int j = 0; j < matrix.GetLength(0); j++)
                        sum[i][j] += sum[i - 1][j] + matrix[j, i];
                }
            }
            else
            {
                throw new ArgumentException("Invalid dimension", "dimension");
            }

            return sum;
        }
        #endregion
    }
}
