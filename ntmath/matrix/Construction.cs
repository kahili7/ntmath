namespace nt.math
{
    using System;
    using System.Collections.Generic;

    public static partial class Matrix
    {
        #region Generic Matrix
        public static T[,] Create<T>(int rows, int cols, T value)
        {
            if (rows < 0)
            {
                throw new ArgumentOutOfRangeException("rows", rows, "Number of rows must be a positive integer.");
            }

            if (cols < 0)
            {
                throw new ArgumentOutOfRangeException("cols", cols, "Number of columns must be a positive integer.");
            }

            T[,] matrix = new T[rows, cols];

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    matrix[i, j] = value;

            return matrix;
        }

        public static T[,] Create<T>(int size, T value)
        {
            if (size < 0)
            {
                throw new ArgumentOutOfRangeException("size", size, "Square matrix's size must be a positive integer.");
            }

            return Create(size, size, value);
        }

        public static T[,] Create<T>(int rows, int cols)
        {
            return Create(rows, cols, default(T));
        }

        public static T[,] Create<T>(int size)
        {
            return Create(size, default(T));
        }

        public static T[][] Jagged<T>(int rows, int cols, T value)
        {
            if (rows < 0)
            {
                throw new ArgumentOutOfRangeException("rows", rows,
                    "Number of rows must be a positive integer.");
            }

            if (cols < 0)
            {
                throw new ArgumentOutOfRangeException("cols", cols,
                    "Number of columns must be a positive integer.");
            }

            T[][] matrix = new T[rows][];

            for (int i = 0; i < rows; i++)
            {
                var row = matrix[i] = new T[cols];

                for (int j = 0; j < row.Length; j++)
                    row[j] = value;
            }

            return matrix;
        }

        public static T[,] Jagged<T>(int size, T value)
        {
            if (size < 0)
            {
                throw new ArgumentOutOfRangeException("size", size, "Square matrix's size must be a positive integer.");
            }

            return Create(size, size, value);
        }

        public static T[][] Jagged<T>(int rows, int cols)
        {
            return Jagged(rows, cols, default(T));
        }

        public static T[,] Jagged<T>(int size)
        {
            return Create(size, default(T));
        }
        #endregion

        #region Create Vector
        public static T[,] RowVector<T>(params T[] values)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            T[,] matrix = new T[1, values.Length];

            for (int i = 0; i < values.Length; i++)
                matrix[0, i] = values[i];

            return matrix;
        }

        public static T[,] ColumnVector<T>(params T[] values)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            T[,] matrix = new T[values.Length, 1];

            for (int i = 0; i < values.Length; i++)
                matrix[i, 0] = values[i];

            return matrix;
        }

        public static T[] Vector<T>(int n, T[] values)
        {
            T[] vector = new T[n];

            if (values != null)
            {
                for (int i = 0; i < values.Length; i++)
                    vector[i] = values[i];
            }

            return vector;
        }

        public static T[] Vector<T>(int n, T value)
        {
            T[] vector = new T[n];

            for (int i = 0; i < n; i++)
                vector[i] = value;

            return vector;
        }

        public static double[] Vector(double a, double b, double increment = 1)
        {
            List<double> list = new List<double>();

            for (double i = a; i < b; i += increment)
                list.Add(i);

            if (list[list.Count - 1] != b)
                list.Add(b);

            return list.ToArray();
        }

        public static int[] Vector(int a, int b, int increment = 1)
        {
            List<int> list = new List<int>();

            for (int i = a; i < b; i += increment)
                list.Add(i);

            if (list[list.Count - 1] != b)
                list.Add(b);

            return list.ToArray();
        }

        public static double[] Vector(double a, double b, int points)
        {
            double[] list = new double[points];

            double increment = (b - a) / points;

            for (int i = 0; i < list.Length; i++)
                list[i] = increment * i;

            return list;
        }
        #endregion

        #region Combine
        public static T[] Concatenate<T>(this T[] a, T[] b)
        {
            T[] r = new T[a.Length + b.Length];

            for (int i = 0; i < a.Length; i++)
                r[i] = a[i];

            for (int i = 0; i < b.Length; i++)
                r[i + a.Length] = b[i];

            return r;
        }

        public static T[] Concatenate<T>(this T[] vector, T element)
        {
            T[] r = new T[vector.Length + 1];
            for (int i = 0; i < vector.Length; i++)
                r[i] = vector[i];

            r[vector.Length] = element;

            return r;
        }

        public static T[] Concatenate<T>(this T element, T[] vector)
        {
            T[] r = new T[vector.Length + 1];

            r[0] = element;

            for (int i = 0; i < vector.Length; i++)
                r[i + 1] = vector[i];

            return r;
        }

        public static T[,] Concatenate<T>(this T[,] matrix, T[] vector)
        {
            return matrix.InsertColumn(vector);
        }

        public static T[,] Concatenate<T>(this T[,] a, T[,] b)
        {
            return Concatenate(new[] { a, b });
        }

        public static T[][] Concatenate<T>(this T[][] a, T[][] b)
        {
            return Concatenate(new[] { a, b });
        }

        public static T[,] Concatenate<T>(params T[][,] matrices)
        {
            int rows = 0;
            int cols = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                cols += matrices[i].GetLength(1);

                if (matrices[i].GetLength(0) > rows)
                    rows = matrices[i].GetLength(0);
            }

            T[,] r = new T[rows, cols];


            int c = 0;

            for (int k = 0; k < matrices.Length; k++)
            {
                int currentRows = matrices[k].GetLength(0);
                int currentCols = matrices[k].GetLength(1);

                for (int j = 0; j < currentCols; j++)
                {
                    for (int i = 0; i < currentRows; i++)
                    {
                        r[i, c] = matrices[k][i, j];
                    }

                    c++;
                }
            }

            return r;
        }

        public static T[][] Concatenate<T>(params T[][][] matrices)
        {
            int rows = 0;
            int cols = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                cols += matrices[i][0].Length;

                if (matrices[i].Length > rows)
                    rows = matrices[i].Length;
            }

            T[][] r = new T[rows][];

            for (int i = 0; i < r.Length; i++)
                r[i] = new T[cols];

            int c = 0;

            for (int k = 0; k < matrices.Length; k++)
            {
                int currentRows = matrices[k].Length;
                int currentCols = matrices[k][0].Length;

                for (int j = 0; j < currentCols; j++)
                {
                    for (int i = 0; i < currentRows; i++)
                    {
                        r[i][c] = matrices[k][i][j];
                    }

                    c++;
                }
            }

            return r;
        }

        public static T[] Concatenate<T>(params T[][] vectors)
        {
            int size = 0;
            for (int i = 0; i < vectors.Length; i++)
                size += vectors[i].Length;

            T[] r = new T[size];

            int c = 0;
            for (int i = 0; i < vectors.Length; i++)
                for (int j = 0; j < vectors[i].Length; j++)
                    r[c++] = vectors[i][j];

            return r;
        }

        public static T[,] Stack<T>(params T[][] vectors)
        {
            return vectors.ToMatrix();
        }

        public static T[,] Stack<T>(params T[] elements)
        {
            return elements.Transpose();
        }

        public static T[,] Stack<T>(T[] vector, T element)
        {
            return vector.Concatenate(element).Transpose();
        }

        public static T[,] Stack<T>(params T[][,] matrices)
        {
            int rows = 0;
            int cols = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                rows += matrices[i].GetLength(0);

                if (matrices[i].GetLength(1) > cols)
                    cols = matrices[i].GetLength(1);
            }

            T[,] r = new T[rows, cols];

            int c = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                for (int j = 0; j < matrices[i].GetLength(0); j++)
                {
                    for (int k = 0; k < matrices[i].GetLength(1); k++)
                        r[c, k] = matrices[i][j, k];

                    c++;
                }
            }

            return r;
        }

        public static T[,] Stack<T>(T[,] matrix, T[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            T[,] r = new T[rows + 1, cols];

            Array.Copy(matrix, r, matrix.Length);

            for (int i = 0; i < vector.Length; i++)
                r[rows, i] = vector[i];

            return r;
        }

        public static T[][] Stack<T>(params T[][][] matrices)
        {
            int rows = 0;
            int cols = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                rows += matrices[i].Length;

                if (matrices[i].Length == 0)
                    continue;

                if (matrices[i][0].Length > cols)
                    cols = matrices[i][0].Length;
            }

            T[][] r = new T[rows][];

            for (int i = 0; i < rows; i++)
                r[i] = new T[cols];

            int c = 0;

            for (int i = 0; i < matrices.Length; i++)
            {
                for (int j = 0; j < matrices[i].Length; j++)
                {
                    for (int k = 0; k < matrices[i][j].Length; k++)
                        r[c][k] = matrices[i][j][k];

                    c++;
                }
            }

            return r;
        }
        #endregion

        #region Diagonal Matrix
        public static T[,] Diagonal<T>(int size, T value)
        {
            if (size < 0)
            {
                throw new ArgumentOutOfRangeException("size", size, "Square matrix's size must be a positive integer.");
            }

            T[,] matrix = new T[size, size];

            for (int i = 0; i < size; i++)
                matrix[i, i] = value;

            return matrix;
        }

        public static T[,] Diagonal<T>(int rows, int cols, T value)
        {
            if (rows < 0)
            {
                throw new ArgumentOutOfRangeException("rows", rows, "Number of rows must be a positive integer.");
            }

            if (cols < 0)
            {
                throw new ArgumentOutOfRangeException("cols", cols, "Number of columns must be a positive integer.");
            }

            T[,] matrix = new T[rows, cols];
            int min = Math.Min(rows, cols);

            for (int i = 0; i < min; i++)
                matrix[i, i] = value;

            return matrix;
        }

        public static T[,] Diagonal<T>(T[] values)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            T[,] matrix = new T[values.Length, values.Length];

            for (int i = 0; i < values.Length; i++)
                matrix[i, i] = values[i];

            return matrix;
        }

        public static T[][] JaggedDiagonal<T>(T[] values)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            T[][] matrix = new T[values.Length][];

            for (int i = 0; i < values.Length; i++)
            {
                matrix[i] = new T[values.Length];
                matrix[i][i] = values[i];
            }

            return matrix;
        }

        public static T[][] JaggedDiagonal<T>(int size, T value)
        {
            if (size < 0)
            {
                throw new ArgumentOutOfRangeException("size", size, "Square matrix's size must be a positive integer.");
            }

            var matrix = new T[size][];

            for (int i = 0; i < matrix.Length; i++)
            {
                matrix[i] = new T[size];
                matrix[i][i] = value;
            }

            return matrix;
        }

        public static T[,] Diagonal<T>(int size, T[] values)
        {
            return Diagonal(size, size, values);
        }

        public static T[,] Diagonal<T>(int rows, int cols, T[] values)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            if (rows < 0)
            {
                throw new ArgumentOutOfRangeException("rows", rows, "Number of rows must be a positive integer.");
            }

            if (cols < 0)
            {
                throw new ArgumentOutOfRangeException("cols", cols, "Number of columns must be a positive integer.");
            }

            T[,] matrix = new T[rows, cols];

            for (int i = 0; i < values.Length; i++)
                matrix[i, i] = values[i];

            return matrix;
        }
        #endregion

        #region Special Matrix
        public static double[,] Identity(int size)
        {
            return Diagonal(size, 1.0);
        }

        public static double[][] JaggedIdentity(int size)
        {
            return JaggedDiagonal(size, 1.0);
        }

        public static double[,] Magic(int size)
        {
            if (size < 3) throw new ArgumentOutOfRangeException("size", size, "The square size must be greater or equal to 3.");

            double[,] matrix = new double[size, size];

            // First algorithm: Odd order
            if ((size % 2) == 1)
            {
                int a = (size + 1) / 2;
                int b = (size + 1);

                for (int j = 0; j < size; j++)
                    for (int i = 0; i < size; i++)
                        matrix[i, j] = size * ((i + j + a) % size) + ((i + 2 * j + b) % size) + 1;
            }
            // Second algorithm: Even order (double)
            else if ((size % 4) == 0)
            {
                for (int j = 0; j < size; j++)
                    for (int i = 0; i < size; i++)
                        if (((i + 1) / 2) % 2 == ((j + 1) / 2) % 2)
                            matrix[i, j] = size * size - size * i - j;
                        else
                            matrix[i, j] = size * i + j + 1;
            }
            // Third algorithm: Even order (single)
            else
            {
                int n = size / 2;
                int p = (size - 2) / 4;
                double t;
                double[,] block = Matrix.Magic(n);

                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        double e = block[i, j];

                        matrix[i, j] = e;
                        matrix[i, j + n] = e + 2 * n * n;
                        matrix[i + n, j] = e + 3 * n * n;
                        matrix[i + n, j + n] = e + n * n;
                    }
                }

                for (int i = 0; i < n; i++)
                {
                    // Swap M[i,j] and M[i+n,j]
                    for (int j = 0; j < p; j++)
                    {
                        t = matrix[i, j];
                        matrix[i, j] = matrix[i + n, j];
                        matrix[i + n, j] = t;
                    }
                    for (int j = size - p + 1; j < size; j++)
                    {
                        t = matrix[i, j];
                        matrix[i, j] = matrix[i + n, j];
                        matrix[i + n, j] = t;
                    }
                }

                // Continue swapping in the boundary
                t = matrix[p, 0];
                matrix[p, 0] = matrix[p + n, 0];
                matrix[p + n, 0] = t;

                t = matrix[p, p];
                matrix[p, p] = matrix[p + n, p];
                matrix[p + n, p] = t;
            }

            return matrix; // return the magic square.
        }

        public static double[,] Centering(int size)
        {
            if (size < 0) throw new ArgumentOutOfRangeException("size", size, "The size of the centering matrix must be a positive integer.");

            double[,] C = Matrix.Create(size, -1.0 / size);

            for (int i = 0; i < size; i++)
                C[i, i] = 1.0 - 1.0 / size;

            return C;
        }
        #endregion

        #region Special Vectors
        public static int[] Indices(int from, int to)
        {
            int[] vector;

            if (to > from)
            {
                vector = new int[to - from];

                for (int i = 0; i < vector.Length; i++)
                    vector[i] = from++;
            }
            else
            {
                vector = new int[from - to];

                for (int i = 0; i < vector.Length; i++)
                    vector[i] = from-- - 1;
            }

            return vector;
        }

        public static int[] Interval(int from, int to)
        {
            int[] vector;

            if (to > from)
            {
                vector = new int[to - from + 1];

                for (int i = 0; i < vector.Length; i++)
                    vector[i] = from++;
            }
            else
            {
                vector = new int[from - to + 1];

                for (int i = 0; i < vector.Length; i++)
                    vector[i] = from--;
            }

            return vector;
        }


        #endregion
    }
}
