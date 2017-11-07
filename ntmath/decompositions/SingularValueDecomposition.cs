namespace nt.math.decompositions
{
    using System;
    using nt.math;

    [Serializable]
    public sealed class SingularValueDecomposition : ICloneable, ISolverMatrixDecomposition<Double>
    {
        private Double[,] u; // left singular vectors
        private Double[,] v; // right singular vectors
        private Double[] s;  // singular values
        private int m;
        private int n;
        private bool swapped;

        private int[] si; // sorting order

        private const Double eps = 2 * Constants.DoubleEpsilon;
        private const Double tiny = Constants.DoubleSmall;

        Double? determinant;
        Double? lndeterminant;
        Double? pseudoDeterminant;
        Double? lnpseudoDeterminant;

        #region Properties
        public Double Condition
        {
            get { return s[0] / s[System.Math.Min(m, n) - 1]; }
        }

        public Double Threshold
        {
            get { return Constants.DoubleEpsilon * System.Math.Max(m, n) * s[0]; }
        }

        public Double TwoNorm
        {
            get { return s[0]; }
        }

        public int Rank
        {
            get
            {
                Double tol = System.Math.Max(m, n) * s[0] * eps;
                int r = 0;

                for (int i = 0; i < s.Length; i++)
                    if (s[i] > tol) r++;

                return r;
            }
        }

        public bool IsSingular
        {
            get { return Rank < Math.Min(m, n); }
        }

        public Double[] Diagonal
        {
            get { return this.s; }
        }

        public Double[,] DiagonalMatrix
        {
            get { return Matrix.Diagonal(s); }
        }

        public Double[,] RightSingularVectors
        {
            get { return v; }
        }

        public Double[,] LeftSingularVectors
        {
            get { return u; }
        }

        public int[] Ordering
        {
            get { return si; }
        }

        public Double AbsoluteDeterminant
        {
            get
            {
                if (!determinant.HasValue)
                {
                    Double det = 1;

                    for (int i = 0; i < s.Length; i++)
                        det *= s[i];

                    determinant = det;
                }

                return determinant.Value;
            }
        }

        public Double LogDeterminant
        {
            get
            {
                if (!lndeterminant.HasValue)
                {
                    double det = 0;

                    for (int i = 0; i < s.Length; i++)
                        det += Math.Log(s[i]);

                    lndeterminant = (Double)det;
                }

                return lndeterminant.Value;
            }
        }

        public Double PseudoDeterminant
        {
            get
            {
                if (!pseudoDeterminant.HasValue)
                {
                    Double det = 1;

                    for (int i = 0; i < s.Length; i++)
                        if (s[i] != 0) det *= s[i];

                    pseudoDeterminant = det;
                }

                return pseudoDeterminant.Value;
            }
        }

        public Double LogPseudoDeterminant
        {
            get
            {
                if (!lnpseudoDeterminant.HasValue)
                {
                    double det = 0;

                    for (int i = 0; i < s.Length; i++)
                        if (s[i] != 0) det += Math.Log(s[i]);

                    lnpseudoDeterminant = (Double)det;
                }

                return lnpseudoDeterminant.Value;
            }
        }
        #endregion

        #region Constructions
        public SingularValueDecomposition(Double[,] value) : this(value, true, true)
        {
        }

        public SingularValueDecomposition(Double[,] value, bool computeLeftSingularVectors, bool computeRightSingularVectors)
            : this(value, computeLeftSingularVectors, computeRightSingularVectors, false)
        {
        }

        public SingularValueDecomposition(Double[,] value, bool computeLeftSingularVectors, bool computeRightSingularVectors, bool autoTranspose)
            : this(value, computeLeftSingularVectors, computeRightSingularVectors, autoTranspose, false)
        {
        }

        public unsafe SingularValueDecomposition(Double[,] value, bool computeLeftSingularVectors, bool computeRightSingularVectors, bool autoTranspose, bool inPlace)
        {
            if (value == null)
                throw new ArgumentNullException("value", "Matrix cannot be null.");

            Double[,] a;
            m = value.GetLength(0); // rows
            n = value.GetLength(1); // cols

            if (m == 0 || n == 0)
                throw new ArgumentException("Matrix does not have any rows or columns.", "value");

            if (m < n) // Check if we are violating JAMA's assumption
            {
                if (!autoTranspose) // Yes, check if we should correct it
                {                   
                    System.Diagnostics.Trace.WriteLine("WARNING: Computing SVD on a matrix with more columns than rows.");

                    // Proceed anyway
                    a = inPlace ? value : (Double[,])value.Clone();
                }
                else
                {
                    // Transposing and swapping
                    a = value.Transpose(inPlace && m == n);
                    m = value.GetLength(1);
                    n = value.GetLength(0);
                    swapped = true;

                    bool aux = computeLeftSingularVectors;
                    computeLeftSingularVectors = computeRightSingularVectors;
                    computeRightSingularVectors = aux;
                }
            }
            else
            {
                // Input matrix is ok
                a = inPlace ? value : (Double[,])value.Clone();
            }

            int nu = System.Math.Min(m, n);
            int ni = System.Math.Min(m + 1, n);

            s = new Double[ni];
            u = new Double[m, nu];
            v = new Double[n, n];

            Double[] e = new Double[n];
            Double[] work = new Double[m];
            bool wantu = computeLeftSingularVectors;
            bool wantv = computeRightSingularVectors;

            fixed (Double* U = u)
            fixed (Double* V = v)
            fixed (Double* A = a)
            {
                // Will store ordered sequence of indices after sorting.
                si = new int[ni]; 
                
                for (int i = 0; i < ni; i++) 
                    si[i] = i;

                // Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
                int nct = System.Math.Min(m - 1, n);
                int nrt = System.Math.Max(0, System.Math.Min(n - 2, m));
                int mrc = System.Math.Max(nct, nrt);

                for (int k = 0; k < mrc; k++)
                {
                    if (k < nct)
                    {
                        // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
                        // Compute 2-norm of k-th column without under/overflow.
                        s[k] = 0;

                        for (int i = k; i < m; i++)
                            s[k] = Tools.Hypotenuse(s[k], a[i, k]);

                        if (s[k] != 0)
                        {
                            if (a[k, k] < 0)
                                s[k] = -s[k];

                            for (int i = k; i < m; i++)
                                a[i, k] /= s[k];

                            a[k, k] += 1;
                        }

                        s[k] = -s[k];
                    }

                    for (int j = k + 1; j < n; j++)
                    {
                        Double* ptr_ak = A + k * n + k; // A[k,k]
                        Double* ptr_aj = A + k * n + j; // A[k,j]

                        if ((k < nct) & (s[k] != 0))
                        {
                            // Apply the transformation.
                            Double t = 0;
                            Double* ak = ptr_ak;
                            Double* aj = ptr_aj;

                            for (int i = k; i < m; i++)
                            {
                                t += (*ak) * (*aj);
                                ak += n; aj += n;
                            }

                            t = -t / *ptr_ak;
                            ak = ptr_ak;
                            aj = ptr_aj;

                            for (int i = k; i < m; i++)
                            {
                                *aj += t * (*ak);
                                ak += n; aj += n;
                            }
                        }

                        // Place the k-th row of A into e for the subsequent calculation of the row transformation.
                        e[j] = *ptr_aj;
                    }

                    if (wantu & (k < nct))
                    {
                        // Place the transformation in U for subsequent back
                        // multiplication.
                        for (int i = k; i < m; i++)
                            u[i, k] = a[i, k];
                    }

                    if (k < nrt)
                    {
                        // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
                        // Compute 2-norm without under/overflow.
                        e[k] = 0;
                        for (int i = k + 1; i < n; i++)
                            e[k] = Tools.Hypotenuse(e[k], e[i]);

                        if (e[k] != 0)
                        {
                            if (e[k + 1] < 0)
                                e[k] = -e[k];

                            for (int i = k + 1; i < n; i++)
                                e[i] /= e[k];

                            e[k + 1] += 1;
                        }

                        e[k] = -e[k];

                        if ((k + 1 < m) & (e[k] != 0))
                        {
                            // Apply the transformation.
                            for (int i = k + 1; i < m; i++)
                                work[i] = 0;

                            int k1 = k + 1;

                            for (int i = k1; i < m; i++)
                            {
                                Double* ai = A + (i * n) + k1;

                                for (int j = k1; j < n; j++, ai++)
                                    work[i] += e[j] * (*ai);
                            }

                            for (int j = k1; j < n; j++)
                            {
                                Double t = -e[j] / e[k1];
                                Double* aj = A + (k1 * n) + j;

                                for (int i = k1; i < m; i++, aj += n)
                                    *aj += t * work[i];
                            }
                        }

                        if (wantv)
                        {
                            // Place the transformation in V for subsequent back multiplication.
                            for (int i = k + 1; i < n; i++)
                                v[i, k] = e[i];
                        }
                    }
                }

                // Set up the final bidiagonal matrix or order p.
                int p = System.Math.Min(n, m + 1);

                if (nct < n) s[nct] = a[nct, nct];
                if (m < p) s[p - 1] = 0;
                if (nrt + 1 < p) e[nrt] = a[nrt, p - 1];

                e[p - 1] = 0;

                // If required, generate U.
                if (wantu)
                {
                    for (int j = nct; j < nu; j++)
                    {
                        for (int i = 0; i < m; i++)
                            u[i, j] = 0;

                        u[j, j] = 1;
                    }

                    for (int k = nct - 1; k >= 0; k--)
                    {
                        if (s[k] != 0)
                        {
                            Double* ptr_uk = U + k * nu + k; // u[k,k]
                            Double* uk, uj;

                            for (int j = k + 1; j < nu; j++)
                            {
                                Double* ptr_uj = U + k * nu + j; // u[k,j]
                                Double t = 0;

                                uk = ptr_uk;
                                uj = ptr_uj;

                                for (int i = k; i < m; i++)
                                {
                                    t += *uk * *uj;
                                    uk += nu; uj += nu;
                                }

                                t = -t / *ptr_uk;
                                uk = ptr_uk;
                                uj = ptr_uj;

                                for (int i = k; i < m; i++)
                                {
                                    *uj += t * (*uk);
                                    uk += nu; uj += nu;
                                }
                            }

                            uk = ptr_uk;

                            for (int i = k; i < m; i++)
                            {
                                *uk = -(*uk);
                                uk += nu;
                            }

                            u[k, k] = 1 + u[k, k];

                            for (int i = 0; i < k - 1; i++)
                                u[i, k] = 0;
                        }
                        else
                        {
                            for (int i = 0; i < m; i++)
                                u[i, k] = 0;
                            u[k, k] = 1;
                        }
                    }
                }

                // If required, generate V.
                if (wantv)
                {
                    for (int k = n - 1; k >= 0; k--)
                    {
                        if ((k < nrt) & (e[k] != 0))
                        {
                            // TODO: The following is a pseudo correction to make SVD
                            //  work on matrices with n > m (less rows than columns).

                            // For the proper correction, compute the decomposition of the
                            //  transpose of A and swap the left and right eigenvectors

                            // Original line:
                            //   for (int j = k + 1; j < nu; j++)
                            // Pseudo correction:
                            //   for (int j = k + 1; j < n; j++)

                            for (int j = k + 1; j < n; j++) // pseudo-correction
                            {
                                Double* ptr_vk = V + (k + 1) * n + k; // v[k + 1, k]
                                Double* ptr_vj = V + (k + 1) * n + j; // v[k + 1, j]

                                Double t = 0;
                                Double* vk = ptr_vk;
                                Double* vj = ptr_vj;

                                for (int i = k + 1; i < n; i++)
                                {
                                    t += *vk * *vj;
                                    vk += n; vj += n;
                                }

                                t = -t / *ptr_vk;
                                vk = ptr_vk;
                                vj = ptr_vj;

                                for (int i = k + 1; i < n; i++)
                                {
                                    *vj += t * (*vk);
                                    vk += n; vj += n;
                                }
                            }
                        }

                        for (int i = 0; i < n; i++)
                            v[i, k] = 0;

                        v[k, k] = 1;
                    }
                }

                // Main iteration loop for the singular values.
                int pp = p - 1;
                int iter = 0;
                while (p > 0)
                {
                    int k, kase;

                    // Here is where a test for too many iterations would go.

                    // This section of the program inspects for
                    // negligible elements in the s and e arrays.  On
                    // completion the variables kase and k are set as follows.

                    // kase = 1     if s(p) and e[k-1] are negligible and k<p
                    // kase = 2     if s(k) is negligible and k<p
                    // kase = 3     if e[k-1] is negligible, k<p, and
                    //              s(k), ..., s(p) are not negligible (qr step).
                    // kase = 4     if e(p-1) is negligible (convergence).

                    for (k = p - 2; k >= -1; k--)
                    {
                        if (k == -1)
                            break;

                        var alpha = tiny + eps * (System.Math.Abs(s[k]) + System.Math.Abs(s[k + 1]));

                        if (System.Math.Abs(e[k]) <= alpha || Double.IsNaN(e[k]))
                        {
                            e[k] = 0;
                            break;
                        }
                    }

                    if (k == p - 2)
                    {
                        kase = 4;
                    }
                    else
                    {
                        int ks;

                        for (ks = p - 1; ks >= k; ks--)
                        {
                            if (ks == k)
                                break;

                            Double t = (ks != p ? System.Math.Abs(e[ks]) : 0) + (ks != k + 1 ? System.Math.Abs(e[ks - 1]) : 0);

                            if (System.Math.Abs(s[ks]) <= tiny + eps * t)
                            {
                                s[ks] = 0;
                                break;
                            }
                        }

                        if (ks == k)
                            kase = 3;
                        else if (ks == p - 1)
                            kase = 1;
                        else
                        {
                            kase = 2;
                            k = ks;
                        }
                    }

                    k++;

                    // Perform the task indicated by kase.
                    switch (kase)
                    {
                        // Deflate negligible s(p).
                        case 1:
                            {
                                Double f = e[p - 2];

                                e[p - 2] = 0;

                                for (int j = p - 2; j >= k; j--)
                                {
                                    Double t = Tools.Hypotenuse(s[j], f);
                                    Double cs = s[j] / t;
                                    Double sn = f / t;

                                    s[j] = t;

                                    if (j != k)
                                    {
                                        f = -sn * e[j - 1];
                                        e[j - 1] = cs * e[j - 1];
                                    }

                                    if (wantv)
                                    {
                                        for (int i = 0; i < n; i++)
                                        {
                                            t = cs * v[i, j] + sn * v[i, p - 1];
                                            v[i, p - 1] = -sn * v[i, j] + cs * v[i, p - 1];
                                            v[i, j] = t;
                                        }
                                    }
                                }
                            }
                            break;

                        // Split at negligible s(k).
                        case 2:
                            {
                                Double f = e[k - 1];

                                e[k - 1] = 0;

                                for (int j = k; j < p; j++)
                                {
                                    Double t = Tools.Hypotenuse(s[j], f);
                                    Double cs = s[j] / t;
                                    Double sn = f / t;

                                    s[j] = t;
                                    f = -sn * e[j];
                                    e[j] = cs * e[j];

                                    if (wantu)
                                    {
                                        for (int i = 0; i < m; i++)
                                        {
                                            t = cs * u[i, j] + sn * u[i, k - 1];
                                            u[i, k - 1] = -sn * u[i, j] + cs * u[i, k - 1];
                                            u[i, j] = t;
                                        }
                                    }
                                }
                            }
                            break;

                        // Perform one qr step.
                        case 3:
                            {
                                // Calculate the shift.
                                Double scale = System.Math.Max(System.Math.Max(System.Math.Max(System.Math.Max(System.Math.Abs(s[p - 1]), System.Math.Abs(s[p - 2])), System.Math.Abs(e[p - 2])), System.Math.Abs(s[k])), System.Math.Abs(e[k]));
                                Double sp = s[p - 1] / scale;
                                Double spm1 = s[p - 2] / scale;
                                Double epm1 = e[p - 2] / scale;
                                Double sk = s[k] / scale;
                                Double ek = e[k] / scale;
                                Double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2;
                                Double c = (sp * epm1) * (sp * epm1);
                                double shift = 0;

                                if ((b != 0) | (c != 0))
                                {
                                    if (b < 0)
                                        shift = -System.Math.Sqrt(b * b + c);
                                    else
                                        shift = System.Math.Sqrt(b * b + c);

                                    shift = c / (b + shift);
                                }

                                Double f = (sk + sp) * (sk - sp) + (Double)shift;
                                Double g = sk * ek;

                                // Chase zeros.
                                for (int j = k; j < p - 1; j++)
                                {
                                    Double t = Tools.Hypotenuse(f, g);
                                    Double cs = f / t;
                                    Double sn = g / t;

                                    if (j != k) e[j - 1] = t;

                                    f = cs * s[j] + sn * e[j];
                                    e[j] = cs * e[j] - sn * s[j];
                                    g = sn * s[j + 1];
                                    s[j + 1] = cs * s[j + 1];

                                    if (wantv)
                                    {
                                        unsafe
                                        {
                                            fixed (Double* ptr_vj = &v[0, j])
                                            {
                                                Double* vj = ptr_vj;
                                                Double* vj1 = ptr_vj + 1;

                                                for (int i = 0; i < n; i++)
                                                {
                                                    /*t = cs * v[i, j] + sn * v[i, j + 1];
                                                    v[i, j + 1] = -sn * v[i, j] + cs * v[i, j + 1];
                                                    v[i, j] = t;*/

                                                    Double vij = *vj;
                                                    Double vij1 = *vj1;

                                                    t = cs * vij + sn * vij1;
                                                    *vj1 = -sn * vij + cs * vij1;
                                                    *vj = t;
                                                    vj += n;
                                                    vj1 += n;
                                                }
                                            }
                                        }
                                    }

                                    t = Tools.Hypotenuse(f, g);
                                    cs = f / t;
                                    sn = g / t;
                                    s[j] = t;
                                    f = cs * e[j] + sn * s[j + 1];
                                    s[j + 1] = -sn * e[j] + cs * s[j + 1];
                                    g = sn * e[j + 1];
                                    e[j + 1] = cs * e[j + 1];

                                    if (wantu && (j < m - 1))
                                    {
                                        fixed (Double* ptr_uj = &u[0, j])
                                        {
                                            Double* uj = ptr_uj;
                                            Double* uj1 = ptr_uj + 1;

                                            for (int i = 0; i < m; i++)
                                            {
                                                /* t = cs * u[i, j] + sn * u[i, j + 1];
                                                 u[i, j + 1] = -sn * u[i, j] + cs * u[i, j + 1];
                                                 u[i, j] = t;*/

                                                Double uij = *uj;
                                                Double uij1 = *uj1;

                                                t = cs * uij + sn * uij1;
                                                *uj1 = -sn * uij + cs * uij1;
                                                *uj = t;
                                                uj += nu; 
                                                uj1 += nu;
                                            }
                                        }
                                    }

                                }

                                e[p - 2] = f;
                                iter = iter + 1;
                            }
                            break;

                        // Convergence.
                        case 4:
                            {
                                // Make the singular values positive.
                                if (s[k] <= 0)
                                {
                                    s[k] = (s[k] < 0 ? -s[k] : 0);
                                    if (wantv)
                                    {
                                        for (int i = 0; i <= pp; i++)
                                            v[i, k] = -v[i, k];
                                    }
                                }

                                // Order the singular values.
                                while (k < pp)
                                {
                                    if (s[k] >= s[k + 1])
                                        break;

                                    Double t = s[k];

                                    s[k] = s[k + 1];
                                    s[k + 1] = t;

                                    int ti = si[k];

                                    si[k] = si[k + 1];
                                    si[k + 1] = ti;

                                    if (wantv && (k < n - 1))
                                    {
                                        for (int i = 0; i < n; i++)
                                        {
                                            t = v[i, k + 1];
                                            v[i, k + 1] = v[i, k];
                                            v[i, k] = t;
                                        }
                                    }

                                    if (wantu && (k < m - 1))
                                    {
                                        for (int i = 0; i < m; i++)
                                        {
                                            t = u[i, k + 1];
                                            u[i, k + 1] = u[i, k];
                                            u[i, k] = t;
                                        }
                                    }

                                    k++;
                                }

                                iter = 0;
                                p--;
                            }
                            break;
                    }
                }
            }

            // If we are violating JAMA's assumption about 
            // the input dimension, we need to swap u and v.
            if (swapped)
            {
                Double[,] temp = this.u;

                this.u = this.v;
                this.v = temp;
            }
        }
        #endregion

        #region Solve
        public Double[,] Solve(Double[,] value)
        {
            Double[,] Y = value;
            Double e = this.Threshold;
            int scols = s.Length;
            var Ls = new Double[scols, scols];

            for (int i = 0; i < s.Length; i++)
            {
                if (System.Math.Abs(s[i]) <= e) Ls[i, i] = 0;
                else Ls[i, i] = 1 / s[i];
            }

            //(V x L*) x Ut x Y
            var VL = v.Multiply(Ls);

            //(V x L* x Ut) x Y
            int vrows = v.GetLength(0);
            int urows = u.GetLength(0);
            var VLU = new Double[vrows, scols];

            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    Double sum = 0;

                    for (int k = 0; k < urows; k++)
                        sum += VL[i, k] * u[j, k];

                    VLU[i, j] = sum;
                }
            }

            //(V x L* x Ut x Y)
            return VLU.Multiply(Y);
        }

        public Double[,] SolveForDiagonal(Double[] value)
        {
            Double[] Y = value;
            Double e = this.Threshold;
            int scols = s.Length;
            var Ls = new Double[scols, scols];

            for (int i = 0; i < s.Length; i++)
            {
                if (System.Math.Abs(s[i]) <= e) Ls[i, i] = 0;
                else Ls[i, i] = 1 / s[i];
            }

            //(V x L*) x Ut x Y
            Double[,] VL = v.Multiply(Ls);

            //(V x L* x Ut) x Y
            int vrows = v.GetLength(0);
            int urows = u.GetLength(0);
            var VLU = new Double[vrows, scols];

            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    Double sum = 0;

                    for (int k = 0; k < urows; k++)
                        sum += VL[i, k] * u[j, k];

                    VLU[i, j] = sum;
                }
            }

            //(V x L* x Ut x Y)
            return VLU.MultiplyByDiagonal(Y);
        }

        public Double[] Solve(Double[] value)
        {
            Double e = this.Threshold;
            var Y = value;
            int scols = s.Length;
            var Ls = new Double[scols, scols];

            for (int i = 0; i < s.Length; i++)
            {
                if (System.Math.Abs(s[i]) <= e) Ls[i, i] = 0;
                else Ls[i, i] = 1 / s[i];
            }

            //(V x L*) x Ut x Y
            var VL = v.Multiply(Ls);

            //(V x L* x Ut) x Y
            int urows = u.GetLength(0);
            int vrows = v.GetLength(0);
            var VLU = new Double[vrows, urows];

            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    Double sum = 0;

                    for (int k = 0; k < scols; k++)
                        sum += VL[i, k] * u[j, k];

                    VLU[i, j] = sum;
                }
            }

            //(V x L* x Ut x Y)
            return VLU.Multiply(Y);
        }
        #endregion

        public Double[,] Inverse()
        {
            Double e = this.Threshold;
            // X = V*S^-1
            int vrows = v.GetLength(0);
            int vcols = v.GetLength(1);
            var X = new Double[vrows, s.Length];

            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < vcols; j++)
                {
                    if (System.Math.Abs(s[j]) > e)
                        X[i, j] = v[i, j] / s[j];
                }
            }

            // Y = X*U'
            int urows = u.GetLength(0);
            int ucols = u.GetLength(1);
            var Y = new Double[vrows, urows];

            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    Double sum = 0;

                    for (int k = 0; k < ucols; k++)
                        sum += X[i, k] * u[j, k];

                    Y[i, j] = sum;
                }
            }

            return Y;
        }

        #region ICloneable Members
        private SingularValueDecomposition()
        {
        }

        public object Clone()
        {
            var svd = new SingularValueDecomposition();

            svd.m = m;
            svd.n = n;
            svd.s = (Double[])s.Clone();
            svd.si = (int[])si.Clone();
            svd.swapped = swapped;

            if (u != null) svd.u = (Double[,])u.Clone();
            if (v != null) svd.v = (Double[,])v.Clone();

            return svd;
        }
        #endregion
    }
}
