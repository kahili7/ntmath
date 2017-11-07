namespace nt.math.statistics.filters
{
    using System;

    public enum KALMAN
    {
        N_MODIF = 1,
        NU_MODIF = 1 << 1,
        NW_MODIF = 1 << 2,
        M_MODIF = 1 << 3,
        NV_MODIF = 1 << 4,
        P_MODIF = 1 << 5,
        LOWMASK = (1 << 8) - 1,
        A_MODIF = 1 << 8,
        W_MODIF = 1 << 9,
        Q_MODIF = 1 << 10,
        MIDMASK = ((1 << 4) - 1) << 8,
        H_MODIF = 1 << 12,
        V_MODIF = 1 << 13,
        R_MODIF = 1 << 14,
        HIGHMASK = ((1 << 4) - 1) << 12,
    }

    public class EKalmanFilter
    {
        private KALMAN _flags;

        private double[,] _P;
        private double[] _x;

        private double[] _a;
        private double[] _d;
        private double[] _v;

        private double[,] _U; //матрица холески факторизации для P
        private double[,] _W;
        private double[,] _Q;
        private double[,] _H;
        private double[,] _R;

        private bool _modified;

        private int _nn;

        private bool _optimizeQ;
        private bool _optimizeVR;

        //вектор состояния системы
        protected double[] x;

        //вектор управляющих воздействий
        protected double[] u;

        //вектор наблюдений
        protected double[] z;

        //вектор инновационный
        protected double[] dz;

        //матрица якоби частных производных относительно x
        protected double[,] A;

        //матрица якоби частных производных относительно w
        protected double[,] W;

        //ковариационная матрица процесса
        protected double[,] Q;

        //матрица наблюдений;
        //матрица якоби частных производных относительно x
        protected double[,] H;

        //матрица якоби частных производных относительно v
        protected double[,] V;

        //ковариационная матрица шума измерений
        protected double[,] R;

        protected int n;
        protected int nu;
        protected int nw;
        protected int m;
        protected int nv;

        public EKalmanFilter()
        {
        }

        public EKalmanFilter(int n, int nu, int nw, int m, int nv)
        {
            SetDim(n, nu, nw, m, nv);
        }

        public void Initialize(double[] x, double[,] P)
        {
            _x.Swap<double>(x);
            _P.Swap<double>(P);
            _flags |= KALMAN.P_MODIF;
        }

        public void TimeUpdateStep(double[] _u)
        {
            SizeUpdate();

            u.Swap<double>(_u);

            MakeCommonProcess();
            MakeAImpl();
            MakeWImpl();
            MakeQImpl();
            MakeProcess();

            if (!_optimizeQ)
            {
                if (_flags == KALMAN.Q_MODIF)
                {
                    _Q = Q;
                    Factor(_Q);
                    UpperInvert(_Q);
                }

                Q.Swap<double>(_Q);

                // W_ = W*U   n.nw = n.nw * nw.nw
                if ((_flags == KALMAN.W_MODIF) || (_flags == KALMAN.Q_MODIF))
                {
                    for (int i = 0; i < n; ++i)
                    {
                        for (int j = 0; j < nw; ++j)
                        {
                            _W[i, j] = W[i, j];

                            for (int k = 0; k < j; ++k)
                                _W[i, j] += W[i, k] * Q[j, k];
                        }

                    }

                }

                W.Swap<double>(_W);
            }

            TimeUpdate();

            if (!_optimizeQ)
            {
                Q.Swap<double>(_Q);
                W.Swap<double>(_W);
            }

            u.Swap<double>(_u);
            _flags &= ~KALMAN.MIDMASK;
        }

        public void MeasureUpdateStep(double[] _z)
        {
            SizeUpdate();

            if (m == 0) return;

            MakeCommonMeasure();
            MakeHImpl();
            MakeVImpl();
            MakeRImpl();
            MakeMeasure();

            for (int i = 0; i < m; ++i)
                dz[i] = _z[i] - z[i];

            MakeDZ();

            if (_optimizeVR)
            {
                if ((_flags == KALMAN.V_MODIF) || (_flags == KALMAN.R_MODIF))
                {
                    for (int i = 0; i < m; ++i)
                        _R[i, i] = V[i, i] * V[i, i] * R[i, i];
                }
            }
            else
            {
                if ((_flags == KALMAN.V_MODIF) || (_flags == KALMAN.R_MODIF))
                {
                    _x = Matrix.Vector<double>(nv, 0.0);

                    // R_ = V*R*V'
                    for (int i = 0; i < m; ++i)
                    {
                        // _x = row i of V*R = (V*R)(i,:)
                        for (int j = 0; j < nv; ++j)
                        {
                            _x[j] = 0.0;

                            for (int k = 0; k < nv; ++k)
                                _x[j] += V[i, k] * R[k, j];

                        }

                        // R_(i,:) = (V*R*V')(i,:) = (V*R)(i,:) * V'
                        for (int j = 0; j < m; ++j)
                        {
                            _R[i, j] = 0.0;

                            for (int k = 0; k < nv; ++k)
                                _R[i, j] += _x[k] * V[j, k];
                        }
                    }

                    // R_ = U*D*U'
                    // diag(R_) = D, upper(R_) = U, lower(R_) = junk
                    Factor(_R);

                    // lower(R_) = (inv(U))'
                    UpperInvert(_R);
                }

                if ((_flags == KALMAN.H_MODIF) || (_flags == KALMAN.V_MODIF) || (_flags == KALMAN.R_MODIF))
                {
                    // H_ = inv(U)*H    m.n = m.m * m.n
                    for (int i = 0; i < m; ++i)
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            _H[i, j] = H[i, j];

                            for (int k = i + 1; k < m; ++k)
                                _H[i, j] += _R[k, i] * H[k, j];
                        }
                    }
                }

                H.Swap<double>(_H);

                // _x = inv(U)*dz    m.1 = m.m * m.1
                _x = Matrix.Vector<double>(m, 0.0);

                for (int i = 0; i < m; ++i)
                {
                    _x[i] = dz[i];

                    for (int k = i + 1; k < m; ++k)
                        _x[i] += _R[k, i] * dz[k];
                }

                dz.Swap<double>(_x);
            }

            _x = Matrix.Vector<double>(n, 0.0); // dx : innovation
            
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                    _a[j] = H[i, j];

                MeasureUpdate(dz[i], _R[i, i]);
            }

            for (int i = 0; i < n; ++i)
                x[i] += _x[i];

            if (!_optimizeVR)
            {
                H.Swap<double>(_H);
            }

            _flags &= ~KALMAN.HIGHMASK;
        }

        public double[] Predict(double[] _u)
        {
            SizeUpdate();
            u.Swap<double>(_u);
            _x = x;

            MakeCommonProcess();
            MakeProcess();

            x.Swap<double>(_x);
            u.Swap<double>(_u);
            return _x;
        }

        public double[] Simulate()
        {
            SizeUpdate();
            _x = z;

            MakeCommonProcess();
            MakeProcess();

            z.Swap<double>(_x);
            return _x;
        }

        public double[,] CalculateP()
        {
            if (!(_flags == KALMAN.P_MODIF))
            {
                _P = Matrix.Create<double>(n);         // keep this resize

                for (int i = 0; i < n; ++i)
                {
                    _P[i, i] = _U[i, i];

                    for (int j = i + 1; j < n; ++j)
                    {
                        _P[i, j] = _U[i, j] * _U[j, j];
                        _P[i, i] += _U[i, j] * _P[i, j];

                        for (int k = j + 1; k < n; ++k)
                        {
                            _P[i, j] += _U[i, k] * _U[j, k] * _U[k, k];
                        }

                        _P[j, i] = _P[i, j];
                    }
                }
            }

            return _P;
        }

        protected void NoModification()
        {
            _modified = false;
        }

        #region Virtual functions
        protected virtual void MakeBaseA()
        {
            NoModification();
        }

        protected virtual void MakeBaseW()
        {
            NoModification();
        }

        protected virtual void MakeBaseQ()
        {
            NoModification();
        }

        protected virtual void MakeBaseH()
        {
            NoModification();
        }

        protected virtual void MakeBaseV()
        {
            NoModification();
        }

        protected virtual void MakeBaseR()
        {
            NoModification();
        }

        protected virtual void MakeCommonProcess()
        {
        }

        protected virtual void MakeCommonMeasure()
        {
        }

        protected virtual void MakeA()
        {
            NoModification();
        }

        protected virtual void MakeW()
        {
            NoModification();
        }

        protected virtual void MakeQ()
        {
            NoModification();
        }

        protected virtual void MakeH()
        {
            NoModification();
        }

        protected virtual void MakeV()
        {
            NoModification();
        }

        protected virtual void MakeR()
        {
            NoModification();
        }

        protected virtual void MakeDZ()
        {
        }

        protected virtual void SizeUpdate()
        {
            if (_flags == 0) return;

            if (_flags == KALMAN.N_MODIF)
            {
                A = Matrix.Create<double>(n);
                MakeBaseAImpl();
            }

            if ((_flags == KALMAN.N_MODIF) || (_flags == KALMAN.NW_MODIF))
            {
                _nn = n + nw;
                _a = Matrix.Vector<double>(_nn, 0.0);
                _v = Matrix.Vector<double>(_nn, 0.0);
                _d = Matrix.Vector<double>(_nn, 0.0);

                if (!_optimizeQ) _W = Matrix.Create<double>(n, nw);

                W = Matrix.Create<double>(n, nw);
                MakeBaseWImpl();
            }

            if (_flags == KALMAN.P_MODIF)
            {
                _U = Matrix.Create<double>(n, _nn);

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        _U[i, j] = _P[i, j];

                Factor(_U);
            }
            else if (_flags == KALMAN.NW_MODIF)
            {
                _P = Matrix.Create<double>(n, _nn);

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        _P[i, j] = _U[i, j];

                Factor(_P);
            }

            if (_flags == KALMAN.NW_MODIF)
            {
                if (!_optimizeQ) _Q = Matrix.Create<double>(nw);

                Q = Matrix.Create<double>(nw);
                MakeBaseQImpl();
            }

            if (m != 0)
            {
                if ((_flags == KALMAN.N_MODIF) || (_flags == KALMAN.M_MODIF))
                {
                    if (!_optimizeVR) _H = Matrix.Create<double>(m, n);

                    H = Matrix.Create<double>(m, n);
                    MakeBaseHImpl();
                }

                if ((_flags == KALMAN.M_MODIF) || (_flags == KALMAN.NV_MODIF))
                {
                    V = Matrix.Create<double>(m, nv);
                    MakeBaseVImpl();
                }

                if (_flags == KALMAN.NV_MODIF)
                {
                    R = Matrix.Create<double>(nv);
                    MakeBaseRImpl();
                }

                if (_flags == KALMAN.M_MODIF)
                {
                    _R = Matrix.Create<double>(m);
                    z = Matrix.Vector<double>(m, 0.0);
                    dz = Matrix.Vector<double>(m, 0.0);
                }
            }

            _flags &= ~KALMAN.LOWMASK;
        }

        protected virtual void MakeProcess() { }

        protected virtual void MakeMeasure() { }
        #endregion

        #region Private functions
        private void TimeUpdate()
        {
            for (int j = n - 1; j > 0; --j)
            {
                for (int i = 0; i <= j; ++i)
                    _d[i] = _U[i, j];

                for (int i = 0; i < n; ++i)
                {
                    _U[i, j] = A[i, j];

                    for (int k = 0; k < j; ++k)
                        _U[i, j] += A[i, k] * _d[k];
                }
            }

            _d[0] = _U[0, 0];

            for (int j = 0; j < n; ++j)
                _U[j, 0] = A[j, 0];

            for (int i = 0; i < nw; ++i)
            {
                _d[i + n] = Q[i, i];

                for (int j = 0; j < n; ++j)
                    _U[j, i + n] = W[j, i];
            }

            int o = -1;

            for (int j = n - 1; j != (uint)o; --j)
            {
                double sigma = 0.0;

                for (int k = 0; k < _nn; ++k)
                {
                    _v[k] = _U[j, k];
                    _a[k] = _d[k] * _v[k];
                    sigma += _v[k] * _a[k];
                }

                _U[j, j] = sigma;

                if (j == 0 || sigma == 0.0) continue;

                double dinv = 1.0 / sigma;

                for (int k = 0; k < j; ++k)
                {
                    sigma = 0.0;

                    for (int i = 0; i < _nn; ++i)
                        sigma += _U[k, i] * _a[i];

                    sigma *= dinv;

                    for (int i = 0; i < _nn; ++i)
                        _U[k, i] -= sigma * _v[i];

                    _U[j, k] = sigma;
                }
            }

            _U = _U.Transpose();
        }

        private void MeasureUpdate(double dz, double r)
        {
            // dz = dz - Hdx
            for (int i = 0; i < n; ++i)
                dz -= _a[i] * _x[i];

            //a = U^T * a
            //d = D * U^T * a
            for (int i = n - 1; i > 0; --i)
            {
                for (int k = 0; k < i; ++k)
                    _a[i] += _U[k, i] * _a[k];

                _d[i] = _U[i, i] * _a[i];
            }

            _d[0] = _U[0, 0] * _a[0];

            double alpha = r + _d[0] * _a[0];
            double gamma = 1.0 / alpha;

            _U[0, 0] = r * gamma * _U[0, 0];

            for (int j = 1; j < n; ++j)
            {
                double beta = alpha;

                alpha += _d[j] * _a[j];

                double lambda = -_a[j] * gamma;

                gamma = 1.0 / alpha;
                _U[j, j] *= beta * gamma;

                for (int i = 0; i < j; ++i)
                {
                    beta = _U[i, j];
                    _U[i, j] = beta + _d[i] * lambda;
                    _d[i] += _d[j] * beta;
                }
            }

            dz *= gamma;

            // dx = dx + K(dz - Hdx)
            for (int j = 0; j < n; ++j)
                _x[j] += _d[j] * dz;
        }

        private void MakeBaseAImpl()
        {
            _modified = true;
            MakeBaseA();

            if (_modified) _flags |= KALMAN.A_MODIF;
        }

        private void MakeBaseWImpl()
        {
            _modified = true;
            MakeBaseW();

            if (_modified) _flags |= KALMAN.W_MODIF;
        }

        private void MakeBaseQImpl()
        {
            _modified = true;
            MakeBaseQ();

            if (_modified) _flags |= KALMAN.Q_MODIF;
        }

        private void MakeBaseHImpl()
        {
            _modified = true;
            MakeBaseH();

            if (_modified) _flags |= KALMAN.H_MODIF;
        }

        private void MakeBaseVImpl()
        {
            _modified = true;
            MakeBaseV();

            if (_modified) _flags |= KALMAN.V_MODIF;
        }

        private void MakeBaseRImpl()
        {
            _modified = true;
            MakeBaseR();

            if (_modified) _flags |= KALMAN.R_MODIF;
        }

        private void MakeAImpl()
        {
            _modified = true;
            MakeA();

            if (_modified) _flags |= KALMAN.A_MODIF;
        }

        private void MakeWImpl()
        {
            _modified = true;
            MakeW();

            if (_modified) _flags |= KALMAN.W_MODIF;
        }

        private void MakeQImpl()
        {
            _modified = true;
            MakeQ();

            if (_modified) _flags |= KALMAN.Q_MODIF;
        }

        private void MakeHImpl()
        {
            _modified = true;
            MakeH();

            if (_modified) _flags |= KALMAN.H_MODIF;
        }

        private void MakeVImpl()
        {
            _modified = true;
            MakeV();

            if (_modified) _flags |= KALMAN.V_MODIF;
        }

        private void MakeRImpl()
        {
            _modified = true;
            MakeR();

            if (_modified) _flags |= KALMAN.R_MODIF;
        }
        #endregion

        #region Static functions
        //верхняя треугольная матрица Холески факторизация UDU
        private static void Factor(double[,] p)
        {
            int n = p.GetLength(0);

            for (int i = n - 1; i > 0; --i)
            {
                double alfa = 1.0 / p[i, i];

                for (int j = 0; j < i; ++j)
                {
                    double beta = p[j, i];

                    p[j, i] = alfa * beta;

                    for (int k = 0; k <= j; ++k)
                        p[k, j] -= beta * p[k, i];
                }
            }
        }

        //верхняя треугольная матрица с единичной диагональю
        private static void UpperInvert(double[,] p)
        {
            int n = p.GetLength(0);
            int o = -1;

            for (int i = n - 2; i != (uint)o; --i)
            {
                for (int k = i + 1; k < n; ++k)
                {

                    double val = p[i, k];

                    for (int j = i + 1; j <= k - 1; ++j)
                        val += p[i, j] * p[k, j];

                    p[k, i] = -val;
                }
            }
        }
        #endregion

        #region Properties
        public int SizeX
        {
            get { return this.n; }
            set
            {
                if (value != this.n)
                {
                    _flags |= KALMAN.N_MODIF;
                    this.n = value;
                }
            }
        }

        public int SizeU
        {
            get { return nu; }
            set
            {
                if (value != this.nu)
                {
                    _flags |= KALMAN.NU_MODIF;
                    this.nu = value;
                }
            }
        }

        public int SizeW
        {
            get { return nw; }
            set
            {
                if (value != this.nw)
                {
                    _flags |= KALMAN.NW_MODIF;
                    this.nw = value;
                }
            }
        }

        public int SizeZ
        {
            get { return m; }
            set
            {
                if (value != this.m)
                {
                    _flags |= KALMAN.M_MODIF;
                    this.m = value;
                }
            }
        }

        public int SizeV
        {
            get { return nv; }
            set
            {
                if (value != this.nv)
                {
                    _flags |= KALMAN.NV_MODIF;
                    this.nv = value;
                }
            }
        }

        public void SetDim(int n, int nu, int nw, int m, int nv)
        {
            SizeX = n;
            SizeU = nu;
            SizeW = nw;
            SizeZ = m;
            SizeV = nv;
        }
        #endregion
    }
}
