namespace nt.math.statistics
{
    using System;
    using nt.math;

    public class ARMAModel : TimeSeriesModel
    {
        private int _p;
        private int _q;
        private double _sigma;
        private double _mean;
        private double _fracdiff;

        private double[] _phi;
        private double[] _theta;
        public static int MaxOrder = 100; 

        // свести вероятность
        private bool _usingWhittle;
        private bool _whittleCache;
        private TimeSeries _whittleTS;

        public ARMAModel(int maxOrder) : base()
        {
            _p = _q = 0;
            _sigma = 1.0;
            _fracdiff = 0.0;

            _phi = Matrix.Vector<double>(maxOrder, 0.0);
            _theta = Matrix.Vector<double>(maxOrder, 0.0);

            Estimation_Mask = new bool[3];
            Estimation_Mask[0] = true;
            Estimation_Mask[1] = false;
            Estimation_Mask[2] = true;

            _usingWhittle = false;
            _whittleCache = false;
        }

        public ARMAModel() : this(ARMAModel.MaxOrder)
        { }

        public double[] Bundle
        {
            get 
            {
                double[] pvec = new double[_p + _q + 3];

                for (int i = 0; i < _p; i++)
                    pvec[i + 3] = _phi[i];

                for (int i = 0; i < _q; i++)
                    pvec[i + _p + 3] = _theta[i];

                pvec[0] = Mean;
                pvec[1] = FracDiff;
                pvec[2] = Sigma;
                return pvec;
            }

            set
            {
                Mean = value[0];
                FracDiff = value[1];
                Sigma2 = (value[2] * value[2]);

                double[] phi = new double[_p];
                double[] theta = new double[_q];

                phi = Matrix.Submatrix<double>(value, 3);
                theta = Matrix.Submatrix<double>(value, _p + 3);
                SetCoeffs(phi, theta);
            }
        }

        public void GetCoeffs(double[] phis, double[] thetas)
        {
            phis = new double[_p];
            thetas = new double[_q];

            for (int i = 0; i < _p; i++)
                phis[i] = _phi[i];

            for (int i = 0; i < _q; i++)
                thetas[i] = _theta[i];
        }

        public void SetCoeffs(double[] phis, double[] thetas)
        {
            for (int i = 0; i < _p; i++)
                _phi[i] = phis[i];

            for (int i = 0; i < _q; i++)
                _theta[i] = thetas[i];
        }

        public bool[] GetMask()
        {
            int i;
            int np = _p + _q + 3;
            bool[] res = new bool[np];

            for (i = 0; i < Math.Min(Estimation_Mask.GetLength(0), np); i++)
                res[i] = Estimation_Mask[i];

            for (; i < np; i++)
            {
                switch (i)
                {
                    case 0:
                    case 2: res[i] = true; break;
                    case 1: res[i] = false; break;
                    default: res[i] = true; break;
                }
            }

            return res;
        }

        #region Properties
        public double Mean
        {
            get { return _mean; }
            set { _mean = value; }
        }

        public int P
        {
            get { return _p; }
            set { _p = value; }
        }

        public int Q
        {
            get { return _q; }
            set { _q = value; }
        }

        public double Sigma
        {
            get { return _sigma; }
            set { _sigma = value; }
        }

        public double Sigma2
        {
            get { return _sigma; }
            set { _sigma = Math.Sqrt(value); }
        }

        public double FracDiff
        {
            get { return _fracdiff; }
            set { _fracdiff = value; }
        }

        public bool Whittle
        {
            get { return _usingWhittle; }
            set { _usingWhittle = value; }
        }
        #endregion
    }
}
