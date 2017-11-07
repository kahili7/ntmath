namespace nt.math.statistics
{
    using System;
    using nt.math;

    public class TimeSeries
    {
        public double[,] _Data;
        public bool[,] _Missing;

        private bool _AnyMissing;

        private int _N;
        private int _NComponents;
        private int _CurrentComponent;

        private int _NAllocated;
        public TimeSeries()
        {
            _AnyMissing = false;
            _N = _CurrentComponent = _NAllocated = 0;
            _NComponents = 1;

            _Data = new double[_NComponents, 1];
            _Missing = new bool[_NComponents, 1];
        }

        public bool HasAnyMissing
        {
            get { return _AnyMissing; }
            set { _AnyMissing = value; }
        }

        public void Append(double x, bool miss)
        {
            if(_NComponents > 1)
                throw new ArgumentOutOfRangeException("_NComponents", _NComponents, " > 1");

            if(_NAllocated < _N + 1)
            {
                _NAllocated += 512;

                _Data = null;
                _Data = new double[_NComponents, _NAllocated];

                _Missing = null;
                _Missing = new bool[_NComponents, _NAllocated];
            }

            _Data[0, _N] = x;
            _Missing[0, _N++] = miss;
            _AnyMissing |= miss;
        }

        public void Append(double[] x, bool[] miss)
        {
            if (_NComponents != x.GetLength(0))
                throw new ArgumentOutOfRangeException("_NComponents", _NComponents, " != x.Count");

            if (_NAllocated < _N + 1)
            {
                _NAllocated += 512;

                _Data = null;
                _Data = new double[_NComponents, _NAllocated];

                _Missing = null;
                _Missing = new bool[_NComponents, _NAllocated];
            }

            for (int i = 0; i < x.GetLength(0); i++)
            {
                _Data[i, _N] = x[i];
                _Missing[i, _N] = miss[i];

                if(miss[i] == true) 
                    _AnyMissing = (_AnyMissing != true) ? true : false;
                else
                    _AnyMissing = (_AnyMissing != false) ? false : true;
            }

            ++_N;
        }

        #region Get && Set
        public double[] GetDataVector()
        {
            double[] tmp = new double[_N];

            for (int i = 0; i < _N; i++)
            {
                tmp[i] = _Data[_CurrentComponent, i];
            }

            return tmp;
        }

        public double[,] GetDataMatrix()
        {
            double[,] tmp = new double[_N, _NComponents];

            for (int i = 0; i < _N; i++)
            {
                for (int j = 0; j < _NComponents; j++)
                    tmp[i, j] = _Data[j, i];
            }

            return tmp;
        }

        public bool[] GetMissingVector()
        {
            bool[] tmp = new bool[_N];

            for (int i = 0; i < _N; i++)
            {
                tmp[i] = _Missing[_CurrentComponent, i];
            }

            return tmp;
        }

        public bool[,] GetMissingMatrix()
        {
            bool[,] tmp = new bool[_N, _NComponents];

            for (int i = 0; i < _N; i++)
            {
                for (int j = 0; j < _NComponents; j++)
                    tmp[i, j] = _Missing[j, i];
            }

            return tmp;
        }

        public void SetMissingVector(bool[] miss)
        {
            _AnyMissing = false;

            for (int i = 0; i < _N; i++)
            {
                _Missing[_CurrentComponent, i] = (miss[i] != false);
                _AnyMissing |= (miss[i] != false);
            }
        }

        public void SetMissingVector(int i, bool miss)
        {
            _Missing[_CurrentComponent, i] = miss;
            _AnyMissing |= miss;
        }
        #endregion

        // среднее значение
        public double SampleMean()
        {
            double total = 0.0;
            int m = 0;

            for(int i = 0; i < _N; i++)
            {
                if(!_Missing[_CurrentComponent, i])
                {
                    total += _Data[_CurrentComponent, i];
                    m++;
                }
            }

            if (m > 0) total /= m;

            return total;
        }

        // минимальное значение
        public double SampleMin()
        {
            double min = Double.MinValue;

            for(int i = 0; i < _N; i++)
            {
                if (!_Missing[_CurrentComponent, i])
                {
                    if (_Data[_CurrentComponent, i] < min)
                    {
                        min = _Data[_CurrentComponent, i];
                    }
                }
            }

            return min;
        }

        // максимальное значение
        public double SampleMax()
        {
            double max = Double.MaxValue;

            for (int i = 0; i < _N; i++)
            {
                if (!_Missing[_CurrentComponent, i])
                {
                    if (_Data[_CurrentComponent, i] > max)
                    {
                        max = _Data[_CurrentComponent, i];
                    }
                }
            }

            return max;
        }

        public void ACFtoPACF(double[] acoff, double[] pacf_storage)
        {
            int maxlag = acoff.GetLength(0) - 1;

            if(maxlag <= 0)
            {
                pacf_storage[0] = 0.0;
                return;
            }

            // рекурсивный алгоритм Левинсона-Дурбина
            double[] phis = new double[maxlag];
            double[] phis2 = new double[maxlag];

            phis[0] = acoff[1] / acoff[0];
            pacf_storage[0] = 1.0;
            pacf_storage[1] = phis[0];

            double vi = acoff[0];
            double phinn = 0.0;

            vi *= (1 - phis[0] * phis[0]);

            for(int i = 2; i <= maxlag; i++)
            {
                for (int j = 0; j < i - 1; j++)
                    phis2[j] = phis2[i - j - 2];

                phinn = acoff[i];

                for (int j = 1; j < i; j++)
                    phinn -= phis[j - 1] * acoff[i - j];

                phinn /= vi;

                for (int j = 0; j < i - 1; j++)
                    phis[j] -= phinn * phis2[j];

                vi *= (1 - phinn * phinn);
                pacf_storage[i] = phis[i - 1] = phinn;
            }
        }

        public unsafe void ComputeSampleACF(double[] acoff, double[] pacf_storage, bool normalize)
        {
            if (_N == 0) return;

            int maxlag = acoff.GetLength(0) - 1;
            double mean = SampleMean();
            
            fixed(double* tdata = &_Data[_CurrentComponent, 0])
            fixed(bool* tmiss = &_Missing[_CurrentComponent, 0])
            // вычисляем автоковариационную функцию
            if(!HasAnyMissing)
            {
                for (int i = 0; i <= maxlag; i++)
                {
                    double total = 0.0;

                    for (int j = i; j < _N; j++)
                        total += (tdata[j] - mean) * (tdata[j - i] - mean);

                    acoff[i] = total / _N;
                }
            }
            else
            {
                double nonmissing = 0.0;

                for (int i = 0; i < _N; i++)
                    nonmissing += Math.Floor(1.5 - (tmiss[i] ? 1.0 : 0.0));

                for (int i = 0; i <= maxlag; i++)
                {
                    double total = 0.0;

                    for (int j = i; j < _N; j++)
                        if ((tmiss[j] == false) && (tmiss[j - i] == false))
                            total += (tdata[j] - mean) * (tdata[j - i] - mean);

                    acoff[i] = total / nonmissing;
                }
            }

            if (pacf_storage != null)
                ACFtoPACF(acoff, pacf_storage);

            if (normalize)
                for (int i = maxlag; i >= 0; i--)
                    acoff[i] /= acoff[0];
        }
    }
}
