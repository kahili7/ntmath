namespace nt.math.statistics
{
    using System;
    using nt.math;

    public class TimeSeriesModel
    {
        protected double LogLikelihood; // логарифм правдоподобия
        protected double AICC;
        protected bool LogLikelihood_Is_Valid;
        protected bool[] Estimation_Mask; // оценочная маска

        public TimeSeriesModel()
        {
            LogLikelihood = 0.0;
            AICC = 0.0;
            LogLikelihood_Is_Valid = false;
            Estimation_Mask = null;
        }

        #region Static functions
        public static double[] ParamMask(double[] parms, double[] mask)
        {
            int sz = (int)(mask.Sum() + 0.5);
            double[] sub = new double[sz];
            int j = 0;

            for (int i = 0; i < parms.GetLength(0); i++)
                if (mask[i] == 1.0)
                    sub[j++] = parms[i];

            return sub;
        }

        public static double[] ParamExpand(double[] x, double[] curparms, double[] mask)
        {
            double[] res = curparms;
            int j = 0;

            for (int i = 0; i < curparms.GetLength(0); i++)
                if (mask[i] == 1.0)
                    res[i] = x[j++];

            return res;
        }
        #endregion

    }
}
