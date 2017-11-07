namespace nt.math.statistics.filters
{
    using System;

    public class KalmanFilter
    {
        private int _r; // количество строк в S
        private int _n; // количество столбцов в y (число наблюдений)
        private int _k; // количество столбцов в A (число эксогенных переменных)
        private int _p; // длина комбинированных возмущений вектора
        private int _T; // количество строк в y (число наблюдений)

        public double[] S;  // остаточная ковариация
        public double[,] P; // предсказание оценки коварияации
        public double[,] F; // Якобиан частные производные df/dx
        public double[,] A; // матрица коэффициентов при эксогенных переменных
        public double[,] H; // Якобиан частные производные dh/dx
        public double[,] Q; // матрица ошибок состояний
        public double[,] R; // матрица ошибок при наблюдении

        public double[,] y; // матрица остаточных измерений
        public double[,] x; // матрица эксогенных переменных
        public double[] m;  // вектор переходов состояний

        public double[,] E; // матрица ошибок прогноза

        //public void Initialize
    }
}
