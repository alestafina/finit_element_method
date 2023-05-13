#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

struct matrix {
   vector<int> ig, jg;
   vector<double> al, au, di;
};

class FEM {
protected:
   matrix A;                        // итоговая матрица
   vector<vector<double>> loc_M;    // локальная матрица масс
   vector<vector<double>> loc_G;    // локальная матрица жесткости
   vector<double> loc_f;            // локальный вектор правой части
   vector<double> F;                // вектор правой части
   vector<double> q;                // вектор весов (ответ)
   vector<double> q_real;           // ожидаемые значения весов

   int hx, hy;
   vector<double> x_nodes, y_nodes; // сетка по x и y
   vector<int> nodes;               // глобальные номера узлов 
   vector<vector<int>> W;           // хранит "адреса" границ каждого КЭ

   vector<int> first_bc, second_bc, thrid_bc; // узлы, в которых выполнены ку

   int N, M, nx, ny;                // кол-во конечных элементов, кол-во узлов
   double betta, lambda, gamma;
public:
   FEM();
   ~FEM();

   void read_data();                // ввод данных
   void making_grid();              // составление стеки
   void build_portrait();           // составление портрета матрицы по сетке
   void glob_matrix();              // сборка глобальной матрицы
   void add_elem(int i, int j, double elem); // добавление элемента в глобальную матрицу
   
   void read_boundary();            // ввод краевых условий
   void thrid_cond();               // учет третьего краевого условия
   void second_cond();              // учет второго краевого условия
   void first_cond();               // учет первого краевого условия

   void print_result();
};

class Functions {
public:
   double func(double x, double y);
   double u_g(double x, double y);
   double u_betta(double x, double y);
   double theta(double x, double y);
   double u_real(double x, double y);
};

class CGM : public FEM {
private:
   vector<double> r, z, buf1, buf2;
   vector<double> al_LU, au_LU, di_LU;
   int maxiter, k;
   double eps;
   double residual, norm;
public: 
   CGM();
   ~CGM();

   void ILU_precond(); // неполная LU факторизация
   
   void set_vector();
   vector<double> matr_vec_mult(vector<double> &x, bool flag);
   double dot_product(vector<double> &a, vector<double> &b);

   vector<double> direct(vector<double> &L, vector<double> &b);
   vector<double> direct(vector<double> &L, vector<double> &D, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &D, vector<double> &b);
   
   void CGM_precond_ILU();
};
